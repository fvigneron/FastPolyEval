//
//  powsr.c
//
//  Authors: Nicolae Mihalache & François Vigneron
//
//  This software is released under BSD licence, with an attribution clause (see Licence file).
//  Please cite the reference below if you use or distribute this software.
//
//  • [1] R. Anton, N. Mihalache & F. Vigneron. Fast evaluation of real and complex polynomials. 2022.
//        https://hal.archives-ouvertes.fr/hal-03820369
//
//  Copyright 2022 Univ. Paris-Est Créteil & Univ. de Reims Champagne-Ardenne.
//

#include <stdlib.h>
#include <string.h>

#include "powsr.h"


#ifdef POWSR_DIRECT_POWERS

// Here is the simple version that directly computes the powers

// MARK: Implementation of public functions

powsr powsr_new(prec_t prec, deg_t size) {
    if(prec < precf) {
        return NULL;
    }
    
    powsr p = malloc(sizeof(powsr_struct));
    p->prec = prec;
    
    mpfr_init2(p->x, prec);
    mpfr_init2(p->res, prec);
    
    return p;
}

bool powsr_free(powsr xn) {
    if(xn == NULL) {
        return false;
    }

    mpfr_clear(xn->x);
    mpfr_clear(xn->res);
    
    free(xn);
    
    return true;
}

bool powsr_set(powsr xn, mpfr_t x) {
    if(xn == NULL || x == NULL) {
        return false;
    }
    
    xn->inited = true;
    mpfr_set(xn->x, x, MPFR_RNDN);
    
    return true;
}

inline mpfr_ptr powsr_pow(powsr xn, deg_t pow) {
    if(pow == 1) {
        return xn->x;
    }
    
    mpfr_pow_ui(xn->res, xn->x, pow, MPFR_RNDN);
    
    return xn->res;
}

inline mpfr_ptr powsr_pow_once(powsr xn, deg_t pow) {
    if(pow == 1) {
        return xn->x;
    }
    
    mpfr_pow_ui(xn->res, xn->x, pow, MPFR_RNDN);
    
    return xn->res;
}

#else

// Here is the full version that computes the powers using a buffered version of the
// repeated squares method

// MARK: Static functions

/// @brief Ensures that @c b has @c size space.
///
/// @param b the powers object
/// @param size the desired minimum size
static inline void resize(powsr b, deg_t size) {
    b->xn = realloc(b->xn, sizeof(__mpfr_struct) * size);
    b->computed = realloc(b->computed, sizeof(bool) * size);
    
    mpfr_ptr a = b->xn;
    bool *c = b->computed;
    
    for (deg_t i = b->size; i < size; i++) {
        mpfr_init2(a + i, b->prec);
        mpfr_set_zero(a + i, 1);
        c[i] = false;
    }
    
    b->size = size;
    
    // compute tps
    byte tp = 0;
    while(size >> tp > 1) {
        tp ++;
    }
    b->tps = tp;
}

/// @brief Retrieves or computes the power @c 2^tp in the buffer @c zn.
///
/// @param zn the powers buffer
/// @param tp the power of two
/// @param save @ref true to save intermediary results in @c zn, @ref false otherwise
///
/// @return the power @c 2^tp of the base of @c zn.
static inline void pows_dyadic(mpfr_t res, powsr xn, byte tp, bool save) {
    deg_t pw = 1U << tp;
    bool *c = xn->computed;
    mpfr_ptr v = xn->xn;
    
    if(c[pw - 1]) {
        mpfr_set(res, v + (pw - 1), MPFR_RNDN);
        
        return;
    }
    
    deg_t ip = pw >> 1;
    while(! c[ip - 1]) {
        ip >>= 1;
    }

    mpfr_set(res, v + (ip - 1), MPFR_RNDN);
    do {
        mpfr_sqr(res, res, MPFR_RNDN);
        ip <<= 1;
        
        if(save) {
            mpfr_set(v + (ip - 1), res, MPFR_RNDN);
            c[ip - 1] = true;
        }
    } while(ip < pw);
}

// MARK: Implementation of public functions

powsr powsr_new(uint prec, deg_t size) {
    if(size > MAX_DEG || size == 0 || prec < precf) {
        return NULL;
    }
    
    powsr p = malloc(sizeof(powsr_struct));
    p->size = size;
    p->prec = prec;
    
    p->xn = malloc(sizeof(__mpfr_struct) * size);
    p->computed = malloc(sizeof(bool) * size);
    mpfr_ptr a = p->xn;
    bool *c = p->computed;
    
    for (deg_t i = 0; i < size; i++) {
        mpfr_init2(a + i, prec);
        mpfr_set_zero(a + i, 1);
        c[i] = false;
    }
    
    // compute tps
    byte tp = 0;
    while(size >> tp > 1) {
        tp ++;
    }
    p->tps = tp;
    
    mpfr_init2(p->pth, prec);
    mpfr_init2(p->res, prec);
    
    return p;
}

bool powsr_free(powsr xn) {
    if(xn == NULL || xn->xn == NULL || xn->computed == NULL) {
        return false;
    }

    mpfr_clear(xn->pth);
    mpfr_clear(xn->res);
    
    free(xn->computed);
    free(xn->xn);
    free(xn);
    
    return true;
}

bool powsr_set(powsr xn, mpfr_t x) {
    if(xn == NULL || x == NULL) {
        return false;
    }
    
    bool *c = xn->computed;
    memset(c + 1, false, zn->size - 1);
    
    mpfr_set(xn->xn, x, MPFR_RNDN);
    c[0] = true;
    
    return true;
}

mpfr_ptr powsr_pow(powsr xn, deg_t pow) {
    if(xn == NULL || pow > MAX_DEG || ! xn->computed[0] || pow == 0) {
        return NULL;
    }
    
    deg_t n = xn->size;
    if(pow > n) {
        long double ns = pow;
        ns *= 1.1;
        ns += 10;
        deg_t size = ns > MAX_DEG ? MAX_DEG : (deg_t) ns;
        
        resize(xn, size);
    }
    
    if(pow <= n && xn->computed[pow - 1]) {
        return xn->xn + (pow - 1);
    }
    
    // compute and store the new value(s)
    byte tl = 1, th = xn->tps + 1, t = (tl + th) / 2;
    deg_t ip = pow & ((1L << t) - 1), lip = pow;
    bool ph = true;
    bool *c = xn->computed;
    mpfr_ptr x = xn->xn;
    
    do {
        if(ip == lip) {
            if(ph) {
                th = t;
            } else {
                tl = t;
            }
            t = (tl + th) / 2;
            lip = ip;
            ip = pow & ((1L << t) - 1);
            
            if(th - tl == 1) {
                break;
            }
            
            continue;
        }
        
        if(c[ip - 1]) {
            ph = false;
            tl = t;
        } else {
            ph = true;
            th = t;
        }
        t = (tl + th) / 2;
        lip = ip;
        ip = pow & ((1L << t) - 1);
    } while(th - tl > 1);
    
    ip = pow & ((1U << tl) - 1);
    // here power ip pre-computed, tl is max with this property and th = tl + 1
    if(ip > 0 && pow - ip <= n && c[pow - ip - 1]) { // cheapest route, just one product
        mpfr_mul(x + (pow - 1), x + (ip - 1), x + (pow - ip - 1), MPFR_RNDN);
        c[pow - 1] = true;
        
        return x + (pow - 1);
    }

    mpfr_ptr pth = xn->pth;
    pows_dyadic(pth, xn, th - 1, true);
    
    mpfr_ptr res = xn->res;
    if(ip > 0) {
        mpfr_set(res, x + (ip - 1), MPFR_RNDN);
    } else {
        mpfr_set_ui(res->x, 1, MPFR_RNDN);
    }
        
    deg_t tth = (deg_t) (1U << (th - 1));
    
    while(pow > ip + tth) {
        if(pow & tth) {
            mpfr_mul(res, res, pth, MPFR_RNDN);
            ip += tth;

            mpfr_set(x + (ip - 1), res, MPFR_RNDN);
            c[ip - 1] = true;
        }
        
        // repeated squares method
        mpfr_sqr(pth, pth, MPFR_RNDN);
        tth *= 2;
        
        if(! c[tth - 1]) {
            mpfr_set(x + (tth - 1), pth, MPFR_RNDN);
            c[tth - 1] = true;
        }
    }
    
    // here pow - ip == tth
    mpfr_mul(res, res, pth, MPFR_RNDN);
    mpfr_set(x + (pow - 1), res, MPFR_RNDN);
    c[pow - 1] = true;
    
    return x + (pow - 1);
}

bool powsr_powOnce(mpfr_t res, powsr xn, deg_t pow) {
    if(res == NULL || xn == NULL || pow > MAX_DEG || ! xn->computed[0] || pow == 0) {
        return false;
    }
    
    deg_t n = xn->size;
    if(pow <= n && xn->computed[pow - 1]) {
        mpfr_set(res, xn->xn + (pow - 1), MPFR_RNDN);
        
        return true;
    }
    
    // compute the new value
    byte tl = 1, th = xn->tps + 1, t = (tl + th) / 2;
    deg_t ip = pow & ((1L << t) - 1), lip = pow;
    bool ph = true;
    bool *c = xn->computed;
    mpfr_ptr x = xn->xn;
    
    do {
        if(ip == lip) {
            if(ph) {
                th = t;
            } else {
                tl = t;
            }
            t = (tl + th) / 2;
            lip = ip;
            ip = pow & ((1L << t) - 1);
            
            if(th - tl == 1) {
                break;
            }
            
            continue;
        }
        
        if(c[ip - 1]) {
            ph = false;
            tl = t;
        } else {
            ph = true;
            th = t;
        }
        t = (tl + th) / 2;
        lip = ip;
        ip = pow & ((1L << t) - 1);
    } while(th - tl > 1);
    
    ip = pow & ((1U << tl) - 1);
    // here power ip pre-computed, tl is max with this property and th = tl + 1
    if(ip > 0 && pow - ip <= n && c[pow - ip - 1]) { // cheapest route, just one product
        mpfr_mul(res, x + (ip - 1), x + (pow - ip - 1), MPFR_RNDN);
        
        return true;
    }
    
    mpfr_ptr pth = xn->pth;
    pows_dyadic(pth, xn, th - 1, true);
    
    if(ip > 0) {
        mpfr_set(res, x + (ip - 1), MPFR_RNDN);
    } else {
        mpfr_set_ui(res->x, 1, MPFR_RNDN);
    }
    
    deg_t tth = (deg_t) (1U << (th - 1));
    
    while(pow > ip + tth) {
        if(pow & tth) {
            mpfr_mul(res, res, pth, MPFR_RNDN);
            ip += tth;
        }
        
        // repeated squares method
        mpfr_sqr(pth, pth, MPFR_RNDN);
        tth *= 2;
    }
    
    // here pow - ip == tth
    mpfr_mul(res, res, pth, MPFR_RNDN);
    
    return true;
}

#endif
