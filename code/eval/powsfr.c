//
//  powsfr.c
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

#include "powsfr.h"

#ifdef POWSFR_DIRECT_POWERS

// Here is the simple version that directly computes the powers

// MARK: Implementation of public functions

powsfr powsfr_new(deg_t size) {
    powsfr p = malloc(sizeof(powsfr_struct));
    p->inited = false;
    
    return p;
}

bool powsfr_free(powsfr xn) {
    if(xn == NULL) {
        return false;
    }
    
    free(xn);
    
    return true;
}

bool powsfr_set(powsfr xn, coeff_t x) {
    if(xn == NULL || isnan(x)) {
        return false;
    }
    
    xn->inited = true;
    xn->x = x;
    
    return true;
}

inline coeff_t powsfr_pow(powsfr xn, deg_t p) {
    return p == 1 ? xn->x : fpow(xn->x, p);
}

inline coeff_t powsfr_pow_once(powsfr xn, deg_t p) {
    return p == 1 ? xn->x : fpow(xn->x, p);
}

#else

// Here is the full version that computes the powers using a buffered version of the
// repeated squares method

// MARK: Static functions

/// @brief Ensures that @c b has @c size space.
///
/// @param b the powers object
/// @param size the desired minimum size
static inline void resize(powsfr b, deg_t size) {
    b->xn = realloc(b->xn, sizeof(coeff_t) * size);
    b->computed = realloc(b->computed, sizeof(bool) * size);
    
    coeff_t *a = b->xn;
    bool *c = b->computed;
    for (deg_t i = b->size; i < size; i++) {
        a[i] = 0;
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

/// @brief Retrieves or computes the power @c 2^tp in the buffer @c xn.
///
/// @param xn the powers buffer
/// @param tp the power of two
/// @param save @ref true to save intermediary results in @c xn, @ref false otherwise
///
/// @return the power @c 2^tp of the base of @c xn.
static inline coeff_t powsfr_dyadic(powsfr xn, byte tp, bool save) {
    deg_t pw = 1U << tp;
    bool *c = xn->computed;
    coeff_t *v = xn->xn;
    
    if(c[pw - 1]) {
        return v[pw - 1];
    }
    
    deg_t ip = pw >> 1;
    while(! c[ip - 1]) {
        ip >>= 1;
    }
    
    coeff_t r = v[ip - 1];
    do {
        r *= r;
        ip <<= 1;
        
        if(save) {
            v[ip - 1] = r;
            c[ip - 1] = true;
        }
    } while(ip < pw);
    
    return r;
}

// MARK: Implementation of public functions

powsfr powsfr_new(deg_t size) {
    if(size > MAX_DEG || size == 0) {
        return NULL;
    }
    
    powsfr p = malloc(sizeof(powsfr_struct));
    p->size = size;
    
    p->xn = malloc(sizeof(coeff_t) * size);
    p->computed = malloc(sizeof(bool) * size);
    coeff_t *a = p->xn;
    bool *c = p->computed;
    for (deg_t i = 0; i < size; i++) {
        a[i] = 0;
        c[i] = false;
    }
    
    // compute tps
    byte tp = 0;
    while(size >> tp > 1) {
        tp ++;
    }
    p->tps = tp;
    
    return p;
}

bool powsfr_free(powsfr xn) {
    if(xn == NULL || xn->xn == NULL || xn->computed == NULL) {
        return false;
    }
    
    free(xn->computed);
    free(xn->xn);
    free(xn);
    
    return true;
}

bool powsfr_set(powsfr xn, coeff_t x) {
    if(xn == NULL || isnan(x)) {
        return false;
    }
    
    bool *c = xn->computed;
    memset(c + 1, false, zn->size - 1);
    
    xn->xn[0] = x;
    c[0] = true;
    
    return true;
}

coeff_t powsfr_pow(powsfr xn, deg_t pow) {
    if(xn == NULL || pow > MAX_DEG || ! xn->computed[0]) {
        return NAN;
    }
    
    if(pow == 0) {
        return 1;
    }
    
    deg_t n = xn->size;
    if(pow > n) {
        long double ns = pow;
        ns *= 1.1;
        ns += 10;
        deg_t size = ns > MAX_DEG ? MAX_DEG : (deg_t) ns;
        
        resize(xn, size);
    }
    
    if(xn->computed[pow - 1]) {
        return xn->xn[pow - 1];
    }
    
    // compute and store the new value(s)
    byte tl = 1, th = xn->tps + 1, t = (tl + th) / 2;
    deg_t ip = pow & ((1L << t) - 1), lip = pow;
    bool ph = true;
    bool *c = xn->computed;
    coeff_t *x = xn->xn;
    
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
    if(ip > 0 && pow - ip <= n && pow - ip <= n && c[pow - ip - 1]) { // cheapest route, just one product
        coeff_t res = x[ip - 1] * x[pow - ip - 1];
        x[pow - 1] = res;
        c[pow - 1] = true;
        
        return res;
    }

    coeff_t pth = powsfr_dyadic(xn, th - 1, true);
    coeff_t res = ip > 0 ? x[ip - 1] : 1;
    deg_t tth = (deg_t) (1U << (th - 1));
    
    while(pow > ip + tth) {
        if(pow & tth) {
            res *= pth;
            ip += tth;
            
            x[ip - 1] = res;
            c[ip - 1] = true;
        }
        
        // repeated squares method
        pth *= pth;
        tth *= 2;
        
        if(! c[tth - 1]) {
            x[tth - 1] = pth;
            c[tth - 1] = true;
        }
    }
    
    // here pow - ip == tth
    res *= pth;
    x[pow - 1] = res;
    c[pow - 1] = true;
    
    return res;
}

coeff_t powsfr_pow_once(powsfr xn, deg_t pow) {
    if(xn == NULL || pow > MAX_DEG || ! xn->computed[0]) {
        return NAN;
    }
    
    if(pow == 0) {
        return 1;
    }
    
    deg_t n = xn->size;
    if(pow <= n && xn->computed[pow - 1]) {
        return xn->xn[pow - 1];
    }
    
    // compute the new value
    byte tl = 1, th = xn->tps + 1, t = (tl + th) / 2;
    deg_t ip = pow & ((1L << t) - 1), lip = pow;
    bool ph = true;
    bool *c = xn->computed;
    coeff_t *x = xn->xn;
    
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
        return x[ip - 1] * x[pow - ip - 1];
    }

    coeff_t pth = powsfr_dyadic(xn, th - 1, false);
    coeff_t res = ip > 0 ? x[ip - 1] : 1;
    deg_t tth = (deg_t) (1U << (th - 1));
    
    while(pow > ip + tth) {
        if(pow & tth) {
            res *= pth;
            ip += tth;
        }
        
        // repeated squares method
        pth *= pth;
        tth *= 2;
    }
    
    // here pow - ip == tth
    res *= pth;
    
    return res;
}

#endif
