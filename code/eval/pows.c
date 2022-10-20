//
//  pows.c
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

#include "pows.h"

// MARK: Static functions

/// @brief Ensures that @c b has @c size space.
///
/// @param b the powers object
/// @param size the desired minimum size
static inline void resize(pows b, deg_t size) {
    b->zn = realloc(b->zn, sizeof(comp_struct) * size);
    b->computed = realloc(b->computed, sizeof(bool) * size);
    
    comp_ptr a = b->zn;
    bool *c = b->computed;
    
    for (deg_t i = b->size; i < size; i++) {
        comp_initz(a + i, b->prec);
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
static inline void pows_dyadic(comp res, pows zn, byte tp, bool save) {
    deg_t pw = 1U << tp;
    bool *c = zn->computed;
    comp_ptr v = zn->zn;
    
    if(c[pw - 1]) {
        comp_set(res, v + (pw - 1));
        
        return;
    }
    
    deg_t ip = pw >> 1;
    while(! c[ip - 1]) {
        ip >>= 1;
    }

    comp_set(res, v + (ip - 1));
    do {
        comp_sqr(res, res, zn->buf1);
        ip <<= 1;
        
        if(save) {
            comp_set(v + (ip - 1), res);
            c[ip - 1] = true;
        }
    } while(ip < pw);
}

// MARK: Implementation of public functions

pows pows_new(prec_t prec, deg_t size) {
    if(size > MAX_DEG || size == 0 || prec < precf) {
        return NULL;
    }
    
    pows p = malloc(sizeof(pows_struct));
    p->size = size;
    p->prec = prec;
    
    p->zn = malloc(sizeof(comp_struct) * size);
    p->computed = malloc(sizeof(bool) * size);
    comp_ptr a = p->zn;
    bool *c = p->computed;
    
    for (deg_t i = 0; i < size; i++) {
        comp_initz(a + i, prec);
        c[i] = false;
    }
    
    // compute tps
    byte tp = 0;
    while(size >> tp > 1) {
        tp ++;
    }
    p->tps = tp;
    
    mpfr_init2(p->buf1, prec);
    mpfr_init2(p->buf2, prec);
    
    comp_init(p->pth, prec);
    comp_init(p->res, prec);
    
    return p;
}

bool pows_free(pows zn) {
    if(zn == NULL || zn->zn == NULL || zn->computed == NULL) {
        return false;
    }
    
    mpfr_clear(zn->buf1);
    mpfr_clear(zn->buf2);
    
    comp_clear(zn->pth);
    comp_clear(zn->res);
    
    free(zn->computed);
    free(zn->zn);
    free(zn);
    
    return true;
}

bool pows_set(pows zn, comp z) {
    if(zn == NULL || z == NULL) {
        return false;
    }
    
    bool *c = zn->computed;
    memset(c + 1, false, zn->size - 1);
    
    comp_set(zn->zn, z);
    c[0] = true;
    
    return true;
}

comp_ptr pows_pow(pows zn, deg_t pow) {
    if(zn == NULL || pow > MAX_DEG || ! zn->computed[0] || pow == 0) {
        return NULL;
    }
    
    deg_t n = zn->size;
    if(pow > n) {
        long double ns = pow;
        ns *= 1.1;
        ns += 10;
        deg_t size = ns > MAX_DEG ? MAX_DEG : (deg_t) ns;
        
        resize(zn, size);
    }
    
    if(pow <= n && zn->computed[pow - 1]) {
        return zn->zn + (pow - 1);
    }
    
    // compute and store the new value(s)
    byte tl = 1, th = zn->tps + 1, t = (tl + th) / 2;
    deg_t ip = pow & ((1L << t) - 1), lip = pow;
    bool ph = true;
    bool *c = zn->computed;
    comp_ptr x = zn->zn;
    
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
        comp_mul(x + (pow - 1), x + (ip - 1), x + (pow - ip - 1), zn->buf1, zn->buf2);
        c[pow - 1] = true;
        
        return x + (pow - 1);
    }

    comp_ptr pth = zn->pth;
    pows_dyadic(pth, zn, th - 1, true);
    
    comp_ptr res = zn->res;
    if(ip > 0) {
        comp_set(res, x + (ip - 1));
    } else {
        mpfr_set_ui(res->x, 1, MPFR_RNDN);
        mpfr_set_zero(res->y, 1);
    }
    
    deg_t tth = (deg_t) (1U << (th - 1));
    
    while(pow > ip + tth) {
        if(pow & tth) {
            comp_mul(res, res, pth, zn->buf1, zn->buf2);
            ip += tth;

            comp_set(x + (ip - 1), res);
            c[ip - 1] = true;
        }
        
        // repeated squares method
        comp_sqr(pth, pth, zn->buf1);
        tth *= 2;
        
        if(! c[tth - 1]) {
            comp_set(x + (tth - 1), pth);
            c[tth - 1] = true;
        }
    }
    
    // here pow - ip == tth
    comp_mul(res, res, pth, zn->buf1, zn->buf2);
    comp_set(x + (pow - 1), res);
    c[pow - 1] = true;
    
    return x + (pow - 1);
}

comp_ptr pows_pow_once(pows zn, deg_t pow) {
    if(zn == NULL || pow > MAX_DEG || ! zn->computed[0] || pow == 0) {
        return NULL;
    }
    
    deg_t n = zn->size;
    if(pow <= n && zn->computed[pow - 1]) {
        return zn->zn + (pow - 1);
    }
    
    // compute the new value
    byte tl = 1, th = zn->tps + 1, t = (tl + th) / 2;
    deg_t ip = pow & ((1L << t) - 1), lip = pow;
    bool ph = true;
    bool *c = zn->computed;
    comp_ptr x = zn->zn;
    comp_ptr res = zn->res;
    
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
        comp_mul(res, x + (ip - 1), x + (pow - ip - 1), zn->buf1, zn->buf2);
        
        return res;
    }
    
    comp_ptr pth = zn->pth;
    pows_dyadic(pth, zn, th - 1, true);
    
    if(ip > 0) {
        comp_set(res, x + (ip - 1));
    } else {
        mpfr_set_ui(res->x, 1, MPFR_RNDN);
        mpfr_set_zero(res->y, 1);
    }
    
    deg_t tth = (deg_t) (1U << (th - 1));
    
    while(pow > ip + tth) {
        if(pow & tth) {
            comp_mul(res, res, pth, zn->buf1, zn->buf2);
            ip += tth;
        }
        
        // repeated squares method
        comp_sqr(pth, pth, zn->buf1);
        tth *= 2;
    }
    
    // here pow - ip == tth
    comp_mul(res, res, pth, zn->buf1, zn->buf2);
    
    return res;
}
