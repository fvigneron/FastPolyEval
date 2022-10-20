//
//  powsf.c
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

#include "powsf.h"

// MARK: Static functions

/// @brief Ensures that @c b has @c size space.
///
/// @param b the powers object
/// @param size the desired minimum size
static inline void resize(powsf b, deg_t size) {
    b->zn = realloc(b->zn, sizeof(compf_struct) * size);
    b->computed = realloc(b->computed, sizeof(bool) * size);
    
    compf_ptr a = b->zn;
    bool *c = b->computed;
    compf z = {{0, 0}};
    for (deg_t i = b->size; i < size; i++) {
        a[i] = *z;
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
static inline void powsf_dyadic(compf res, powsf zn, byte tp, bool save) {
    deg_t pw = 1U << tp;
    bool *c = zn->computed;
    compf_ptr v = zn->zn;
    
    if(c[pw - 1]) {
        *res = v[pw - 1];
        
        return;
    }
    
    deg_t ip = pw >> 1;
    while(! c[ip - 1]) {
        ip >>= 1;
    }
    
    *res = v[ip - 1];
    do {
        compf_sqr(res, res);
        ip <<= 1;
        
        if(save) {
            v[ip - 1] = *res;
            c[ip - 1] = true;
        }
    } while(ip < pw);
}

// MARK: Implementation of public functions

powsf powsf_new(deg_t size) {
    if(size > MAX_DEG || size == 0) {
        return NULL;
    }
    
    powsf p = malloc(sizeof(powsf_struct));
    p->size = size;
    
    p->zn = malloc(sizeof(compf_struct) * size);
    p->computed = malloc(sizeof(bool) * size);
    compf_ptr a = p->zn;
    bool *c = p->computed;
    compf z = {{0, 0}};
    
    for (deg_t i = 0; i < size; i++) {
        a[i] = *z;
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

bool powsf_free(powsf zn) {
    if(zn == NULL || zn->zn == NULL || zn->computed == NULL) {
        return false;
    }
    
    free(zn->computed);
    free(zn->zn);
    free(zn);
    
    return true;
}

bool powsf_set(powsf zn, compf z) {
    if(zn == NULL || z == NULL) {
        return false;
    }
    
    bool *c = zn->computed;
    memset(c + 1, false, zn->size - 1);
    
    zn->zn[0] = *z;
    c[0] = true;
    
    return true;
}

compf_ptr powsf_pow(powsf zn, deg_t pow) {
    if(zn == NULL || pow > MAX_DEG || ! zn->computed[0] || pow == 0) {
        return NULL;
    }
    
    if(pow > zn->size) {
        long double ns = pow;
        ns *= 1.1;
        ns += 10;
        deg_t size = ns > MAX_DEG ? MAX_DEG : (deg_t) ns;
        
        resize(zn, size);
    }
    
    if(zn->computed[pow - 1]) {
        return zn->zn + (pow - 1);
    }
    
    // compute and store the new value(s)
    byte tl = 1, th = zn->tps + 1, t = (tl + th) / 2;
    deg_t ip = pow & ((1L << t) - 1), lip = pow;
    bool ph = true;
    bool *c = zn->computed;
    compf_ptr x = zn->zn;
    
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
    // here if ip > 0, then power ip pre-computed, tl is max with this property and th = tl + 1
    if(ip > 0 && c[pow - ip - 1]) { // cheapest route, just one product
        compf_mul(x + (pow - 1), x + (ip - 1), x + (pow - ip - 1));
        c[pow - 1] = true;
        
        return x + (pow - 1);
    }

    compf pth;
    powsf_dyadic(pth, zn, th - 1, true);
    compf res;
    if(ip == 0) {
        res->x = 1;
        res->y = 0;
    } else {
        *res = x[ip - 1];
    }
    deg_t tth = (deg_t) (1U << (th - 1));
    
    while(pow > ip + tth) {
        if(pow & tth) {
            compf_mul(res, res, pth);
            ip += tth;
            
            x[ip - 1] = *res;
            c[ip - 1] = true;
        }
        
        // repeated squares method
        compf_sqr(pth, pth);
        tth *= 2;
        
        if(! c[tth - 1]) {
            x[tth - 1] = *pth;
            c[tth - 1] = true;
        }
    }
    
    // here pow - ip == tth
    compf_mul(res, res, pth);
    x[pow - 1] = *res;
    c[pow - 1] = true;
    
    return x + (pow - 1);
}

compf_ptr powsf_pow_once(powsf zn, deg_t pow) {
    if(zn == NULL || pow == 0 || pow > MAX_DEG || ! zn->computed[0]) {
        return NULL;
    }
    
    deg_t n = zn->size;
    if(pow <= n && zn->computed[pow - 1]) {
        return zn->zn  + (pow - 1);
    }
    
    // compute the new value
    byte tl = 1, th = zn->tps + 1, t = (tl + th) / 2;
    deg_t ip = pow & ((1L << t) - 1), lip = pow;
    bool ph = true;
    bool *c = zn->computed;
    compf_ptr x = zn->zn;
    compf_ptr res = zn->res;
    
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
    // here if ip > 0, then power ip pre-computed, tl is max with this property and th = tl + 1
    if(ip > 0 && pow - ip <= n && c[pow - ip - 1]) { // cheapest route, just one product
        compf_mul(res, x + (ip - 1), x + (pow - ip - 1));
        
        return res;
    }
    
    compf pth;
    powsf_dyadic(pth, zn, th - 1, false);
    if(ip == 0) {
        res->x = 1;
        res->y = 0;
    } else {
        *res = x[ip - 1];
    }
    deg_t tth = (deg_t) (1U << (th - 1));
    
    while(pow > ip + tth) {
        if(pow & tth) {
            compf_mul(res, res, pth);
            ip += tth;
        }
        
        // repeated squares method
        compf_sqr(pth, pth);
        tth *= 2;
    }
    
    // here pow - ip == tth
    compf_mul(res, res, pth);
    
    return res;
}
