//
//  polyfr.c
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

#include "polyfr.h"


polyfr polyfr_new(deg_t degree) {
    if(degree > MAX_DEG) {
        return NULL;
    }
    
    polyfr P = malloc(sizeof(polyfr_struct));
    P->degree = degree;
    P->modified = true;
    
    P->a = malloc(sizeof(coeff_t) * (degree + 1));
    coeff_t *a = P->a;
    for (deg_t i = 0; i <= degree; i++) {
        a[i] = 0;
    }
    
    return P;
}

bool polyfr_free(polyfr P) {
    if(P == NULL || P->a == NULL) {
        return false;
    }
    
    free(P->a);
    free(P);
    
    return true;
}

bool polyfr_set(polyfr P, coeff_t coeff, deg_t ind) {
    if(P == NULL || isnan(coeff) || ind > P->degree) {
        return false;
    }
    
    P->a[ind] = coeff;
    P->modified = true;
    
    return true;
}

coeff_t polyfr_get(polyfr P, deg_t ind) {
    if(P == NULL || ind > P->degree) {
        return NAN;
    }
    
    return P->a[ind];
}

bool polyfr_eval_c(compf res, polyfr P, compf z) {
    if(P == NULL || z == NULL || res == NULL) {
        return false;
    }
    
    if(z->y == 0) {
        res->x = polyfr_eval(P, z->x);
        res->y = 0;
        
        return true;
    }
    
    coeff_t *c = P->a;
    deg_t deg = P->degree;
    res->x = c[deg];
    res->y = 0;
    
    for (deg_t i = deg - 1; i <= MAX_DEG; i--) {
        compf_mul(res, res, z);
        res->x += c[i];
    }
    
    return true;
}

coeff_t polyfr_eval(polyfr P, coeff_t x) {
    if(P == NULL || isnan(x)) {
        return NAN;
    }
    
    coeff_t *c = P->a;
    deg_t deg = P->degree;
    coeff_t v = c[deg];
    
    for (deg_t i = deg - 1; i <= MAX_DEG; i--) {
        v *= x;
        v += c[i];
    }
    
    return v;
}

polyf polyfr_comp(polyfr P) {
    if(P == NULL) {
        return NULL;
    }
    
    deg_t dp = P->degree;
    polyf S = polyf_new(dp);
    
    coeff_t *a = P->a;
    compf_ptr c = S->a;
    
    for (deg_t i = 0; i <= dp; i++) {
        c[i].x = a[i];
        c[i].y = 0;
    }
    
    return S;
}

polyfr polyfr_derivative(polyfr P) {
    if(P == NULL) {
        return NULL;
    }
    
    deg_t deg = P->degree == 0 ? 0 : P->degree - 1;
    polyfr D = polyfr_new(deg);
    D->a[0] = 0;
    
    coeff_t *c = P->a;
    coeff_t *b = D->a;
    
    deg_t j;
    for (deg_t i = P->degree; i > 0; i--) {
        j = i - 1;
        b[j] = c[i] * i;
    }
    
    return D;
}

polyfr polyfr_sum(polyfr P, polyfr Q) {
    if(P == NULL || Q == NULL) {
        return NULL;
    }
    
    deg_t dp = P->degree;
    deg_t dq = Q->degree;
    deg_t ds = dp > dq ? dp : dq;
    deg_t dm = dp + dq - ds;
    
    polyfr S = polyfr_new(ds);
    
    coeff_t *a = P->a;
    coeff_t *b = Q->a;
    coeff_t *c = S->a;
    
    if(dp > dq) {
        for (deg_t i = dm + 1; i <= ds; i++) {
            c[i] = a[i];
        }
    } else {
        for (deg_t i = dm + 1; i <= ds; i++) {
            c[i] = b[i];
        }
    }
    
    for (deg_t i = 0; i <= dm; i++) {
        c[i] = a[i] + b[i];
    }
    
    // check if there are dominant coefficients that are null
    while(c[ds] == 0 && ds > 0) {
        ds --;
    }
    
    if(ds < S->degree) {
        S->a = realloc(S->a, (ds + 1) * sizeof(coeff_t));
        S->degree = ds;
    }
    
    return S;
}

polyfr polyfr_diff(polyfr P, polyfr Q) {
    if(P == NULL || Q == NULL) {
        return NULL;
    }
    
    deg_t dp = P->degree;
    deg_t dq = Q->degree;
    deg_t ds = dp > dq ? dp : dq;
    deg_t dm = dp + dq - ds;
    
    polyfr S = polyfr_new(ds);
    
    coeff_t *a = P->a;
    coeff_t *b = Q->a;
    coeff_t *c = S->a;
    
    if(dp > dq) {
        for (deg_t i = dm + 1; i <= ds; i++) {
            c[i] = a[i];
        }
    } else {
        for (deg_t i = dm + 1; i <= ds; i++) {
            c[i] = -b[i];
        }
    }
    
    for (deg_t i = 0; i <= dm; i++) {
        c[i] = a[i] - b[i];
    }
    
    // check if there are dominant coefficients that are null
    while(c[ds] == 0 && ds > 0) {
        ds --;
    }
    
    if(ds < S->degree) {
        S->a = realloc(S->a, (ds + 1) * sizeof(coeff_t));
        S->degree = ds;
    }
    
    return S;
}

polyfr polyfr_prod(polyfr P, polyfr Q) {
    if(P == NULL || Q == NULL) {
        return NULL;
    }
    
    deg_t dp = P->degree;
    deg_t dq = Q->degree;
    
    // check for arithmetic overflow
    long double ddp = dp;
    long double ddq = dq;
    long double dds = ddp + ddq;
    if(dds > MAX_DEG) {
        return NULL;
    }
    
    deg_t ds = dp + dq;
    
    polyfr S = polyfr_new(ds);
    
    coeff_t *a = P->a;
    coeff_t *b = Q->a;
    coeff_t *c = S->a;
    
    // compute the coefficients
    for (deg_t i = 0; i <= ds; i++) {
        coeff_t p = 0;
        deg_t min = i <= dq ? 0 : i - dq;
        deg_t max = dp <= i ? dp : i;
        
        for (deg_t j = min; j <= max; j++) {
            p += a[j] * b[i - j];
        }
        
        c[i] = p;
    }
    
    return S;
}

polyfr polyfr_sqr(polyfr P) {
    if(P == NULL) {
        return NULL;
    }
    
    deg_t dp = P->degree;
    
    // check for arithmetic overflow
    long double ddp = dp;
    long double dds = 2 * ddp;
    if(dds > MAX_DEG) {
        return NULL;
    }
    
    deg_t ds = dp << 1;
    
    polyfr S = polyfr_new(ds);
    
    coeff_t *a = P->a;
    coeff_t *c = S->a;
    
    // compute the coefficients
    for (deg_t i = 0; i <= ds; i++) {
        coeff_t p = 0;
        deg_t min = i <= dp ? 0 : i - dp;
        deg_t max = dp <= i ? dp : i;
        deg_t h = (min + max + 1) / 2;
        
        for (deg_t j = min; j < h; j++) {
            p += a[j] * a[i - j];
        }
        
        c[i] = 2 * p + ((max - min) & 1 ? 0 : a[h] * a[h]);
    }
    
    return S;
}

polyfr polyf_hyp(int n) {
    if(n < 0 || n > 35) {
        return NULL;
    }
    
    polyfr pn = polyfr_new(1);
    polyfr_set(pn, 0, 0);
    polyfr_set(pn, 1, 1);
    
    for (int i = 1; i < n; i++) {
        polyfr p = polyfr_sqr(pn);
        polyfr_set(p, polyfr_get(p, 1) + 1, 1);
        
        polyfr_free(pn);
        pn = p;
    }
    
    return pn;
}

// returns (a x P(x) - b P(x) - c Q(x)) / d
static polyfr polyfr_comb(polyfr P, polyfr Q, long a, long b, long c, long d) {
    deg_t deg = P->degree + 1;
    deg = Q->degree > deg ? Q->degree : deg;
    polyfr R = polyfr_new(deg);
    
    if(a != 0) {
        for (deg_t i = 0; i <= P->degree; i++) {
            R->a[i + 1] = a * P->a[i];
        }
    }
    
    if(b != 0) {
        for (deg_t i = 0; i <= P->degree; i++) {
            R->a[i] -= b * P->a[i];
        }
    }
    
    if(c != 0) {
        for (deg_t i = 0; i <= Q->degree; i++) {
            R->a[i] -= c * Q->a[i];
        }
    }
    
    if(d != 1) {
        for (deg_t i = 0; i <= deg; i++) {
            R->a[i] /= d;
        }
    }
    
    return R;
}

polyfr polyf_cheb(int n) {
    if(n < 0 || n > 1000000000) {
        return NULL;
    }
    
    polyfr pnm1 = polyfr_new(0);
    polyfr_set(pnm1, 1, 0);
    if(n == 0) {
        return pnm1;
    }
    
    polyfr pn = polyfr_new(1);
    polyfr_set(pn, 1, 1);
    if(n == 1) {
        polyfr_free(pnm1);
        
        return pn;
    }
        
    for (int i = 2; i <= n; i++) {
        polyfr p = polyfr_comb(pn, pnm1, 2, 0, 1, 1);
        polyfr_free(pnm1);
        
        pnm1 = pn;
        pn = p;
    }
    
    polyfr_free(pnm1);
    
    return pn;
}

polyfr polyf_leg(int n) {
    if(n < 0 || n > 1000000000) {
        return NULL;
    }
    
    polyfr pnm1 = polyfr_new(0);
    polyfr_set(pnm1, 1, 0);
    if(n == 0) {
        return pnm1;
    }
    
    polyfr pn = polyfr_new(1);
    polyfr_set(pn, 1, 1);
    if(n == 1) {
        polyfr_free(pnm1);
        
        return pn;
    }
        
    for (int i = 2; i <= n; i++) {
        polyfr p = polyfr_comb(pn, pnm1, 2 * i - 1, 0, i - 1, i);
        polyfr_free(pnm1);
        
        pnm1 = pn;
        pn = p;
    }
    
    polyfr_free(pnm1);
    
    return pn;
}

polyfr polyf_her(int n) {
    if(n < 0 || n > 1000000000) {
        return NULL;
    }
    
    polyfr pnm1 = polyfr_new(0);
    polyfr_set(pnm1, 1, 0);
    if(n == 0) {
        return pnm1;
    }
    
    polyfr pn = polyfr_new(1);
    polyfr_set(pn, 2, 1);
    if(n == 1) {
        polyfr_free(pnm1);
        
        return pn;
    }
        
    for (int i = 2; i <= n; i++) {
        polyfr p = polyfr_comb(pn, pnm1, 2, 0, 2 * i - 2, 1);
        polyfr_free(pnm1);
        
        pnm1 = pn;
        pn = p;
    }
    
    polyfr_free(pnm1);
    
    return pn;
}

polyfr polyf_lag(int n) {
    if(n < 0 || n > 1000000000) {
        return NULL;
    }
    
    polyfr pnm1 = polyfr_new(0);
    polyfr_set(pnm1, 1, 0);
    if(n == 0) {
        return pnm1;
    }
    
    polyfr pn = polyfr_new(1);
    polyfr_set(pn, 1, 0);
    polyfr_set(pn, -1, 1);
    if(n == 1) {
        polyfr_free(pnm1);
        
        return pn;
    }
        
    for (int i = 2; i <= n; i++) {
        polyfr p = polyfr_comb(pn, pnm1, -1, 1 - 2 * i, i - 1, i);
        polyfr_free(pnm1);
        
        pnm1 = pn;
        pn = p;
    }
    
    polyfr_free(pnm1);
    
    return pn;
}
