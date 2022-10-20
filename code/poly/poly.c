//
//  poly.c
//
//  Authors: Nicolae Mihalache & François Vigneron
//
//  This software is released under BSD licence, with an attribution clause (see Licence file).
//  Please cite the reference below if you use or distribute this software.
//
//  • [1] R. Anton, N. Mihalache & F. Vigneron. Fast evaluation of real and complex polynomials. 2022.
//        https://hal.archives-ouvertes.fr/hal-03820369w
//
//  Copyright 2022 Univ. Paris-Est Créteil & Univ. de Reims Champagne-Ardenne.
//

#include <stdlib.h>

#include "poly.h"


poly poly_new(deg_t degree, prec_t prec) {
    if(degree > MAX_DEG || prec < precf) {
        return NULL;
    }
    
    poly P = malloc(sizeof(poly_struct));
    P->degree = degree;
    P->prec = prec;
    P->modified = true;
    
    P->a = malloc(sizeof(comp_struct) * (degree + 1));
    comp_ptr a = P->a;
    for (deg_t i = 0; i <= degree; i++) {
        comp_initz(a + i, prec);
    }
    
    mpfr_init2(P->buf1, prec);
    mpfr_init2(P->buf2, prec);
    
    return P;
}

poly poly_from_roots(comp_ptr roots, deg_t degree, prec_t prec) {
    if(degree <= 2) {
        poly P = poly_new(degree, prec);
        mpfr_set_si(P->a[degree].x, 1, MPFR_RNDN);
        mpfr_set_zero(P->a[degree].y, 1);
        
        switch(degree) {
            case 1:
                mpfr_neg(P->a[0].x, roots[0].x, MPFR_RNDN);
                mpfr_neg(P->a[0].y, roots[0].y, MPFR_RNDN);
                break;
              
            case 2:
                mpfr_add(P->a[1].x, roots[0].x, roots[1].x, MPFR_RNDN);
                mpfr_neg(P->a[1].x, P->a[1].x, MPFR_RNDN);
                
                mpfr_add(P->a[1].y, roots[0].y, roots[1].y, MPFR_RNDN);
                mpfr_neg(P->a[1].y, P->a[1].y, MPFR_RNDN);
                
                comp_mul(P->a, roots, roots + 1, P->buf1, P->buf2);
                break;
        }
        
        return P;
    }
    
    deg_t d2 = degree >> 1;
    poly P = poly_from_roots(roots, d2, prec);
    poly Q = poly_from_roots(roots + d2, degree - d2, prec);
    
    poly R = poly_prod(P, Q);
    poly_free(P);
    poly_free(Q);
    
    return R;
}

bool poly_free(poly P) {
    if(P == NULL || P->a == NULL) {
        return false;
    }
    
    mpfr_clear(P->buf1);
    mpfr_clear(P->buf2);
    
    deg_t degree = P->degree;
    comp_ptr a = P->a;
    for (deg_t i = 0; i <= degree; i++) {
        comp_clear(a + i);
    }
    
    free(a);
    free(P);
    
    return true;
}

bool poly_set(poly P, comp coeff, deg_t ind) {
    if(P == NULL || coeff == NULL || ind > P->degree) {
        return false;
    }
    
    comp_set(P->a + ind, coeff);
    P->modified = true;
    
    return true;
}

bool poly_eval(comp res, poly P, comp z) {
    if(P == NULL || res == NULL || z == NULL) {
        return false;
    }
    
    comp_ptr c = P->a;
    deg_t deg = P->degree;
    comp_set(res, c + deg);
    
    mpfr_ptr b1 = P->buf1, b2 = P->buf2;
    
    for (deg_t i = deg - 1; i <= MAX_DEG; i--) {
        comp_mul(res, res, z, b1, b2);
        comp_add(res, res, c + i);
    }
    
    return true;
}

bool poly_eval_r(comp res, poly P, mpfr_t x) {
    if(P == NULL || res == NULL || x == NULL || mpfr_nan_p(x)) {
        return false;
    }
    
    comp_ptr c = P->a;
    deg_t deg = P->degree;
    comp_set(res, c + deg);
    
    for (deg_t i = deg - 1; i <= MAX_DEG; i--) {
        comp_mulr(res, res, x);
        comp_add(res, res, c + i);
    }
    
    return true;
}

poly poly_derivative(poly P) {
    if(P == NULL) {
        return NULL;
    }
    
    deg_t deg = P->degree == 0 ? 0 : P->degree - 1;
    poly D = poly_new(deg, P->prec);
    mpfr_set_zero(D->a[0].x, 1);
    mpfr_set_zero(D->a[0].y, 1);
    
    comp_ptr c = P->a;
    comp_ptr b = D->a;
    
    deg_t j;
    for (deg_t i = P->degree; i > 0; i--) {
        j = i - 1;
        comp_muli(b + j, c + i, i);
    }
    
    return D;
}

poly poly_sum(poly P, poly Q) {
    if(P == NULL || Q == NULL) {
        return NULL;
    }
    
    prec_t prec = P->prec > Q->prec ? P->prec : Q->prec;
    
    deg_t dp = P->degree;
    deg_t dq = Q->degree;
    deg_t ds = dp > dq ? dp : dq;
    deg_t dm = dp + dq - ds;
    
    poly S = poly_new(ds, prec);
    
    comp_ptr a = P->a;
    comp_ptr b = Q->a;
    comp_ptr c = S->a;
    
    if(dp > dq) {
        for (deg_t i = dm + 1; i <= ds; i++) {
            comp_set(c + i, a + i);
        }
    } else {
        for (deg_t i = dm + 1; i <= ds; i++) {
            comp_set(c + i, b + i);
        }
    }
    
    for (deg_t i = 0; i <= dm; i++) {
        comp_add(c + i, a + i, b + i);
    }
    
    // check if there are dominant coefficients that are null
    while(comp_zero(c + ds) && ds > 0) {
        ds --;
    }
    
    if(ds < S->degree) {
        S->a = realloc(S->a, (ds + 1) * sizeof(comp_struct));
        S->degree = ds;
    }
    
    return S;
}

poly poly_diff(poly P, poly Q) {
    if(P == NULL || Q == NULL) {
        return NULL;
    }
    
    prec_t prec = P->prec > Q->prec ? P->prec : Q->prec;
    
    deg_t dp = P->degree;
    deg_t dq = Q->degree;
    deg_t ds = dp > dq ? dp : dq;
    deg_t dm = dp + dq - ds;
    
    poly S = poly_new(ds, prec);
    
    comp_ptr a = P->a;
    comp_ptr b = Q->a;
    comp_ptr c = S->a;
    
    if(dp > dq) {
        for (deg_t i = dm + 1; i <= ds; i++) {
            comp_set(c + i, a + i);
        }
    } else {
        for (deg_t i = dm + 1; i <= ds; i++) {
            comp_neg(c + i, b + i);
        }
    }
    
    for (deg_t i = 0; i <= dm; i++) {
        comp_sub(c + i, a + i, b + i);
    }
    
    // check if there are dominant coefficients that are null
    while(comp_zero(c + ds) && ds > 0) {
        ds --;
    }
    
    if(ds < S->degree) {
        S->a = realloc(S->a, (ds + 1) * sizeof(comp_struct));
        S->degree = ds;
    }
    
    return S;
}

poly poly_prod(poly P, poly Q) {
    if(P == NULL || Q == NULL) {
        return NULL;
    }
    
    prec_t prec = P->prec > Q->prec ? P->prec : Q->prec;
    
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
    
    poly S = poly_new(ds, prec);
    
    comp_ptr a = P->a;
    comp_ptr b = Q->a;
    comp_ptr c = S->a;
    
    mpfr_ptr b1 = P->buf1, b2 = P->buf2;
    
    // compute the coefficients
    for (deg_t i = 0; i <= ds; i++) {
        comp p, buf;
        comp_initz(p, prec);
        comp_init(buf, prec);
        
        deg_t min = i <= dq ? 0 : i - dq;
        deg_t max = dp <= i ? dp : i;
        
        for (deg_t j = min; j <= max; j++) {
            comp_mul(buf, a + j, b + (i - j), b1, b2);
            comp_add(p, p, buf);
        }
        
        comp_set(c + i, p);
    }
    
    return S;
}

poly poly_sqr(poly P) {
    if(P == NULL) {
        return NULL;
    }
    
    deg_t dp = P->degree;
    prec_t prec = P->prec;
    
    // check for arithmetic overflow
    long double ddp = dp;
    long double dds = 2 * ddp;
    if(dds > MAX_DEG) {
        return NULL;
    }
    
    deg_t ds = dp << 1;
    
    poly S = poly_new(ds, prec);
    
    comp_ptr a = P->a;
    comp_ptr c = S->a;
    
    mpfr_ptr b1 = P->buf1, b2 = P->buf2;
    
    // compute the coefficients
    for (deg_t i = 0; i <= ds; i++) {
        comp p, buf;
        comp_initz(p, prec);
        comp_init(buf, prec);
        
        deg_t min = i <= dp ? 0 : i - dp;
        deg_t max = dp <= i ? dp : i;
        deg_t h = (min + max + 1) / 2;
        
        for (deg_t j = min; j < h; j++) {
            comp_mul(buf, a + j, a + (i - j), b1, b2);
            comp_add(p, p, buf);
        }
        
        comp_muli(p, p, 2);
        if(((max - min) & 1) == 0) {
            comp_sqr(buf, a + h, b1);
            comp_add(p, p, buf);
        }
        
        comp_set(c + i, p);
    }
    
    return S;
}
