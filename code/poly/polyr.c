//
//  polyr.c
//
//  Authors: Nicolae Mihalache & François Vigneron
//
//  This software is released under BSD licence, with an attribution clause (see Licence file).
//  Please cite the reference below if you use or distribute this software.
//
//  • [1] R. Anton, N. Mihalache & F. Vigneron. Fast evaluation of real and complex polynomials. 2022.
//        https://hal.archives-ouvertes.fr/hal-03820369
//
//  Copyright 2022 Univ. Paris-Est Créteil & Univ. de Reims Champagne-Ardenne.w
//

#include <stdlib.h>

#include "polyr.h"


polyr polyr_new(deg_t degree, prec_t prec) {
    if(degree > MAX_DEG || prec < precf) {
        return NULL;
    }
    
    polyr P = malloc(sizeof(polyr_struct));
    P->degree = degree;
    P->prec = prec;
    P->modified = true;
    
    P->a = malloc(sizeof(__mpfr_struct) * (degree + 1));
    mpfr_ptr a = P->a;
    for (deg_t i = 0; i <= degree; i++) {
        mpfr_init2(a + i, prec);
        mpfr_set_zero(a + i, 1);
    }
    
    mpfr_init2(P->buf1, prec);
    mpfr_init2(P->buf2, prec);
    
    return P;
}

bool polyr_free(polyr P) {
    if(P == NULL || P->a == NULL) {
        return false;
    }
    
    mpfr_clear(P->buf1);
    mpfr_clear(P->buf2);
    
    deg_t degree = P->degree;
    mpfr_ptr a = P->a;
    for (deg_t i = 0; i <= degree; i++) {
        mpfr_clear(a + i);
    }
    
    free(a);
    free(P);
    
    P->a = NULL;
    
    return true;
}

bool polyr_set(polyr P, mpfr_t coeff, deg_t ind) {
    if(P == NULL || coeff == NULL || ind > P->degree) {
        return false;
    }
    
    mpfr_set(P->a + ind, coeff, MPFR_RNDN);
    P->modified = true;
    
    return true;
}

bool polyr_seti(polyr P, long coeff, deg_t ind) {
    if(P == NULL || ind > P->degree) {
        return false;
    }
    
    mpfr_set_si(P->a + ind, coeff, MPFR_RNDN);
    P->modified = true;
    
    return true;
}

bool polyr_eval_c(comp res, polyr P, comp z) {
    if(P == NULL || res == NULL || z == NULL) {
        return false;
    }
    
    mpfr_ptr c = P->a;
    deg_t deg = P->degree;
    mpfr_set(res->x, c + deg, MPFR_RNDN);
    mpfr_set_zero(res->y, 1);
    
    mpfr_ptr b1 = P->buf1, b2 = P->buf2;
    
    for (deg_t i = deg - 1; i <= MAX_DEG; i--) {
        comp_mul(res, res, z, b1, b2);
        mpfr_add(res->x, res->x, c + i, MPFR_RNDN);
    }
    
    return true;
}

bool polyr_eval(mpfr_t res, polyr P, mpfr_t x) {
    if(P == NULL || res == NULL || x == NULL || mpfr_nan_p(x)) {
        return false;
    }
    
    mpfr_ptr c = P->a;
    deg_t deg = P->degree;
    mpfr_set(res, c + deg, MPFR_RNDN);
    
    for (deg_t i = deg - 1; i <= MAX_DEG; i--) {
        mpfr_mul(res, res, x, MPFR_RNDN);
        mpfr_add(res, res, c + i, MPFR_RNDN);
    }
    
    return true;
}

polyr polyr_derivative(polyr P) {
    if(P == NULL) {
        return NULL;
    }
    
    deg_t deg = P->degree == 0 ? 0 : P->degree - 1;
    polyr D = polyr_new(deg, P->prec);
    mpfr_set_zero(D->a, 1);
    
    mpfr_ptr c = P->a;
    mpfr_ptr b = D->a;
    
    deg_t j;
    for (deg_t i = P->degree; i > 0; i--) {
        j = i - 1;
        mpfr_mul_si(b + j, c + i, i, MPFR_RNDN);
    }
    
    return D;
}

polyr polyr_sum(polyr P, polyr Q) {
    if(P == NULL || Q == NULL) {
        return NULL;
    }
    
    prec_t prec = P->prec > Q->prec ? P->prec : Q->prec;
    
    deg_t dp = P->degree;
    deg_t dq = Q->degree;
    deg_t ds = dp > dq ? dp : dq;
    deg_t dm = dp + dq - ds;
    
    polyr S = polyr_new(ds, prec);
    
    mpfr_ptr a = P->a;
    mpfr_ptr b = Q->a;
    mpfr_ptr c = S->a;
    
    if(dp > dq) {
        for (deg_t i = dm + 1; i <= ds; i++) {
            mpfr_set(c + i, a + i, MPFR_RNDN);
        }
    } else {
        for (deg_t i = dm + 1; i <= ds; i++) {
            mpfr_set(c + i, b + i, MPFR_RNDN);
        }
    }
    
    for (deg_t i = 0; i <= dm; i++) {
        mpfr_add(c + i, a + i, b + i, MPFR_RNDN);
    }
    
    // check if there are dominant coefficients that are null
    while(mpfr_zero_p(c + ds) && ds > 0) {
        ds --;
    }
    
    if(ds < S->degree) {
        S->a = realloc(S->a, (ds + 1) * sizeof(__mpfr_struct));
        S->degree = ds;
    }
    
    return S;
}

polyr polyr_diff(polyr P, polyr Q) {
    if(P == NULL || Q == NULL) {
        return NULL;
    }
    
    prec_t prec = P->prec > Q->prec ? P->prec : Q->prec;
    
    deg_t dp = P->degree;
    deg_t dq = Q->degree;
    deg_t ds = dp > dq ? dp : dq;
    deg_t dm = dp + dq - ds;
    
    polyr S = polyr_new(ds, prec);
    
    mpfr_ptr a = P->a;
    mpfr_ptr b = Q->a;
    mpfr_ptr c = S->a;
    
    if(dp > dq) {
        for (deg_t i = dm + 1; i <= ds; i++) {
            mpfr_set(c + i, a + i, MPFR_RNDN);
        }
    } else {
        for (deg_t i = dm + 1; i <= ds; i++) {
            mpfr_neg(c + i, b + i, MPFR_RNDN);
        }
    }
    
    for (deg_t i = 0; i <= dm; i++) {
        mpfr_sub(c + i, a + i, b + i, MPFR_RNDN);
    }
    
    // check if there are dominant coefficients that are null
    while(mpfr_zero_p(c + ds) && ds > 0) {
        ds --;
    }
    
    if(ds < S->degree) {
        S->a = realloc(S->a, (ds + 1) * sizeof(__mpfr_struct));
        S->degree = ds;
    }
    
    return S;
}

polyr polyr_prod(polyr P, polyr Q) {
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
    
    polyr S = polyr_new(ds, prec);
    
    mpfr_ptr a = P->a;
    mpfr_ptr b = Q->a;
    mpfr_ptr c = S->a;
    
    mpfr_ptr b1 = P->buf1, b2 = P->buf2;
    
    // compute the coefficients
    for (deg_t i = 0; i <= ds; i++) {
        mpfr_set_zero(b1, 1);
        
        deg_t min = i <= dq ? 0 : i - dq;
        deg_t max = dp <= i ? dp : i;
        
        for (deg_t j = min; j <= max; j++) {
            mpfr_mul(b2, a + j, b + (i - j), MPFR_RNDN);
            mpfr_add(b1, b1, b2, MPFR_RNDN);
        }
        
        mpfr_set(c + i, b1, MPFR_RNDN);
    }
    
    return S;
}

polyr polyr_sqr(polyr P) {
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
    
    polyr S = polyr_new(ds, prec);
    
    mpfr_ptr a = P->a;
    mpfr_ptr c = S->a;
    
    mpfr_ptr b1 = P->buf1, b2 = P->buf2;
    
    // compute the coefficients
    for (deg_t i = 0; i <= ds; i++) {
        mpfr_set_zero(b1, 1);
        
        deg_t min = i <= dp ? 0 : i - dp;
        deg_t max = dp <= i ? dp : i;
        deg_t h = (min + max + 1) / 2;
        
        for (deg_t j = min; j < h; j++) {
            mpfr_mul(b2, a + j, a + (i - j), MPFR_RNDN);
            mpfr_add(b1, b1, b2, MPFR_RNDN);
        }
        
        mpfr_mul_2si(b1, b1, 1, MPFR_RNDN);
        if(((max - min) & 1) == 0) {
            mpfr_sqr(b2, a + h, MPFR_RNDN);
            mpfr_add(b1, b1, b2, MPFR_RNDN);
        }
        
        mpfr_set(c + i, b1, MPFR_RNDN);
    }
    
    return S;
}

polyr poly_hyp(int n, prec_t prec) {
    if(n < 0 || n > 35) {
        return NULL;
    }
    
    mpfr_t a;
    mpfr_init2(a, prec);
    
    polyr pn = polyr_new(1, prec);
    
    mpfr_set_zero(a, 1);
    polyr_set(pn, a, 0);
    
    mpfr_set_ui(a, 1, MPFR_RNDN);
    polyr_set(pn, a, 1);
    
    for (int i = 1; i < n; i++) {
        polyr p = polyr_sqr(pn);
        mpfr_add(p->a + 1, p->a + 1, a, MPFR_RNDN);
        
        polyr_free(pn);
        pn = p;
    }
    
    mpfr_clear(a);
    
    return pn;
}

// returns (a x P(x) - b P(x) - c Q(x)) / d
static polyr polyr_comb(polyr P, polyr Q, long a, long b, long c, long d) {
    prec_t prec = P->prec;
    
    deg_t deg = P->degree + 1;
    deg = Q->degree > deg ? Q->degree : deg;
    polyr R = polyr_new(deg, prec);
    
    mpfr_t x;
    mpfr_init2(x, prec);
    
    if(a != 0) {
        for (deg_t i = 0; i <= P->degree; i++) {
            mpfr_mul_si(R->a + (i + 1), P->a + i, a, MPFR_RNDN);
        }
    }
    
    if(b != 0) {
        for (deg_t i = 0; i <= P->degree; i++) {
            mpfr_mul_si(x, P->a + i, b, MPFR_RNDN);
            mpfr_sub(R->a + i, R->a + i, x, MPFR_RNDN);
        }
    }
    
    if(c != 0) {
        for (deg_t i = 0; i <= Q->degree; i++) {
            mpfr_mul_si(x, Q->a + i, c, MPFR_RNDN);
            mpfr_sub(R->a + i, R->a + i, x, MPFR_RNDN);
        }
    }
    
    if(d != 1) {
        for (deg_t i = 0; i <= deg; i++) {
            mpfr_div_si(R->a + i, R->a + i, d, MPFR_RNDN);
        }
    }
    
    mpfr_clear(x);
    
    return R;
}

polyr poly_cheb(int n, prec_t prec) {
    if(n < 0 || n > 1000000000) {
        return NULL;
    }
    
    polyr pnm1 = polyr_new(0, prec);
    polyr_seti(pnm1, 1, 0);
    if(n == 0) {
        return pnm1;
    }
    
    polyr pn = polyr_new(1, prec);
    polyr_seti(pn, 1, 1);
    if(n == 1) {
        polyr_free(pnm1);
        
        return pn;
    }
        
    for (int i = 2; i <= n; i++) {
        polyr p = polyr_comb(pn, pnm1, 2, 0, 1, 1);
        polyr_free(pnm1);
        
        pnm1 = pn;
        pn = p;
    }
    
    polyr_free(pnm1);
    
    return pn;
}

polyr poly_leg(int n, prec_t prec) {
    if(n < 0 || n > 1000000000) {
        return NULL;
    }
    
    polyr pnm1 = polyr_new(0, prec);
    polyr_seti(pnm1, 1, 0);
    if(n == 0) {
        return pnm1;
    }
    
    polyr pn = polyr_new(1, prec);
    polyr_seti(pn, 1, 1);
    if(n == 1) {
        polyr_free(pnm1);
        
        return pn;
    }
        
    for (int i = 2; i <= n; i++) {
        polyr p = polyr_comb(pn, pnm1, 2 * i - 1, 0, i - 1, i);
        polyr_free(pnm1);
        
        pnm1 = pn;
        pn = p;
    }
    
    polyr_free(pnm1);
    
    return pn;
}

polyr poly_her(int n, prec_t prec) {
    if(n < 0 || n > 1000000000) {
        return NULL;
    }
    
    polyr pnm1 = polyr_new(0, prec);
    polyr_seti(pnm1, 1, 0);
    if(n == 0) {
        return pnm1;
    }
    
    polyr pn = polyr_new(1, prec);
    polyr_seti(pn, 2, 1);
    if(n == 1) {
        polyr_free(pnm1);
        
        return pn;
    }
        
    for (int i = 2; i <= n; i++) {
        polyr p = polyr_comb(pn, pnm1, 2, 0, 2 * i - 2, 1);
        polyr_free(pnm1);
        
        pnm1 = pn;
        pn = p;
    }
    
    polyr_free(pnm1);
    
    return pn;
}

polyr poly_lag(int n, prec_t prec) {
    if(n < 0 || n > 1000000000) {
        return NULL;
    }
    
    polyr pnm1 = polyr_new(0, prec);
    polyr_seti(pnm1, 1, 0);
    if(n == 0) {
        return pnm1;
    }
    
    polyr pn = polyr_new(1, prec);
    polyr_seti(pn, 1, 0);
    polyr_seti(pn, -1, 1);
    if(n == 1) {
        polyr_free(pnm1);
        
        return pn;
    }
        
    for (int i = 2; i <= n; i++) {
        polyr p = polyr_comb(pn, pnm1, -1, 1 - 2 * i, i - 1, i);
        polyr_free(pnm1);
        
        pnm1 = pn;
        pn = p;
    }
    
    polyr_free(pnm1);
    
    return pn;
}
