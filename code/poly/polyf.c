//
//  polyf.c
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

#include "polyf.h"


polyf polyf_new(deg_t degree) {
    if(degree > MAX_DEG) {
        return NULL;
    }
    
    polyf P = malloc(sizeof(polyf_struct));
    P->degree = degree;
    P->modified = true;
    
    P->a = malloc(sizeof(compf_struct) * (degree + 1));
    compf_ptr a = P->a;
    compf z = {{0, 0}};
    for (deg_t i = 0; i <= degree; i++) {
        a[i] = *z;
    }
    
    return P;
}

polyf polyf_from_roots(compf_ptr roots, deg_t degree) {
    if(degree <= 2) {
        polyf P = polyf_new(degree);
        P->a[degree].x = 1;
        P->a[degree].y = 0;
        
        switch(degree) {                
            case 1:
                P->a[0].x = -roots[0].x;
                P->a[0].y = -roots[0].y;
                break;
              
            case 2:
                P->a[1].x = -roots[0].x - roots[1].x;
                P->a[1].y = -roots[0].y - roots[1].y;
                
                compf_mul(P->a, roots, roots + 1);
                break;
        }
        
        return P;
    }
    
    deg_t d2 = degree >> 1;
    polyf P = polyf_from_roots(roots, d2);
    polyf Q = polyf_from_roots(roots + d2, degree - d2);
    
    polyf R = polyf_prod(P, Q);
    polyf_free(P);
    polyf_free(Q);
    
    return R;
}

bool polyf_free(polyf P) {
    if(P == NULL || P->a == NULL) {
        return false;
    }
    
    free(P->a);
    free(P);
    
    return true;
}

bool polyf_set(polyf P, compf coeff, deg_t ind) {
    if(P == NULL || coeff == NULL || ind > P->degree) {
        return false;
    }
    
    P->a[ind] = *coeff;
    P->modified = true;
    
    return true;
}

bool polyf_eval(compf res, polyf P, compf z) {
    if(P == NULL || res == NULL || z == NULL) {
        return false;
    }
    
    compf_ptr c = P->a;
    deg_t deg = P->degree;
    *res = c[deg];
    
    for (deg_t i = deg - 1; i <= MAX_DEG; i--) {
        compf_mul(res, res, z);
        compf_add(res, res, c + i);
    }
    
    return true;
}

bool polyf_eval_r(compf res, polyf P, coeff_t x) {
    if(P == NULL || res == NULL || isnan(x)) {
        return false;
    }
    
    compf_ptr c = P->a;
    deg_t deg = P->degree;
    *res = c[deg];
    
    for (deg_t i = deg - 1; i <= MAX_DEG; i--) {
        compf_mulr(res, res, x);
        compf_add(res, res, c + i);
    }
    
    return true;
}

polyf polyf_derivative(polyf P) {
    if(P == NULL) {
        return NULL;
    }
    
    deg_t deg = P->degree == 0 ? 0 : P->degree - 1;
    polyf D = polyf_new(deg);
    D->a[0].x = 0;
    D->a[0].y = 0;
    
    compf_ptr c = P->a;
    compf_ptr b = D->a;
    
    deg_t j;
    for (deg_t i = P->degree; i > 0; i--) {
        j = i - 1;
        compf_mulr(b + j, c + i, i);
    }
    
    return D;
}

polyf polyf_sum(polyf P, polyf Q) {
    if(P == NULL || Q == NULL) {
        return NULL;
    }
    
    deg_t dp = P->degree;
    deg_t dq = Q->degree;
    deg_t ds = dp > dq ? dp : dq;
    deg_t dm = dp + dq - ds;
    
    polyf S = polyf_new(ds);
    
    compf_ptr a = P->a;
    compf_ptr b = Q->a;
    compf_ptr c = S->a;
    
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
        compf_add(c + i, a + i, b + i);
    }
    
    // check if there are dominant coefficients that are null
    while(c[ds].x == 0 && c[ds].y == 0 && ds > 0) {
        ds --;
    }
    
    if(ds < S->degree) {
        S->a = realloc(S->a, (ds + 1) * sizeof(compf_struct));
        S->degree = ds;
    }
    
    return S;
}

polyf polyf_diff(polyf P, polyf Q) {
    if(P == NULL || Q == NULL) {
        return NULL;
    }
    
    deg_t dp = P->degree;
    deg_t dq = Q->degree;
    deg_t ds = dp > dq ? dp : dq;
    deg_t dm = dp + dq - ds;
    
    polyf S = polyf_new(ds);
    
    compf_ptr a = P->a;
    compf_ptr b = Q->a;
    compf_ptr c = S->a;
    
    if(dp > dq) {
        for (deg_t i = dm + 1; i <= ds; i++) {
            c[i] = a[i];
        }
    } else {
        for (deg_t i = dm + 1; i <= ds; i++) {
            compf_neg(c + i, b + i);
        }
    }
    
    for (deg_t i = 0; i <= dm; i++) {
        compf_sub(c + i, a + i, b + i);
    }
    
    // check if there are dominant coefficients that are null
    while(c[ds].x == 0 && c[ds].y == 0 && ds > 0) {
        ds --;
    }
    
    if(ds < S->degree) {
        S->a = realloc(S->a, (ds + 1) * sizeof(compf_struct));
        S->degree = ds;
    }
    
    return S;
}

polyf polyf_prod(polyf P, polyf Q) {
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
    
    polyf S = polyf_new(ds);
    
    compf_ptr a = P->a;
    compf_ptr b = Q->a;
    compf_ptr c = S->a;
    
    // compute the coefficients
    for (deg_t i = 0; i <= ds; i++) {
        compf p = {{0, 0}}, buf;
        deg_t min = i <= dq ? 0 : i - dq;
        deg_t max = dp <= i ? dp : i;
        
        for (deg_t j = min; j <= max; j++) {
            compf_mul(buf, a + j, b + (i - j));
            compf_add(p, p, buf);
        }
        
        c[i] = *p;
    }
    
    return S;
}

polyf polyf_sqr(polyf P) {
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
    
    polyf S = polyf_new(ds);
    
    compf_ptr a = P->a;
    compf_ptr c = S->a;
    
    // compute the coefficients
    for (deg_t i = 0; i <= ds; i++) {
        compf p = {{0, 0}}, buf;
        deg_t min = i <= dp ? 0 : i - dp;
        deg_t max = dp <= i ? dp : i;
        deg_t h = (min + max + 1) / 2;
        
        for (deg_t j = min; j < h; j++) {
            compf_mul(buf, a + j, a + (i - j));
            compf_add(p, p, buf);
        }
        
        compf_mulr(p, p, 2);
        if(((max - min) & 1) == 0) {
            compf_sqr(buf, a + h);
            compf_add(p, p, buf);
        }
        
        c[i] = *p;
    }
    
    return S;
}
