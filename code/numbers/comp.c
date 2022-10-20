//
//  comp.c
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

#include "comp.h"

inline real_t real_log2(mpfr_t x) {
    if(x == NULL || mpfr_nan_p(x)) {
        return NAN;
    }
    
    if(mpfr_inf_p(x)) {
        return INF_P;
    }
    
    if(mpfr_zero_p(x)) {
        return INF_M;
    }
    
    // x is regular
    mpfr_exp_t ex = x->_mpfr_exp;
    x->_mpfr_exp = 0;
        
    real_t l = mpfr_getp(x, MPFR_RNDN);
    x->_mpfr_exp = ex;
    l = plog2(l < 0 ? -l : l);
    
    return l + ex;
}

inline real_t mpfr_s(mpfr_t x) {
    return pfloor(real_log2(x)) + 1;
}

real_t comp_log2(comp z) {
    if(z == NULL) {
        return NAN;
    }
    
    if(mpfr_inf_p(z->x) || mpfr_inf_p(z->y)) {
        return INF_P;
    }
    
    if(mpfr_nan_p(z->x) || mpfr_nan_p(z->y)) {
        return NAN;
    }
    
    if(mpfr_zero_p(z->x)) {
        if(mpfr_zero_p(z->y)) {
            return INF_M;
        }
        
        return real_log2(z->y);
    }
           
    if(mpfr_zero_p(z->y)) {
        return real_log2(z->x);
    }
    
    // both coordinates are regular here
    mpfr_exp_t ex = z->x->_mpfr_exp;
    mpfr_exp_t ey = z->y->_mpfr_exp;
    
    if(ex > ey + precf / 2 + 1) {
        return real_log2(z->x);
    }
    
    if(ey > ex + precf / 2 + 1) {
        return real_log2(z->y);
    }
    
    z->x->_mpfr_exp = 0;
    z->y->_mpfr_exp -= ex;
    real_t x = mpfr_getp(z->x, MPFR_RNDN);
    real_t y = mpfr_getp(z->y, MPFR_RNDN);
    z->x->_mpfr_exp = ex;
    z->y->_mpfr_exp += ex;
    
    real_t l = plog2(phypot(x, y));
    
    return l + ex;
}

real_t comp_s(comp z) {
    return pfloor(comp_log2(z)) + 1;
}
