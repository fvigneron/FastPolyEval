//
//  ntypes.c
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

#include <mpfr.h>

#include "ntypes.h"

bool ntypes_check(void) {
    if(sizeof(byte) != 1 || sizeof(word) != 2 ||
       sizeof(uint) != 4 || sizeof(ulong) != 8) {
        return false;
    }
    
#ifdef MACHINE_EXTRA_PREC
    if(sizeof(coeff_t) < 10) {
        return false;
    }
#else
#ifndef MACHINE_LOW_PREC
    if(sizeof(coeff_t) != 8) {
        return false;
    }
#else
    if(sizeof(coeff_t) != 4) {
        return false;
    }
#endif
#endif
    
    bool mpexp = mpfr_set_emin(-MAX_MP_EXP) == 0;
    mpexp = mpexp && mpfr_set_emax(MAX_MP_EXP) == 0;
    
    return mpexp && (1.0L + 1E-18L) != 1.0L;
}

real_t bits_sum(real_t b1, real_t b2) {
    if(b1 >= b2 + precf) {
        return b1 + pldexp(b1, 1 - precf); // ad ulp(b1), as log(1+x) ~ x when 0 < x << 1
    }
    
    if(b2 >= b1 + precf) {
        return b2 + pldexp(b2, 1 - precf); // ad ulp(b2)
    }
    
    real_t b = pfloor(b1);
    
    return plog2(pexp2(b1 - b) + pexp2(b2 - b)) + b;
}

real_t nt_err(real_t vb, real_t evb, real_t db, real_t edb) {
    real_t v = evb - vb, d = edb - db;
    real_t m = bits_sum(v, d);
    real_t terr = bits_sum(m, v + d);
    real_t absErr = vb - db + terr;
    
    return absErr;
}
