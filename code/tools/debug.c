//
//  debug.c
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

#include "debug.h"

char *pm(mpfr_t x, int dig) {
    char fmt[30];
    snprintf(fmt, 30, "%%.%dRg", dig);
    
    int len = 2 * dig + 10;
    char *str = malloc(len);
    mpfr_snprintf(str, len, fmt, x);
    
    return str;
}

char *pc(comp z, int dig) {
    char fmt[40];
    snprintf(fmt, 40, "%%.%dRg, %%.%dRg", dig, dig);
    
    int len = 4 * dig + 25;
    char *str = malloc(len);
    mpfr_snprintf(str, len, fmt, z->x, z->y);
    
    return str;
}

char *pcf(compf z) {
    int len = 100;
    char *str = malloc(len);
    mpfr_snprintf(str, len, "%.19Lg, %.19Lg", z->x, z->y);
    
    return str;    
}
