//
//  debug.h
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

/**
 \file debug.h
 \brief A few basic tool for easy debug when usinf MPFR.
*/

#ifndef debug_h
#define debug_h

#include <stdio.h>
#include <stdlib.h>
#include <mpfr.h>

#include "ntypes.h"
#include "compf.h"
#include "comp.h"

/// @brief Return a string containing the value of @c x with @c dig digits.
///
/// @warning The user is responsible to free the memory containg the results.
///
/// @param x a real number, MPFR
/// @param dig the number of digits to print
///
/// @return the string containing the value of @c z, @c NULL if some error occurred
char *pm(mpfr_t x, int dig);

/// @brief Return a string containing the value of @c x with @c dig digits.
///
/// The format is that of a CSV file line, the real and imaginary components are separated by a comma.
///
/// @warning The user is responsible to free the memory containg the results.
///
/// @param z a complex number, with MPFR real and imaginary parts
/// @param dig the number of digits to print
///
/// @return the string containing the value of @c z, @c NULL if some error occurred
char *pc(comp z, int dig);

/// @brief Return a string containing the value of @c z.
///
/// The format is that of a CSV file line, the real and imaginary components are separated by a comma.
///
/// @warning The user is responsible to free the memory containg the results.
///
/// @param z a complex number, with machine numbers real and imaginary parts
///
/// @return the string containing the value of @c z, @c NULL if some error occurred
char *pcf(compf z);

#endif /* debug_h */
