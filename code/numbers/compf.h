//
//  compf.h
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
 \file compf.h
 \brief Definition of machine complex numbers.
*/

#ifndef compf_h
#define compf_h

#include "ntypes.h"

// MARK: data types

/// @struct compf_struct
/// @brief Machine complex numbers.
typedef struct {
    coeff_t x;  ///< the real part of the complex number
    coeff_t y;  ///< the imaginary part of the complex number
} compf_struct;

/// \brief Practical wrapper for @c compf_struct
///
/// To avoid the constant use @c * and @c & the type  @c compf is a pointer.
///
/// Example of use:
/// @code
///    compf c;
///    polyf cx = c->x;
/// @endcode
typedef compf_struct compf[1];

/// Convenience pointer to @c compf_struct
typedef compf_struct *compf_ptr;

// MARK: macros & functions prototypes

/// Sets @c d to @c a.
#define compf_set(d, a)         (d)->x = (a)->x; \
                                (d)->y = (a)->y;

/// Sets @c d to the real value @c a.
#define compf_setr(d, a)        (d)->x = (a); \
                                (d)->y = 0;

/// Sets @c d to @c -a.
#define compf_neg(d, a)         (d)->x = -(a)->x; \
                                (d)->y = -(a)->y;

/// Adds @c a to @c b and stores the result in @c d, all of type @ref compf or @ref compf_ptr.
#define compf_add(d, a, b)      (d)->x = (a)->x + (b)->x; \
                                (d)->y = (a)->y + (b)->y;

/// Adds the complex number @c a to the real number @c b and stores the result in @c d.
#define compf_addr(d, a, r)     (d)->x = (a)->x + (r); \
                                (d)->y = (a)->y;

/// Subtracts @c b from @c a and stores the result in @c d, all of type @ref compf or @ref compf_ptr.
#define compf_sub(d, a, b)      (d)->x = (a)->x - (b)->x; \
                                (d)->y = (a)->y - (b)->y;

/// Subtracts the real number @c b from the complex number @c a and stores the result in @c d.
#define compf_subr(d, a, r)     (d)->x = (a)->x - (r); \
                                (d)->y = (a)->y;

#define compf_mul(d, a, b)      coeff_t px = (a)->x * (b)->x - (a)->y * (b)->y; \
                                (d)->y = (a)->x * (b)->y + (a)->y * (b)->x; \
                                (d)->x = px;

/// Multiplies the complex number @c a to the real number @c r and stores the result in @c d.
#define compf_mulr(d, a, r)     (d)->x = (a)->x * (r); \
                                (d)->y = (a)->y * (r);

/// Adds the complex number @c a*r to @c d, where @c d,a are complex and @c r is real.
#define compf_amr(d, a, r)      (d)->x += (a)->x * (r); \
                                (d)->y += (a)->y * (r);

/// Squares @c a and stores the result in @c d, all of type @ref compf or @ref compf_ptr.
#define compf_sqr(d, a)         coeff_t px = (a)->x * (a)->x - (a)->y * (a)->y; \
                                (d)->y = 2 * (a)->x * (a)->y; \
                                (d)->x = px;

/// Divides @c a by @c b and stores the result in @c d, all of type @ref compf or @ref compf_ptr.
#define compf_div(d, a, b)      { \
                                coeff_t m2 = compf_mod2(b); \
                                coeff_t px = ((a)->x * (b)->x + (a)->y * (b)->y) / m2; \
                                (d)->y = ((a)->y * (b)->x - (a)->x * (b)->y) / m2; \
                                (d)->x = px;\
                                }

/// Computes the modulus of the complex number @c a.
#define compf_mod(a)            fhypot((a)->x, (a)->y)

/// Computes the distance betwwen the complex numbera @c a and @c b.
#define compf_dist(a, b)        fhypot((a)->x - (b)->x, (a)->y - (b)->y)

/// Computes the square of the modulus of the complex number @c a.
/// \warning There is a danger of overflow, better use the slower @c compf_mod() instead.
#define compf_mod2(a)           ((a)->x * (a)->x + (a)->y * (a)->y)

/// Computes the @c log_2 of the modulus of the complex number @c a.
#define compf_log2(a)           (plog2(compf_mod(a)))

/// Computes the @c log_2 of the absolute value of the real number @c a.
#define coeff_log2(a)           (plog2((a) < 0 ? -(a) : (a)))

/// Computes the scale of the complex number @c a, see [1].
#define compf_s(a)              (pfloor(plog2(compf_mod(a))) + 1)

/// Computes the scale of the real number @c a, see [1].
#define coeff_s(a)              (pfloor(coeff_log2(a)) + 1)

#endif /* compf_h */
