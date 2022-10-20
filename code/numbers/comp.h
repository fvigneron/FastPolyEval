//
//  comp.h
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
 \file comp.h
 \brief Definition of MPFR complex numbers.
*/

#ifndef comp_h
#define comp_h

#include <mpfr.h>

#include "ntypes.h"

// MARK: data types

/// @struct comp_struct
/// @brief Multi-precision floating point complex numbers.
typedef struct {
    mpfr_t x;  ///< the real part of the complex number
    mpfr_t y;  ///< the imaginary part of the complex number
} comp_struct;

/// \brief Practical wrapper for @c comp_struct
///
/// To avoid the constant use @c * and @c & the type  @c compf is a pointer.
typedef comp_struct comp[1];

/// Convenience pointer to @c comp_struct
typedef comp_struct *comp_ptr;

// MARK: macros & functions prototypes

/// Allocates memory for the digits of @c d.
#define comp_init(d, prec)      mpfr_init2((d)->x, prec); \
                                mpfr_init2((d)->y, prec);

/// Allocates memory for the digits of @c d and sets its value to @c 0.
#define comp_initz(d, prec)     mpfr_init2((d)->x, prec); \
                                mpfr_init2((d)->y, prec); \
                                mpfr_set_zero((d)->x, 1); \
                                mpfr_set_zero((d)->y, 1);

/// De-allocates the memory used by the digits of @c d.
#define comp_clear(d)           mpfr_clear((d)->x); \
                                mpfr_clear((d)->y);

/// Tests if @c d==0.
#define comp_zero(d)            (mpfr_zero_p((d)->x) && mpfr_zero_p((d)->y))

/// Adds @c a to @c b and stores the result in @c d, all of type @ref comp or @ref comp_ptr.
#define comp_add(d, a, b)       mpfr_add((d)->x, (a)->x, (b)->x, MPFR_RNDN); \
                                mpfr_add((d)->y, (a)->y, (b)->y, MPFR_RNDN);

/// Sets @c d to @c a.
#define comp_set(d, a)          mpfr_set((d)->x, (a)->x, MPFR_RNDN); \
                                mpfr_set((d)->y, (a)->y, MPFR_RNDN);

/// Sets @c d to the real value @c a.
#define comp_setr(d, a)         mpfr_set((d)->x, (a), MPFR_RNDN); \
                                mpfr_set_zero((d)->y, 1);

/// Sets @c d to @c -a.
#define comp_neg(d, a)          mpfr_neg((d)->x, (a)->x, MPFR_RNDN); \
                                mpfr_neg((d)->y, (a)->y, MPFR_RNDN);

/// Adds the complex number @c a to the real number @c b and stores the result in @c d.
#define comp_addr(d, a, r)      mpfr_add((d)->x, (a)->x, (r), MPFR_RNDN); \
                                mpfr_set((d)->y, (a)->y, MPFR_RNDN);

/// Adds the complex number @c a*r to @c d, where @c d,a are complex and @c r is unsigned integer.
#define comp_amu(d, a, r, buf)  mpfr_mul_ui((buf), (a)->x, (r), MPFR_RNDN); \
                                mpfr_add((d)->x, (d)->x, (buf), MPFR_RNDN); \
                                mpfr_mul_ui((buf), (a)->y, (r), MPFR_RNDN); \
                                mpfr_add((d)->y, (d)->y, (buf), MPFR_RNDN);

/// Subtracts @c b from @c a and stores the result in @c d, all of type @ref comp or @ref comp_ptr.
#define comp_sub(d, a, b)       mpfr_sub((d)->x, (a)->x, (b)->x, MPFR_RNDN); \
                                mpfr_sub((d)->y, (a)->y, (b)->y, MPFR_RNDN);

/// Subtracts the real number @c b from the complex number @c a and stores the result in @c d.
#define comp_subr(d, a, r)      mpfr_sub((d)->x, (a)->x, (r), MPFR_RNDN); \
                                mpfr_set((d)->y, (a)->y, MPFR_RNDN);

/// Multiplies @c a to @c b and stores the result in @c d, all of type @ref comp or @ref comp_ptr.
#define comp_mul(d, a, b, buf1, buf2)   mpfr_mul((buf1), (a)->x, (b)->x, MPFR_RNDN); \
                                        mpfr_mul((buf2), (a)->y, (b)->y, MPFR_RNDN); \
                                        mpfr_sub((buf1), (buf1), (buf2), MPFR_RNDN); \
                                        mpfr_mul((buf2), (a)->x, (b)->y, MPFR_RNDN); \
                                        mpfr_mul((d)->y, (a)->y, (b)->x, MPFR_RNDN); \
                                        mpfr_add((d)->y, (d)->y, (buf2), MPFR_RNDN); \
                                        mpfr_set((d)->x, (buf1), MPFR_RNDN);

/// Divides @c a by @c b and stores the result in @c d, all of type @ref comp or @ref comp_ptr.
#define comp_div(d, a, b, b1, b2, b3)   mpfr_sqr((b3), (b)->x, MPFR_RNDN); \
                                        mpfr_sqr((b2), (b)->y, MPFR_RNDN); \
                                        mpfr_add((b3), (b3), (b2), MPFR_RNDN); \
                                        mpfr_mul((b1), (a)->x, (b)->x, MPFR_RNDN); \
                                        mpfr_mul((b2), (a)->y, (b)->y, MPFR_RNDN); \
                                        mpfr_add((b1), (b1), (b2), MPFR_RNDN); \
                                        mpfr_mul((b2), (a)->x, (b)->y, MPFR_RNDN); \
                                        mpfr_mul((d)->y, (a)->y, (b)->x, MPFR_RNDN); \
                                        mpfr_sub((d)->y, (d)->y, (b2), MPFR_RNDN); \
                                        mpfr_div((d)->y, (d)->y, (b3), MPFR_RNDN); \
                                        mpfr_div((d)->x, (b1), (b3), MPFR_RNDN);

/// Computes the modulus of @c a and stores the result in @c m, @c a of type @ref comp and @c m of type mpfr_t.
#define comp_mod(m, a)   mpfr_hypot((m), (a)->x, (a)->y, MPFR_RNDN);

/// Multiplies the complex number @c a to the real number @c r and stores the result in @c d.
#define comp_mulr(d, a, r)      mpfr_mul((d)->x, (a)->x, (r), MPFR_RNDN); \
                                mpfr_mul((d)->y, (a)->y, (r), MPFR_RNDN);

/// Multiplies the complex number @c a to the integer @c i and stores the result in @c d.
#define comp_muli(d, a, i)      mpfr_mul_si((d)->x, (a)->x, (i), MPFR_RNDN); \
                                mpfr_mul_si((d)->y, (a)->y, (i), MPFR_RNDN);

/// Multiplies the complex number @c a to the unsigned integer @c i and stores the result in @c d.
#define comp_mulu(d, a, i)      mpfr_mul_ui((d)->x, (a)->x, (i), MPFR_RNDN); \
                                mpfr_mul_ui((d)->y, (a)->y, (i), MPFR_RNDN);

/// Squares @c a and stores the result in @c d, all of type @ref comp or @ref comp_ptr.
#define comp_sqr(d, a, buf)     mpfr_sqr((buf), (a)->y, MPFR_RNDN); \
                                mpfr_mul((d)->y, (a)->x, (a)->y, MPFR_RNDN); \
                                mpfr_mul_2si((d)->y, (d)->y, 1, MPFR_RNDN); \
                                mpfr_sqr((d)->x, (a)->x, MPFR_RNDN); \
                                mpfr_sub((d)->x, (d)->x, (buf), MPFR_RNDN);

/// @brief Computes the base @c 2 @c log of @c |z|.
///
/// @param z the complex number
///
/// @return @c log_2(|z|)
real_t comp_log2(comp z);

/// @brief Computes the base @c 2 @c log of @c |x|.
///
/// @param x the real number
///
/// @return @c log_2(|x|)
real_t real_log2(mpfr_t x);

/// @brief Computes @c s(z).
///
/// @see [1]
///
/// @param z the complex number
///
/// @return @c [log_2(|z|)]+1
real_t comp_s(comp z);

/// @brief Computes @c s(x).
///
/// @see [1]
///
/// @param x the real number
///
/// @return @c [log_2(|z|)]+1
real_t mpfr_s(mpfr_t x);

#endif /* comp_h */
