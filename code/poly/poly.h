//
//  poly.h
//
//  Authors: Nicolae Mihalache & François Vigneron
//
//  This software is released under BSD licence, with an attribution clause (see Licence file).
//  Please cite the reference below if you use or distribute this software.
//
//  • [1] R. Anton, N. Mihalache & F. Vigneron. Fast evaluation of real and complex polynomials. 2022.
//        https://hal.archives-ouvertes.fr/hal-03820369
//
//  Copyright 2022 Univ. Paris-Est Créteil & Univ. de Reims Champagne-Ardenne.//

/**
 \file poly.h
 \brief Definition of complex polynomials with arbitary precision coefficients.
*/

#ifndef poly_h
#define poly_h

#include <mpfr.h>

#include "ntypes.h"
#include "comp.h"

// MARK: data types

/// @struct poly_struct
/// @brief Polynomial with multi-precision floating point complex coefficients.
typedef struct {
    deg_t degree;   ///< the degree of the polynomial
    prec_t prec;    ///< the precision of the coefficients, in bits
    comp_ptr a;     ///< the coefficients
    bool modified;  ///< the status of the coefficients
    mpfr_t buf1;    ///< a buffer
    mpfr_t buf2;    ///< another buffer
} poly_struct;

/// @brief Practical wrapper for @c poly_struct
///
/// To avoid the constant use @c * and @c & the type  @c poly_t is a pointer.
typedef poly_struct poly_t[1];

/// Convenience pointer to @c poly_struct
typedef poly_struct *poly;

// MARK: functions prototypes

/// @brief Returns a new complex polynomial of given @c degree, with coefficients of precision @c prec.
///
/// @warning @c prec must be at least @c precf
///
/// @param degree the degree of the polynomial
/// @param prec the precision of the coefficients, in bits
///
/// @return the new polynomial, @c NULL if the degree is larger than @c MAX_DEG.
poly poly_new(deg_t degree, prec_t prec);

/// @brief Returns a new complex polynomial given the list of its @c roots, with coefficients of precision @c prec.
///
/// @warning @c prec must be at least @c precf
///
/// @param roots the roots of the polynomial
/// @param degree the degree of the polynomial
/// @param prec the precision of the coefficients, in bits 
///
/// @return the new polynomial, @c NULL if some error occurred.
poly poly_from_roots(comp_ptr roots, deg_t degree, prec_t prec);

/// @brief Frees all the memory used by the polynomial @c P, assuming the struct has been allocated with
/// @c malloc(), for example with @c poly_new().
///
/// @param P the polynomial
///
/// @return @ref true if successfull, @ref false otherwise.
bool poly_free(poly P);

/// @brief Sets the coefficient of the polynomial @c P corresponding to the power @c ind to @c coeff.
///
/// @param P the polynomial
/// @param coeff the coefficient
/// @param ind the index
///
/// @return @ref true if successfull, @ref false otherwise.
bool poly_set(poly P, comp coeff, deg_t ind);

/// @brief Evaluates @c P(z) using Horner's method.
///
/// @param res the result
/// @param P the polynomial
/// @param z the argument
///
/// @return @ref true if successfull, @ref false otherwise.
bool poly_eval(comp res, poly P, comp z);

/// @brief Evaluates @c P(x) using Horner's method.
///
/// @param res the result
/// @param P the polynomial
/// @param x the real argument
///
/// @return @ref true if successfull, @ref false otherwise.
bool poly_eval_r(comp res, poly P, mpfr_t x);

/// @brief Computes the derivative of @c P.
///
/// @param P the polynomial
///
/// @return the resulted polynomial, @c NULL if some error occurred.
poly poly_derivative(poly P);

/// @brief Computes @c P+Q.
///
/// @param P a polynomial
/// @param Q another polynomial
///
/// @return the resulted polynomial, @c NULL if some error occurred.
poly poly_sum(poly P, poly Q);

/// @brief Computes @c P-Q.
///
/// @param P a polynomial
/// @param Q another polynomial
///
/// @return the resulted polynomial, @c NULL if some error occurred.
poly poly_diff(poly P, poly Q);

/// @brief Computes @c P*Q.
///
/// @param P a polynomial
/// @param Q another polynomial
///
/// @return the resulted polynomial, @c NULL if some error occurred.
poly poly_prod(poly P, poly Q);

/// @brief Computes the square of @c P.
///
/// @param P the polynomial
///
/// @return the resulted polynomial, @c NULL if some error occurred.
poly poly_sqr(poly P);

#endif /* poly_h */
