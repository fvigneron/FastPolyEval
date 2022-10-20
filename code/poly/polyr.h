//
//  polyr.h
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

/**
 \file polyr.h
 \brief Definition of real polynomials with arbitary precision coefficients.
*/

#ifndef polyr_h
#define polyr_h

#include <mpfr.h>

#include "ntypes.h"
#include "comp.h"

// MARK: data types

/// @struct polyr_struct
/// @brief Polynomial with multi-precision floating point complex coefficients.
typedef struct {
    deg_t degree;   ///< the degree of the polynomial
    prec_t prec;    ///< the precision of the coefficients, in bits
    mpfr_ptr a;     ///< the coefficients
    bool modified;  ///< the status of the coefficients
    mpfr_t buf1;    ///< a buffer
    mpfr_t buf2;    ///< another buffer
} polyr_struct;

/// @brief Practical wrapper for @c polyr_struct
///
/// To avoid the constant use @c * and @c & the type  @c polyr_t is a pointer.
typedef polyr_struct polyr_t[1];

/// Convenience pointer to @c polyr_struct
typedef polyr_struct *polyr;

// MARK: functions prototypes

/// @brief Returns a new real polynomial of given @c degree, with coefficients of precision @c prec.
///
/// @warning @c prec must be at least @c precf and the degree at most @c MAX_DEG.
///
/// @param degree the degree of the polynomial
/// @param prec the precision of the coefficients, in bits
///
/// @return the new buffer, @c NULL if some error occurred.
polyr polyr_new(deg_t degree, prec_t prec);

/// @brief Frees all the memory used by the polynomial @c P, assuming the struct has been allocated with
/// @c malloc(), for example with @c polyr_new().
///
/// @param P the polynomial
///
/// @return @ref true if successfull, @ref false otherwise.
bool polyr_free(polyr P);

/// @brief Sets the coefficient of the polynomial @c P corresponding to the power @c ind to @c coeff.
///
/// @param P the polynomial
/// @param coeff the coefficient
/// @param ind the index
///
/// @return @ref true if successfull, @ref false otherwise.
bool polyr_set(polyr P, mpfr_t coeff, deg_t ind);

/// @brief Sets the coefficient of the polynomial @c P corresponding to the power @c ind to @c coeff.
///
/// @param P the polynomial
/// @param coeff the coefficient
/// @param ind the index
///
/// @return @ref true if successfull, @ref false otherwise.
bool polyr_seti(polyr P, long coeff, deg_t ind);

/// @brief Evaluates @c P(z) using Horner's method.
///
/// @param res the result
/// @param P the polynomial
/// @param z the argument
///
/// @return @ref true if successfull, @ref false otherwise.
bool polyr_eval_c(comp res, polyr P, comp z);

/// @brief Evaluates @c P(x) using Horner's method.
///
/// @param res the result
/// @param P the polynomial
/// @param x the real argument
///
/// @return @ref true if successfull, @ref false otherwise.
bool polyr_eval(mpfr_t res, polyr P, mpfr_t x);

/// @brief Computes the derivative of @c P.
///
/// @param P the polynomial
///
/// @return the resulted polynomial, @c NULL if some error occurred.
polyr polyr_derivative(polyr P);

/// @brief Computes @c P+Q.
///
/// @param P a polynomial
/// @param Q another polynomial
///
/// @return the resulted polynomial, @c NULL if some error occurred.
polyr polyr_sum(polyr P, polyr Q);

/// @brief Computes @c P-Q.
///
/// @param P a polynomial
/// @param Q another polynomial
///
/// @return the resulted polynomial, @c NULL if some error occurred.
polyr polyr_diff(polyr P, polyr Q);

/// @brief Computes @c P*Q.
///
/// @param P a polynomial
/// @param Q another polynomial
///
/// @return the resulted polynomial, @c NULL if some error occurred.
polyr polyr_prod(polyr P, polyr Q);

/// @brief Computes the square of @c P.
///
/// @param P the polynomial
///
/// @return the resulted polynomial, @c NULL if some error occurred.
polyr polyr_sqr(polyr P);

/// @brief Computes the @c n-th hyperbolic polynomial, the @c n-th image of @c 0
/// under the iteration of @c z->z^2+c. It is a polynomial of degree @c 2^{n-1} in @c c.
///
/// @param n the order of the polynomial
/// @param prec the precision of the coefficients
///
/// @return the hyperbolic polynomial, @c NULL if some error occurred.
polyr poly_hyp(int n, prec_t prec);

/// @brief Computes the Chebyshev polynomial of degree @c n.
///
/// @param n the degree of the polynomial
/// @param prec the precision of the coefficients
///
/// @return the Chebyshev polynomial, @c NULL if some error occurred.
polyr poly_cheb(int n, prec_t prec);

/// @brief Computes the Legendre polynomial of degree @c n.
///
/// @param n the degree of the polynomial
/// @param prec the precision of the coefficients
///
/// @return the Legendre polynomial, @c NULL if some error occurred.
polyr poly_leg(int n, prec_t prec);

/// @brief Computes the Hermite polynomial of degree @c n.
///
/// @param n the degree of the polynomial
/// @param prec the precision of the coefficients
///
/// @return the Hermite polynomial, @c NULL if some error occurred.
polyr poly_her(int n, prec_t prec);

/// @brief Computes the Laguerre polynomial of degree @c n.
///
/// @param n the degree of the polynomial
/// @param prec the precision of the coefficients
///
/// @return the Laguerre polynomial, @c NULL if some error occurred.
polyr poly_lag(int n, prec_t prec);

#endif /* polyr_h */
