//
//  polyfr.h
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
 \file polyfr.h
 \brief Definition of real polynomials with machine floating point coefficients.
*/

#ifndef polyfr_h
#define polyfr_h

#include "ntypes.h"
#include "compf.h"
#include "polyf.h"

// MARK: data types

/// @struct polyfr_struct
/// @brief Polynomial with machine floating point real coefficients.
typedef struct {
    deg_t degree;   ///< the degree of the polynomial
    coeff_t *a;     ///< the coefficients
    bool modified;  ///< the status of the coefficients
} polyfr_struct;

/// @brief Practical wrapper for @c polyfr_struct
///
/// To avoid the constant use @c * and @c & the type  @c polyfr_t is a pointer.
typedef polyfr_struct polyfr_t[1];

/// Convenience pointer to @c polyfr_struct
typedef polyfr_struct *polyfr;

// MARK: functions prototypes

/// @brief Returns a new real polynomial of given @c degree, with machine floating point coefficients.
///
/// @param degree the degree of the polynomial
///
/// @return the new polynomial, @c NULL if the degree is larger than @c MAX_DEG.
polyfr polyfr_new(deg_t degree);

/// @brief Frees all the memory used by the polynomial @c P, assuming the struct has been allocated with
/// @c malloc(), for example with @c polyfr_new().
///
/// @param P the polynomial
///
/// @return @ref true if successfull, @ref false otherwise.
bool polyfr_free(polyfr P);

/// @brief Sets the coefficient corresponding to the power @c ind to @c coeff.
///
/// @param P the polynomial
/// @param coeff the coefficient
/// @param ind the index
///
/// @return @ref true if successfull, @ref false otherwise.
bool polyfr_set(polyfr P, coeff_t coeff, deg_t ind);

/// @brief Returns the coefficient corresponding to the power @c ind.
///
/// @param P the polynomial
/// @param ind the index
///
/// @return the coefficient corresponding to the power @c ind.
coeff_t polyfr_get(polyfr P, deg_t ind);

/// @brief Evaluates @c P(z) using Horner's method.
///
/// @param res the result
/// @param P the polynomial
/// @param z the argument
///
/// @return @ref true if successfull, @ref false otherwise.
bool polyfr_eval_c(compf res, polyfr P, compf z);

/// @brief Evaluates @c P(x) using Horner's method.
///
/// @param P the polynomial
/// @param x the argument
///
/// @return the result, @c NaN if some error occurred.
coeff_t polyfr_eval(polyfr P, coeff_t x);

/// @brief Return a complex version of the real polynomial @c P.
///
/// @param P a polynomial
///
/// @return the resulted polynomial, @c NULL if some error occurred.
polyf polyfr_comp(polyfr P);

/// @brief Computes the derivative of @c P.
///
/// @param P the polynomial
///
/// @return the resulted polynomial, @c NULL if some error occurred.
polyfr polyfr_derivative(polyfr P);

/// @brief Computes @c P+Q.
///
/// @param P a polynomial
/// @param Q another polynomial
///
/// @return the resulted polynomial, @c NULL if some error occurred.
polyfr polyfr_sum(polyfr P, polyfr Q);

/// @brief Computes @c P-Q.
///
/// @param P a polynomial
/// @param Q another polynomial
///
/// @return the resulted polynomial, @c NULL if some error occurred.
polyfr polyfr_diff(polyfr P, polyfr Q);

/// @brief Computes @c P*Q.
///
/// @param P a polynomial
/// @param Q another polynomial
///
/// @return the resulted polynomial, @c NULL if some error occurred.
polyfr polyfr_prod(polyfr P, polyfr Q);

/// @brief Computes the square of @c P.
///
/// @param P the polynomial
///
/// @return the resulted polynomial, @c NULL if some error occurred.
polyfr polyfr_sqr(polyfr P);

/// @brief Computes the @c n-th hyperbolic polynomial, the @c n-th image of @c 0
/// under the iteration of @c z->z^2+c. It is a polynomial of degree @c 2^{n-1} in @c c.
///
/// @param n the order of the polynomial
///
/// @return the hyperbolic polynomial, @c NULL if some error occurred.
polyfr polyf_hyp(int n);

/// @brief Computes the Chebyshev polynomial of degree @c n.
///
/// @param n the degree of the polynomial
///
/// @return the Chebyshev polynomial, @c NULL if some error occurred.
polyfr polyf_cheb(int n);

/// @brief Computes the Legendre polynomial of degree @c n.
///
/// @param n the degree of the polynomial
///
/// @return the Legendre polynomial, @c NULL if some error occurred.
polyfr polyf_leg(int n);

/// @brief Computes the Hermite polynomial of degree @c n.
///
/// @param n the degree of the polynomial
///
/// @return the Hermite polynomial, @c NULL if some error occurred.
polyfr polyf_her(int n);

/// @brief Computes the Laguerre polynomial of degree @c n.
///
/// @param n the degree of the polynomial
///
/// @return the Laguerre polynomial, @c NULL if some error occurred.
polyfr polyf_lag(int n);

#endif /* polyfr_h */
