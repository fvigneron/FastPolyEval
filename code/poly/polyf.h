//
//  polyf.h
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
 \file polyf.h
 \brief Definition of complex polynomials with machine floating point coefficients.
*/

#ifndef polyf_h
#define polyf_h

#include "ntypes.h"
#include "compf.h"

// MARK: data types

/// @struct polyf_struct
/// @brief Polynomial with machine floating point complex coefficients.
typedef struct {
    deg_t degree;   ///< the degree of the polynomial
    compf_ptr a;    ///< the coefficients
    bool modified;  ///< the status of the coefficients
} polyf_struct;

/// @brief Practical wrapper for @c polyf_struct
///
/// To avoid the constant use @c * and @c & the type  @c polyf_t is a pointer.
typedef polyf_struct polyf_t[1];

/// Convenience pointer to @c polyf_struct
typedef polyf_struct *polyf;

// MARK: functions prototypes

/// @brief Returns a new complex polynomial of given @c degree, with machine floating point coefficients.
///
/// @param degree the degree of the polynomial
///
/// @return the new polynomial, @c NULL if the degree is larger than @c MAX_DEG.
polyf polyf_new(deg_t degree);

/// @brief Returns a new complex polynomial given the list of its @c roots, with machine floating point coefficients.
///
/// @param roots the roots of the polynomial
/// @param degree the degree of the polynomial
///
/// @return the new polynomial, @c NULL if some error occurred.
polyf polyf_from_roots(compf_ptr roots, deg_t degree);

/// @brief Frees all the memory used by the polynomial @c P, assuming the struct has been allocated with
/// @c malloc(), for example with @c polyr_new().
///
/// @param P the polynomial
///
/// @return @ref true if successfull, @ref false otherwise.
bool polyf_free(polyf P);

/// @brief Sets the coefficient of the polynomial @c P corresponding to the power @c ind to @c coeff.
///
/// @param P the polynomial
/// @param coeff the coefficient
/// @param ind the index
///
/// @return @ref true if successfull, @ref false otherwise.
bool polyf_set(polyf P, compf coeff, deg_t ind);

/// @brief Evaluates @c P(z) using Horner's method.
///
/// @param res the result
/// @param P the polynomial
/// @param z the argument
///
/// @return @ref true if successfull, @ref false otherwise.
bool polyf_eval(compf res, polyf P, compf z);

/// @brief Evaluates @c P(x) using Horner's method.
///
/// @param res the result
/// @param P the polynomial
/// @param x the real argument
///
/// @return @ref true if successfull, @ref false otherwise.
bool polyf_eval_r(compf res, polyf P, coeff_t x);

/// @brief Computes the derivative of @c P.
///
/// @param P the polynomial
///
/// @return the resulted polynomial, @c NULL if some error occurred.
polyf polyf_derivative(polyf P);

/// @brief Computes @c P+Q.
///
/// @param P a polynomial
/// @param Q another polynomial
///
/// @return the resulted polynomial, @c NULL if some error occurred.
polyf polyf_sum(polyf P, polyf Q);

/// @brief Computes @c P-Q.
///
/// @param P a polynomial
/// @param Q another polynomial
///
/// @return the resulted polynomial, @c NULL if some error occurred.
polyf polyf_diff(polyf P, polyf Q);

/// @brief Computes @c P*Q.
///
/// @param P a polynomial
/// @param Q another polynomial
///
/// @return the resulted polynomial, @c NULL if some error occurred.
polyf polyf_prod(polyf P, polyf Q);

/// @brief Computes the square of @c P.
///
/// @param P the polynomial
///
/// @return the resulted polynomial, @c NULL if some error occurred.
polyf polyf_sqr(polyf P);

#endif /* polyf_h */
