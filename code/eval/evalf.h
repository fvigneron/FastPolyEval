//
//  evalf.h
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
 \file evalf.h
 \brief Definition of polynomials evaluator with machine floating point coefficients.
*/

#ifndef evalf_h
#define evalf_h

#include "ntypes.h"
#include "compf.h"
#include "concave.h"
#include "polyf.h"
#include "polyfr.h"
#include "powsf.h"
#include "powsfr.h"

// MARK: data types

/// @struct evalf_struct
/// @brief Evaluator of polynomials with machine floating point coefficients.
typedef struct {
    polyf P;        ///< the polynomial, with complex coefficients
    polyfr Q;       ///< the polynomial, with real coefficients
    bool real;      ///< the type of polynomial to evaluate
    concave f;      ///< concave cover of the magnitude of coefficients
    powsf zn;       ///< powers of a complex argument
    powsfr xn;      ///< powers of a real argument
    real_t valErr;  ///< the [approximative] upper bound for the absolute error of the last evaluation, in bits
    real_t derErr;  ///< the [approximative] upper bound for the absolute error of the last derivative evaluation, in bits
    real_t ntErr;   ///< the [approximative] upper bound for the absolute error of the last Newton term evaluation, in bits
    deg_t terms;    ///< the number of polynomial terms computed by the last operation
} evalf_struct;

/// @brief Practical wrapper for @c evalf_struct
///
/// To avoid the constant use @c * and @c & the type  @c evalf_t is a pointer.
typedef evalf_struct evalf_t[1];

/// Convenience pointer to @c evalf_struct
typedef evalf_struct *evalf;

// MARK: macros & functions prototypes

/// The [approximative] upper bound for the error of the last evaluation, in bits. @b Not reported by the Newton term computation.
#define evalf_lastValError(ev) (ev->valErr)

/// The [approximative] upper bound for the error of the last derivative evaluation, in bits.  @b Not reported by the Newton term computation.
#define evalf_lastDerError(ev) (ev->derErr)

/// @brief Returns a new evaluator of the complex polynomial @c P.
///
/// @param P the complex polynomial
///
/// @return the new evaluator, @c NULL if some error occurred.
evalf evalf_new(polyf P);

/// @brief Returns a new evaluator of the real polynomial @c Q.
///
/// @param Q the real polynomial
///
/// @return the new evaluator, @c NULL if some error occurred.
evalf evalf_new_r(polyfr Q);

/// @brief Frees all the memory used by the evaluator @c ev, assuming the struct has been allocated with
/// @c malloc(), for example with @c evalf_new() or with @c evalf_new_r().
///
/// @param ev the evaluator
///
/// @return @ref true if successfull, @ref false otherwise.
bool evalf_free(evalf ev);

/// @brief Evaluates @c ev->P(z) (or @c ev->Q(z)) using the method described in [1].
///
/// @note This function chooses the quickest variant to compute, depending if the polynomial of @c ev
/// and / or @c z are real.
///
/// @param v the result
/// @param ev the evaluator
/// @param z the complex argument
///
/// @return @ref true if successfull, @ref false otherwise.
bool evalf_val(compf v, evalf ev, compf z);

/// @brief Evaluates @c ev->P'(z) (or @c ev->Q'(z)) using the method described in [1].
///
/// @note This function chooses the quickest variant to compute, depending if the polynomial of @c ev
/// and / or @c z are real.
///
/// @param d the result
/// @param ev the evaluator
/// @param z the complex argument
///
/// @return @ref true if successfull, @ref false otherwise.
bool evalf_der(compf d, evalf ev, compf z);

/// @brief Evaluates @c ev->P(z) and @c ev->P'(z) (or @c ev->Q(z) and @c ev->Q'(z)) using
/// the method described in [1].
///
/// @note This function chooses the quickest variant to compute, depending if the polynomial of @c ev
/// and / or @c z are real.
///
/// @param v the value
/// @param d the derivative
/// @param ev the evaluator
/// @param z the complex argument
///
/// @return @ref true if successfull, @ref false otherwise.
bool evalf_val_der(compf v, compf d, evalf ev, compf z);

/// @brief Computes the Newthon method step of @c ev->P (or @c ev->Q) using
/// the method described in [1].
///
/// @note This function chooses the quickest variant to compute, depending if the polynomial of @c ev
/// and / or @c z are real.
///
/// @param nt the Newton term
/// @param ev the evaluator
/// @param z the complex argument
///
/// @return @ref true if successfull, @ref false otherwise.
bool evalf_newton(compf nt, evalf ev, compf z);

// MARK: Bare-metal variants

/// @brief Analyses the complex polynomial @c ev->P, after some of its coefficients have been changed.
///
/// @note The general functions @c evalf_val(), @c evalf_der(), @c evalf_val_der() and
/// @c evalf_newton() automatically analyse the appropriate polynomial. Use this function only with
/// the optimized versions like @c evalf_val_cx().
///
/// @param ev the evaluator
///
/// @return @ref true if @c ev is ready to use, @ref false otherwise.
bool evalf_analyse(evalf ev);

/// @brief Analyses the real polynomial @c ev->Q, after some of its coefficients have been changed.
///
/// @note The general functions @c evalf_val(), @c evalf_der(), @c evalf_val_der() and
/// @c evalf_newton() automatically analyse the appropriate polynomial. Use this function only with
/// the optimized versions like @c evalf_val_rx().
///
/// @param ev the evaluator
///
/// @return @ref true if @c ev is ready to use, @ref false otherwise.
bool evalf_analyse_r(evalf ev);

/// @brief Evaluates @c ev->P(x) using the method described in [1].
///
/// @warning For maximum speed, no checks are performed on the parameters, on the type
/// of polynomial represented by @c ev nor if the polynomial has been pre-conditionned after
/// the last coefficient update.
///
/// @param v the result
/// @param ev the evaluator
/// @param z the complex argument
///
/// @return @ref true if successfull, @ref false otherwise.
bool evalf_val_cc(compf v, evalf ev, compf z);

/// @brief Evaluates @c ev->P'(x) using the method described in [1].
///
/// @warning For maximum speed, no checks are performed on the parameters, on the type
/// of polynomial represented by @c ev nor if the polynomial has been pre-conditionned after
/// the last coefficient update.
///
/// @param d the result
/// @param ev the evaluator
/// @param z the complex argument
///
/// @return @ref true if successfull, @ref false otherwise.
bool evalf_der_cc(compf d, evalf ev, compf z);

/// @brief Evaluates @c ev->P(x) and @c ev->P'(x) using
/// the method described in [1].
///
/// @warning For maximum speed, no checks are performed on the parameters, on the type
/// of polynomial represented by @c ev nor if the polynomial has been pre-conditionned after
/// the last coefficient update.
///
/// @param v the value
/// @param d the derivative
/// @param ev the evaluator
/// @param z the complex argument
///
/// @return @ref true if successfull, @ref false otherwise.
bool evalf_val_der_cc(compf v, compf d, evalf ev, compf z);

/// @brief Computes the Newthon method step of @c ev->P using
/// the method described in [1].
///
/// @warning For maximum speed, no checks are performed on the parameters, on the type
/// of polynomial represented by @c ev nor if the polynomial has been pre-conditionned after
/// the last coefficient update.
///
/// @param nt the Newton term
/// @param ev the evaluator
/// @param z the complex argument
///
/// @return @ref true if successfull, @ref false otherwise.
bool evalf_newton_cc(compf nt, evalf ev, compf z);

/// @brief Evaluates @c ev->P(x) using the method described in [1].
///
/// @warning For maximum speed, no checks are performed on the parameters, on the type
/// of polynomial represented by @c ev nor if the polynomial has been pre-conditionned after
/// the last coefficient update.
///
/// @param v the result
/// @param ev the evaluator
/// @param x the real argument
///
/// @return @ref true if successfull, @ref false otherwise.
bool evalf_val_cr(compf v, evalf ev, coeff_t x);

/// @brief Evaluates @c ev->P'(x) using the method described in [1].
///
/// @warning For maximum speed, no checks are performed on the parameters, on the type
/// of polynomial represented by @c ev nor if the polynomial has been pre-conditionned after
/// the last coefficient update.
///
/// @param d the result
/// @param ev the evaluator
/// @param x the real argument
///
/// @return @ref true if successfull, @ref false otherwise.
bool evalf_der_cr(compf d, evalf ev, coeff_t x);

/// @brief Evaluates @c ev->P(x) and @c ev->P'(x) using
/// the method described in [1].
///
/// @warning For maximum speed, no checks are performed on the parameters, on the type
/// of polynomial represented by @c ev nor if the polynomial has been pre-conditionned after
/// the last coefficient update.
///
/// @param v the value
/// @param d the derivative
/// @param ev the evaluator
/// @param x the real argument
///
/// @return @ref true if successfull, @ref false otherwise.
bool evalf_val_der_cr(compf v, compf d, evalf ev, coeff_t x);

/// @brief Computes the Newthon method step of @c ev->P using
/// the method described in [1].
///
/// @warning For maximum speed, no checks are performed on the parameters, on the type
/// of polynomial represented by @c ev nor if the polynomial has been pre-conditionned after
/// the last coefficient update.
///
/// @param nt the Newton term
/// @param ev the evaluator
/// @param x the real argument
///
/// @return @ref true if successfull, @ref false otherwise.
bool evalf_newton_cr(compf nt, evalf ev, coeff_t x);

/// @brief Evaluates @c ev->Q(x) using the method described in [1].
///
/// @warning For maximum speed, no checks are performed on the parameters, on the type
/// of polynomial represented by @c ev nor if the polynomial has been pre-conditionned after
/// the last coefficient update.
///
/// @param v the result
/// @param ev the evaluator
/// @param z the complex argument
///
/// @return @ref true if successfull, @ref false otherwise.
bool evalf_val_rc(compf v, evalf ev, compf z);

/// @brief Evaluates @c ev->Q'(x) using the method described in [1].
///
/// @warning For maximum speed, no checks are performed on the parameters, on the type
/// of polynomial represented by @c ev nor if the polynomial has been pre-conditionned after
/// the last coefficient update.
///
/// @param d the result
/// @param ev the evaluator
/// @param z the complex argument
///
/// @return @ref true if successfull, @ref false otherwise.
bool evalf_der_rc(compf d, evalf ev, compf z);

/// @brief Evaluates @c ev->Q(x) and @c ev->Q'(x) using
/// the method described in [1].
///
/// @warning For maximum speed, no checks are performed on the parameters, on the type
/// of polynomial represented by @c ev nor if the polynomial has been pre-conditionned after
/// the last coefficient update.
///
/// @param v the value
/// @param d the derivative
/// @param ev the evaluator
/// @param z the complex argument
///
/// @return @ref true if successfull, @ref false otherwise.
bool evalf_val_der_rc(compf v, compf d, evalf ev, compf z);

/// @brief Computes the Newthon method step of @c ev->Q using
/// the method described in [1].
///
/// @warning For maximum speed, no checks are performed on the parameters, on the type
/// of polynomial represented by @c ev nor if the polynomial has been pre-conditionned after
/// the last coefficient update.
///
/// @param nt the Newton term
/// @param ev the evaluator
/// @param z the complex argument
///
/// @return @ref true if successfull, @ref false otherwise.
bool evalf_newton_rc(compf nt, evalf ev, compf z);

/// @brief Evaluates the real polynomial @c ev->Q(x) using the method described in [1].
///
/// @warning For maximum speed, no checks are performed on the parameters, on the type
/// of polynomial represented by @c ev nor if the polynomial has been pre-conditionned after
/// the last coefficient update.
///
/// @param ev the evaluator
/// @param x the real argument
///
/// @return the result, @c NaN is some error occurred.
coeff_t evalf_val_rr(evalf ev, coeff_t x);

/// @brief Evaluates the derivative of real polynomial @c ev->Q'(x) using the method described in [1].
///
/// @warning For maximum speed, no checks are performed on the parameters, on the type
/// of polynomial represented by @c ev nor if the polynomial has been pre-conditionned after
/// the last coefficient update.
///
/// @param ev the evaluator
/// @param x the real argument
///
/// @return the result, @c NaN is some error occurred.
coeff_t evalf_der_rr(evalf ev, coeff_t x);

/// @brief Evaluates @c ev->Q(x) and @c ev->Q'(x) using
/// the method described in [1].
///
/// @warning For maximum speed, no checks are performed on the parameters, on the type
/// of polynomial represented by @c ev nor if the polynomial has been pre-conditionned after
/// the last coefficient update.
///
/// @param v the value
/// @param d the derivative
/// @param ev the evaluator
/// @param x the real argument
///
/// @return @ref true if successfull, @ref false otherwise.
bool evalf_val_der_rr(coeff_t *v, coeff_t *d, evalf ev, coeff_t x);

/// @brief Computes the Newthon method step of the real polynomial @c ev->Q using
/// the method described in [1].
///
/// @warning For maximum speed, no checks are performed on the parameters, on the type
/// of polynomial represented by @c ev nor if the polynomial has been pre-conditionned after
/// the last coefficient update. 
///
/// @param ev the evaluator
/// @param x the real argument
///
/// @return the result, @c NaN is some error occurred.
coeff_t evalf_newton_rr(evalf ev, coeff_t x);

#endif /* evalf_h */
