//
//  eval.h
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
 \file eval.h
 \brief Definition of polynomial evaluator with arbitary precision coefficients.
*/

#ifndef eval_h
#define eval_h

#include <mpfr.h>

#include "ntypes.h"
#include "comp.h"
#include "concave.h"
#include "poly.h"
#include "polyr.h"
#include "pows.h"
#include "powsr.h"

// MARK: constants &data types

/// Minimum prcision for evaluators.
#define MIN_EVAL_PREC      8
#define MIN_EVAL_PREC_STR "8"

/// @struct eval_struct
/// @brief Evaluator of polynomials with multi-precision floating point coefficients.
typedef struct {
    poly P;         ///< the polynomial, with complex coefficients
    polyr Q;        ///< the polynomial, with real coefficients
    bool real;      ///< the type of polynomial to evaluate
    prec_t prec;    ///< the precision of evaluations, in bits
    concave f;      ///< concave cover of the scales of coefficients
    pows zn;        ///< powers of a complex argument
    powsr xn;       ///< powers of a real argument
    real_t valErr;  ///< the [approximative] upper bound for the absolute error of the last evaluation, in bits
    real_t derErr;  ///< the [approximative] upper bound for the absolute error of the last derivative evaluation, in bits
    real_t ntErr;   ///< the [approximative] upper bound for the absolute error of the last Newton term evaluation, in bits
    deg_t terms;    ///< the number of polynomial terms computed by the last operation
    comp buf;       ///< a buffer for internal computations
    mpfr_t br;      ///< a real buffer for internal computations
} eval_struct;

/// @brief Practical wrapper for @c eval_struct
///
/// To avoid the constant use @c * and @c & the type  @c eval_t is a pointer.
typedef eval_struct eval_t[1];

/// Convenience pointer to @c eval_struct
typedef eval_struct *eval;

// MARK: macros & functions prototypes

/// The [approximative] upper bound for the error of the last evaluation, in bits. @b Not reported by the Newton term computation.
#define eval_lastValError(ev) (ev->valErr)

/// The [approximative] upper bound for the error of the last derivative evaluation, in bits.  @b Not reported by the Newton term computation.
#define eval_lastDerError(ev) (ev->derErr)

/// @brief Returns a new evaluator of the complex polynomial @c P.
///
/// @param P the complex polynomial
/// @param prec the precision of evaluations, in bits
///
/// @return the new evaluator, @c NULL if some error occurred.
eval eval_new(poly P, prec_t prec);

/// @brief Returns a new evaluator of the real polynomial @c Q.
///
/// @param Q the real polynomial
/// @param prec the precision of evaluations, in bits
///
/// @return the new evaluator, @c NULL if some error occurred.
eval eval_new_r(polyr Q, prec_t prec);

/// @brief Frees all the memory used by the evaluator @c ev, assuming the struct has been allocated with
/// @c malloc(), for example with @c eval_new() or with @c eval_new_r().
///
/// @param ev the evaluator
///
/// @return @ref true if successfull, @ref false otherwise.
bool eval_free(eval ev);

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
bool eval_val(comp v, eval ev, comp z);

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
bool eval_der(comp d, eval ev, comp z);

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
bool eval_val_der(comp v, comp d, eval ev, comp z);

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
bool eval_newton(comp nt, eval ev, comp z);

// MARK: Bare-metal variants

/// @brief Analyses the complex polynomial @c ev->P, after some of its coefficients have been changed.
///
/// @note The general functions @c eval_val(), @c eval_der(), @c eval_val_der() and
/// @c eval_newton() automatically analyse the appropriate polynomial. Use this function only with
/// the optimized versions like @c eval_val_cx().
///
/// @param ev the evaluator
///
/// @return @ref true if @c ev is ready to use, @ref false otherwise.
bool eval_analyse(eval ev);

/// @brief Analyses the real polynomial @c ev->Q, after some of its coefficients have been changed.
///
/// @note The general functions @c eval_val(), @c eval_der(), @c eval_val_der() and
/// @c eval_newton() automatically analyse the appropriate polynomial. Use this function only with
/// the optimized versions like @c eval_val_rx().
///
/// @param ev the evaluator
///
/// @return @ref true if @c ev is ready to use, @ref false otherwise.
bool eval_analyse_r(eval ev);

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
bool eval_val_cc(comp v, eval ev, comp z);

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
bool eval_der_cc(comp d, eval ev, comp z);

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
bool eval_val_der_cc(comp v, comp d, eval ev, comp z);

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
bool eval_newton_cc(comp nt, eval ev, comp z);

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
bool eval_val_cr(comp v, eval ev, mpfr_t x);

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
bool eval_der_cr(comp d, eval ev, mpfr_t x);

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
bool eval_val_der_cr(comp v, comp d, eval ev, mpfr_t x);

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
bool eval_newton_cr(comp nt, eval ev, mpfr_t x);

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
bool eval_val_rc(comp v, eval ev, comp z);

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
bool eval_der_rc(comp d, eval ev, comp z);

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
bool eval_val_der_rc(comp v, comp d, eval ev, comp z);

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
bool eval_newton_rc(comp nt, eval ev, comp z);

/// @brief Evaluates the real polynomial @c ev->Q(x) using the method described in [1].
///
/// @warning For maximum speed, no checks are performed on the parameters, on the type
/// of polynomial represented by @c ev nor if the polynomial has been pre-conditionned after
/// the last coefficient update.
///
/// @param v the value
/// @param ev the evaluator
/// @param x the real argument
///
/// @return @ref true if successfull, @ref false otherwise.
bool eval_val_rr(mpfr_t v, eval ev, mpfr_t x);

/// @brief Evaluates the derivative of real polynomial @c ev->Q'(x) using the method described in [1].
///
/// @warning For maximum speed, no checks are performed on the parameters, on the type
/// of polynomial represented by @c ev nor if the polynomial has been pre-conditionned after
/// the last coefficient update.
///
/// @param d the derivative
/// @param ev the evaluator
/// @param x the real argument
///
/// @return @ref true if successfull, @ref false otherwise.
bool eval_der_rr(mpfr_t d, eval ev, mpfr_t x);

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
bool eval_val_der_rr(mpfr_t v, mpfr_t d, eval ev, mpfr_t x);

/// @brief Computes the Newthon method step of the real polynomial @c ev->Q using
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
bool eval_newton_rr(mpfr_t nt, eval ev, mpfr_t x);

#endif /* eval_h */
