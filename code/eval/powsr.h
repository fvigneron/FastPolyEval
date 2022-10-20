//
//  powsr.h
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
 \file powsr.h
 \brief Definition of a buffer for pre-computed powers of a real number with arbitary precision.
*/

#ifndef powsr_h
#define powsr_h

#include <mpfr.h>

#include "ntypes.h"

// MARK: data types

// comment this line to use repeated squares method and buffers to compute powers of @c x
// leave uncommented tu use a direct computation with @c mpfr_pow_ui()
#define POWSR_DIRECT_POWERS

#ifdef POWSR_DIRECT_POWERS

/// @struct powsr_struct
/// @brief The powers of the real number @c x using multi-precision floating point numbers.
typedef struct {
    prec_t prec;    ///< the precision of the powers of @c x, in bits
    bool inited;    ///< the status of the value @c x
    mpfr_t x;       ///< the real number
    mpfr_t res;     ///< a buffer
} powsr_struct;

#else

/// @struct powsr_struct
/// @brief The powers of the real number @c x using multi-precision floating point numbers.
typedef struct {
    prec_t prec;    ///< the precision of the powers of @c x, in bits
    deg_t size;     ///< the memory size allocated
    byte tps;       ///< the largest non-negative integer such that @c 2^tps<=size
    bool *computed; ///< the status of powers
    mpfr_ptr xn;    ///< the powers of @c x
    mpfr_t pth;     ///< a buffer
    mpfr_t res;     ///< another buffer
} powsr_struct;

#endif

/// @brief Practical wrapper for @c powsr_struct
///
/// To avoid the constant use @c * and @c & the type  @c powsr_t is a pointer.
typedef powsr_struct powsr_t[1];

/// Convenience pointer to @c eval_struct
typedef powsr_struct *powsr;

// MARK: macros & functions prototypes

/// @brief Returns a new buffer of powers of real numbers of precision @c prec,
/// with initial storage space for @c size powers.
///
/// @warning @c prec must be at least @c precf and @c size at most @c MAX_DEG.
///
/// @param prec the precision of the coefficients, in bits
/// @param size thessize of the buffer
///
/// @return the new buffer, @c NULL if some error occurred.
powsr powsr_new(prec_t prec, deg_t size);

/// @brief Frees all the memory used by the buffer @c xn, assuming the struct has been allocated with
/// @c malloc(), for example with @c powsr_new().
///
/// @param xn the powers buffer
///
/// @return @ref true if successfull, @ref false otherwise.
bool powsr_free(powsr xn);

/// @brief Sets the real number of which the powers will be computed by the buffer @c xn.
///
/// @param xn the powers buffer
/// @param x the complex number
///
/// @return @ref true if successfull, @ref false otherwise.
bool powsr_set(powsr xn, mpfr_t x);

/// @brief Computes @c x^pow using repeated squares method and the cache of previously computed powers.
///
/// It caches intermediary powers of @c x that have been computed to accelerate later calls of this function. Also, if
/// @c pow is larger than the size of the buffer @c xn, it is automatically increased to store the result.
///
/// @param xn the powers buffer
/// @param pow the power to compute
///
/// @return the result @c x^pow, @c NULL if some error occurred.
mpfr_ptr powsr_pow(powsr xn, deg_t pow);

/// @brief Computes @c x^pow using repeated squares method and the cache of previously computed powers.
///
/// It does @b NOT cache intermediary powers of @c x into @c xn.
///
/// @param xn the powers buffer
/// @param pow the power to compute
///
/// @return the result @c x^pow, @c NULL if some error occurred.
mpfr_ptr powsr_pow_once(powsr xn, deg_t pow);

#endif /* powsr_h */
