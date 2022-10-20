//
//  powsfr.h
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
 \file powsfr.h
 \brief Definition of a buffer for pre-computed powers of a machine real number.
*/

#ifndef powsfr_h
#define powsfr_h

#include "ntypes.h"

// MARK: data types

// comment this line to use repeated squares method and buffers to compute powers of @c x
// leave uncommented tu use a direct computation with @c pow() or @c powl()
#define POWSFR_DIRECT_POWERS

#ifdef POWSFR_DIRECT_POWERS

/// @struct powsfr_struct
/// @brief The powers of the real number @c x using machine floating point numbers.
typedef struct {
    bool inited;  ///< the status of the value @c x
    coeff_t x;    ///< the real number
} powsfr_struct;

#else

/// @struct powsfr_struct
/// @brief The powers of the real number @c x using machine floating point numbers.
typedef struct {
    deg_t size;     ///< the memory size allocated
    byte tps;       ///< the largest non-negative integer such that @c 2^tps<=size
    bool *computed; ///< the status of powers
    coeff_t *xn;    ///< the powers of @c x
} powsfr_struct;

#endif

/// @brief Practical wrapper for @c powsfr_struct
///
/// To avoid the constant use @c * and @c & the type  @c powsfr_t is a pointer.
typedef powsfr_struct powsfr_t[1];

/// Convenience pointer to @c eval_struct
typedef powsfr_struct *powsfr;

// MARK: macros & functions prototypes

/// @brief Returns a new buffer of powers of machine real numbers, with initial storage space for @c size powers.
///
/// @warning @c size must be at most @c MAX_DEG.
///
/// @param size the size of the buffer
///
/// @return the new buffer, @c NULL if some error occurred.
powsfr powsfr_new(deg_t size);

/// @brief Frees all the memory used by the buffer @c xn, assuming the struct has been allocated with
/// @c malloc(), for example with @c powsfr_new().
///
/// @param xn the powers buffer
///
/// @return @ref true if successfull, @ref false otherwise.
bool powsfr_free(powsfr xn);

/// @brief Sets the real number of which the powers will be computed by the buffer @c xn.
///
/// @param xn the powers buffer
/// @param x the real number
///
/// @return @ref true if successfull, @ref false otherwise.
bool powsfr_set(powsfr xn, coeff_t x);

/// @brief Computes @c x^pow using repeated squares method and the cache of previously computed powers.
///
/// It caches intermediary powers of @c x that have been computed to accelerate later calls of this function. Also, if
/// @c pow is larger than the size of the buffer @c xn, it is automatically increased to store the result.
///
/// @param xn the powers buffer
/// @param pow the power to compute
///
/// @return the result @c x^pow, @c NaN if some error occurred.
coeff_t powsfr_pow(powsfr xn, deg_t pow);

/// @brief Computes @c x^pow using repeated squares method and the cache of previously computed powers.
///
/// It does @b NOT cache intermediary powers of @c x.
///
/// @param xn the powers buffer
/// @param pow the power to compute
///
/// @return the result @c x^pow, @c NULL if some error occurred.
coeff_t powsfr_pow_once(powsfr xn, deg_t pow);

#endif /* powsfr_h */
