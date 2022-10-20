//
//  powsf.h
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
 \file powsf.h
 \brief Definition of a buffer for pre-computed powers of a machine complex number.
*/

#ifndef powsf_h
#define powsf_h

#include "ntypes.h"
#include "compf.h"

// MARK: data types

/// @struct powsf_struct
/// @brief The powers of the complex number @c z using machine floating point numbers.
typedef struct {
    deg_t size;     ///< the memory size allocated
    byte tps;       ///< the largest non-negative integer such that @c 2^tps<=size
    bool *computed; ///< the status of powers
    compf_ptr zn;   ///< the powers of @c z
    compf res;      ///< a buffer for results
} powsf_struct;

/// @brief Practical wrapper for @c powsf_struct
///
/// To avoid the constant use @c * and @c & the type  @c powsf_t is a pointer.
typedef powsf_struct powsf_t[1];

/// Convenience pointer to @c eval_struct
typedef powsf_struct *powsf;

// MARK: macros & functions prototypes

/// @brief Returns a new buffer of powers of machine complex numbers, with initial storage space for @c size powers.
///
/// @warning @c size must be at most @c MAX_DEG.
///
/// @param size the size of the buffer
///
/// @return the new buffer, @c NULL if some error occurred.
powsf powsf_new(deg_t size);

/// @brief Frees all the memory used by the buffer @c zn, assuming the struct has been allocated with
/// @c malloc(), for example with @c powsf_new().
///
/// @param zn the powers buffer
///
/// @return @ref true if successfull, @ref false otherwise.
bool powsf_free(powsf zn);

/// @brief Sets the complex number of which the powers will be computed by the buffer @c zn.
///
/// @param zn the powers buffer
/// @param z the complex number
///
/// @return @ref true if successfull, @ref false otherwise.
bool powsf_set(powsf zn, compf z);

/// @brief Computes @c z^pow using repeated squares method and the cache of previously computed powers.
///
/// It caches intermediary powers of @c z that have been computed to accelerate later calls of this function. Also, if
/// @c pow is larger than the size of the buffer @c zn, it is automatically increased to store the result.
///
/// @warning fails if pow==0
///
/// @param zn the powers buffer
/// @param pow the power to compute
///
/// @return the result @c z^pow, @c NULL if some error occurred.
compf_ptr powsf_pow(powsf zn, deg_t pow);

/// @brief Computes @c z^pow using repeated squares method and the cache of previously computed powers.
///
/// It does @b NOT cache intermediary powers of @c z into @c zn.
///
/// @param zn the powers buffer
/// @param pow the power to compute
///
/// @return the result @c z^pow, @c NULL if some error occurred.
compf_ptr powsf_pow_once(powsf zn, deg_t pow);

#endif /* powsf_h */
