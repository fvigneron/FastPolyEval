//
//  pows.h
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
 \file pows.h
 \brief Definition of a buffer for pre-computed powers of a complex number with arbitary precision.
*/

#ifndef pows_h
#define pows_h

#include "ntypes.h"
#include "comp.h"

// MARK: data types

/// @struct pows_struct
/// @brief The powers of the complex number @c z using multi-precision floating point numbers.
typedef struct {
    prec_t prec;    ///< the precision of the powers of @c z, in bits
    deg_t size;     ///< the memory size allocated
    byte tps;       ///< the largest non-negative integer such that @c 2^tps<=size
    bool *computed; ///< the status of powers
    comp_ptr zn;    ///< the powers of @c z
    mpfr_t buf1;    ///< a buffer
    mpfr_t buf2;    ///< another buffer
    comp pth;       ///< a buffer
    comp res;       ///< another buffer
} pows_struct;

/// @brief Practical wrapper for @c pows_struct
///
/// To avoid the constant use @c * and @c & the type  @c pows_t is a pointer.
typedef pows_struct pows_t[1];

/// Convenience pointer to @c eval_struct
typedef pows_struct *pows;

// MARK: functions prototypes

/// @brief Returns a new buffer of powers of complex numbers of precision @c prec,
/// with initial storage space for @c size powers.
///
/// @warning @c prec must be at least @c precf and @c size at most @c MAX_DEG.
///
/// @param prec the precision of the coefficients, in bits
/// @param size the size of the buffer
///
/// @return the new buffer, @c NULL if some error occurred.
pows pows_new(prec_t prec, deg_t size);

/// @brief Frees all the memory used by the buffer @c zn, assuming the struct has been allocated with
/// @c malloc(), for example with @c pows_new().
///
/// @param zn the powers buffer
///
/// @return @ref true if successfull, @ref false otherwise.
bool pows_free(pows zn);

/// @brief Sets the complex number of which the powers will be computed by the buffer @c zn.
///
/// @param zn the powers buffer
/// @param z the complex number
///
/// @return @ref true if successfull, @ref false otherwise.
bool pows_set(pows zn, comp z);

/// @brief Computes @c z^pow using repeated squares method and the cache of previously computed powers.
///
/// It caches intermediary powers of @c z that have been computed to accelerate later calls of this function. Also, if
/// @c pow is larger than the size of the buffer @c zn, it is automatically increased to store the result.
///
/// @param zn the powers buffer
/// @param pow the power to compute
///
/// @return the result @c z^pow, @c NULL if some error occurred.
comp_ptr pows_pow(pows zn, deg_t pow);

/// @brief Computes @c z^pow using repeated squares method and the cache of previously computed powers.
///
/// It does @b NOT cache intermediary powers of @c z into @c zn.
///
/// @param zn the powers buffer
/// @param pow the power to compute
///
/// @return the result @c z^pow, @c NULL if some error occurred.
comp_ptr pows_pow_once(pows zn, deg_t pow);

#endif /* pows_h */
