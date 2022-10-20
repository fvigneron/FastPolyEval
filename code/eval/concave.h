//
//  concave.h
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
 \file concave.h
 \brief Definition of a concave function computed from the coefficients of some polynomial.
 
  Admissible powers are for any of the original polynomial @c P or its derivative @c P'.
*/

#ifndef concave_h
#define concave_h

#include "ntypes.h"
#include "list.h"

// MARK: data types

/// The extra bits that are needed as explained in [1].
#define BITS_GUARD 6

/// @struct concave_struct
/// @brief Description of a concave function computed from the coefficients of some polynomial.
typedef struct {
    prec_t extraBits; ///< extra bits for guarding, depending on the degree of the polynomial
    prec_t prec;      ///< the precision that will be used to evaluate the polynomial, excluding @c extraBits
    list_t def;       ///< the definition of the concave map
    list_t all;       ///< all terms of the polynomial that may be used for evaluation
    deg_t start;      ///< the position in @c allPow of the largest power to evaluate (for given slope)
    deg_t mid;        ///< the position in @c allPow of the power which gives the maximum modulus
    deg_t end;        ///< the position in @c allPow of the least power to evaluate (for given slope)
} concave_struct;

/// Convenience pointer to @c concave_struct
typedef concave_struct *concave;

// MARK: macros & functions prototypes

/// @brief Computes the concave cover of the graph of a list @c l sorted in decreasing order.
///
/// @note Retains the list of non-zero coefficients in @c all.
///
/// @param l the sorted list
/// @param prec the precision that will be used to evaluate the polynomial
///
/// @return the concave function above the graph of @c l, @c NULL if some error occurred.
concave conc_new(list l, prec_t prec);

/// @brief Frees all the memory used by the concave map @c f, assuming the struct has been allocated with
/// @c malloc(), for example with @c conc_new().
///
/// @param f the map
///
/// @return @ref true if successfull, @ref false otherwise.
bool conc_free(concave f);

/// @brief Computes the range of indexes in @c f->allPow to use for evaluating the original polynomial at @c z with
/// @c la=log_2|z|.
///
/// Stores the result in the @c f->start and respectively @c f->end.
///
/// @note @c f->start>=f->end, otherwise the value of the polynomial is @c 0.
///
/// @param f the map 
/// @param la the slope, or @c log_2|z|
///
/// @return @ref true if successfull, @ref false otherwise.
bool conc_range(concave f, real_t la);

/// @brief Computes the range of indexes in @c f->allPow to use for evaluating the derivative of the
/// original polynomial at @c z with @c la=log_2|z|.
///
/// Stores the result in the @c f->start and respectively @c f->end.
///
/// @note @c f->start>=f->end, otherwise the value of the derivative is @c 0.
///
/// @param f the map
/// @param la the slope, or @c log_2|z|
///
/// @return @ref true if successfull, @ref false otherwise.
bool conc_range_der(concave f, real_t la);

#endif /* conc_h */
