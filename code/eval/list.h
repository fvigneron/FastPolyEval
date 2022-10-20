//
//  list.h
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
 \file list.h
 \brief Definition of a list that can quickly sort real machine floating-point numbers and keep track of their permutation.
*/

#ifndef list_h
#define list_h

#include "ntypes.h"

/// The multiplicative factor by which to grow the lists.
#define LIST_FACTOR 1.2

/// The additive term by which to grow the lists.
#define LIST_TERM 100

 // MARK: data types

 /// @struct list_struct
 /// @brief A list of real numbers that can be sorted, while keeping track of original order.
 typedef struct {
     deg_t size;     ///< the memory size allocated
     deg_t count;    ///< the number of elements
     deg_t *k;       ///< indexes, powers, the permutation of the numbers when sorted
     real_t *s;      ///< the list of numbers
     bool sorted;    ///< the state of the list
 } list_struct;

 /// @brief Practical wrapper for @c list_struct
 ///
 /// To avoid the constant use @c * and @c & the type  @c list_t is a pointer.
 typedef list_struct list_t[1];

 /// Convenience pointer to @c list_struct
 typedef list_struct *list;

 // MARK: macros & functions prototypes

/// @brief Returns a new list of machine real numbers, with storage space for @c size numbers.
///
/// @param size the size of the list
///
/// @return the new list, @c NULL if some error occurred.
list list_new(deg_t size);

/// @brief Initializes an existing list @c l  with storage space for @c size numbers.
///
/// @param l the list
/// @param size the size of the list
///
/// @return @ref true if successfull, @ref false otherwise.
bool list_init(list l, deg_t size);

/// @brief Frees all the memory used by the list @c l, assuming the struct has been allocated with
/// @c malloc(), for example with @c list_new().
///
/// @param l the list
///
/// @return @ref true if successfull, @ref false otherwise.
bool list_free(list l);

/// @brief Frees all the memory used by the list @c l, but not the list @c l itself.
///
/// @param l the list
///
/// @return @ref true if successfull, @ref false otherwise.
bool list_clear(list l);

/// @brief Returns a copy of the list @c l.
///
/// @param l the list
///
/// @return the new list, @c NULL if some error occurred.
list list_clone(list l);

/// @brief Adds the couple @c (k,s) at a new position at the end of the list @c l.
///
/// @param l the list
/// @param k the index
/// @param s the real number
///
/// @return @ref true if successfull, @ref false otherwise.
bool list_add(list l, deg_t k, real_t s);

/// @brief Trims the list @c l to minimal size to contain all its elements.
///
/// @param l the list
///
/// @return @ref true if successfull, @ref false otherwise.
bool list_trim(list l);

/// @brief Sorts the list in descending order while preserving couples @c (k[i],s[i]) for all @c i.
///
/// @warning The memory usage will double during the operation, all cleared up before returning. In case of failure,
/// the original list is destroyed.
///
/// @param l the list
///
/// @return @ref true if successfull, @ref false otherwise.
bool list_sort(list l);

#endif /* list_h */
