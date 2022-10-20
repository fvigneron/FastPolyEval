//
//  arrayf.h
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
 \file arrayf.h
 \brief Definition of a variable length array of machine floating point complex numbers.
*/

#ifndef arrayf_h
#define arrayf_h

#include <stdio.h>
#include <stdlib.h>

#include "ntypes.h"
#include "compf.h"
#include "polyf.h"
#include "polyfr.h"

// MARK: constants & data types

#define LISTF_MIN_SIZE       100
#define LISTF_SIZE_INCREASE  1.25

/// @struct arrayf_struct
/// @brief A variable length list of machine floating point complex numbers.
typedef struct {
    ulong len;      ///< the number of complex numbers in the list
    ulong size;     ///< the present capacity of the list (memory size in number of complex numbers)
    compf_ptr zi;   ///< the elements of the list
} arrayf_struct;

/// @brief Practical wrapper for @c arrayf_struct
///
/// To avoid the constant use @c * and @c & the type  @c arrayf_t is a pointer.
typedef arrayf_struct arrayf_t[1];

/// Convenience pointer to @c arrayf_struct
typedef arrayf_struct *arrayf;

// MARK: functions prototypes

/// @brief Returns a new empty list with memory size at least @c size.
///
/// @param size the number of elements to allocate memory for
///
/// @return the list, @c NULL if not enougn memory is available
arrayf arrayf_new(ulong size);

/// @brief Returns the list of coefficients of the real polynomial @c P.
///
/// @param P the real polynomial
///
/// @return the list, @c NULL if some error occurred.
arrayf arrayf_new_polyfr(polyfr P);

/// @brief Returns the list of coefficients of the polynomial @c P.
///
/// @param P the polynomial
///
/// @return the list, @c NULL if some error occurred.
arrayf arrayf_new_polyf(polyf P);

/// @brief Frees all the memory used by the list @c l, assuming the struct has been allocated with
/// @c malloc(), for example with @c list_new().
///
/// @param l the list
///
/// @return @ref true if successfull, @ref false otherwise.
bool arrayf_free(arrayf l);

/// @brief Adds the complex number @c z at the end of the list. The size of the list is automatically
/// increased, if needed.
///
/// @param l the list
/// @param z the new element to add to the list
///
/// @return @ref true if successfull, @ref false otherwise.
bool arrayf_add(arrayf l, compf z);

/// @brief Returns the pointer of the element on position @c pos in the list @c l.
///
/// @param l the list
/// @param pos the position of the element
///
/// @return the pointer to the element on position @c pos in @c l, @c NULL if some error occurred.
compf_ptr arrayf_get(arrayf l, ulong pos);

/// @brief Checks if all complex numbers in the list @c l are real.
///
/// @param l the list
///
/// @return @ref true if all elements in the list are real, @ref false otherwise.
bool arrayf_is_real(arrayf l);

/// @brief Checks if there are infinite values in the list @c l.
///
/// @param l the list
///
/// @return the position of the first infinite value in @c l, @c -1 if there is none, LONG_MAX if some error occurred.
long arrayf_first_inf(arrayf l);

/// @brief Checks if all complex numbers in the list @c l are well defined.
///
/// @param l the list
///
/// @return the position of the first undefined value in @c l, @c -1 if there is none, LONG_MAX if some error occurred.
long arrayf_first_nan(arrayf l);

/// @brief Returns the real polynomial with the coefficients the real parts of the complex numbers in the list @c l.
///
/// @param l the list
///
/// @return the polynomial, @c NULL if some coefficients are not real or some other error occurred.
polyfr arrayf_polyfr(arrayf l);

/// @brief Returns the polynomial with the coefficients in the list @c l.
///
/// @param l the list
///
/// @return the polynomial, @c NULL if some error occurred.
polyf arrayf_polyf(arrayf l);

// MARK: IO operations

/// @brief Writes the list @c l to a CSV file with the given path @c fileName.
///
/// @param l the list
/// @param fileName the path of the CSV file
/// @param verbose @ref true to print information about the error, if any
///
/// @return @ref true if successfull, @ref false otherwise.
bool arrayf_write(arrayf l, char *fileName, bool verbose);

/// @brief Writes the list @c l to the end of a CSV file with the given path @c fileName.
///
/// @param l the list
/// @param fileName the path of the CSV file
/// @param verbose @ref true to print information about the error, if any
///
/// @return @ref true if successfull, @ref false otherwise.
bool arrayf_append(arrayf l, char *fileName, bool verbose);

/// @brief Reads a list of complex numbers from the CSV file with path @c fileName.
///
/// @param fileName the path of the CSV file
/// @param verbose @ref true to print information about the file and the error, if any
///
/// @return the list of complex numbers in the CSV file, @c NULL if some error occurred.
arrayf arrayf_read(char *fileName, bool verbose);

#endif /* arrayf_h */
