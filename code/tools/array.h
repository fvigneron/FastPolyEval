//
//  array.h
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
 \file array.h
 \brief Definition of a variable length array of arbitary precision complex numbers that are based on [mpfr] (https://www.mpfr.org).
*/

#ifndef array_h
#define array_h

#include <stdio.h>
#include <stdlib.h>
#include <mpfr.h>

#include "ntypes.h"
#include "comp.h"
#include "poly.h"
#include "polyr.h"

#define ARRAY_MIN_SIZE       100
#define ARRAY_SIZE_INCREASE  1.25

/// @struct array_struct
/// @brief A variable length array of machine floating point complex numbers.
typedef struct {
    ulong len;      ///< the number of complex numbers in the array
    ulong size;     ///< the present capacity of the array (memory size in number of complex numbers)
    comp_ptr zi;    ///< the elements of the array
} array_struct;

/// @brief Practical wrapper for @c array_struct
///
/// To avoid the constant use @c * and @c & the type  @c array_t is a pointer.
typedef array_struct array_t[1];

/// Convenience pointer to @c array_struct
typedef array_struct *array;

// MARK: functions prototypes

/// @brief Returns a new empty array with memory size at least @c size.
///
/// @param size the number of elements to allocate memory for
///
/// @return the array, @c NULL if not enougn memory is available
array array_new(ulong size);

/// @brief Returns the array of coefficients of the real polynomial @c P.
///
/// @param P the real polynomial
/// @param prec the precision of the numbers in the array
///
/// @return the array, @c NULL if some error occurred.
array array_new_polyr(polyr P, prec_t prec);

/// @brief Returns the array of coefficients of the polynomial @c P.
///
/// @param P the polynomial
/// @param prec the precision of the numbers in the array
///
/// @return the array, @c NULL if some error occurred.
array array_new_poly(poly P, prec_t prec);

/// @brief Frees all the memory used by the array @c l, assuming the struct has been allocated with
/// @c malloc(), for example with @c array_new().
///
/// @param l the array
///
/// @return @ref true if successfull, @ref false otherwise.
bool array_free(array l);

/// @brief Adds the complex number @c z at the end of the array. The size of the array is automatically
/// increased, if needed.
///
/// @param l the array
/// @param z the new element to add to the array
/// @param prec the precision of the number in the array
///
/// @return @ref true if successfull, @ref false otherwise.
bool array_add(array l, comp z, prec_t prec);

/// @brief Returns the pointer of the element on position @c pos in the array @c l.
///
/// @param l the array
/// @param pos the position of the element
///
/// @return the pointer to the element on position @c pos in @c l, @c NULL if some error occurred.
comp_ptr array_get(array l, ulong pos);

/// @brief Checks if all complex numbers in the array @c l are real.
///
/// @param l the array
///
/// @return @ref true if all elements in the array are real, @ref false otherwise.
bool array_is_real(array l);

/// @brief Checks if there are infinite values in the list @c l.
///
/// @param l the list
///
/// @return the position of the first infinite value in @c l, @c -1 if there is none, LONG_MAX if some error occurred.
long array_first_inf(array l);

/// @brief Checks if all complex numbers in the list @c l are well defined.
///
/// @param l the list
///
/// @return the position of the first undefined value in @c l, @c -1 if there is none, LONG_MAX if some error occurred.
long array_first_nan(array l);

/// @brief Returns the real polynomial with the coefficients the real parts of the complex numbers in the array @c l.
///
/// @param l the array
/// @param prec the precision of the coefiicients of the polynomial
///
/// @return the polynomial, @c NULL if some coefficients are not real or some other error occurred.
polyr array_polyr(array l, prec_t prec);

/// @brief Returns the polynomial with the coefficients in the array @c l.
///
/// @param l the array
/// @param prec the precision of the coefiicients of the polynomial
///
/// @return the polynomial, @c NULL if some error occurred.
poly array_poly(array l, prec_t prec);

// MARK: IO operations

/// @brief Writes the array @c l to a CSV file with the given path @c fileName.
///
/// @param l the array
/// @param fileName the path of the CSV file
/// @param digits the number of digits after the decimal point of the numbers in the file
/// @param verbose @ref true to print information about the error, if any
///
/// @return @ref true if successfull, @ref false otherwise.
bool array_write(array l, char *fileName, int digits, bool verbose);

/// @brief Writes the array @c l to the end of a CSV file with the given path @c fileName.
///
/// @param l the array
/// @param fileName the path of the CSV file
/// @param digits the number of digits after the decimal point of the numbers in the file
/// @param verbose @ref true to print information about the error, if any
///
/// @return @ref true if successfull, @ref false otherwise.
bool array_append(array l, char *fileName, int digits, bool verbose);

/// @brief Reads a array of complex numbers of precision @c prec from the CSV file with path @c fileName.
///
/// @param fileName the path of the CSV file
/// @param prec the precision
/// @param verbose @ref true to print information about the file and the error, if any
///
/// @return the array of complex numbers in the CSV file, @c NULL if some error occurred.
array array_read(char *fileName, prec_t prec, bool verbose);

#endif /* array_h */
