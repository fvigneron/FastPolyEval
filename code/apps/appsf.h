//
//  appsf.h
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
 \file appsf.h
 \brief The implementation of the mini-apps with machine numbers.
*/

#ifndef appsf_h
#define appsf_h

#include <stdio.h>
#include "arrayf.h"

/// @brief Joins the real parts of the sequences @c z1 and @c z2 as the real parts and imaginary
/// parts of a new list.
///
/// @param z1 the first list
/// @param z2 the second list
///
/// @return the joint list, @c NULL if some error occurred.
arrayf appf_join(arrayf z1, arrayf z2);

/// @brief Computes the set product of the real parts of two sequences @c z1 and @c z2.
///
/// @param z1 the first list
/// @param z2 the second list
///
/// @return the set product, @c NULL if some error occurred.
arrayf appf_grid(arrayf z1, arrayf z2);

/// @brief Computes the complex exponential of a sequence @c z.
///
/// @param z the list of points
///
/// @return the image by the complex exponential, @c NULL if some error occurred.
arrayf appf_exp(arrayf z);

/// @brief Computes @c a*exp(ib) for each term @c a+ib of a sequence @c z of complex numbers.
///
/// @param z the list of points
///
/// @return the image by @c a+ib->a*exp(ib), @c NULL if some error occurred.
arrayf appf_rot(arrayf z);

/// @brief Computes the complex numbers given by a sequence @c z of their polar coordinates on the sphere.
///
/// @param z the list of polar coordinates
///
/// @return the complex numbers given by their polar coordinates, @c NULL if some error occurred.
arrayf appf_polar(arrayf z);

/// @brief Returns the list of @c n real numbers that split the interval @c [st,en] in @c n-1 equal intervals.
///
/// @param n the number of points, at least @c 2
/// @param st the start point
/// @param en the end point
///
/// @return the list of real numbers (in complex format), @c NULL if some error occurred.
arrayf appf_unif(ulong n, coeff_t st, coeff_t en);

/// @brief Compares two list of complex values and reports the eventual errors.
///
/// @param z1 the fist list
/// @param z2 the second list
/// @param output a file name for the output file, @c NULL if not used
///
/// @return @ref true if the comparison is complete, @ref false if some error occurred.
bool appf_compare(arrayf z1, arrayf z2, char *output);

/// @brief Returns the real parts of the numbers in the sequence @c z.
///
/// @param z the list of points
///
/// @return the real parts in @c z, @c NULL if some error occurred.
arrayf appf_re(arrayf z);

/// @brief Returns the imaginary parts of the numbers in the sequence @c z.
///
/// @param z the list of points
///
/// @return the imaginary parts in @c z, @c NULL if some error occurred.
arrayf appf_im(arrayf z);

/// @brief Returns the conjugates of the numbers in the sequence @c z.
///
/// @param z the list of points
///
/// @return the conjugates of @c z, @c NULL if some error occurred.
arrayf appf_conj(arrayf z);

/// @brief Computes tensor (or piontwise) product of the sequences @c z1 and @c z2.
///
/// @param z1 the first list
/// @param z2 the second list
///
/// @return the tensor product, @c NULL if some error occurred.
arrayf appf_tensor(arrayf z1, arrayf z2);

/// @brief Evaluates the polynomial @c P in all points in the list @c z, using the @b FPE algorithm.
///
/// @param P the polynomial
/// @param z the list of evaluation points
/// @param outFile the results output file name
/// @param outError the errors list file name
/// @param count the number of repetitions of the computing task
/// @param outHorner the results output file name, computed by Horner's method
/// @param countHorner the number of repetitions of the computing task by Horner's method
/// @param inFiles a vector with the names of input files, for printing stats only
///
/// @return @ref true if the evaluaion is complete, @ref false if some error occurred.
bool appf_eval_p(polyf P, arrayf z, char *outFile, char *outError,
                 long count, char *outHorner, long countHorner, const char **inFiles);

/// @brief Evaluates the derivative of the polynomial @c P in all points in the list @c z, using the @b FPE algorithm.
///
/// @param P the polynomial
/// @param z the list of evaluation points
/// @param outFile the results output file name
/// @param outError the errors list file name
/// @param count the number of repetitions of the computing task
/// @param outHorner the results output file name, computed by Horner's method
/// @param countHorner the number of repetitions of the computing task by Horner's method
/// @param inFiles a vector with the names of input files, for printing stats only
///
/// @return @ref true if the evaluaion is complete, @ref false if some error occurred.
bool appf_eval_d(polyf P, arrayf z, char *outFile, char *outError,
                 long count, char *outHorner, long countHorner, const char **inFiles);

/// @brief Evaluates the Newton terms of the polynomial @c P in all points in the list @c z, using the @b FPE algorithm.
///
/// @param P the polynomial
/// @param z the list of evaluation points
/// @param outFile the results output file name
/// @param outError the errors list file name
/// @param count the number of repetitions of the computing task
/// @param outHorner the results output file name, computed by Horner's method
/// @param countHorner the number of repetitions of the computing task by Horner's method
/// @param inFiles a vector with the names of input files, for printing stats only
///
/// @return @ref true if the evaluaion is complete, @ref false if some error occurred.
bool appf_eval_n(polyf P, arrayf z, char *outFile, char *outError,
                 long count, char *outHorner, long countHorner, const char **inFiles);

/// @brief Evaluates the Newton iterates of the polynomial @c P in all points in the list @c z,
/// using the @b FPE algorithm.
///
/// This is a method to approximate the roots of @c P, but there are no guarantees for the convergence. If the iterates
/// escape far from the origin, the alogrith stops (for those starting points).
///
/// @param P the polynomial
/// @param z the list of evaluation points
/// @param maxIter the maximum iterates for each starting point
/// @param outFile the results output file name
/// @param outError the errors list file name
/// @param count the number of repetitions of the computing task
/// @param outHorner the results output file name, computed by Horner's method
/// @param countHorner the number of repetitions of the computing task by Horner's method
/// @param inFiles a vector with the names of input files, for printing stats only
///
/// @return @ref true if the evaluaion is complete, @ref false if some error occurred.
bool appf_eval_r(polyf P, arrayf z, ulong maxIter, char *outFile, char *outError,
                 long count, char *outHorner, long countHorner, const char **inFiles);

#endif /* appsf_h */
