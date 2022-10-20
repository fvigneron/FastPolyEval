//
//  apps.h
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
 \file apps.h
 \brief The implementation of the mini-apps with MPFR numbers.
*/

#ifndef apps_h
#define apps_h

#include <stdio.h>

#include "array.h"

/// @brief Joins the real parts of the sequences @c z1 and @c z2 as the real parts and imaginary
/// parts of a new list.
///
/// @param prec the precision of the result
/// @param z1 the first list
/// @param z2 the second list
///
/// @return the joint list, @c NULL if some error occurred.
array app_join(long prec, array z1, array z2);

/// @brief Computes the set product of the real parts of two sequences @c z1 and @c z2.
///
/// @param prec the precision of the result
/// @param z1 the first list
/// @param z2 the second list
///
/// @return the set product, @c NULL if some error occurred.
array app_grid(long prec, array z1, array z2);

/// @brief Computes the complex exponential of a sequence @c z.
///
/// @param prec the precision of the result
/// @param z the list of points
///
/// @return the image by the complex exponential, @c NULL if some error occurred.
array app_exp(long prec, array z);

/// @brief Computes @c a*exp(ib) for each term @c a+ib of a sequence @c z of complex numbers.
///
/// @param prec the precision of the result
/// @param z the list of points
///
/// @return the image by @c a+ib->a*exp(ib), @c NULL if some error occurred.
array app_rot(long prec, array z);

/// @brief Computes the complex numbers given by a sequence @c z of their polar coordinates on the sphere.
///
/// @param prec the precision of the result
/// @param z the list of polar coordinates
///
/// @return the complex numbers given by their polar coordinates, @c NULL if some error occurred.
array app_polar(long prec, array z);

/// @brief Returns the list of @c n real numbers that split the interval @c [st,en] in @c n-1 equal intervals.
///
/// @param prec the precision of the result
/// @param n the number of points, at least @c 2
/// @param st the start point
/// @param en the end point
///
/// @return the list of real numbers (in complex format), @c NULL if some error occurred.
array app_unif(long prec, ulong n, mpfr_t st, mpfr_t en);

/// @brief Returns the list of @c n random real numbers in the interval @c [st,en], uniformly distributed.
///
/// @param prec the precision of the result
/// @param n the number of points, at least @c 2
/// @param st the start point
/// @param en the end point
///
/// @return the list of random real numbers (in complex format), @c NULL if some error occurred.
array app_rand(long prec, ulong n, mpfr_t st, mpfr_t en);

/// @brief Returns the list of @c n random real numbers, normally distributed with center @c cent and variance @c var.
///
/// @param prec the precision of the result
/// @param n the number of points, at least @c 2
/// @param cent the center of the normal distribution
/// @param var the variance of the normal distribution
///
/// @return the list of random real numbers (in complex format), @c NULL if some error occurred.
array app_normal(long prec, ulong n, mpfr_t cent, mpfr_t var);

/// @brief Returns a list of @c n compelx numbers, uniformly distributed on the Riemann sphere.
///
/// @param prec the precision of the result
/// @param n the number of points, at least @c 2
///
/// @return the list of points uniformly distributed on the sphere, @c NULL if some error occurred.
array app_sphere(long prec, ulong n);

/// @brief Compares two list of complex values and reports the eventual errors.
///
/// @param prec the precision of the result
/// @param z1 the fist list
/// @param z2 the second list
/// @param output a file name for the output file, @c NULL if not used
///
/// @return @ref true if the comparison is complete, @ref false if some error occurred.
bool  app_compare(long prec, array z1, array z2, char *output);

/// @brief Returns the real parts of the numbers in the sequence @c z.
///
/// @param prec the precision of the result
/// @param z the list of points
///
/// @return the real parts in @c z, @c NULL if some error occurred.
array app_re(long prec, array z);

/// @brief Returns the imaginary parts of the numbers in the sequence @c z.
///
/// @param prec the precision of the result
/// @param z the list of points
///
/// @return the imaginary parts in @c z, @c NULL if some error occurred.
array app_im(long prec, array z);

/// @brief Returns the conjugates of the numbers in the sequence @c z.
///
/// @param prec the precision of the result
/// @param z the list of points
///
/// @return the conjugates of @c z, @c NULL if some error occurred.
array app_conj(long prec, array z);

/// @brief Computes tensor (or piontwise) product of the sequences @c z1 and @c z2.
///
/// @param prec the precision of the result
/// @param z1 the first list
/// @param z2 the second list
///
/// @return the tensor product, @c NULL if some error occurred.
array app_tensor(long prec, array z1, array z2);

/// @brief Evaluates the polynomial @c P in all points in the list @c z, using the @b FPE algorithm.
///
/// @param prec the precision of intermediary computations and of the result
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
bool app_eval_p(long prec, poly P, array z, char *outFile, char *outError,
                 long count, char *outHorner, long countHorner, const char **inFiles);

/// @brief Evaluates the derivative of the polynomial @c P in all points in the list @c z, using the @b FPE algorithm.
///
/// @param prec the precision of intermediary computations and of the result
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
bool app_eval_d(long prec, poly P, array z, char *outFile, char *outError,
                 long count, char *outHorner, long countHorner, const char **inFiles);

/// @brief Evaluates the Newton terms of the polynomial @c P in all points in the list @c z, using the @b FPE algorithm.
///
/// @param prec the precision of intermediary computations and of the result
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
bool app_eval_n(long prec, poly P, array z, char *outFile, char *outError,
                 long count, char *outHorner, long countHorner, const char **inFiles);

/// @brief Evaluates the Newton iterates of the polynomial @c P in all points in the list @c z,
/// using the @b FPE algorithm.
///
/// This is a method to approximate the roots of @c P, but there are no guarantees for the convergence. If the iterates
/// escape far from the origin, the alogrith stops (for those starting points).
///
/// @param prec the precision of intermediary computations and of the result
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
bool app_eval_r(long prec, poly P, array z, ulong maxIter, char *outFile, char *outError,
                 long count, char *outHorner, long countHorner, const char **inFiles);

/// @brief Analyses the polynomail @c P and outputs its concave cover, ignored coefficients and (optionally) the
/// modules of complex numbers for which the reduced polynomial changes.
///
///These analysis is performed for the evaluations of the @c P with @c prec bits of precision. The concave cover
///corresponds to the pre-processing phase of the @b FPE algortihm. The map @c k->log_2|a_k| is considered.
///
/// @param prec the precision of evaluations
/// @param P the polynomial
/// @param outFile the output file containing the concave cover and the ognored coefficients
/// @param outChange the output file containing the intervals for @c log_2|z| and the corresponding reduced polynomials of @c P
///
/// @return @ref true if the analysis is complete, @ref false if some error occurred.
bool app_analyse(long prec, poly P, char *outFile, char *outChange);

#endif /* apps_h */
