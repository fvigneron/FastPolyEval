//
//  help.h
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
 \file help.h
 \brief A basic help system for the mini-apps.
*/

#ifndef help_h
#define help_h

#include <stdio.h>

#include "ntypes.h"

// MARK: constants & data types

/// Indexes for help pages.
#define HELP_MAIN       0
#define HELP_HYP        1
#define HELP_CHEBYshev  2
#define HELP_LEGENDRE   3
#define HELP_HERMITE    4
#define HELP_LAGUERRE   5
#define HELP_SUM        6
#define HELP_DIFF       7
#define HELP_PROD       8
#define HELP_CONCAT     9
#define HELP_JOIN      10
#define HELP_GRID      11
#define HELP_EXP       12
#define HELP_ROT       13
#define HELP_POLAR     14
#define HELP_ROOTS     15
#define HELP_DER       16
#define HELP_UNIF      17
#define HELP_RAND      18
#define HELP_NORM      19
#define HELP_SPHERE    20
#define HELP_EVAL      21
#define HELP_EVALD     22
#define HELP_EVALN     23
#define HELP_NEWTON    24
#define HELP_COMP      25
#define HELP_RE        26
#define HELP_IM        27
#define HELP_CONJ      28
#define HELP_TENSOR    29
#define HELP_ANALYSE   30

/// The number of help pages.
#define HELP_COUNT     31

/// The number of mini-apps.
#define APP_COUNT      (HELP_COUNT - 1)

/// @struct help_struct
/// @brief Description of a basic help screen.
///
/// Each help screen contains a brief summary (@c before), some @c lines of messages (commands) and a final discussion (@c after).
/// Lines are structures as an array of variable size, with descriptive headers.
typedef struct {
    char *before;        ///< the first line to display
    int columnCount;     ///< the number of columns
    int *columnWidths;   ///< the widths of all but the last column (array length should be 1 less than columnCount)
    char **headers;      ///< the column headers
    int linesCount;      ///< the number of lines / commands
    char ***lines;       ///< the double array of messages
    char *after;         ///< the last line to display
} help_struct;

/// @brief Practical wrapper for @c arrayf_struct
///
/// To avoid the constant use @c * and @c & the type  @c help_t is a pointer.
typedef help_struct *help;

// MARK: functions prototypes

/// @brief Prints the help screen described by @c h.
///
/// @param h the help screen to display
///
/// @return Boolean value expressing wether help could indeed be provided.
bool fpe_help_print(help h);

/// @brief Return the help system in the list with index @c ind.
///
/// @param ind the index, use constants HELP_XXX
///
/// @return the help sytem with the given index.
/// Returns the list, @c NULL if @c ind is an invalid request.
help fpe_help_get(int ind);

/// @brief Prints computation stats in a standard format for automatic processing.
///
/// @param mes a short description of the computing task
/// @param prec the precision used, in bits
/// @param deg the degree of the polynomial
/// @param numType the number type
/// @param real @ref true if it is a real polynomial, @ref false otherwise
/// @param polyFile the polynomial file name
/// @param ptsFile the points list file name
/// @param outFile the results output file name
/// @param outError the errors list file name
/// @param timePP pre-processing time in ns
/// @param count the number of repetitions of the computing task
/// @param timeEV evaluation time in ns
/// @param tpts the total number of evaluation points (taking into acount the repetitions)
/// @param outHorner the results output file name, computed by Horner's method
/// @param timeHo the computing time in ns, by Horner's method
/// @param tptsho the total number of evaluation by Horner's method
/// @param tit total number of iterates, only for the iterative Newton method
/// @param titho total number of iterates, only for the iterative Newton method, computed by Horner's method
void stats_print(char *mes, long prec, deg_t deg, char *numType, bool real,
                 char *polyFile, char *ptsFile, char *outFile, char *outError,
                 double timePP, long count, double timeEV, ulong tpts,
                 char *outHorner, double timeHo, ulong tptsho, ulong tit, ulong titho);

#endif /* help_h */
