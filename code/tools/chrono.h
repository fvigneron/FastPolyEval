//
//  chrono.h
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
 \file chrono.h
 \brief Tools for precise time measuring.
*/


#ifndef chrono_h
#define chrono_h

#include <stdio.h>

#include <sys/timeb.h>
#include <time.h>

#include "ntypes.h"

/// One billion (nanoseconds in one second).
#define BILLION    1000000000L

/// One billion (miliseconds in one second).
#define MILLION       1000000L

/// Miliseconds in one minute.
#define MINUTE          60000L

/// Miliseconds in one hour.
#define HOUR              (60 * MINUTE)

/// Miliseconds in one day.
#define DAY               (24 * HOUR)

/// The ID of the system clock to use.
#define USED_CLOCK CLOCK_MONOTONIC

/// Short type name for precise time.
typedef struct timespec ptime;

/// @brief Formats the time duration in miliseconds @c ms into human readable string @c str.
///
/// @warning The string @c str should be at least @c 90 bytes long.
///
/// @param ms the time lapse in miliseconds
/// @param str the human readable print
void millis(long ms, char *str);

/// @brief Formats the time duration in nanoseconds @c ns into human readable string @c str.
///
/// @warning The string @c str should be at least @c 90 bytes long.
///
/// @param ns the time lapse in miliseconds
/// @param str the human readable print
void nanos(long ns, char *str);

/// @brief Computes and prints the time lapse since @c ts.
///
/// Pretty prints the elapsed time into the string @c str and updates @c ts to the current time.
/// Does not fail if @c str==NULL, can also be used as a starter for the chronometer.
///
/// @warning The string @c str should be at least @c 90 bytes long.
///
/// @param ts old time stamp
/// @param str human readable string for the time increment
///
/// @return the total time in nanos
ulong lap(ptime *ts, char *str);

/// @brief Prints the timer resolution.
void print_time_res(void);

#endif /* chrono_h */
