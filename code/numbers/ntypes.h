//
//  ntypes.h
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
 \file ntypes.h
 \brief Definition of basic types.
*/

#ifndef ntypes_h
#define ntypes_h

#include <inttypes.h>
#include <stdint.h>
#include <math.h>

/// @c byte is @c uint8
typedef uint8_t byte;

/// @c word is @c uint16
typedef uint16_t word;

/// @c uint is @c uint32
typedef uint32_t uint;

/// @c ulong is @c uint64
typedef uint64_t ulong;

/// Logic type @c bool can take values \ref true or \ref false
typedef byte bool;

/// Boolean value true.
#define true   1

/// Boolean value false.
#define false  0

/// \f$ \pi \f$ with fp80 precision
#define PI (3.14159265358979323846264338327950288L)

// comment this line to use the fastest precise machine floating point numbers
// for polynomial evaluation (@c fp64 numbers)
// leave the line uncommented to use @c fp80 numbers on x86 machines, that occupy
// 16 bytes of memory, but only use 10 bytes for storing the numbers
// #define MACHINE_EXTRA_PREC

// Comment the previous line and uncomment this line to use 32-bit floats for coefficients
// #define MACHINE_LOW_PREC

/// The machine number type to use for polynomial coefficients and evaluation.
#ifdef MACHINE_EXTRA_PREC
  typedef long double coeff_t;
  #define flog2 log2l
  #define fhypot hypotl
  #define ffloor floorl
  #define fpow powl
  #define fexp expl
  #define fsin sinl
  #define fcos cosl
  #define ftan tanl
  #define mpfr_getf mpfr_get_ld
  #define precf 64
  #define PRECF_STR "64"
  #define TYPEFP_STR "FP80"
  #define TYPEF_STR "long double"
  #define INF_M (-HUGE_VALL)
  #define INF_P HUGE_VALL
  #define FMT_COEFF "L"
#else
#ifndef MACHINE_LOW_PREC
  typedef double coeff_t;
  #define flog2 log2
  #define fhypot hypot
  #define ffloor floor
  #define fpow pow
  #define fexp exp
  #define fsin sin
  #define fcos cos
  #define ftan tan
  #define mpfr_getf mpfr_get_d
  #define precf 53
  #define PRECF_STR "53"
  #define TYPEFP_STR "FP64"
  #define TYPEF_STR "double"
  #define INF_M (-HUGE_VAL)
  #define INF_P HUGE_VAL
  #define FMT_COEFF "l"
#else
  typedef float coeff_t;
  #define flog2 log2f
  #define fhypot hypotf
  #define ffloor floorf
  #define fpow powf
  #define fexp expf
  #define fsin sinf
  #define fcos cosf
  #define ftan tanf
  #define mpfr_getf mpfr_get_f
  #define precf 24
  #define PRECF_STR "24"
  #define TYPEFP_STR "FP32"
  #define TYPEF_STR "float"
  #define INF_M (-HUGE_VALF)
  #define INF_P HUGE_VALF
  #define FMT_COEFF ""
#endif
#endif

#define MAX_EXP (sizeof(coeff_t) == 8 ? 960 : 16320)


// Comment this line to limit the degree of the polynomials to 2^32-2
// leave this line uncommented to use degrees up to 2^64-2
#define HUGE_DEGREES

/// The integer number type to use for polynomial degrees and indexes.
#ifdef HUGE_DEGREES
  typedef ulong deg_t;
  #define MAX_DEG (UINT64_MAX - 1)
  #define FMT_DEG PRIu64
#else
  typedef uint deg_t;
  #define MAX_DEG (UINT32_MAX - 1)
  #define FMT_DEG PRIu32
#endif


// Comment this line to use the fastest precise machine floating point numbers
// for polynomial pre-condtionning; largely enough in 99+% usage cases
// leave the line uncommented to use @c fp80 numbers on x86 machines, that occupy
// 16 bytes of memory, but only use 10 bytes for storing the numbers
//#define INTERNAL_EXTRA_PREC

/// The machine number type to use for polynomial analysis and preconditionning.
#ifdef INTERNAL_EXTRA_PREC
  typedef long double real_t;
  #define plog2 log2l
  #define pexp2 exp2l
  #define phypot hypotl
  #define pfloor floorl
  #define pceil ceill
  #define mpfr_getp mpfr_get_ld
  #define PEPS 1E-3000L
  #define PDEL 1E-15L
  #define pexp expl
  #define pcos cosl
  #define psin sinl
  #define ptan tanl
  #define pldexp ldexpl
  #define pfrexp frexpl
  #define FMT_REAL "L"
#else
  typedef double real_t;
  #define plog2 log2
  #define pexp2 exp2
  #define phypot hypot
  #define pfloor floor
  #define pceil ceil
  #define mpfr_getp mpfr_get_d
  #define PEPS 1E-200
  #define PDEL 1E-12
  #define pexp exp
  #define pcos cos
  #define psin sin
  #define ptan tan
  #define pldexp ldexp
  #define pfrexp frexp
  #define FMT_REAL "l"
#endif


// Comment this line to limit the precision of mpfr numbers to to 4 * 10^9
// leave this line uncommented to use up to 16 * 10^18 bits for each number
#define HUGE_PREC

/// The integer number type to use for polynomial degrees and indexes.
#ifdef HUGE_PREC
  typedef ulong prec_t;
  #define MAX_PREC (16000000000000000000UL)
#else
  typedef uint prec_t;
  #define MAX_PREC (4000000000UL)
#endif

// Comment this line to limit the exponents of mpfr numbers to to +/- 2 * 10^9
// leave this line uncommented to use exponents in the range +/- 4 * 10^18
#define HUGE_MP_EXP

/// The integer number type to use for polynomial degrees and indexes.
#ifdef HUGE_MP_EXP
  #define MAX_MP_EXP (4000000000000000000L)
#else
  #define MAX_MP_EXP (2000000000L)
#endif

#define NEWTON_CONV_BITS     5
#define NEWTON_ESCAPE_BITS   4

/// @brief Use this function to check the compatibility of the machine with the settings above.
bool ntypes_check(void);

/// @brief Computes @c log_2(2^b1+2^b2), even if @c b1 or @c b2 are outside the exponent range of @ref real_t.
///
/// @param b1 first operand, in bits
/// @param b2 second operand, in bits
///
/// @return the bits sum of the operands, that is the base @c 2  @c log of their base @c 2 exponentials
real_t bits_sum(real_t b1, real_t b2);

/// @brief Computes @c log_2(ntErr), where @c ntErr is an upper bound for the error of the
/// Newton term.
///
/// @param vb @c log_2(|value|)
/// @param evb @c log_2(valErr)
/// @param db @c log_2(|derivative|)
/// @param edb @c log_2(derErr)
///
/// @return absolute error of the Newton term, as a power of @c 2
real_t nt_err(real_t vb, real_t evb, real_t db, real_t edb);

#endif /* ntypes_h */
