//
//  main.c
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

#include <stdio.h>
#include <string.h>

#include "main.h"
#include "help.h"
#include "ntypes.h"
#include "polyfr.h"
#include "arrayf.h"
#include "array.h"
#include "evalf.h"
#include "eval.h"
#include "chrono.h"
#include "debug.h"

#include "apps.h"
#include "appsf.h"

// MARK: apps with command line parameters

static bool write_polyr(polyr P, long prec, const char *fileName, char *comm[], int commCount) {
    array l = array_new_polyr(P, prec);
    int dig = (int) (prec / 3);
    
    if(l == NULL || array_first_inf(l) >= 0 || array_first_nan(l) >= 0) {
        array_write(l, (char *) fileName, dig, true); // only for the error messages
        array_free(l);
        
        return false;
    }
    
    FILE *f = fopen(fileName, "w");
    if(f == NULL) {
        printf("\n**** Could not open \"%s\" for writing ! ****\n\n", fileName);
        array_free(l);
        
        return false;
    }
    
    for (int i = 0; i < commCount; i++) {
        fprintf(f, "# %s\n", comm[i]);
    }
    
    fclose(f);
    
    bool ok = array_append(l, (char *) fileName, dig, true);
    array_free(l);
    
    return ok;
}

static bool write_polyfr(polyfr P, const char *fileName, char *comm[], int commCount) {
    arrayf l = arrayf_new_polyfr(P);
    
    if(l == NULL || arrayf_first_inf(l) >= 0 || arrayf_first_nan(l) >= 0) {
        arrayf_write(l, (char *) fileName, true); // only for the error messages
        arrayf_free(l);
        
        return false;
    }
    
    FILE *f = fopen(fileName, "w");
    if(f == NULL) {
        printf("\n**** Could not open \"%s\" for writing ! ****\n\n", fileName);
        arrayf_free(l);
        
        return false;
    }
    
    for (int i = 0; i < commCount; i++) {
        fprintf(f, "# %s\n", comm[i]);
    }
    
    fclose(f);
    
    bool ok = arrayf_append(l, (char *) fileName, true);
    arrayf_free(l);
    
    return ok;
}

static int hyp(long prec, int argc, const char * argv[]) {
    int per;
    
    if(argc < 3 || sscanf(argv[1], "%d", &per) < 1 || per < 1) {
        return APP_PARAMS;
    }
        
    polyr P = poly_hyp(per, prec);
    
    if(P == NULL) {
        printf("Could not generate the hyperbolic polynomial of order %d with precision %ld ! \n",
               per, prec);
        
        return APP_ERROR;
    }
    
    char comm[200], *pcomm[] = {comm};
    snprintf(comm, 200, "The hyperbolic polynomial of order %d with precision %ld (MPFR)",
             per, prec);
    
    if(! write_polyr(P, prec, argv[2], pcomm, 1)) {
        printf("Could not write the hyperbolic polynomial to %s ! \n", argv[2]);
        polyr_free(P);
        
        return APP_ERROR;
    }
    
    printf("The hyperbolic polynomial of order %d is written to %s.\n", per, argv[2]);
    polyr_free(P);
    
    return APP_OK;
}

static int hypf(int argc, const char * argv[]) {
    int per;
    if(argc < 3 || sscanf(argv[1], "%d", &per) < 1 || per < 1) {
        return APP_PARAMS;
    }
    
    polyfr P = polyf_hyp(per);
    
    if(P == NULL) {
        printf("Could not generate the hyperbolic polynomial of order %d with machine"
               " floating point numbers ! \n", per);
        
        return APP_ERROR;
    }
    
    char comm[200], *pcomm[] = {comm};
    snprintf(comm, 200, "The hyperbolic polynomial of order %d with precision "
             PRECF_STR " (" TYPEFP_STR ")", per);
    
    if(! write_polyfr(P, argv[2], pcomm, 1)) {
        printf("Could not write the hyperbolic polynomial to %s ! \n", argv[2]);
        polyfr_free(P);
        
        return APP_ERROR;
    }
    
    printf("The hyperbolic polynomial of order %d is written to %s.\n", per, argv[2]);
    polyfr_free(P);
    
    return APP_OK;
}

static int cheb(long prec, int argc, const char * argv[]) {
    int d;
    if(argc < 3 || sscanf(argv[1], "%d", &d) < 1 || d < 0) {
        return APP_PARAMS;
    }
    
    polyr P = poly_cheb(d, prec);
    
    if(P == NULL) {
        printf("Could not generate the Chebyshev polynomial of degree %d with precision %ld ! \n",
               d, prec);
        
        return APP_ERROR;
    }
    
    char comm[200], *pcomm[] = {comm};
    snprintf(comm, 200, "The Chebyshev polynomial of degree %d with precision %ld (MPFR)",
             d, prec);
    
    if(! write_polyr(P, prec, argv[2], pcomm, 1)) {
        printf("Could not write the Chebyshev polynomial to %s ! \n", argv[2]);
        polyr_free(P);
        
        return APP_ERROR;
    }
    
    printf("The Chebyshev polynomial of degree %d is written to %s.\n", d, argv[2]);
    polyr_free(P);
    
    return APP_OK;
}

static int chebf(int argc, const char * argv[]) {
    int d;
    if(argc < 3 || sscanf(argv[1], "%d", &d) < 1 || d < 0) {
        return APP_PARAMS;
    }
    
    polyfr P = polyf_cheb(d);
    
    if(P == NULL) {
        printf("Could not generate the Chebyshev polynomial of degree %d with machine"
               " floating point numbers ! \n", d);
        
        return APP_ERROR;
    }
    
    char comm[200], *pcomm[] = {comm};
    snprintf(comm, 200, "The Chebyshev polynomial of degree %d with precision "
             PRECF_STR " (" TYPEFP_STR ")", d);
    
    if(! write_polyfr(P, argv[2], pcomm, 1)) {
        printf("Could not write the Chebyshev polynomial to %s ! \n", argv[2]);
        polyfr_free(P);
        
        return APP_ERROR;
    }
    
    printf("The Chebyshev polynomial of degree %d is written to %s.\n", d, argv[2]);
    polyfr_free(P);
    
    return APP_OK;
}

static int leg(long prec, int argc, const char * argv[]) {
    int d;
    if(argc < 3 || sscanf(argv[1], "%d", &d) < 1 || d < 0) {
        return APP_PARAMS;
    }
    
    polyr P = poly_leg(d, prec);
    
    if(P == NULL) {
        printf("Could not generate the Legendre polynomial of degree %d with precision %ld ! \n",
               d, prec);
        
        return APP_ERROR;
    }
    
    char comm[200], *pcomm[] = {comm};
    snprintf(comm, 200, "The Legendre polynomial of degree %d with precision %ld (MPFR)",
             d, prec);
    
    if(! write_polyr(P, prec, argv[2], pcomm, 1)) {
        printf("Could not write the Legendre polynomial to %s ! \n", argv[2]);
        polyr_free(P);
        
        return APP_ERROR;
    }
    
    printf("The Legendre polynomial of degree %d is written to %s.\n", d, argv[2]);
    polyr_free(P);
    
    return APP_OK;
}

static int legf(int argc, const char * argv[]) {
    int d;
    if(argc < 3 || sscanf(argv[1], "%d", &d) < 1 || d < 0) {
        return APP_PARAMS;
    }
    
    polyfr P = polyf_leg(d);
    
    if(P == NULL) {
        printf("Could not generate the Legendre polynomial of degree %d with machine"
               " floating point numbers ! \n", d);
        
        return APP_ERROR;
    }
    
    char comm[200], *pcomm[] = {comm};
    snprintf(comm, 200, "The Legendre polynomial of degree %d with precision "
             PRECF_STR " (" TYPEFP_STR ")", d);
    
    if(! write_polyfr(P, argv[2], pcomm, 1)) {
        printf("Could not write the Legendre polynomial to %s ! \n", argv[2]);
        polyfr_free(P);
        
        return APP_ERROR;
    }
    
    printf("The Legendre polynomial of degree %d is written to %s.\n", d, argv[2]);
    polyfr_free(P);
    
    return APP_OK;
}

static int her(long prec, int argc, const char * argv[]) {
    int d;
    if(argc < 3 || sscanf(argv[1], "%d", &d) < 1 || d < 0) {
        return APP_PARAMS;
    }
        
    polyr P = poly_her(d, prec);
    
    if(P == NULL) {
        printf("Could not generate the Hermite polynomial of degree %d with precision %ld ! \n",
               d, prec);
        
        return APP_ERROR;
    }
    
    char comm[200], *pcomm[] = {comm};
    snprintf(comm, 200, "The Hermite polynomial of degree %d with precision %ld (MPFR)",
             d, prec);
    
    if(! write_polyr(P, prec, argv[2], pcomm, 1)) {
        printf("Could not write the Hermite polynomial to %s ! \n", argv[2]);
        polyr_free(P);
        
        return APP_ERROR;
    }
    
    printf("The Hermite polynomial of degree %d is written to %s.\n", d, argv[2]);
    polyr_free(P);
    
    return APP_OK;
}

static int herf(int argc, const char * argv[]) {
    int d;
    if(argc < 3 || sscanf(argv[1], "%d", &d) < 1 || d < 0) {
        return APP_PARAMS;
    }
    
    polyfr P = polyf_her(d);
    
    if(P == NULL) {
        printf("Could not generate the Hermite polynomial of degree %d with machine"
               " floating point numbers ! \n", d);
        
        return APP_ERROR;
    }
    
    char comm[200], *pcomm[] = {comm};
    snprintf(comm, 200, "The Hermite polynomial of degree %d with precision "
             PRECF_STR " (" TYPEFP_STR ")", d);
    
    if(! write_polyfr(P, argv[2], pcomm, 1)) {
        printf("Could not write the Hermite polynomial to %s ! \n", argv[2]);
        polyfr_free(P);
        
        return APP_ERROR;
    }
    
    printf("The Hermite polynomial of degree %d is written to %s.\n", d, argv[2]);
    polyfr_free(P);
    
    return APP_OK;
}

static int lag(long prec, int argc, const char * argv[]) {
    int d;
    if(argc < 3 || sscanf(argv[1], "%d", &d) < 1 || d < 0) {
        return APP_PARAMS;
    }
        
    polyr P = poly_lag(d, prec);
    
    if(P == NULL) {
        printf("Could not generate the Laguerre polynomial of degree %d with precision %ld ! \n",
               d, prec);
        
        return APP_ERROR;
    }
    
    char comm[200], *pcomm[] = {comm};
    snprintf(comm, 200, "The Laguerre polynomial of degree %d with precision %ld (MPFR)",
             d, prec);
    
    if(! write_polyr(P, prec, argv[2], pcomm, 1)) {
        printf("Could not write the Laguerre polynomial to %s ! \n", argv[2]);
        polyr_free(P);
        
        return APP_ERROR;
    }
    
    printf("The Laguerre polynomial of degree %d is written to %s.\n", d, argv[2]);
    polyr_free(P);
    
    return APP_OK;
}

static int lagf(int argc, const char * argv[]) {
    int d;
    if(argc < 3 || sscanf(argv[1], "%d", &d) < 1 || d < 0) {
        return APP_PARAMS;
    }
    
    polyfr P = polyf_lag(d);
    
    if(P == NULL) {
        printf("Could not generate the Laguerre polynomial of degree %d with machine"
               " floating point numbers ! \n", d);
        
        return APP_ERROR;
    }
    
    char comm[200], *pcomm[] = {comm};
    snprintf(comm, 200, "The Laguerre polynomial of degree %d with precision "
             PRECF_STR " (" TYPEFP_STR ")", d);
    
    if(! write_polyfr(P, argv[2], pcomm, 1)) {
        printf("Could not write the Laguerre polynomial to %s ! \n", argv[2]);
        polyfr_free(P);
        
        return APP_ERROR;
    }
    
    printf("The Laguerre polynomial of degree %d is written to %s.\n", d, argv[2]);
    polyfr_free(P);
    
    return APP_OK;
}

static polyf read_polyf(const char *fileName) {
    arrayf cf1 = arrayf_read((char *) fileName, true);
    
    if(cf1 == NULL) {
        printf("Could not read coefficients from %s ! \n", fileName);
        
        return NULL;
    }
    
    polyf P = arrayf_polyf(cf1);
    arrayf_free(cf1);
    
    if(P != NULL) {
        printf("Read a polynomial of degree %" FMT_DEG " from %s.\n", P->degree, fileName);
    }
    
    return P;
}

static bool write_polyf(polyf S, const char *fileName) {
    arrayf l = arrayf_new_polyf(S);
    if(l != NULL && (arrayf_first_inf(l) >= 0 || arrayf_first_nan(l) >= 0)) {
        printf("The polynomial has invalid coefficients !\n");
    }
    
    if(l == NULL || ! arrayf_write(l, (char *) fileName, true)) {
        printf("Could not write the result to %s !\n", fileName);
        arrayf_free(l);
        
        return false;
    }
    
    arrayf_free(l);
    
    return true;
}

static poly read_poly(const char *fileName, prec_t prec) {
    array cf1 = array_read((char *) fileName, prec, true);
    
    if(cf1 == NULL) {
        printf("Could not read coefficients from %s ! \n", fileName);
        
        return NULL;
    }
    
    poly P = array_poly(cf1, prec);
    array_free(cf1);
    
    if(P != NULL) {
        printf("Read a polynomial of degree %" FMT_DEG " from %s.\n", P->degree, fileName);
    }
    
    return P;
}

static bool write_poly(poly S, const char *fileName, prec_t prec) {
    array l = array_new_poly(S, prec);
    if(l != NULL && (array_first_inf(l) >= 0 || array_first_nan(l) >= 0)) {
        printf("The polynomial has invalid coefficients !\n");
    }
    
    if(l == NULL || ! array_write(l, (char *) fileName, (int) (prec / 3), true)) {
        printf("Could not write the result to %s ! \n", fileName);
        array_free(l);
        
        return false;
    }
    
    array_free(l);
    
    return true;
}

static int sum(long prec, int argc, const char * argv[]) {
    if(argc < 4) {
        return APP_PARAMS;
    }
    
    poly P = read_poly(argv[1], prec);
    poly Q = read_poly(argv[2], prec);
    
    poly S = poly_sum(P, Q);
    poly_free(P);
    poly_free(Q);
    
    bool ok = write_poly(S, argv[3], prec);
    if(ok) {
        printf("The sum is written to %s.\n", argv[3]);
    }
    
    poly_free(S);
    
    return ok ? APP_OK : APP_ERROR;
}

static int sumf(int argc, const char * argv[]) {
    if(argc < 4) {
        return APP_PARAMS;
    }
    
    polyf P = read_polyf(argv[1]);
    polyf Q = read_polyf(argv[2]);
    
    polyf S = polyf_sum(P, Q);
    polyf_free(P);
    polyf_free(Q);
    
    bool ok = write_polyf(S, argv[3]);
    if(ok) {
        printf("The sum is written to %s.\n", argv[3]);
    }
    
    polyf_free(S);
    
    return ok ? APP_OK : APP_ERROR;
}

static int diff(long prec, int argc, const char * argv[]) {
    if(argc < 4) {
        return APP_PARAMS;
    }
    
    poly P = read_poly(argv[1], prec);
    poly Q = read_poly(argv[2], prec);
    
    poly S = poly_diff(P, Q);
    poly_free(P);
    poly_free(Q);
    
    bool ok = write_poly(S, argv[3], prec);
    if(ok) {
        printf("The difference is written to %s.\n", argv[3]);
    }
    
    poly_free(S);
    
    return ok ? APP_OK : APP_ERROR;
}

static int difff(int argc, const char * argv[]) {
    if(argc < 4) {
        return APP_PARAMS;
    }
    
    polyf P = read_polyf(argv[1]);
    polyf Q = read_polyf(argv[2]);
    
    polyf S = polyf_diff(P, Q);
    polyf_free(P);
    polyf_free(Q);
    
    bool ok = write_polyf(S, argv[3]);
    if(ok) {
        printf("The difference is written to %s.\n", argv[3]);
    }
    
    polyf_free(S);
    
    return ok ? APP_OK : APP_ERROR;
}

static int prod(long prec, int argc, const char * argv[]) {
    if(argc < 4) {
        return APP_PARAMS;
    }
    
    poly P = read_poly(argv[1], prec);
    poly Q = read_poly(argv[2], prec);
    
    poly S = poly_prod(P, Q);
    poly_free(P);
    poly_free(Q);
    
    bool ok = write_poly(S, argv[3], prec);
    if(ok) {
        printf("The product is written to %s.\n", argv[3]);
    }
    
    poly_free(S);
    
    return ok ? APP_OK : APP_ERROR;
}

static int prodf(int argc, const char * argv[]) {
    if(argc < 4) {
        return APP_PARAMS;
    }
    
    polyf P = read_polyf(argv[1]);
    polyf Q = read_polyf(argv[2]);
    
    polyf S = polyf_prod(P, Q);
    polyf_free(P);
    polyf_free(Q);
    
    bool ok = write_polyf(S, argv[3]);
    if(ok) {
        printf("The product is written to %s.\n", argv[3]);
    }
    
    polyf_free(S);
    
    return ok ? APP_OK : APP_ERROR;
}

static int concat(long prec, int argc, const char * argv[]) {
    if(argc < 4) {
        return APP_PARAMS;
    }
    
    array z1 = array_read((char *) argv[1], prec, true);
    if(z1 == NULL) {
        return APP_ERROR;
    }
    
    array z2 = array_read((char *) argv[2], prec, true);
    if(z2 == NULL) {
        array_free(z1);
        
        return APP_ERROR;
    }
    
    bool ok = array_write(z1, (char *) argv[3], 1 + prec / 3.2, true);
    ok = ok && array_append(z2, (char *) argv[3], 1 + prec / 3.2, true);
    if(ok) {
        printf("The concatenation with %" PRIu64 " points is written to %s.\n",
               z1->len + z2->len, argv[3]);
    } else {
        printf("IO error when writing the sequences to %s.\n", argv[3]);
    }
    
    array_free(z1);
    array_free(z2);
    
    return ok ? APP_OK : APP_ERROR;
}

static int concatf(int argc, const char * argv[]) {
    if(argc < 4) {
        return APP_PARAMS;
    }
    
    arrayf z1 = arrayf_read((char *) argv[1], true);
    if(z1 == NULL) {
        return APP_ERROR;
    }
    
    arrayf z2 = arrayf_read((char *) argv[2], true);
    if(z2 == NULL) {
        return APP_ERROR;
    }
    
    bool ok = arrayf_write(z1, (char *) argv[3], true);
    ok = ok && arrayf_append(z2, (char *) argv[3], true);
    if(ok) {
        printf("The concatenation with %" PRIu64 " points is written to %s.\n",
               z1->len + z2->len, argv[3]);
    } else {
        printf("IO error when writing the sequences to %s.\n", argv[3]);
    }
    
    arrayf_free(z1);
    arrayf_free(z2);
    
    return ok ? APP_OK : APP_ERROR;
}

static int join(long prec, int argc, const char * argv[]) {
    if(argc < 4) {
        return APP_PARAMS;
    }
    
    array z1 = array_read((char *) argv[1], prec, true);
    if(z1 == NULL) {
        return APP_ERROR;
    }
    
    array z2 = array_read((char *) argv[2], prec, true);
    if(z2 == NULL) {
        array_free(z1);
        
        return APP_ERROR;
    }
    
    array r = app_join(prec, z1, z2);
    
    array_free(z1);
    array_free(z2);
    
    if(r == NULL) {
        printf("Could not join the two lists !\n");
        
        return APP_ERROR;
    }
    
    bool ok = array_write(r, (char *) argv[3], 1 + prec / 3.2, true);
    if(ok) {
        printf("The joint sequence of length %" PRIu64 " is written to %s.\n", r->len, argv[3]);
    } else {
        printf("Could not write to %s !\n", argv[3]);
    }
    
    array_free(r);
    
    return ok ? APP_OK : APP_ERROR;
}

static int joinf(int argc, const char * argv[]) {
    if(argc < 4) {
        return APP_PARAMS;
    }
    
    arrayf z1 = arrayf_read((char *) argv[1], true);
    if(z1 == NULL) {
        return APP_ERROR;
    }
    
    arrayf z2 = arrayf_read((char *) argv[2], true);
    if(z2 == NULL) {
        arrayf_free(z1);
        
        return APP_ERROR;
    }
    
    arrayf r = appf_join(z1, z2);
    
    arrayf_free(z1);
    arrayf_free(z2);
    
    bool ok = arrayf_write(r, (char *) argv[3], true);
    if(ok) {
        printf("The joint sequence of length %" PRIu64 " is written to %s.\n", r->len, argv[3]);
    }
    
    arrayf_free(r);
    
    return ok ? APP_OK : APP_ERROR;
}

static int grid(long prec, int argc, const char * argv[]) {
    if(argc < 4) {
        return APP_PARAMS;
    }
    
    array z1 = array_read((char *) argv[1], prec, true);
    if(z1 == NULL) {
        return APP_ERROR;
    }
    
    array z2 = array_read((char *) argv[2], prec, true);
    if(z2 == NULL) {
        array_free(z1);
        
        return APP_ERROR;
    }
    
    array r = app_grid(prec, z1, z2);
        
    array_free(z1);
    array_free(z2);
    
    if(r == NULL) {
        printf("Could not create a grid !\n");
        
        return APP_ERROR;
    }
    
    bool ok = array_write(r, (char *) argv[3], 1 + prec / 3.2, true);
    if(ok) {
        printf("The grid with %" PRIu64 " points is written to %s.\n", r->len, argv[3]);
    }
    
    array_free(r);
    
    return ok ? APP_OK : APP_ERROR;
}

static int gridf(int argc, const char * argv[]) {
    if(argc < 4) {
        return APP_PARAMS;
    }
    
    arrayf z1 = arrayf_read((char *) argv[1], true);
    if(z1 == NULL) {
        return APP_ERROR;
    }
    
    arrayf z2 = arrayf_read((char *) argv[2], true);
    if(z2 == NULL) {
        return APP_ERROR;
    }
    
    arrayf r = appf_grid(z1, z2);
    
    arrayf_free(z1);
    arrayf_free(z2);
    
    bool ok = arrayf_write(r, (char *) argv[3], true);
    if(ok) {
        printf("The grid with %" PRIu64 " points is written to %s.\n", r->len, argv[3]);
    }
    
    arrayf_free(r);
    
    return ok ? APP_OK : APP_ERROR;
}

static int expz(long prec, int argc, const char * argv[]) {
    if(argc < 3) {
        return APP_PARAMS;
    }

    array z = array_read((char *) argv[1], prec, true);
    if(z == NULL) {
        return APP_ERROR;
    }
    
    array ez = app_exp(prec, z);
    bool ok = ez != NULL && ez->len == z->len
                && array_write(ez, (char *) argv[2], 1 + prec / 3.2, true);
    
    if(! ok) {
        printf("Could not write the result to %s.\n", argv[2]);
    } else {
        printf("The exponential is written to %s.\n", argv[2]);
    }

    array_free(z);
    array_free(ez);
    
    return ok ? APP_OK : APP_ERROR;
}

static int expzf(int argc, const char * argv[]) {
    if(argc < 3) {
        return APP_PARAMS;
    }
    
    arrayf z = arrayf_read((char *) argv[1], true);
    if(z == NULL) {
        return APP_ERROR;
    }
    
    arrayf ez = appf_exp(z);
    bool ok = ez != NULL && ez->len == z->len && arrayf_write(ez, (char *) argv[2], true);
    
    if(! ok) {
        printf("Could not write the result to %s.\n", argv[2]);
    } else {
        printf("The exponential is written to %s.\n", argv[2]);
    }
    
    arrayf_free(z);
    arrayf_free(ez);
    
    return ok ? APP_OK : APP_ERROR;
}

static int rot(long prec, int argc, const char * argv[]) {
    if(argc < 3) {
        return APP_PARAMS;
    }

    array z = array_read((char *) argv[1], prec, true);
    
    if(z == NULL) {
        return APP_ERROR;
    }
    
    array ez = app_rot(prec, z);
    bool ok = ez != NULL && ez->len == z->len
                && array_write(ez, (char *) argv[2], 1 + prec / 3.2, true);
    
    if(! ok) {
        printf("Could not write the result to %s.\n", argv[2]);
    } else {
        printf("The rotation is written to %s.\n", argv[2]);
    }

    array_free(z);
    array_free(ez);
    
    return ok ? APP_OK : APP_ERROR;
}

static int rotf(int argc, const char * argv[]) {
    if(argc < 3) {
        return APP_PARAMS;
    }
    
    arrayf z = arrayf_read((char *) argv[1], true);
    if(z == NULL) {
        return APP_ERROR;
    }
    
    arrayf ez = appf_rot(z);
    bool ok = ez != NULL && ez->len == z->len && arrayf_write(ez, (char *) argv[2], true);
    
    if(! ok) {
        printf("Could not write the result to %s.\n", argv[2]);
    } else {
        printf("The rotation is written to %s.\n", argv[2]);
    }
    
    arrayf_free(z);
    arrayf_free(ez);
    
    return ok ? APP_OK : APP_ERROR;
}

static int polar(long prec, int argc, const char * argv[]) {
    if(argc < 3) {
        return APP_PARAMS;
    }

    array z = array_read((char *) argv[1], prec, true);
    if(z == NULL) {
        return APP_ERROR;
    }
    
    array ez = app_polar(prec, z);
    bool ok = ez != NULL && ez->len == z->len
                && array_write(ez, (char *) argv[2], 1 + prec / 3.2, true);
    
    if(! ok) {
        printf("Could not write the result to %s.\n", argv[2]);
    } else {
        printf("The projected points have been written to %s.\n", argv[2]);
    }
    
    array_free(z);
    array_free(ez);
    
    return ok ? APP_OK : APP_ERROR;
}

static int polarf(int argc, const char * argv[]) {
    if(argc < 3) {
        return APP_PARAMS;
    }
    
    arrayf z = arrayf_read((char *) argv[1], true);
    if(z == NULL) {
        return APP_ERROR;
    }
    
    arrayf ez = appf_polar(z);
    bool ok = ez != NULL && ez->len == z->len && arrayf_write(ez, (char *) argv[2], true);
    
    if(! ok) {
        printf("Could not write the result to %s.\n", argv[2]);
    } else {
        printf("The projected points have been written to %s.\n", argv[2]);
    }
    
    arrayf_free(z);
    arrayf_free(ez);
    
    return ok ? APP_OK : APP_ERROR;
}

static int roots(long prec, int argc, const char * argv[]) {
    if(argc < 3) {
        return APP_PARAMS;
    }
    
    array z = array_read((char *) argv[1], prec, true);
    if(z == NULL || z->len == 0) {
        return APP_ERROR;
    }
    
    poly S = poly_from_roots(z->zi, z->len, prec);
    
    bool ok = write_poly(S, argv[2], prec);
    if(ok) {
        printf("The polynomial is written to %s.\n", argv[2]);
    } else {
        printf("Could not write the polynomial!\n");
    }
    
    poly_free(S);
    
    return ok ? APP_OK : APP_ERROR;
}

static int rootsf(int argc, const char * argv[]) {
    if(argc < 3) {
        return APP_PARAMS;
    }
    
    arrayf z = arrayf_read((char *) argv[1], true);
    if(z == NULL || z->len == 0) {
        return APP_ERROR;
    }
    
    polyf S = polyf_from_roots(z->zi, z->len);
    
    bool ok = write_polyf(S, argv[2]);
    if(ok) {
        printf("The polynomial is written to %s.\n", argv[2]);
    } else {
        printf("Could not write the polynomial!\n");
    }
    
    polyf_free(S);
    
    return ok ? APP_OK : APP_ERROR;
}

static int der(long prec, int argc, const char * argv[]) {
    if(argc < 3) {
        return APP_PARAMS;
    }
    
    poly P = read_poly(argv[1], prec);
    if(P == NULL) {
        printf("Could not read the coefficients of the polynomial from %s.\n", argv[1]);
        
        return APP_ERROR;
    }
    
    poly S = poly_derivative(P);
    poly_free(P);
    
    bool ok = write_poly(S, argv[2], prec);
    if(ok) {
        printf("The derivative is written to %s.\n", argv[2]);
    } else {
        printf("Could not write the derivative!\n");
    }
    
    poly_free(S);
    
    return ok ? APP_OK : APP_ERROR;
}

static int derf(int argc, const char * argv[]) {
    if(argc < 3) {
        return APP_PARAMS;
    }
    
    polyf P = read_polyf(argv[1]);
    if(P == NULL) {
        printf("Could not read the coefficients of the polynomial from %s.\n", argv[1]);
        
        return APP_ERROR;
    }
    
    polyf S = polyf_derivative(P);
    polyf_free(P);
    
    bool ok = write_polyf(S, argv[2]);
    if(ok) {
        printf("The derivative is written to %s.\n", argv[2]);
    } else {
        printf("Could not write the derivative!\n");
    }
    
    polyf_free(S);
    
    return ok ? APP_OK : APP_ERROR;
}

static int unif(long prec, int argc, const char * argv[]) {
    long n;

    if(argc < 3 || sscanf(argv[1], "%ld", &n) < 1 || n < 2) {
        return APP_PARAMS;
    }

    mpfr_t st, en;
    mpfr_init2(st, prec);
    mpfr_set_zero(st, 1);
    
    mpfr_init2(en, prec);
    mpfr_set_ui(en, 1, MPFR_RNDN);
    
    if(argc >= 4) {
        if(mpfr_set_str(st, argv[3], 10, MPFR_RNDN) != 0) {
            return APP_PARAMS;
        }
    }
    
    if(argc >= 5) {
        if(mpfr_set_str(en, argv[4], 10, MPFR_RNDN) != 0) {
            return APP_PARAMS;
        }
    }
        
    array ez = app_unif(prec, n, st, en);
    bool ok = ez != NULL && ez->len == n
                && array_write(ez, (char *) argv[2], 1 + prec / 3.2, true);
    
    if(! ok) {
        printf("Could not write the result to %s.\n", argv[2]);
    } else {
        printf("The points are written to %s.\n", argv[2]);
    }

    mpfr_clear(st);
    mpfr_clear(en);
    
    array_free(ez);
    
    return ok ? APP_OK : APP_ERROR;
}

static int uniff(int argc, const char * argv[]) {
    long n;
    
    if(argc < 3 || sscanf(argv[1], "%ld", &n) < 1 || n < 2) {
        return APP_PARAMS;
    }
    
    coeff_t st = 0, en = 1;
    
    if(argc >= 4) {
        if(sscanf(argv[3], "%" FMT_COEFF "g", &st) < 1) {
            return APP_PARAMS;
        }
    }
    
    if(argc >= 5) {
        if(sscanf(argv[4], "%" FMT_COEFF "g", &en) < 1) {
            return APP_PARAMS;
        }
    }
    
    arrayf z = appf_unif(n, st, en);
    bool ok = z != NULL && z->len == n && arrayf_write(z, (char *) argv[2], true);
    
    if(! ok) {
        printf("Could not write the result to %s.\n", argv[2]);
    } else {
        printf("The points are written to %s.\n", argv[2]);
    }
    
    arrayf_free(z);
    
    return ok ? APP_OK : APP_ERROR;
}

static int ev_rand(long prec, int argc, const char * argv[]) {
    long n;

    if(argc < 3 || sscanf(argv[1], "%ld", &n) < 1 || n < 1) {
        return APP_PARAMS;
    }

    mpfr_t st, en;
    mpfr_init2(st, prec);
    mpfr_set_zero(st, 1);
    
    mpfr_init2(en, prec);
    mpfr_set_ui(en, 1, MPFR_RNDN);
    
    if(argc >= 4) {
        if(mpfr_set_str(st, argv[3], 10, MPFR_RNDN) != 0) {
            return APP_PARAMS;
        }
    }
    
    if(argc >= 5) {
        if(mpfr_set_str(en, argv[4], 10, MPFR_RNDN) != 0) {
            return APP_PARAMS;
        }
    }
    
    array ez = app_rand(prec, n, st, en);
    bool ok = ez != NULL && ez->len == n
                && array_write(ez, (char *) argv[2], 1 + prec / 3.2, true);
    
    if(! ok) {
        printf("Could not write the result to %s.\n", argv[2]);
    } else {
        printf("The random points are written to %s.\n", argv[2]);
    }
    
    mpfr_clear(st);
    mpfr_clear(en);

    array_free(ez);
    
    return ok ? APP_OK : APP_ERROR;
}

static int ev_randf(int argc, const char * argv[]) {
    long prec;

    if(argc < 1 || sscanf(argv[0], "%ld", &prec) < 1) {
        return APP_PARAMS;
    }

    return ev_rand(prec, argc, argv);
}

static int normal(long prec, int argc, const char * argv[]) {
    long n;

    if(argc < 3 || sscanf(argv[1], "%ld", &n) < 1 || n < 1) {
        return APP_PARAMS;
    }

    mpfr_t cent, var;
    mpfr_init2(cent, prec);
    mpfr_set_zero(cent, 1);
    
    mpfr_init2(var, prec);
    mpfr_set_ui(var, 1, MPFR_RNDN);
    
    if(argc >= 4) {
        if(mpfr_set_str(cent, argv[3], 10, MPFR_RNDN) != 0) {
            return APP_PARAMS;
        }
    }
    
    if(argc >= 5) {
        if(mpfr_set_str(var, argv[4], 10, MPFR_RNDN) != 0) {
            return APP_PARAMS;
        }
    }
        
    array ez = app_normal(prec, n, cent, var);
    bool ok = ez != NULL && ez->len == n
                && array_write(ez, (char *) argv[2], 1 + prec / 3.2, true);
    
    if(! ok) {
        printf("Could not write the result to %s.\n", argv[2]);
    } else {
        printf("The normally distributed points are written to %s.\n", argv[2]);
    }

    mpfr_clear(cent);
    mpfr_clear(var);
    
    array_free(ez);
    
    return ok ? APP_OK : APP_ERROR;
}

static int normalf(int argc, const char * argv[]) {
    long prec;

    if(argc < 1 || sscanf(argv[0], "%ld", &prec) < 1) {
        return APP_PARAMS;
    }

    return normal(prec, argc, argv);
}

static int sphere(long prec, int argc, const char * argv[]) {
    long n;

    if(argc < 3 || sscanf(argv[1], "%ld", &n) < 1 || n < 1) {
        return APP_PARAMS;
    }
    
    array z = app_sphere(prec, n);
    bool ok = array_write(z, (char *) argv[2], 1 + prec / 3.2, true);
    
    if(! ok) {
        printf("Could not write the points to %s.\n", argv[2]);
    } else {
        printf("%" PRIu64 " polar coordinates are written to %s.\n", z->len, argv[2]);
    }

    array_free(z);
    
    return ok ? APP_OK : APP_ERROR;
}

static int spheref(int argc, const char * argv[]) {
    long prec;

    if(argc < 1 || sscanf(argv[0], "%ld", &prec) < 1) {
        return APP_PARAMS;
    }
    
    return sphere(prec, argc, argv);
}

static int eval_pf(int argc, const char * argv[]) {
    if(argc < 4) {
        return APP_PARAMS;
    }
    
    ulong count = 1;
    if(argc > 5) {
        if(sscanf(argv[5], "%" PRIu64, &count) < 1 || count < 1) {
            return APP_PARAMS;
        }
    }
    
    // machine floating numbers
    polyf P = read_polyf(argv[1]);
    if(P == NULL) {
        printf("Could not read the polynomial from %s.\n", argv[1]);
        
        return APP_ERROR;
    }
    
    arrayf z = arrayf_read((char *) argv[2], true);
    if(z == NULL || z->len == 0) {
        return APP_ERROR;
    }
    
    char *outError = argc > 4 ? (char *) argv[4] : NULL;
    char *outHorner = argc > 6 ? (char *) argv[6] : NULL;
    ulong c2 = 1;
    if(argc > 7) {
        sscanf(argv[7], "%" PRIu64, &c2);
    }
    
    bool ok = appf_eval_p(P, z, (char *) argv[3], outError, count, outHorner, c2, argv + 1);
    
    return ok ? APP_OK : APP_ERROR;
}

static int eval_p(long prec, int argc, const char * argv[]) {
    ulong count = 1;
    if(argc >= 6) {
        if(sscanf(argv[5], "%" PRIu64, &count) < 1 || count < 1) {
            return APP_PARAMS;
        }
    }
    
    // MPFR numbers with precision prec
    poly P = read_poly(argv[1], prec);
    if(P == NULL) {
        printf("Could not read the polynomial from %s.\n", argv[1]);
        
        return APP_ERROR;
    }
    
    array z = array_read((char *) argv[2], prec, true);
    if(z == NULL || z->len == 0) {
        return APP_ERROR;
    }
    
    char *outError = argc > 4 ? (char *) argv[4] : NULL;
    char *outHorner = argc > 6 ? (char *) argv[6] : NULL;
    ulong c2 = 1;
    if(argc > 7) {
        sscanf(argv[7], "%" PRIu64, &c2);
    }
    
    bool ok = app_eval_p(prec, P, z, (char *) argv[3], outError, count, outHorner, c2, argv + 1);
    
    return ok ? APP_OK : APP_ERROR;
}

static int eval_df(int argc, const char * argv[]) {
    ulong count = 1;
    if(argc >= 6) {
        if(sscanf(argv[5], "%" PRIu64, &count) < 1 || count < 1) {
            return APP_PARAMS;
        }
    }
    
    // machine floating numbers
    polyf P = read_polyf(argv[1]);
    if(P == NULL) {
        printf("Could not read the polynomial from %s.\n", argv[1]);
        
        return APP_ERROR;
    }
    
    arrayf z = arrayf_read((char *) argv[2], true);
    if(z == NULL || z->len == 0) {
        return APP_ERROR;
    }
    
    char *outError = argc > 4 ? (char *) argv[4] : NULL;
    char *outHorner = argc > 6 ? (char *) argv[6] : NULL;
    ulong c2 = 1;
    if(argc > 7) {
        sscanf(argv[7], "%" PRIu64, &c2);
    }
    
    bool ok = appf_eval_d(P, z, (char *) argv[3], outError, count, outHorner, c2, argv + 1);
    
    return ok ? APP_OK : APP_ERROR;
}

static int eval_d(long prec, int argc, const char * argv[]) {
    ulong count = 1;
    if(argc >= 6) {
        if(sscanf(argv[5], "%" PRIu64, &count) < 1 || count < 1) {
            return APP_PARAMS;
        }
    }
    
    // MPFR numbers with precision prec
    poly P = read_poly(argv[1], prec);
    if(P == NULL) {
        printf("Could not read the polynomial from %s.\n", argv[1]);
        
        return APP_ERROR;
    }
    
    array z = array_read((char *) argv[2], prec, true);
    if(z == NULL || z->len == 0) {
        return APP_ERROR;
    }
    
    char *outError = argc > 4 ? (char *) argv[4] : NULL;
    char *outHorner = argc > 6 ? (char *) argv[6] : NULL;
    ulong c2 = 1;
    if(argc > 7) {
        sscanf(argv[7], "%" PRIu64, &c2);
    }
    
    bool ok = app_eval_d(prec, P, z, (char *) argv[3], outError, count, outHorner, c2, argv + 1);
    
    return ok ? APP_OK : APP_ERROR;
}

static int eval_nf(int argc, const char * argv[]) {
    ulong count = 1;
    if(argc >= 6) {
        if(sscanf(argv[5], "%" PRIu64, &count) < 1 || count < 1) {
            return APP_PARAMS;
        }
    }
    
    // machine floating numbers
    polyf P = read_polyf(argv[1]);
    if(P == NULL) {
        printf("Could not read the polynomial from %s.\n", argv[1]);
        
        return APP_ERROR;
    }
    
    arrayf z = arrayf_read((char *) argv[2], true);
    if(z == NULL || z->len == 0) {
        return APP_ERROR;
    }
    
    char *outError = argc > 4 ? (char *) argv[4] : NULL;
    char *outHorner = argc > 6 ? (char *) argv[6] : NULL;
    ulong c2 = 1;
    if(argc > 7) {
        sscanf(argv[7], "%" PRIu64, &c2);
    }
    
    bool ok = appf_eval_n(P, z, (char *) argv[3], outError, count, outHorner, c2, argv + 1);
    
    return ok ? APP_OK : APP_ERROR;
}

static int eval_n(long prec, int argc, const char * argv[]) {
    ulong count = 1;
    if(argc >= 6) {
        if(sscanf(argv[5], "%" PRIu64, &count) < 1 || count < 1) {
            return APP_PARAMS;
        }
    }
    
    // MPFR numbers with precision prec
    poly P = read_poly(argv[1], prec);
    if(P == NULL) {
        printf("Could not read the polynomial from %s.\n", argv[1]);
        
        return APP_ERROR;
    }
    
    array z = array_read((char *) argv[2], prec, true);
    if(z == NULL || z->len == 0) {
        return APP_ERROR;
    }
    
    char *outError = argc > 4 ? (char *) argv[4] : NULL;
    char *outHorner = argc > 6 ? (char *) argv[6] : NULL;
    ulong c2 = 1;
    if(argc > 7) {
        sscanf(argv[7], "%" PRIu64, &c2);
    }
    
    bool ok = app_eval_n(prec, P, z, (char *) argv[3], outError, count, outHorner, c2, argv + 1);
    
    return ok ? APP_OK : APP_ERROR;
}

static int eval_rf(int argc, const char * argv[]) {
    if (argc < 5) {
        return APP_PARAMS;
    }
    
    ulong maxIter = 0;
    if(sscanf(argv[3], "%" PRIu64, &maxIter) < 1 || maxIter < 1) {
        return APP_PARAMS;
    }
    
    ulong count = 1;
    if(argc >= 7) {
        if(sscanf(argv[6], "%" PRIu64, &count) < 1 || count < 1) {
            return APP_PARAMS;
        }
    }
    
    // machine floating numbers
    polyf P = read_polyf(argv[1]);
    if(P == NULL) {
        printf("Could not read the polynomial from %s.\n", argv[1]);
        
        return APP_ERROR;
    }
    
    arrayf z = arrayf_read((char *) argv[2], true);
    if(z == NULL || z->len == 0) {
        return APP_ERROR;
    }
    
    char *outError = argc > 5 ? (char *) argv[5] : NULL;
    char *outHorner = argc > 7 ? (char *) argv[7] : NULL;
    ulong c2 = 1;
    if(argc > 8) {
        sscanf(argv[8], "%" PRIu64, &c2);
    }
    
    bool ok = appf_eval_r(P, z, maxIter, (char *) argv[4], outError, count, outHorner, c2, argv + 1);
    
    return ok ? APP_OK : APP_ERROR;
}

static int eval_r(long prec, int argc, const char * argv[]) {
    ulong maxIter = 0;
    if(sscanf(argv[3], "%" PRIu64, &maxIter) < 1 || maxIter < 1) {
        return APP_PARAMS;
    }
    
    ulong count = 1;
    if(argc >= 7) {
        if(sscanf(argv[6], "%" PRIu64, &count) < 1 || count < 1) {
            return APP_PARAMS;
        }
    }
    
    // MPFR numbers with precision prec
    poly P = read_poly(argv[1], prec);
    if(P == NULL) {
        printf("Could not read the polynomial from %s.\n", argv[1]);
        
        return APP_ERROR;
    }
    
    array z = array_read((char *) argv[2], prec, true);
    if(z == NULL || z->len == 0) {
        return APP_ERROR;
    }
    
    char *outError = argc > 5 ? (char *) argv[5] : NULL;
    char *outHorner = argc > 7 ? (char *) argv[7] : NULL;
    ulong c2 = 1;
    if(argc > 8) {
        sscanf(argv[8], "%" PRIu64, &c2);
    }
    
    bool ok = app_eval_r(prec, P, z, maxIter, (char *) argv[4],
                         outError, count, outHorner, c2, argv + 1);
    
    return ok ? APP_OK : APP_ERROR;
}

static int comparef(int argc, const char * argv[]) {
    arrayf z1 = arrayf_read((char *) argv[1], true);
    if(z1 == NULL) {
        return APP_ERROR;
    }
    
    arrayf z2 = arrayf_read((char *) argv[2], true);
    if(z2 == NULL) {
        arrayf_free(z1);
        
        return APP_ERROR;
    }
    
    bool ok = appf_compare(z1, z2, argc > 3 ? (char *) argv[3] : NULL);
    
    arrayf_free(z1);
    arrayf_free(z2);
    
    return ok ? APP_OK : APP_ERROR;
}

static int compare(long prec, int argc, const char * argv[]) {
    array z1 = array_read((char *) argv[1], prec, true);
    if(z1 == NULL) {
        return APP_ERROR;
    }
    
    array z2 = array_read((char *) argv[2], prec, true);
    if(z2 == NULL) {
        array_free(z1);
        
        return APP_ERROR;
    }
    
    bool ok = app_compare(prec, z1, z2, argc > 3 ? (char *) argv[3] : NULL);
    
    array_free(z1);
    array_free(z2);
    
    return ok ? APP_OK : APP_ERROR;
}

static int ref(int argc, const char * argv[]) {
    arrayf z = arrayf_read((char *) argv[1], true);
    if(z == NULL) {
        return APP_ERROR;
    }
    
    arrayf ez = appf_re(z);
    
    bool ok = ez != NULL && ez->len == z->len && arrayf_write(ez, (char *) argv[2], true);
    
    if(! ok) {
        printf("Could not write the result to %s.\n", argv[2]);
    } else {
        printf("The real parts are written to %s.\n", argv[2]);
    }
    
    arrayf_free(z);
    arrayf_free(ez);
    
    return ok ? APP_OK : APP_ERROR;
}

static int re(long prec, int argc, const char * argv[]) {
    array z = array_read((char *) argv[1], prec, true);
    if(z == NULL) {
        return APP_ERROR;
    }
    
    array ez = app_re(prec, z);
    bool ok = ez != NULL && ez->len == z->len
                && array_write(ez, (char *) argv[2], 1 + prec / 3.2, true);
    
    if(! ok) {
        printf("Could not write the result to %s.\n", argv[2]);
    } else {
        printf("The real parts are written to %s.\n", argv[2]);
    }
    
    array_free(z);
    array_free(ez);
    
    return ok ? APP_OK : APP_ERROR;
}

static int imf(int argc, const char * argv[]) {
    arrayf z = arrayf_read((char *) argv[1], true);
    if(z == NULL) {
        return APP_ERROR;
    }
    
    arrayf ez = appf_im(z);
    bool ok = ez != NULL && ez->len == z->len && arrayf_write(ez, (char *) argv[2], true);
    
    if(! ok) {
        printf("Could not write the result to %s.\n", argv[2]);
    } else {
        printf("The real parts are written to %s.\n", argv[2]);
    }
    
    arrayf_free(z);
    arrayf_free(ez);
    
    return ok ? APP_OK : APP_ERROR;
}

static int im(long prec, int argc, const char * argv[]) {
    array z = array_read((char *) argv[1], prec, true);
    if(z == NULL) {
        return APP_ERROR;
    }
    
    array ez = app_im(prec, z);
    bool ok = ez != NULL && ez->len == z->len
                && array_write(ez, (char *) argv[2], 1 + prec / 3.2, true);
    
    if(! ok) {
        printf("Could not write the result to %s.\n", argv[2]);
    } else {
        printf("The real parts are written to %s.\n", argv[2]);
    }
    
    array_free(z);
    array_free(ez);
    
    return ok ? APP_OK : APP_ERROR;
}

static int conjugatef(int argc, const char * argv[]) {
    arrayf z = arrayf_read((char *) argv[1], true);
    if(z == NULL) {
        return APP_ERROR;
    }
    
    arrayf ez = appf_conj(z);
    bool ok = ez != NULL && ez->len == z->len && arrayf_write(ez, (char *) argv[2], true);
    
    if(! ok) {
        printf("Could not write the result to %s.\n", argv[2]);
    } else {
        printf("The conjugates are written to %s.\n", argv[2]);
    }
    
    arrayf_free(z);
    arrayf_free(ez);
    
    return ok ? APP_OK : APP_ERROR;
}

static int conjugate(long prec, int argc, const char * argv[]) {
    array z = array_read((char *) argv[1], prec, true);
    if(z == NULL) {
        return APP_ERROR;
    }
    
    array ez = app_conj(prec, z);
    bool ok = ez != NULL && ez->len == z->len
                && array_write(ez, (char *) argv[2], 1 + prec / 3.2, true);
    
    if(! ok) {
        printf("Could not write the result to %s.\n", argv[2]);
    } else {
        printf("The conjugates are written to %s.\n", argv[2]);
    }
    
    array_free(z);
    array_free(ez);
    
    return ok ? APP_OK : APP_ERROR;
}

static int tensor(long prec, int argc, const char * argv[]) {
    array z1 = array_read((char *) argv[1], prec, true);
    if(z1 == NULL) {
        return APP_ERROR;
    }
    
    array z2 = array_read((char *) argv[2], prec, true);
    if(z2 == NULL) {
        array_free(z1);
        
        return APP_ERROR;
    }
    
    array r = app_tensor(prec, z1, z2);
    
    array_free(z1);
    array_free(z2);
    
    if(r == NULL) {
        printf("Out of memory !\n");
        
        return APP_ERROR;
    }
    
    bool ok = array_write(r, (char *) argv[3], 1 + prec / 3.2, true);
    if(ok) {
        printf("The tensor product sequence of length %" PRIu64 " is written to %s.\n", r->len, argv[3]);
    }
    
    array_free(r);
    
    return ok ? APP_OK : APP_ERROR;
}

static int tensorf(int argc, const char * argv[]) {
    arrayf z1 = arrayf_read((char *) argv[1], true);
    if(z1 == NULL) {
        return APP_ERROR;
    }
    
    arrayf z2 = arrayf_read((char *) argv[2], true);
    if(z2 == NULL) {
        arrayf_free(z1);
        
        return APP_ERROR;
    }
    
    arrayf r = appf_tensor(z1, z2);
    
    arrayf_free(z1);
    arrayf_free(z2);
    
    bool ok = arrayf_write(r, (char *) argv[3], true);
    if(ok) {
        printf("The tensor product sequence of length %" PRIu64 " is written to %s.\n",
               r->len, argv[3]);
    }
    
    arrayf_free(r);
    
    return ok ? APP_OK : APP_ERROR;
}

static int analyse(long prec, int argc, const char * argv[]) {
    if(argc < 3) {
        return APP_PARAMS;
    }
    
    if(prec < MIN_EVAL_PREC) {
        printf("Please use at least " MIN_EVAL_PREC_STR " bits of precision");
        
        return APP_PARAMS;
    }
    
    poly P = read_poly(argv[1], 120);
    if(P == NULL) {
        printf("Could not read the polynomial from %s.\n", argv[1]);
        
        return APP_ERROR;
    }
    
    char *outChange = argc < 4 ? NULL : (char *) argv[3];
    bool ok = app_analyse(prec, P, (char *) argv[2], outChange);
    
    poly_free(P);
    
    return ok ? APP_OK : APP_ERROR;
}

static int analysef(int argc, const char * argv[]) {
    long prec;
    if(argc < 3 || sscanf(argv[0], "%ld", &prec) < 1) {
        return APP_PARAMS;
    }
    
    return analyse(prec, argc, argv);
}

// MARK: the list of available utilities

static mainFunc programs[APP_COUNT] = {
    &hyp,
    &cheb,
    &leg,
    &her,
    &lag,
    &sum,
    &diff,
    &prod,
    &concat,
    &join,
    &grid,
    &expz,
    &rot,
    &polar,
    &roots,
    &der,
    &unif,
    &ev_rand,
    &normal,
    &sphere,
    &eval_p,
    &eval_d,
    &eval_n,
    &eval_r,
    &compare,
    &re,
    &im,
    &conjugate,
    &tensor,
    &analyse
};

static mainFuncf programsf[APP_COUNT] = {
    &hypf,
    &chebf,
    &legf,
    &herf,
    &lagf,
    &sumf,
    &difff,
    &prodf,
    &concatf,
    &joinf,
    &gridf,
    &expzf,
    &rotf,
    &polarf,
    &rootsf,
    &derf,
    &uniff,
    &ev_randf,
    &normalf,
    &spheref,
    &eval_pf,
    &eval_df,
    &eval_nf,
    &eval_rf,
    &comparef,
    &ref,
    &imf,
    &conjugatef,
    &tensorf,
    &analysef
};

/// The entry point in the app. depending on the first agument, it launches one of the incorporated
/// utilities or prints a @c help screen.
///
/// @param argc the number of command line argument
/// @param argv the command line arguments
/// @return @c 0 if successfull, @c 1 otherwise
int main(int argc, const char * argv[]) {
    if(! ntypes_check()) {
        printf("Sorry, your system is not comptible with the choices from ntypes.h.\n"
               "Please check and update.\n");
        
        return 1;
    }
    
    help hm = fpe_help_get(HELP_MAIN);
    
    if(argc < 2) {
        fpe_help_print(hm);
        
        return 1;
    }
    
    const char *command = argv[1];
    for(int i = 0; i < hm->linesCount; i++) {
        if(! strncmp(command, hm->lines[0][i], 30)) {
            if(argc > 2 && ! strncmp("-help", argv[2], 10)) {
                fpe_help_print(fpe_help_get(i + 1));
                
                return 1;
            } else {
                long prec = -1;
                int err = -1;
                if(argc > 2 && sscanf(argv[2], "%ld", &prec) == 1) {
                    prec = prec < precf ? precf : prec;
                    bool machine = prec <= precf;
                    
                    if(machine) {
                        printf("Working with machine numbers with %ld bits (" TYPEFP_STR
                               ", " TYPEF_STR ").\n", prec);
                        
                        err = programsf[i](argc - 2, argv + 2);
                    } else {
                        printf("Working with MPFR with precision %ld bits.\n", prec);
                        
                        err = programs[i](prec, argc - 2, argv + 2);
                    }
                } else {
                    printf("Could not read the precision !\n\n");
                    
                    err = APP_PARAMS;
                }
                
                switch(err) {
                    case APP_OK:
                        printf("\n");
                        return 0;
                        
                    case APP_PARAMS:
                        fpe_help_print(fpe_help_get(i + 1));
                        printf("\n");
                        return 1;
                        
                    case APP_ERROR:
                        printf("\n\nThe operation did not succeed !!!\n\n");
                        return 2;
                }
                
                return err;
            }
        }
    }
    
    fpe_help_print(hm);
    
    return 1;
}
