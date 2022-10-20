//
//  array.c
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

#include "array.h"

array array_new(ulong size) {
    ulong max = UINT64_MAX / sizeof(comp_struct);
    if(size > max) {
        return NULL;
    }
    
    array l = malloc(sizeof(array_struct));
    if(l == NULL) {
        return NULL;
    }
    
    l->len = 0;
    l->size = size < ARRAY_MIN_SIZE ? ARRAY_MIN_SIZE : size;
    
    l->zi = malloc(sizeof(comp_struct) * l->size);
    if(l->zi == NULL) {
        free(l);
        
        return NULL;
    }
    
    return l;
}

array array_new_polyr(polyr P, prec_t prec) {
    if(P == NULL || P->degree >= UINT64_MAX) {
        return NULL;
    }
    
    ulong n = P->degree + 1;
    array l = array_new(n);
    if(l == NULL) {
        return NULL;
    }
    
    for (ulong i = 0; i < n; i++) {
        comp_init(l->zi + i, prec);
        mpfr_set(l->zi[i].x, P->a + i, MPFR_RNDN);
        mpfr_set_zero(l->zi[i].y, 1);
    }
    
    l->len = n;
    
    return l;
}

array array_new_poly(poly P, prec_t prec) {
    if(P == NULL || P->degree >= UINT64_MAX) {
        return NULL;
    }
    
    ulong n = P->degree + 1;
    array l = array_new(n);
    if(l == NULL) {
        return NULL;
    }
    
    for (ulong i = 0; i < n; i++) {
        comp_init(l->zi + i, prec);
        comp_set(l->zi + i, P->a + i);
    }
    
    l->len = n;
    
    return l;
}

bool array_free(array l) {
    if(l == NULL || l->zi == NULL) {
        return false;
    }
    
    for (ulong i = 0; i < l->len; i++) {
        comp_clear(l->zi + i);
    }
    
    free(l->zi);
    free(l);
    
    return true;
}

bool array_add(array l, comp z, prec_t prec) {
    if(l == NULL || z == NULL) {
        return false;
    }
    
    if(l->len == l->size) {
        l->size *= ARRAY_SIZE_INCREASE;
        
        l->zi = realloc(l->zi, sizeof(comp_struct) * l->size);
        
        if(l->zi == NULL) {
            return false;
        }
    }
    
    comp_init(l->zi + l->len, prec);
    comp_set(l->zi + l->len, z);
    l->len ++;
    
    return true;
}

comp_ptr array_get(array l, ulong pos) {
    if(l == NULL || pos >= l->len) {
        return NULL;
    }
    
    return l->zi + pos;
}

bool array_is_real(array l) {
    if(l == NULL) {
        return false;
    }
    
    bool real = true;
    comp_ptr zi = l->zi;
    for (ulong i = 0; i < l->len && real; i++) {
        real &= mpfr_zero_p(zi[i].y);
    }
    
    return real;
}

long array_first_inf(array l) {
    if(l == NULL) {
        return LONG_MAX;
    }
    
    comp_ptr zi = l->zi;
    for (ulong i = 0; i < l->len; i++) {
        if(mpfr_inf_p(zi[i].x) || mpfr_inf_p(zi[i].y)) {
            return i;
        }
    }
    
    return -1;
}

long array_first_nan(array l) {
    if(l == NULL) {
        return LONG_MAX;
    }
    
    comp_ptr zi = l->zi;
    for (ulong i = 0; i < l->len; i++) {
        if(mpfr_nan_p(zi[i].x) || mpfr_nan_p(zi[i].y)) {
            return i;
        }
    }
    
    return -1;
}

polyr array_polyr(array l, prec_t prec) {
    if(l == NULL || l->len == 0) {
        return false;
    }
    
    polyr P = polyr_new(l->len - 1, prec);
    if(P == NULL) {
        return NULL;
    }
    
    bool ok = true;
    for (ulong i = 0; i < l->len && ok; i++) {
        ok &= polyr_set(P, l->zi[i].x, i);
    }
    
    if(! ok) {
        polyr_free(P);
        
        return NULL;
    }
    
    return P;
}

poly array_poly(array l, prec_t prec) {
    if(l == NULL || l->len == 0) {
        return false;
    }
    
    poly P = poly_new(l->len - 1, prec);
    if(P == NULL) {
        return NULL;
    }
    
    bool ok = true;
    for (ulong i = 0; i < l->len && ok; i++) {
        ok &= poly_set(P, l->zi + i, i);
    }
    
    if(! ok) {
        poly_free(P);
        
        return NULL;
    }
    
    return P;
}

bool array_write(array l, char *fileName, int digits, bool verbose) {
    if(fileName == NULL) {
        if(verbose) {
            printf("\n**** The file name for writing is NULL ! ****\n\n");
        }
        
        return false;
    }
    
    if(l == NULL) {
        if(verbose) {
            printf("\n**** The list of points to write is NULL ! ****\n\n");
        }
        
        return false;
    }
    
    if(l->zi == NULL) {
        if(verbose) {
            printf("\n**** The list of points to write is invalid ! ****\n\n");
        }
        
        return false;
    }
    
    if(digits < 2) {
        if(verbose) {
            printf("\n**** Cannot write numbers with %d digits ! ****\n\n", digits);
        }
        
        return false;
    }
    
    long fnan = array_first_nan(l);
    if(fnan >= 0) {
        if(verbose) {
            printf("\n**** The list to write to \"%s\" contains undefined numbers ! ****\n", fileName);
            printf("**** The first undefined number is on position %ld ! ****\n\n", fnan);
        }
        
        return false;
    }
    
    long finf = array_first_inf(l);
    if(finf >= 0) {
        if(verbose) {
            printf("\n**** The list to write to \"%s\" contains infinite numbers ! ****\n", fileName);
            printf("**** The first infinite number is on position %ld ! ****\n\n", finf);
        }
        
        return false;
    }
    
    FILE *f = fopen(fileName, "w");
    if(f == NULL) {
        if(verbose) {
            printf("\n**** Could not open \"%s\" for writing ! ****\n\n", fileName);
        }
        
        return false;
    }
    
    bool ok = true;
    comp_ptr zi = l->zi;
    char format[100];
    snprintf(format, 100, "%%.%dRg, %%.%dRg\n", digits, digits);
    
    for (ulong i = 0; i < l->len && ok; i++) {
        ok = mpfr_fprintf(f, format, zi[i].x, zi[i].y) > 3;
    }
    
    fclose(f);
    
    if(! ok && verbose) {
        printf("\n**** \"%s\" contains some conversion errors ! ****\n\n", fileName);
    }
    
    return ok;
}

bool array_append(array l, char *fileName, int digits, bool verbose) {
    if(fileName == NULL) {
        if(verbose) {
            printf("\n**** The file name for writing is NULL ! ****\n\n");
        }
        
        return false;
    }
    
    if(l == NULL) {
        if(verbose) {
            printf("\n**** The list of points to write is NULL ! ****\n\n");
        }
        
        return false;
    }
    
    if(l->zi == NULL) {
        if(verbose) {
            printf("\n**** The list of points to write is invalid ! ****\n\n");
        }
        
        return false;
    }
    
    if(digits < 2) {
        if(verbose) {
            printf("\n**** Cannot write numbers with %d digits ! ****\n\n", digits);
        }
        
        return false;
    }
    
    long fnan = array_first_nan(l);
    if(fnan >= 0) {
        if(verbose) {
            printf("\n**** The list to write to \"%s\" contains undefined numbers ! ****\n", fileName);
            printf("**** The first undefined number is on position %ld ! ****\n\n", fnan);
        }
        
        return false;
    }
    
    long finf = array_first_inf(l);
    if(finf >= 0) {
        if(verbose) {
            printf("\n**** The list to write to \"%s\" contains infinite numbers ! ****\n", fileName);
            printf("**** The first infinite number is on position %ld ! ****\n\n", finf);
        }
        
        return false;
    }
    
    FILE *f = fopen(fileName, "a");
    if(f == NULL) {
        if(verbose) {
            printf("\n**** Could not open \"%s\" for appending ! ****\n\n", fileName);
        }
        
        return false;
    }
    
    bool ok = true;
    comp_ptr zi = l->zi;
    char format[100];
    snprintf(format, 100, "%%.%dRg, %%.%dRg\n", digits, digits);
    
    for (ulong i = 0; i < l->len && ok; i++) {
        ok = mpfr_fprintf(f, format, zi[i].x, zi[i].y) > 3;
    }
    
    fclose(f);
    
    if(! ok && verbose) {
        printf("\n**** \"%s\" contains some conversion errors ! ****\n\n", fileName);
    }
    
    return ok;
}

array array_read(char *fileName, prec_t prec, bool verbose) {
    if(fileName == NULL) {
        if(verbose) {
            printf("\n**** The file name is NULL ! ****\n\n");
        }
        
        return NULL;
    }
    
    if(prec < precf) {
        if(verbose) {
            printf("\n**** The precision should be at least %d ! ****\n\n", precf);
        }
        
        return NULL;
    }
    
    FILE *f = fopen(fileName, "r");
    if(f == NULL) {
        if(verbose) {
            printf("\n**** Could not open \"%s\" for reading ! ****\n\n", fileName);
        }
        
        return NULL;
    }
    
    array l = array_new(0);
    char line[10000], *test;
    line[9999] = 0;
    comp z;
    comp_init(z, prec);
    
    while(! feof(f)) {
        test = fgets(line, 10000, f);
        if(test == NULL && feof(f)) {
            break;
        }
        
        if(test == NULL || line[9999] != 0) { // error or line too long, precision > 15000 !!!
            array_free(l);
            comp_clear(z);
            
            if(verbose) {
                printf("\n**** There is an invalid line in \"%s\" ! ****\n\n", fileName);
            }
            
            return NULL;
        }
        
        if(line[0] == '#' || line[0] == ';' || line[0] == '/') {
            continue;
        }
        
        int pv = -1;
        for (int i = 0; line[i] != 0 && pv < 0; i++) {
            pv = line[i] == ',' ? i : -1;
        }
        int end = -1;
        for (int i = pv + 1; line[i] != 0 && end < 0; i++) {
            end = line[i] == '\n' ? i : -1;
        }
        
        if(pv < 0) { // not a complex number
            continue;
        }
        
        line[pv] = 0; // split the string
        if(end > pv) {
            line[end] = 0;
        }

        bool ok = mpfr_set_str(z->x, line, 10, MPFR_RNDN) == 0;
        ok = ok && mpfr_set_str(z->y, line + (pv + 1), 10, MPFR_RNDN) == 0;
        
        if(! ok || ! array_add(l, z, prec)) {
            array_free(l);
            comp_clear(z);
            
            if(verbose) {
                printf("\n**** Could not add a number to the list, not enough memory ? ****\n\n");
            }
            
            return NULL;
        }
    }
    
    fclose(f);
    comp_clear(z);
    
    long fnan = array_first_nan(l);
    if(fnan >= 0) {
        if(verbose) {
            printf("\n**** The list to write to \"%s\" contains undefined numbers ! ****\n", fileName);
            printf("**** The first undefined number is on position %ld ! ****\n\n", fnan);
        }
        
        return NULL;
    }
    
    long finf = array_first_inf(l);
    if(finf >= 0) {
        if(verbose) {
            printf("\n**** The list to write to \"%s\" contains infinite numbers ! ****\n", fileName);
            printf("**** The first infinite number is on position %ld ! ****\n\n", finf);
        }
        
        return NULL;
    }
    
    if(verbose) {
        printf("Read %" PRIu64 " points from \"%s\".\n", l->len, fileName);
    }
    
    return l;
}

