//
//  arrayf.c
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

#include <mpfr.h>
#include <inttypes.h>

#include "arrayf.h"

arrayf arrayf_new(ulong size) {
    ulong max = UINT64_MAX / sizeof(compf_struct);
    if(size > max) {
        return NULL;
    }
    
    arrayf l = malloc(sizeof(arrayf_struct));
    if(l == NULL) {
        return NULL;
    }
    
    l->len = 0;
    l->size = size < LISTF_MIN_SIZE ? LISTF_MIN_SIZE : size;
    
    l->zi = malloc(sizeof(compf_struct) * l->size);
    if(l->zi == NULL) {
        free(l);
        
        return NULL;
    }
    
    return l;
}

arrayf arrayf_new_polyfr(polyfr P) {
    if(P == NULL || P->degree >= UINT64_MAX) {
        return NULL;
    }
    
    ulong n = P->degree + 1;
    arrayf l = arrayf_new(n);
    if(l == NULL) {
        return NULL;
    }
    
    for (ulong i = 0; i < n; i++) {
        l->zi[i].x = P->a[i];
        l->zi[i].y = 0;
    }
    
    l->len = n;
    
    return l;
}

arrayf arrayf_new_polyf(polyf P) {
    if(P == NULL || P->degree >= UINT64_MAX) {
        return NULL;
    }
    
    ulong n = P->degree + 1;
    arrayf l = arrayf_new(n);
    if(l == NULL) {
        return NULL;
    }
    
    for (ulong i = 0; i < n; i++) {
        l->zi[i] = P->a[i];
    }
    
    l->len = n;
    
    return l;
}

bool arrayf_free(arrayf l) {
    if(l == NULL || l->zi == NULL) {
        return false;
    }
    
    free(l->zi);
    l->zi = NULL;
    
    free(l);
    
    return true;
}

bool arrayf_add(arrayf l, compf z) {
    if(l == NULL || z == NULL) {
        return false;
    }
    
    if(l->len == l->size) {
        l->size *= LISTF_SIZE_INCREASE;
        
        ulong size = sizeof(compf_struct) * l->size;
        l->zi = realloc(l->zi, size);
        
        if(l->zi == NULL) {
            return false;
        }
    }
    
    l->zi[l->len ++] = *z;
    
    return true;
}

compf_ptr arrayf_get(arrayf l, ulong pos) {
    if(l == NULL || pos >= l->len) {
        return NULL;
    }
    
    return l->zi + pos;
}

bool arrayf_is_real(arrayf l) {
    if(l == NULL) {
        return false;
    }
    
    bool real = true;
    compf_ptr zi = l->zi;
    for (ulong i = 0; i < l->len && real; i++) {
        real &= zi[i].y == 0;
    }
    
    return real;
}

long arrayf_first_inf(arrayf l) {
    if(l == NULL) {
        return LONG_MAX;
    }
    
    compf_ptr zi = l->zi;
    for (ulong i = 0; i < l->len; i++) {
        if(isinf(zi[i].x) || isinf(zi[i].y)) {
            return i;
        }
    }
    
    return -1;
}

long arrayf_first_nan(arrayf l) {
    if(l == NULL) {
        return LONG_MAX;
    }
    
    compf_ptr zi = l->zi;
    for (ulong i = 0; i < l->len; i++) {
        if(isnan(zi[i].x) || isnan(zi[i].y)) {
            return i;
        }
    }
    
    return -1;
}

polyfr arrayf_polyfr(arrayf l) {
    if(l == NULL || l->len == 0) {
        return false;
    }
    
    polyfr P = polyfr_new(l->len - 1);
    if(P == NULL) {
        return NULL;
    }
    
    for (ulong i = 0; i < l->len; i++) {
        P->a[i] = l->zi[i].x;
    }
    
    return P;
}

polyf arrayf_polyf(arrayf l) {
    if(l == NULL || l->len == 0) {
        return false;
    }
    
    polyf P = polyf_new(l->len - 1);
    if(P == NULL) {
        return NULL;
    }
    
    for (ulong i = 0; i < l->len; i++) {
        P->a[i] = l->zi[i];
    }
    
    return P;
}

bool arrayf_write(arrayf l, char *fileName, bool verbose) {
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
    
    long fnan = arrayf_first_nan(l);
    if(fnan >= 0) {
        if(verbose) {
            printf("\n**** The list to write to \"%s\" contains undefined numbers ! ****\n", fileName);
            printf("**** The first undefined number is on position %ld ! ****\n", fnan);
            printf("**** Please try again using precision at least %d ! ****\n\n", precf + 1);
        }
        
        return false;
    }
    
    long finf = arrayf_first_inf(l);
    if(finf >= 0) {
        if(verbose) {
            printf("\n**** The list to write to \"%s\" contains infinite numbers ! ****\n", fileName);
            printf("**** The first infinite number is on position %ld ! ****\n", finf);
            printf("**** Please try again using precision at least %d ! ****\n\n", precf + 1);
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
    compf_ptr zi = l->zi;
    for (ulong i = 0; i < l->len && ok; i++) {
        ok = fprintf(f, "%.19" FMT_COEFF "g, %.19" FMT_COEFF "g\n", zi[i].x, zi[i].y) > 3;
    }
    
    fclose(f);
    
    if(! ok && verbose) {
        printf("\n**** \"%s\" contains some conversion errors ! ****\n\n", fileName);
    }
    
    return ok;
}

bool arrayf_append(arrayf l, char *fileName, bool verbose) {
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
    
    long fnan = arrayf_first_nan(l);
    if(fnan >= 0) {
        if(verbose) {
            printf("\n**** The list to write to \"%s\" contains undefined numbers ! ****\n", fileName);
            printf("**** The first undefined number is on position %ld ! ****\n", fnan);
            printf("**** Please try again using precision at least %d ! ****\n\n", precf + 1);
        }
        
        return false;
    }
    
    long finf = arrayf_first_inf(l);
    if(finf >= 0) {
        if(verbose) {
            printf("\n**** The list to write to \"%s\" contains infinite numbers ! ****\n", fileName);
            printf("**** The first infinite number is on position %ld ! ****\n", finf);
            printf("**** Please try again using precision at least %d ! ****\n\n", precf + 1);
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
    compf_ptr zi = l->zi;
    for (ulong i = 0; i < l->len && ok; i++) {
        ok = fprintf(f, "%.19" FMT_COEFF "g, %.19" FMT_COEFF "g\n", zi[i].x, zi[i].y) > 3;
    }
    
    fclose(f);
    
    if(! ok && verbose) {
        printf("\n**** \"%s\" contains some conversion errors ! ****\n\n", fileName);
    }
    
    return ok;
}

arrayf arrayf_read(char *fileName, bool verbose) {
    if(fileName == NULL) {
        if(verbose) {
            printf("\n**** The file name is NULL ! ****\n\n");
        }
        
        return NULL;
    }
    
    FILE *f = fopen(fileName, "r");
    if(f == NULL) {
        if(verbose) {
            printf("\n**** Could not open \"%s\" ! ****\n\n", fileName);
        }
        
        return NULL;
    }
    
    arrayf l = arrayf_new(0);
    char line[10000], *test;
    line[9999] = 0;
    compf z;
    coeff_t x = 0, y = 0;
    
    while(! feof(f)) {
        test = fgets(line, 10000, f);
        if(test == NULL && feof(f)) {
            break;
        }
        
        if(test == NULL || line[9999] != 0) { // error or line too long, precision > 15000 !!!
            arrayf_free(l);
            
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
        
        if(pv < 0) { // not a complex number
            continue;
        }
        
        line[pv] = ' '; // split the string

        if(sscanf(line, "%" FMT_COEFF "g%" FMT_COEFF "g", &x, &y) < 2) {
            arrayf_free(l);
            
            if(verbose) {
                printf("\n**** There is an invalid line in \"%s\" ! ****\n\n", fileName);
            }

            return NULL;
        }
        
        z->x = x;
        z->y = y;
        
        if(! arrayf_add(l, z)) {
            arrayf_free(l);
            
            if(verbose) {
                printf("\n**** Could not add a number to the list, not enough memory ? ****\n\n");
            }
                       
            return NULL;
        }
    }
    
    fclose(f);
    
    long fnan = arrayf_first_nan(l);
    if(fnan >= 0) {
        if(verbose) {
            printf("\n**** The list to write to \"%s\" contains undefined numbers ! ****\n", fileName);
            printf("**** The first undefined number is on position %ld ! ****\n\n", fnan);
        }
        
        return NULL;
    }
    
    long finf = arrayf_first_inf(l);
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
