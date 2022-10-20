//
//  appsf.c
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


#include "appsf.h"
#include "evalf.h"
#include "chrono.h"
#include "help.h"

arrayf appf_join(arrayf z1, arrayf z2) {
    ulong len = z1->len < z2->len ? z1->len : z2->len;
    arrayf r = arrayf_new(len);
    if(r == NULL) {
        printf("Out of memory !\n");
        
        return NULL;
    }
    
    compf z;
    for (ulong i = 0; i < len; i++) {
        z->x = arrayf_get(z1, i)->x;
        z->y = arrayf_get(z2, i)->x;
        arrayf_add(r, z);
    }
    
    return r;
}

arrayf appf_grid(arrayf z1, arrayf z2) {
    ulong n = z1->len;
    ulong len = n * z2->len;
    arrayf r = arrayf_new(len);
    if(r == NULL) {
        printf("Out of memory !\n");
        
        return NULL;
    }
    
    compf z;
    for (ulong i = 0; i < len; i++) {
        z->x = arrayf_get(z1, i % n)->x;
        z->y = arrayf_get(z2, i / n)->x;
        arrayf_add(r, z);
    }
    
    return r;
}

arrayf appf_exp(arrayf z) {
    arrayf ez = arrayf_new(z->size);
    compf_ptr p;
    compf ep;
    coeff_t a;
    bool ok = true;
    for (ulong i = 0; i < z->len; i++) {
        p = arrayf_get(z, i);
        a = fexp(p->x);
        ep->x = a * fcos(p->y);
        ep->y = a * fsin(p->y);
        
        ok &= arrayf_add(ez, ep);
    }
    
    if(! ok) {
        arrayf_free(ez);
        
        return NULL;
    }
    
    return ez;
}

arrayf appf_rot(arrayf z) {
    arrayf ez = arrayf_new(z->size);
    compf_ptr p;
    compf ep;
    coeff_t a;
    bool ok = true;
    for (ulong i = 0; i < z->len; i++) {
        p = arrayf_get(z, i);
        a = p->x;
        ep->x = a * fcos(p->y);
        ep->y = a * fsin(p->y);
        
        ok &= arrayf_add(ez, ep);
    }
    
    if(! ok) {
        arrayf_free(ez);
        
        return NULL;
    }
    
    return ez;
}

arrayf appf_polar(arrayf z) {
    arrayf ez = arrayf_new(z->size);
    compf_ptr p;
    compf ep;
    coeff_t a, hpi = PI / 2;
    bool ok = true;
    for (ulong i = 0; i < z->len; i++) {
        p = arrayf_get(z, i);
        a = ftan((hpi + p->x) / 2);
        ep->x = a * fcos(p->y);
        ep->y = a * fsin(p->y);
        
        ok &= arrayf_add(ez, ep);
    }
    
    if(! ok) {
        arrayf_free(ez);
        
        return NULL;
    }
    
    return ez;
}

arrayf appf_unif(ulong n, coeff_t st, coeff_t en) {
    arrayf z = arrayf_new(n);
    if(z == NULL) {
        printf("Out of memory !\n");
        
        return NULL;
    }
    
    compf ep = {{0, 0}};
    coeff_t d = en - st;
    bool ok = true;
    ulong nm1 = n - 1;
    for (ulong i = 0; i < n && ok; i++) {
        ep->x = st + i * d / nm1;
        
        ok &= arrayf_add(z, ep);
    }
    
    if(! ok) {
        arrayf_free(z);
        
        return NULL;
    }
    
    return z;
}

bool appf_compare(arrayf z1, arrayf z2, char *output) {
    FILE *out = NULL;
    if(output != NULL) {
        out = fopen(output, "w");
    }
    
    if(z1->len != z2->len) {
        ulong min = z1->len < z2->len ? z1->len : z2->len;
        
        printf("Comparing the first %" PRIu64 " points of the two lists.\n\n", min);
    }
    
    // compare values
    long pa = -1, pr = -1;
    coeff_t da = 0, dr = 0, ea, er, ma, mb;
    compf e;
    
    long i = 0;
    for (; i < z2->len && i < z1->len; i++) {
        compf_sub(e, z1->zi + i, z2->zi + i);
        ea = compf_mod(e);
        
        if(ea > da) {
            da = ea;
            pa = i;
        }
        
        ma = compf_mod(z1->zi + i);
        mb = compf_mod(z2->zi + i);
        ma = ma > mb ? ma : mb;
        
        if(ma == 0) {
            if(out != NULL) {
                fprintf(out, "0\n");
            }
            
            continue;
        }
        
        er = ea / ma;
        if(er > dr) {
            dr = er;
            pr = i;
        }
        
        if(out != NULL) {
            fprintf(out, "%.10" FMT_COEFF "g\n", er);
        }
    }
    
    if(out != NULL) {
        fclose(out);
        
        printf("Written %ld relative errors to %s.\n\n", i, output);
    }
    
    if(pa < 0) {
        printf("The two list are identical, using machine numbers.\n");
    } else {
        printf("The largest absolute error found on line %ld:\n", pa + 1);
        printf("p: %.19" FMT_COEFF "g, %.19" FMT_COEFF "g\n", z1->zi[pa].x, z1->zi[pa].y);
        printf("q: %.19" FMT_COEFF "g, %.19" FMT_COEFF "g\n", z2->zi[pa].x, z2->zi[pa].y);
        printf("|p - q| = %.10" FMT_COEFF "g\n\n", da);
        
        if(pr >= 0) {
            printf("The largest relative error found on line %ld:\n", pr + 1);
            printf("p: %.19" FMT_COEFF "g, %.19" FMT_COEFF "g\n", z1->zi[pr].x, z1->zi[pr].y);
            printf("q: %.19" FMT_COEFF "g, %.19" FMT_COEFF "g\n", z2->zi[pr].x, z2->zi[pr].y);
            printf("|p - q| / max(|p|, |q|) = %.10" FMT_COEFF "g\n\n", dr);
        }
    }
    
    return true;
}

arrayf appf_re(arrayf z) {
    arrayf ez = arrayf_new(z->size);
    if(ez == NULL) {
        printf("Out of memory !\n");
        
        return NULL;
    }
    
    compf_ptr p;
    compf ep = {{0, 0}};
    bool ok = true;
    for (ulong i = 0; i < z->len; i++) {
        p = arrayf_get(z, i);
        ep->x = p->x;
        
        ok &= arrayf_add(ez, ep);
    }
    
    if(! ok) {
        arrayf_free(ez);
        
        return NULL;
    }
    
    return ez;
}

arrayf appf_im(arrayf z) {
    arrayf ez = arrayf_new(z->size);
    if(ez == NULL) {
        printf("Out of memory !\n");
        
        return NULL;
    }
    
    compf_ptr p;
    compf ep = {{0, 0}};
    bool ok = true;
    for (ulong i = 0; i < z->len; i++) {
        p = arrayf_get(z, i);
        ep->x = p->y;
        
        ok &= arrayf_add(ez, ep);
    }
    
    if(! ok) {
        arrayf_free(ez);
        
        return NULL;
    }
    
    return ez;
}

arrayf appf_conj(arrayf z) {
    arrayf ez = arrayf_new(z->size);
    if(ez == NULL) {
        printf("Out of memory !\n");
        
        return NULL;
    }
    
    compf_ptr p;
    compf ep;
    bool ok = true;
    for (ulong i = 0; i < z->len; i++) {
        p = arrayf_get(z, i);
        ep->x = p->x;
        ep->y = -p->y;
        
        ok &= arrayf_add(ez, ep);
    }
    
    if(! ok) {
        arrayf_free(ez);
        
        return NULL;
    }
    
    return ez;
}

arrayf appf_tensor(arrayf z1, arrayf z2) {
    ulong len = z1->len < z2->len ? z1->len : z2->len;
    arrayf r = arrayf_new(len);
    if(r == NULL) {
        printf("Out of memory !\n");
        
        return NULL;
    }
    
    compf z;
    compf_ptr a, b;
    bool ok = true;
    for (ulong i = 0; i < len && ok; i++) {
        a = arrayf_get(z1, i);
        b = arrayf_get(z2, i);
        compf_mul(z, a, b);
        
        ok = ok && arrayf_add(r, z);
    }
    
    if(! ok) {
        arrayf_free(r);
        
        return NULL;
    }
    
    return r;
}

static polyfr is_realf(polyf P) {
    if(P == NULL) {
        return NULL;
    }
    
    for (deg_t i = 0; i <= P->degree; i++) {
        if(P->a[i].y != 0) {
            return NULL;
        }
    }
    
    polyfr R = polyfr_new(P->degree);
    if(R == NULL) { // out of memory
        return NULL;
    }
    
    for (ulong i = 0; i <= R->degree; i++) {
        R->a[i] = P->a[i].x;
    }
    
    return R;
}

static bool write_errf(real_t v[], arrayf pz, deg_t terms[], const char *fileName) {
    if(fileName == NULL || v == NULL || pz == NULL || pz->len == 0 || terms == NULL) {
        return false;
    }
    
    FILE *f = fopen(fileName, "w");
    if(f == NULL) {
        return false;
    }
    
    bool ok = fprintf(f, "Correct bits, log_2(abs error), # evaluated terms\n") > 1;
    for (long i = 0; i < pz->len && ok; i++) {
        real_t cb = compf_log2(pz->zi + i) - v[i]; // correct bits
        cb = cb < 1 ? 0 : cb;
        ok = fprintf(f, "%.5" FMT_REAL "g, %.5" FMT_REAL "g, %" FMT_DEG "\n",
                     cb, v[i], terms[i]) > 1;
    }
    
    fclose(f);
    
    return ok;
}

static void statsf(char *mes, deg_t deg, bool real, char *polyFile, char *ptsFile, char *outFile,
                   char *outError, double timePP, long count, double timeEV, ulong tpts,
                   char *outHorner, double timeHo, ulong tptsho, ulong tit, ulong titho) {
    
    stats_print(mes, precf, deg, TYPEFP_STR, real, polyFile, ptsFile, outFile, outError, timePP,
                count, timeEV, tpts, outHorner, timeHo, tptsho, tit, titho);
}

bool appf_eval_p(polyf P, arrayf z, char *outFile, char *outError,
                 long count, char *outHorner, long countHorner, const char **inFiles) {
    print_time_res();
    ptime ts;
    char time[100];
    
    deg_t deg, pts;
    double timePP, timeEV, timeHo = 0;
    
    deg = P->degree;
    
    polyfr R = is_realf(P);
    evalf eval;
    bool real = R != NULL;
    
    if(real) {
        printf("Read the real polynomial of degree %" FMT_DEG ".\n", R->degree);
        
        lap(&ts, NULL);
        if(count == 1) {
            eval = evalf_new_r(R);
            if(eval == NULL) {
                printf("Could not create an evaluator for this polynomial.\n");
                
                return false;
            }
            
            timePP = lap(&ts, time);
            printf("Evaluator created in %s.\n", time);
        } else {
            for (ulong i = 1; i < count; i++) {
                eval = evalf_new_r(R);
                
                evalf_free(eval);
            }
            
            eval = evalf_new_r(R);
            if(eval == NULL) {
                printf("Could not create an evaluator for this polynomial.\n");
                
                return false;
            }
            
            timePP = lap(&ts, time);
            printf("Evaluator created %ld times in %s.\n", count, time);
        }
    } else {
        printf("Read the complex polynomial of degree %" FMT_DEG ".\n", P->degree);
        
        lap(&ts, NULL);
        if(count == 1) {
            eval = evalf_new(P);
            if(eval == NULL) {
                printf("Could not create an evaluator for this polynomial.\n");
                
                return false;
            }
            
            timePP = lap(&ts, time);
            printf("Evaluator created in %s.\n", time);
        } else {
            for (ulong i = 1; i < count; i++) {
                eval = evalf_new(P);
                
                evalf_free(eval);
            }
            
            eval = evalf_new(P);
            if(eval == NULL) {
                printf("Could not create an evaluator for this polynomial.\n");
                
                return false;
            }
            
            timePP = lap(&ts, time);
            printf("Evaluator created %ld times in %s.\n", count, time);
        }
    }
    
    pts = z->len;
    
    arrayf pz = arrayf_new(pts);
    if(pz == NULL) {
        printf("Out of memory !\n");
        
        return false;
    }
    
    real_t *err = NULL;
    deg_t *terms = NULL;
    if(outError != NULL) {
        err = malloc(sizeof(real_t) * pts);
        terms = malloc(sizeof(deg_t) * pts);
    }
    
    if(count == 1) {
        printf("Evaluating the polynomial in %" PRIu64 " points ... ", pts);
    } else {
        printf("Evaluating %ld times the polynomial in %" PRIu64 " points ... ", count, pts);
    }
    
    lap(&ts, NULL);
    
    compf v;
    compf_ptr p;
    
    // order of iterations avoids eventual optimizations by caches
    for (ulong j = 0; j < count; j++) {
        for (ulong i = 0; i < pts; i++) {
            p = z->zi + i;
            
            evalf_val(v, eval, p);
            
            if(j == 0) {
                arrayf_add(pz, v);
                
                if(err != NULL) {
                    err[i] = eval->valErr;
                    terms[i] = eval->terms;
                }
            }
        }
    }
    
    timeEV = lap(&ts, time);
    printf("done in %s\n", time);
    
    bool ok = arrayf_write(pz, outFile, true);
    if(ok) {
        printf("The values are written to %s.\n", outFile);
    }
    
    if(ok && err != NULL) {
        if(write_errf(err, pz, terms, outError)) {
            printf("The error bounds in bits are written to %s.\n", outError);
        } else {
            printf("Could not write the error bounds to %s.\n", outError);
        }
    }
    
    if(ok && outHorner != NULL) {
        arrayf cz = arrayf_new(z->len);
        if(cz == NULL) {
            printf("Out of memory !\n");
            
            return false;
        }
        
        if(countHorner == 1) {
            printf("Evaluating the polynomial by Horner's method in %" PRIu64 " points ... ", z->len);
        } else {
            printf("Evaluating %ld times the polynomial by Horner's method in %"
                   FMT_DEG " points ... ", countHorner, z->len);
        }
        lap(&ts, NULL);
        
        for (ulong j = 0; j < countHorner; j++) {
            for (ulong i = 0; i < z->len; i++) {
                polyf_eval(v, P, z->zi + i);
                
                if(j == 0) {
                    arrayf_add(cz, v);
                    
                    if(err != NULL) {
                        err[i] = eval->valErr;
                    }
                }
            }
        }
        
        timeHo = lap(&ts, time);
        printf("done in %s\n", time);
        
        ok = arrayf_write(cz, outHorner, true);
        if(ok) {
            printf("The values by Horner's method are written to %s.\n", outHorner);
        }
        
        arrayf_free(cz);
    }
    
    if(err != NULL) {
        free(err);
        free(terms);
    }
    
    arrayf_free(z);
    arrayf_free(pz);
    evalf_free(eval);
    polyfr_free(R);
    polyf_free(P);
    
    if (ok) {
        statsf("Stats", deg, real, (char *) inFiles[0], (char *) inFiles[1], outFile, outError,
               timePP, count, timeEV, pts * count, outHorner, timeHo, pts * countHorner, 0, 0);
    }
    
    return ok;
}

bool appf_eval_d(polyf P, arrayf z, char *outFile, char *outError,
                 long count, char *outHorner, long countHorner, const char **inFiles) {
    print_time_res();
    ptime ts;
    char time[100];
        
    deg_t deg, pts;
    double timePP, timeEV, timeHo = 0;
    
    deg = P->degree;
    
    polyfr R = is_realf(P);
    evalf eval;
    bool real = R != NULL;
    
    if(real) {
        printf("Read the real polynomial of degree %" FMT_DEG ".\n", R->degree);
        
        lap(&ts, NULL);
        if(count == 1) {
            eval = evalf_new_r(R);
            if(eval == NULL) {
                printf("Could not create an evaluator for this polynomial.\n");
                
                return false;
            }
            
            timePP = lap(&ts, time);
            printf("Evaluator created in %s.\n", time);
        } else {
            for (ulong i = 0; i < count - 1; i++) {
                eval = evalf_new_r(R);
                
                evalf_free(eval);
            }
            
            eval = evalf_new_r(R);
            if(eval == NULL) {
                printf("Could not create an evaluator for this polynomial.\n");
                
                return false;
            }
            
            timePP = lap(&ts, time);
            printf("Evaluator created %ld times in %s.\n", count, time);
        }
    } else {
        printf("Read the complex polynomial of degree %" FMT_DEG ".\n", P->degree);
        
        lap(&ts, NULL);
        if(count == 1) {
            eval = evalf_new(P);
            if(eval == NULL) {
                printf("Could not create an evaluator for this polynomial.\n");
                
                return false;
            }
            
            timePP = lap(&ts, time);
            printf("Evaluator created in %s.\n", time);
        } else {
            for (ulong i = 0; i < count - 1; i++) {
                eval = evalf_new(P);
                
                evalf_free(eval);
            }
            
            eval = evalf_new(P);
            if(eval == NULL) {
                printf("Could not create an evaluator for this polynomial.\n");
                
                return false;
            }
            
            timePP = lap(&ts, time);
            printf("Evaluator created %ld times in %s.\n", count, time);
        }
    }
    
    pts = z->len;
    
    arrayf pz = arrayf_new(z->len);
    if(pz == NULL) {
        printf("Out of memory !\n");
        
        return false;
    }
    
    real_t *err = NULL;
    deg_t *terms = NULL;
    if(outError != NULL) {
        err = malloc(sizeof(real_t) * pts);
        terms = malloc(sizeof(deg_t) * pts);
    }
    
    if(count == 1) {
        printf("Evaluating its derivative in %" PRIu64 " points ... ", z->len);
    } else {
        printf("Evaluating %ld times its derivative in %" PRIu64 " points ... ", count, z->len);
    }
    
    lap(&ts, NULL);
    
    compf v;
    compf_ptr p;
    // order of iterations avoids eventual optimizations by caches
    for (ulong j = 0; j < count; j++) {
        for (ulong i = 0; i < z->len; i++) {
            p = z->zi + i;
            
            evalf_der(v, eval, p);
            
            if(j == 0) {
                arrayf_add(pz, v);
                
                if(err != NULL) {
                    err[i] = eval->derErr;
                    terms[i] = eval->terms;
                }
            }
        }
    }
    
    timeEV = lap(&ts, time);
    printf("done in %s\n", time);
    
    bool ok = arrayf_write(pz, outFile, true);
    if(ok) {
        printf("The values are written to %s.\n", outFile);
    }
    
    if(err != NULL) {
        if(write_errf(err, pz, terms, outError)) {
            printf("The error bounds in bits are written to %s.\n", outError);
        } else {
            printf("Could not write the error bounds to %s.\n", outError);
        }
        
        free(err);
        free(terms);
    }
    
    if(ok && outHorner != NULL) {
        arrayf cz = arrayf_new(z->len);
        if(cz == NULL) {
            printf("Out of memory !\n");
            
            return false;
        }
        
        if(countHorner == 1) {
            printf("Evaluating the polynomial by Horner's method in %" PRIu64 " points ... ", z->len);
        } else {
            printf("Evaluating %ld times the polynomial by Horner's method in %"
                   PRIu64 " points ... ", countHorner, z->len);
        }
        lap(&ts, NULL);
        
        polyf D = polyf_derivative(P);
        
        for (ulong j = 0; j < countHorner; j++) {
            for (ulong i = 0; i < z->len; i++) {
                polyf_eval(v, D, z->zi + i);
                
                if(j == 0) {
                    arrayf_add(cz, v);
                }
            }
        }
        
        timeHo = lap(&ts, time);
        printf("done in %s\n", time);
        
        ok = arrayf_write(cz, outHorner, true);
        if(ok) {
            printf("The values by Horner's method are written to %s.\n", outHorner);
        }
        
        polyf_free(D);
        arrayf_free(cz);
    }
    
    arrayf_free(z);
    arrayf_free(pz);
    evalf_free(eval);
    polyfr_free(R);
    polyf_free(P);
    
    if (ok) {
        statsf("StatsD", deg, real, (char *) inFiles[0], (char *) inFiles[1], outFile, outError,
               timePP, count, timeEV, pts * count, outHorner, timeHo, pts * countHorner, 0, 0);
    }
    
    return ok;
}

bool appf_eval_n(polyf P, arrayf z, char *outFile, char *outError,
                 long count, char *outHorner, long countHorner, const char **inFiles) {
    print_time_res();
    ptime ts;
    char time[100];
    
    deg_t deg, pts;
    double timePP, timeEV, timeHo = 0;
        
    deg = P->degree;
    
    polyfr R = is_realf(P);
    evalf eval;
    bool real = R != NULL;
    
    if(real) {
        printf("Read the real polynomial of degree %" FMT_DEG ".\n", R->degree);
        
        lap(&ts, NULL);
        if(count == 1) {
            eval = evalf_new_r(R);
            if(eval == NULL) {
                printf("Could not create an evaluator for this polynomial.\n");
                
                return false;
            }
            
            timePP = lap(&ts, time);
            printf("Evaluator created in %s.\n", time);
        } else {
            for (ulong i = 0; i < count - 1; i++) {
                eval = evalf_new_r(R);
                
                evalf_free(eval);
            }
            
            eval = evalf_new_r(R);
            if(eval == NULL) {
                printf("Could not create an evaluator for this polynomial.\n");
                
                return false;
            }
            
            timePP = lap(&ts, time);
            printf("Evaluator created %lu times in %s.\n", count, time);
        }
    } else {
        printf("Read the complex polynomial of degree %" FMT_DEG ".\n", P->degree);
        
        lap(&ts, NULL);
        if(count == 1) {
            eval = evalf_new(P);
            if(eval == NULL) {
                printf("Could not create an evaluator for this polynomial.\n");
                
                return false;
            }
            
            timePP = lap(&ts, time);
            printf("Evaluator created in %s.\n", time);
        } else {
            for (ulong i = 0; i < count - 1; i++) {
                eval = evalf_new(P);
                
                evalf_free(eval);
            }
            
            eval = evalf_new(P);
            if(eval == NULL) {
                printf("Could not create an evaluator for this polynomial.\n");
                
                return false;
            }
            
            timePP = lap(&ts, time);
            printf("Evaluator created %lu times in %s.\n", count, time);
        }
    }
    
    pts = z->len;
    
    arrayf pz = arrayf_new(z->len);
    if(pz == NULL) {
        printf("Out of memory !\n");
        
        return false;
    }
    
    real_t *err = NULL;
    deg_t *terms = NULL;
    if(outError != NULL) {
        err = malloc(sizeof(real_t) * pts);
        terms = malloc(sizeof(deg_t) * pts);
    }
    
    if(count == 1) {
        printf("Evaluating its Newton term in %" PRIu64 " points ... ", z->len);
    } else {
        printf("Evaluating %lu times its Newton term in %" PRIu64 " points ... ", count, z->len);
    }
    
    lap(&ts, NULL);
    
    compf v;
    compf_ptr p;
    // order of iterations avoids eventual optimizations by caches
    for (ulong j = 0; j < count; j++) {
        for (ulong i = 0; i < z->len; i++) {
            p = z->zi + i;
            
            evalf_newton(v, eval, p);
            
            if(j == 0) {
                arrayf_add(pz, v);
                
                if(err != NULL) {
                    err[i] = eval->ntErr;
                    terms[i] = eval->terms;
                }
            }
        }
    }
    
    timeEV = lap(&ts, time);
    printf("done in %s\n", time);
    
    bool ok = arrayf_write(pz, outFile, true);
    if(ok) {
        printf("The values are written to %s.\n", outFile);
    }
    
    if(err != NULL) {
        if(write_errf(err, pz, terms, outError)) {
            printf("The error bounds in bits are written to %s.\n", outError);
        } else {
            printf("Could not write the error bounds to %s.\n", outError);
        }
        
        free(err);
        free(terms);
    }
    
    if(ok && outHorner != NULL) {
        arrayf cz = arrayf_new(z->len);
        if(cz == NULL) {
            printf("Out of memory !\n");
            
            return false;
        }
        
        if(countHorner == 1) {
            printf("Evaluating the polynomial by Horner's method in %" PRIu64 " points ... ", z->len);
        } else {
            printf("Evaluating %ld times the polynomial by Horner's method in %"
                   PRIu64 " points ... ", countHorner, z->len);
        }
        lap(&ts, NULL);
        
        polyf D = polyf_derivative(P);
        compf d, nt;
        
        for (ulong j = 0; j < countHorner; j++) {
            for (ulong i = 0; i < z->len; i++) {
                polyf_eval(v, P, z->zi + i);
                polyf_eval(d, D, z->zi + i);
                compf_div(nt, v, d);
                
                if(j == 0) {
                    arrayf_add(cz, nt);
                }
            }
        }
        
        timeHo = lap(&ts, time);
        printf("done in %s\n", time);
        
        ok = arrayf_write(cz, outHorner, true);
        if(ok) {
            printf("The values by Horner's method are written to %s.\n", outHorner);
        }
        
        polyf_free(D);
        arrayf_free(cz);
    }
    
    arrayf_free(z);
    arrayf_free(pz);
    evalf_free(eval);
    polyfr_free(R);
    polyf_free(P);
    
    if (ok) {
        statsf("StatsN", deg, real, (char *) inFiles[0], (char *) inFiles[1], outFile, outError,
               timePP, count, timeEV, pts * count, outHorner, timeHo, pts * countHorner, 0, 0);
    }
    
    return ok;
}

static bool write_errfn(real_t v[], arrayf pz, real_t terms[], long iters[], const char *fileName) {
    if(fileName == NULL || v == NULL || pz == NULL || pz->len == 0 || terms == NULL) {
        return false;
    }
    
    FILE *f = fopen(fileName, "w");
    if(f == NULL) {
        return false;
    }
    
    bool ok = fprintf(f, "Correct bits, log_2(abs error), avg # terms, iterates\n") > 1;
    for (long i = 0; i < pz->len && ok; i++) {
        real_t cb = v[i] == INF_M ? precf : compf_log2(pz->zi + i) - v[i]; // correct bits
        cb = cb < 1 ? 0 : cb;
        ok = fprintf(f, "%.5" FMT_REAL "g, %.5" FMT_REAL "g, %.2" FMT_REAL "f, %ld\n",
                     cb, v[i], terms[i], iters[i]) > 1;
    }
    
    fclose(f);
    
    return ok;
}

bool appf_eval_r(polyf P, arrayf z, ulong maxIter, char *outFile, char *outError,
                 long count, char *outHorner, long countHorner, const char **inFiles) {
    print_time_res();
    ptime ts;
    char time[100];
    
    deg_t deg, pts;
    double timePP, timeEV, timeHo = 0;
        
    deg = P->degree;
    
    polyfr R = is_realf(P);
    evalf eval;
    bool real = R != NULL;
    
    if(real) {
        printf("Read the real polynomial of degree %" FMT_DEG ".\n", R->degree);
        
        lap(&ts, NULL);
        if(count == 1) {
            eval = evalf_new_r(R);
            if(eval == NULL) {
                printf("Could not create an evaluator for this polynomial.\n");
                
                return false;
            }
            
            timePP = lap(&ts, time);
            printf("Evaluator created in %s.\n", time);
        } else {
            for (ulong i = 0; i < count - 1; i++) {
                eval = evalf_new_r(R);
                
                evalf_free(eval);
            }
            
            eval = evalf_new_r(R);
            if(eval == NULL) {
                printf("Could not create an evaluator for this polynomial.\n");
                
                return false;
            }
            
            timePP = lap(&ts, time);
            printf("Evaluator created %lu times in %s.\n", count, time);
        }
    } else {
        printf("Read the complex polynomial of degree %" FMT_DEG ".\n", P->degree);
        
        lap(&ts, NULL);
        if(count == 1) {
            eval = evalf_new(P);
            if(eval == NULL) {
                printf("Could not create an evaluator for this polynomial.\n");
                
                return false;
            }
            
            timePP = lap(&ts, time);
            printf("Evaluator created in %s.\n", time);
        } else {
            for (ulong i = 0; i < count - 1; i++) {
                eval = evalf_new(P);
                
                evalf_free(eval);
            }
            
            eval = evalf_new(P);
            if(eval == NULL) {
                printf("Could not create an evaluator for this polynomial.\n");
                
                return false;
            }
            
            timePP = lap(&ts, time);
            printf("Evaluator created %lu times in %s.\n", count, time);
        }
    }
    
    pts = z->len;
    
    arrayf pz = arrayf_new(z->len);
    if(pz == NULL) {
        printf("Out of memory !\n");
        
        return false;
    }
    
    real_t *err = NULL;
    real_t *terms = NULL;
    long *iters = NULL;
    if(outError != NULL) {
        err = malloc(sizeof(real_t) * pts);
        terms = malloc(sizeof(real_t) * pts);
        iters = malloc(sizeof(long) * pts);
    }
    
    if(count == 1) {
        printf("Searching for roots from %" PRIu64 " starting points ... ", z->len);
    } else {
        printf("Searching %lu times for roots from %" PRIu64 " starting points ... ", count, z->len);
    }
    
    lap(&ts, NULL);
    
    compf v, p;
    ulong tit = 0;
    // order of iterations avoids eventual optimizations by caches
    for (ulong j = 0; j < count; j++) {
        for (ulong i = 0; i < z->len; i++) {
            bool conv = false, out = false;
            *p = z->zi[i];
            real_t vb = 0, sb = compf_log2(p), pb = sb, nb;
            
            long k = 0;
            real_t tr = 0;
            for (; k < maxIter && ! conv && ! out; k++) {
                evalf_newton(v, eval, p);
                vb = compf_log2(v);
                nb = eval->ntErr;
                // convergence or loss of precision
                conv = vb == INF_M || nb == INF_M || vb - nb < NEWTON_CONV_BITS;
                out = pb - sb >= NEWTON_ESCAPE_BITS; // divergence
                compf_sub(p, p, v);
                
                pb = compf_log2(p);
                
                tr += eval->terms;
            }
            
            if(j == 0) {
                arrayf_add(pz, p);
                
                if(err != NULL) {
                    err[i] = vb > eval->ntErr ? vb : eval->ntErr;
                    terms[i] = tr / k;
                    iters[i] = out ? -k : k;
                    tit += k;
                }
            }
        }
    }
    
    timeEV = lap(&ts, time);
    printf("done in %s\n", time);
    
    bool ok = arrayf_write(pz, outFile, true);
    if(ok) {
        printf("The values are written to %s.\n", outFile);
    }
    
    if(err != NULL) {
        if(write_errfn(err, pz, terms, iters, outError)) {
            printf("The error bounds in bits are written to %s.\n", outError);
        } else {
            printf("Could not write the error bounds to %s.\n", outError);
        }
        
        free(err);
        free(terms);
        free(iters);
    }
    
    ulong tith = 0;
    if(ok && outHorner) {
        arrayf cz = arrayf_new(z->len);
        if(cz == NULL) {
            printf("Out of memory !\n");
            
            return false;
        }
        
        if(countHorner == 1) {
            printf("Searching with Horner's method for roots from %" PRIu64 " starting points ... ",
                   z->len);
        } else {
            printf("Searching %ld times with Horner's method for roots from %"
                   PRIu64 " starting points ... ", countHorner, z->len);
        }
        lap(&ts, NULL);
        
        polyf D = polyf_derivative(P);
        compf s, d, nt;
        
        real_t lb = plog2(z->len);
        real_t convb = 3 * precf / 4 - lb - NEWTON_CONV_BITS;
        for (ulong j = 0; j < countHorner; j++) {
            for (ulong i = 0; i < z->len; i++) {
                *s = z->zi[i];
                bool conv = false, out = false;
                real_t vb = 0, sb = compf_log2(s), pb = sb;
                
                ulong k = 0;
                for (; k < maxIter && ! conv && ! out; k++) {
                    polyf_eval(v, P, s);
                    polyf_eval(d, D, s);
                    compf_div(nt, v, d);
                    
                    vb = compf_log2(nt);
                    conv = pb - vb >= convb || vb < -convb;
                    out = vb - sb >= NEWTON_ESCAPE_BITS;
                    compf_sub(s, s, nt);
                    
                    pb = compf_log2(s);
                }
                
                if(j == 0) {
                    arrayf_add(cz, s);
                    tith += k;
                }
            }
        }
        
        timeHo = lap(&ts, time);
        printf("done in %s\n", time);
        
        ok = arrayf_write(cz, outHorner, true);
        if(ok) {
            printf("The roots by Horner's method are written to %s.\n", outHorner);
        }
        
        polyf_free(D);
        arrayf_free(cz);
    }
    
    arrayf_free(z);
    arrayf_free(pz);
    evalf_free(eval);
    polyfr_free(R);
    polyf_free(P);
    
    if (ok) {
        statsf("StatsR", deg, real, (char *) inFiles[0], (char *) inFiles[1], outFile, outError,
               timePP, count, timeEV, pts * count, outHorner, timeHo, pts * countHorner, tit, tith);
    }
    
    return ok;
}
