//
//  apps.c
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

#include "ntypes.h"
#include "apps.h"
#include "eval.h"
#include "chrono.h"
#include "help.h"

array app_join(long prec, array z1, array z2) {
    ulong len = z1->len < z2->len ? z1->len : z2->len;
    array r = array_new(len);
    if(r == NULL) {
        printf("Out of memory !\n");
        
        return NULL;
    }
    
    comp z;
    comp_init(z, prec);
    bool ok = true;
    for (ulong i = 0; i < len && ok; i++) {
        mpfr_set(z->x, array_get(z1, i)->x, MPFR_RNDN);
        mpfr_set(z->y, array_get(z2, i)->x, MPFR_RNDN);
        
        ok &= array_add(r, z, prec);
    }
    
    comp_clear(z);
    
    if(! ok) {
        array_free(r);
        
        return NULL;
    }
    
    return r;
}

array app_grid(long prec, array z1, array z2) {
    ulong n = z1->len;
    ulong len = n * z2->len;
    array r = array_new(len);
    if(r == NULL) {
        printf("Out of memory !\n");
        
        return NULL;
    }
    
    comp z;
    comp_init(z, prec);
    bool ok = true;
    for (ulong i = 0; i < len && ok; i++) {
        mpfr_set(z->x, array_get(z1, i % n)->x, MPFR_RNDN);
        mpfr_set(z->y, array_get(z2, i / n)->x, MPFR_RNDN);
        
        ok &= array_add(r, z, prec);
    }
    
    comp_clear(z);
    
    if(! ok) {
        array_free(r);
        
        return NULL;
    }
    
    return r;
}

array app_exp(long prec, array z) {
    array ez = array_new(z->size);
    comp_ptr p;
    comp ep;
    comp_init(ep, prec);
    
    mpfr_t a;
    mpfr_init2(a, prec);
    bool ok = true;
    for (ulong i = 0; i < z->len; i++) {
        p = array_get(z, i);
        mpfr_exp(a, p->x, MPFR_RNDN);
        
        mpfr_cos(ep->x, p->y, MPFR_RNDN);
        mpfr_mul(ep->x, ep->x, a, MPFR_RNDN);
        
        mpfr_sin(ep->y, p->y, MPFR_RNDN);
        mpfr_mul(ep->y, ep->y, a, MPFR_RNDN);
        
        ok &= array_add(ez, ep, prec);
    }
    
    mpfr_clear(a);
    comp_clear(ep);
    
    if(! ok) {
        array_free(ez);
        
        return NULL;
    }
    
    return ez;
}

array app_rot(long prec, array z) {
    array ez = array_new(z->size);
    comp_ptr p;
    comp ep;
    comp_init(ep, prec);
    
    bool ok = true;
    for (ulong i = 0; i < z->len; i++) {
        p = array_get(z, i);
        
        mpfr_cos(ep->x, p->y, MPFR_RNDN);
        mpfr_mul(ep->x, ep->x, p->x, MPFR_RNDN);
        
        mpfr_sin(ep->y, p->y, MPFR_RNDN);
        mpfr_mul(ep->y, ep->y, p->x, MPFR_RNDN);
        
        ok &= array_add(ez, ep, prec);
    }
    
    comp_clear(ep);
    
    if(! ok) {
        array_free(ez);
        
        return NULL;
    }
    
    return ez;
}

array app_polar(long prec, array z) {
    array ez = array_new(z->size);
    comp_ptr p;
    comp ep;
    comp_init(ep, prec);
    
    mpfr_t a, hpi;
    mpfr_init2(a, prec);
    mpfr_init2(hpi, prec);
    mpfr_const_pi(hpi, MPFR_RNDN);
    mpfr_mul_2si(hpi, hpi, -1, MPFR_RNDN);
    
    bool ok = true;
    for (ulong i = 0; i < z->len; i++) {
        p = array_get(z, i);
        mpfr_add(a, hpi, p->x, MPFR_RNDN);
        mpfr_mul_2si(a, a, -1, MPFR_RNDN);
        mpfr_tan(a, a, MPFR_RNDN);
        
        mpfr_cos(ep->x, p->y, MPFR_RNDN);
        mpfr_mul(ep->x, ep->x, a, MPFR_RNDN);
        
        mpfr_sin(ep->y, p->y, MPFR_RNDN);
        mpfr_mul(ep->y, ep->y, a, MPFR_RNDN);
        
        ok &= array_add(ez, ep, prec);
    }
    
    mpfr_clear(a);
    comp_clear(ep);
    
    if(! ok) {
        array_free(ez);
        
        return NULL;
    }
    
    return ez;
}

array app_unif(long prec, ulong n, mpfr_t st, mpfr_t en) {    
    array ez = array_new(n);
    if(ez == NULL) {
        printf("Out of memory !\n");
        
        return NULL;
    }
    
    comp ep;
    comp_init(ep, prec);
    mpfr_set_zero(ep->y, 1);
    
    mpfr_t d;
    mpfr_init2(d, prec);
    mpfr_sub(d, en, st, MPFR_RNDN);
    
    bool ok = true;
    ulong nm1 = n - 1;
    for (ulong i = 0; i < n && ok; i++) {
        mpfr_mul_ui(ep->x, d, i, MPFR_RNDN);
        mpfr_div_ui(ep->x, ep->x, nm1, MPFR_RNDN);
        mpfr_add(ep->x, ep->x, st, MPFR_RNDN);
        
        ok &= array_add(ez, ep, prec);
    }
    
    mpfr_clear(d);
    comp_clear(ep);
    
    if(! ok) {
        array_free(ez);
        
        return NULL;
    }
    
    return ez;
}

array app_rand(long prec, ulong n, mpfr_t st, mpfr_t en) {
    array ez = array_new(n);
    if(ez == NULL) {
        printf("Out of memory !\n");
        
        return NULL;
    }
    
    comp ep;
    comp_init(ep, prec);
    mpfr_set_zero(ep->y, 1);
    
    mpfr_t d;
    mpfr_init2(d, prec);
    mpfr_sub(d, en, st, MPFR_RNDN);
    
    gmp_randstate_t rnd;
    gmp_randinit_default(rnd);
    
    ptime now;
    clock_gettime(USED_CLOCK, &now);
    gmp_randseed_ui(rnd, now.tv_nsec);
    
    bool ok = true;
    for (ulong i = 0; i < n && ok; i++) {
        mpfr_urandom(ep->x, rnd, MPFR_RNDN);
        mpfr_mul(ep->x, ep->x, d, MPFR_RNDN);
        mpfr_add(ep->x, ep->x, st, MPFR_RNDN);
        
        ok &= array_add(ez, ep, prec);
    }
    
    gmp_randclear(rnd);
    
    comp_clear(ep);
    mpfr_clear(d);
    
    if(! ok) {
        array_free(ez);
        
        return NULL;
    }
    
    return ez;
}

array app_normal(long prec, ulong n, mpfr_t cent, mpfr_t var) {
    array ez = array_new(n);
    if(ez == NULL) {
        printf("Out of memory !\n");
        
        return NULL;
    }
    
    comp ep;
    comp_init(ep, prec);
    mpfr_set_zero(ep->y, 1);
    
    gmp_randstate_t rnd;
    gmp_randinit_default(rnd);
    
    ptime now;
    clock_gettime(USED_CLOCK, &now);
    gmp_randseed_ui(rnd, now.tv_nsec);
    
    bool ok = true;
    for (ulong i = 0; i < n && ok; i++) {
        mpfr_nrandom(ep->x, rnd, MPFR_RNDN);
        mpfr_mul(ep->x, ep->x, var, MPFR_RNDN);
        mpfr_add(ep->x, ep->x, cent, MPFR_RNDN);
        
        ok &= array_add(ez, ep, prec);
    }
    
    gmp_randclear(rnd);
    
    comp_clear(ep);
    
    if(! ok) {
        array_free(ez);
        
        return NULL;
    }
    
    return ez;
}

array app_sphere(long prec, ulong n) {
    // at least the machine float precision
    prec = prec < precf ? precf : prec;

    // pi, pi/2, 2pi, 2pi/n, pi/n, parallel angle
    mpfr_t pi, hpi, tpi, ang, hang, pang;
    mpfr_init2(pi, prec);
    mpfr_const_pi(pi, MPFR_RNDN);
    
    mpfr_init2(hpi, prec);
    mpfr_mul_2si(hpi, pi, -1, MPFR_RNDN);
    
    mpfr_init2(tpi, prec);
    mpfr_mul_2si(tpi, pi, 1, MPFR_RNDN);
    
    mpfr_init2(ang, prec);
    mpfr_mul_2si(ang, pi, 1, MPFR_RNDN);
    mpfr_div_ui(ang, ang, n, MPFR_RNDN);
    
    mpfr_init2(hang, prec);
    mpfr_mul_2si(hang, ang, -1, MPFR_RNDN);
    
    mpfr_init2(pang, prec);
    
    array z = array_new(n);
    if(z == NULL) {
        printf("Out of memory !\n");
        
        return NULL;
    }
        
    comp ep;
    comp_init(ep, prec);
    bool ok = true;
    
    // equator
    mpfr_set_zero(ep->x, 1);
    for (ulong i = 0; i < n && ok; i++) {
        mpfr_mul_ui(ep->y, ang, i, MPFR_RNDN);
        
        ok &= array_add(z, ep, prec);
    }
    
    // other parallels
    ulong vert = (n - 1) / 4;
    for (ulong j = 1; j <= vert && ok; j++) {
        mpfr_mul_ui(ep->x, ang, j, MPFR_RNDN);
        
        mpfr_cos(pang, ep->x, MPFR_RNDN);
        mpfr_div(pang, ang, pang, MPFR_RNDN); // enlarged angle for the j-th parallel
        
        mpfr_div(pang, tpi, pang, MPFR_RNDN);
        ulong count = mpfr_get_ui(pang, MPFR_RNDN); // the number of arcs on the parallel
        mpfr_div_ui(pang, tpi, count, MPFR_RNDN);   // the actual angle on the parallel
        
        for (ulong i = 0; i < count; i++) {
            mpfr_mul_ui(ep->y, pang, i, MPFR_RNDN);
            
            // north
            ok &= array_add(z, ep, prec);
            mpfr_neg(ep->x, ep->x, MPFR_RNDN);
            
            // south
            ok &= array_add(z, ep, prec);
            mpfr_neg(ep->x, ep->x, MPFR_RNDN);
        }
    }
    
    comp_clear(ep);
    mpfr_clears(pang, hang, ang, hpi, tpi, pi, NULL);
    
    if(! ok) {
        array_free(z);
        
        return NULL;
    }
    
    return z;
}

bool app_compare(long prec, array z1, array z2, char *output) {
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
    mpfr_t da, dr, ea, er, ma, mb;
    mpfr_init2(da, prec);
    mpfr_set_zero(da, 1);
    mpfr_init2(dr, prec);
    mpfr_set_zero(dr, 1);
    mpfr_init2(ea, prec);
    mpfr_init2(er, prec);
    mpfr_init2(ma, prec);
    mpfr_init2(mb, prec);
    
    comp e;
    comp_init(e, prec);
    
    long i = 0;
    for (; i < z2->len && i < z1->len; i++) {
        comp_sub(e, z1->zi + i, z2->zi + i);
        comp_mod(ea, e);
        
        if(mpfr_cmp(ea, da) > 0) {
            mpfr_set(da, ea, MPFR_RNDN);
            pa = i;
        }
        
        comp_mod(ma, z1->zi + i);
        comp_mod(mb, z2->zi + i);
        mpfr_max(ma, ma, mb, MPFR_RNDN);
        
        if(mpfr_zero_p(ma)) {
            if(out != NULL) {
                fprintf(out, "0\n");
            }
            
            continue;
        }
        
        mpfr_div(er, ea, ma, MPFR_RNDN);
        if(mpfr_cmp(er, dr) > 0) {
            mpfr_set(dr, er, MPFR_RNDN);
            pr = i;
        }
        
        if(out != NULL) {
            mpfr_fprintf(out, "%.10Rg\n", er);
        }
    }
    
    if(out != NULL) {
        fclose(out);
        
        printf("Written %ld relative errors to %s.\n\n", i, output);
    }
    
    if(pa < 0) {
        printf("The two list are identical up to %ld bits of precision.\n", prec);
    } else {
        int dig = 1 + prec / 3.2;
        char fp[40], fq[40];
        snprintf(fp, 40, "p: %%.%dRg, %%.%dRg\n", dig, dig);
        snprintf(fq, 40, "q: %%.%dRg, %%.%dRg\n", dig, dig);
        
        printf("The largest absolute error found on line %ld:\n", pa + 1);
        mpfr_printf(fp, z1->zi[pa].x, z1->zi[pa].y);
        mpfr_printf(fq, z2->zi[pa].x, z2->zi[pa].y);
        mpfr_printf("|p - q| = %.10Rg\n\n", da);
        
        if(pr >= 0) {
            printf("The largest relative error found on line %ld:\n", pr + 1);
            mpfr_printf(fp, z1->zi[pr].x, z1->zi[pr].y);
            mpfr_printf(fq, z2->zi[pr].x, z2->zi[pr].y);
            mpfr_printf("|p - q| / max(|p|, |q|) = %.10Rg\n\n", dr);
        }
    }
    
    comp_clear(e);
    
    mpfr_clears(da, dr, ea, er, ma, mb, NULL);
    
    return true;
}

array app_re(long prec, array z) {
    array ez = array_new(z->size);
    if(ez == NULL) {
        printf("Out of memory !\n");
        
        return NULL;
    }
    
    comp_ptr p;
    comp ep;
    comp_initz(ep, prec);
    
    bool ok = true;
    for (ulong i = 0; i < z->len; i++) {
        p = array_get(z, i);
        mpfr_set(ep->x, p->x, MPFR_RNDN);
        
        ok &= array_add(ez, ep, prec);
    }
    
    comp_clear(ep);
    
    if(! ok) {
        array_free(ez);
        
        return NULL;
    }
    
    return ez;
}

array app_im(long prec, array z) {
    array ez = array_new(z->size);
    if(ez == NULL) {
        printf("Out of memory !\n");
        
        return NULL;
    }
    
    comp_ptr p;
    comp ep;
    comp_initz(ep, prec);
    
    bool ok = true;
    for (ulong i = 0; i < z->len; i++) {
        p = array_get(z, i);
        mpfr_set(ep->x, p->y, MPFR_RNDN);
        
        ok &= array_add(ez, ep, prec);
    }
    
    comp_clear(ep);
    
    if(! ok) {
        array_free(ez);
        
        return NULL;
    }
    
    return ez;
}

array app_conj(long prec, array z) {
    array ez = array_new(z->size);
    if(ez == NULL) {
        printf("Out of memory !\n");
        
        return NULL;
    }
    
    comp_ptr p;
    comp ep;
    comp_init(ep, prec);
    
    bool ok = true;
    for (ulong i = 0; i < z->len; i++) {
        p = array_get(z, i);
        mpfr_set(ep->x, p->x, MPFR_RNDN);
        mpfr_set(ep->y, p->y, MPFR_RNDN);
        mpfr_neg(ep->y, ep->y,  MPFR_RNDN);
        
        ok &= array_add(ez, ep, prec);
    }
    
    comp_clear(ep);
    
    if(! ok) {
        array_free(ez);
        
        return NULL;
    }
    
    return ez;
}

array app_tensor(long prec, array z1, array z2) {
    ulong len = z1->len < z2->len ? z1->len : z2->len;
    array r = array_new(len);
    if(r == NULL) {
        printf("Out of memory !\n");
        
        return NULL;
    }
    
    mpfr_t b1, b2;
    mpfr_init2(b1, prec);
    mpfr_init2(b2, prec);
    
    comp z;
    comp_init(z, prec);
    comp_ptr a ,b;
    bool ok = true;
    for (ulong i = 0; i < len && ok; i++) {
        a = array_get(z1, i);
        b = array_get(z2, i);
        comp_mul(z, a, b, b1, b2);
        
        ok &= array_add(r, z, prec);
    }
    
    comp_clear(z);
    mpfr_clears(b1, b2, NULL);
    
    if(! ok) {
        array_free(r);
        
        return NULL;
    }
    
    return r;
}

static polyr is_real(poly P) {
    if(P == NULL) {
        return NULL;
    }
    
    for (deg_t i = 0; i <= P->degree; i++) {
        if(! mpfr_zero_p(P->a[i].y)) {
            return NULL;
        }
    }
    
    polyr R = polyr_new(P->degree, P->prec);
    if(R == NULL) { // out of memory
        return NULL;
    }
    
    for (ulong i = 0; i <= R->degree; i++) {
        polyr_set(R, P->a[i].x, i);
    }
    
    return R;
}

static bool write_err(real_t v[], array pz, deg_t terms[], const char *fileName) {
    if(fileName == NULL || v == NULL || pz == NULL || pz->len == 0 || terms == NULL) {
        return false;
    }
    
    FILE *f = fopen(fileName, "w");
    if(f == NULL) {
        return false;
    }
    
    bool ok = fprintf(f, "Correct bits, log_2(abs error), # evaluated terms\n") > 1;
    for (long i = 0; i < pz->len && ok; i++) {
        real_t cb = comp_log2(pz->zi + i) - v[i]; // correct bits
        cb = cb < 1 ? 0 : cb;
        ok = fprintf(f, "%.5" FMT_REAL "g, %.5" FMT_REAL "g, %" FMT_DEG "\n",
                     cb, v[i], terms[i]) > 1;
    }
    
    fclose(f);
    
    return ok;
}

static void statsmp(char *mes, long prec, deg_t deg, bool real, char *polyFile, char *ptsFile, char *outFile,
                   char *outError, double timePP, long count, double timeEV, ulong tpts,
                   char *outHorner, double timeHo, ulong tptsho, ulong tit, ulong titho) {
    
    stats_print(mes, prec, deg, "mpfr", real, polyFile, ptsFile, outFile, outError, timePP,
                count, timeEV, tpts, outHorner, timeHo, tptsho, tit, titho);
}

bool app_eval_p(long prec, poly P, array z, char *outFile, char *outError,
                long count, char *outHorner, long countHorner, const char **inFiles) {
    print_time_res();
    ptime ts;
    char time[100];
    
    deg_t deg, pts;
    double timePP, timeEV, timeHo = 0;
        
    deg = P->degree;
    
    polyr R = is_real(P);
    eval eval;
    bool real = R != NULL;
    
    if(real) {
        printf("Read the real polynomial of degree %" FMT_DEG ".\n", R->degree);
        
        lap(&ts, NULL);
        if(count == 1) {
            eval = eval_new_r(R, prec);
            if(eval == NULL) {
                printf("Could not create an evaluator for this polynomial.\n");
                
                return false;
            }
            
            timePP = lap(&ts, time);
            printf("Evaluator created in %s.\n", time);
        } else {
            for (ulong i = 0; i < count - 1; i++) {
                eval = eval_new_r(R, prec);
                
                eval_free(eval);
            }
            
            eval = eval_new_r(R, prec);
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
            eval = eval_new(P, prec);
            if(eval == NULL) {
                printf("Could not create an evaluator for this polynomial.\n");
                
                return false;
            }
            
            timePP = lap(&ts, time);
            printf("Evaluator created in %s.\n", time);
        } else {
            for (ulong i = 0; i < count - 1; i++) {
                eval = eval_new(P, prec);
                
                eval_free(eval);
            }
            
            eval = eval_new(P, prec);
            if(eval == NULL) {
                printf("Could not create an evaluator for this polynomial.\n");
                
                return false;
            }
            
            timePP = lap(&ts, time);
            printf("Evaluator created %ld times in %s.\n", count, time);
        }
    }
    
    pts = z->len;
    
    array pz = array_new(pts);
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
    
    lap(&ts, NULL);
    if(count == 1) {
        printf("Evaluating the polynomial in %" PRIu64 " points ... ", pts);
    } else {
        printf("Evaluating %ld times the polynomial in %" PRIu64 " points ... ", count, pts);
    }
    
    comp v;
    comp_init(v, prec);
    
    comp_ptr p;
    // order of iterations avoids eventual optimizations by caches
    for (ulong j = 0; j < count; j++) {
        for (ulong i = 0; i < pts; i++) {
            p = z->zi + i;
            
            eval_val(v, eval, p);
            
            if(j == 0) {
                array_add(pz, v, prec);
                
                if(err != NULL) {
                    err[i] = eval->valErr;
                    terms[i] = eval->terms;
                }
            }
        }
    }

    timeEV = lap(&ts, time);
    printf("done in %s\n", time);
    
    bool ok = array_write(pz, outFile, 1 + prec / 3.2, true);
    if(ok) {
        printf("The values are written to %s.\n", outFile);
    }
    
    if(ok && err != NULL) {
        if(write_err(err, pz, terms, outError)) {
            printf("The error bounds in bits are written to %s.\n", outError);
        } else {
            printf("Could not write the error bounds to %s.\n", outError);
        }
        
        free(err);
        free(terms);
    }
    
    if(ok && outHorner != NULL) {
        array cz = array_new(z->len);
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
        
        for (ulong j = 0; j < countHorner; j++) {
            for (ulong i = 0; i < z->len; i++) {
                poly_eval(v, P, z->zi + i);
                
                if(j == 0) {
                    array_add(cz, v, prec);
                }
            }
        }
        
        timeHo = lap(&ts, time);
        printf("done in %s\n", time);
        
        ok = array_write(cz, outHorner, 1 + prec / 3.2, true);
        if(ok) {
            printf("The values by Horner's method are written to %s.\n", outHorner);
        }
        
        array_free(cz);
    }
    
    array_free(z);
    array_free(pz);
    eval_free(eval);
    polyr_free(R);
    poly_free(P);
    
    comp_clear(v);
    
    if (ok) {
        statsmp("Stats", prec, deg, real, (char *) inFiles[0], (char *) inFiles[1], outFile, outError,
               timePP, count, timeEV, pts * count, outHorner, timeHo, pts * countHorner, 0, 0);
    }
    
    return ok;
}

bool app_eval_d(long prec, poly P, array z, char *outFile, char *outError,
                long count, char *outHorner, long countHorner, const char **inFiles) {
    print_time_res();
    ptime ts;
    char time[100];
    
    deg_t deg, pts;
    double timePP, timeEV, timeHo = 0;
        
    deg = P->degree;
    
    polyr R = is_real(P);
    eval eval;
    bool real = R != NULL;
    
    if(real) {
        printf("Read the real polynomial of degree %" FMT_DEG ".\n", R->degree);
        
        lap(&ts, NULL);
        if(count == 1) {
            eval = eval_new_r(R, prec);
            if(eval == NULL) {
                printf("Could not create an evaluator for this polynomial.\n");
                
                return false;
            }
            
            timePP = lap(&ts, time);
            printf("Evaluator created in %s.\n", time);
        } else {
            for (ulong i = 0; i < count - 1; i++) {
                eval = eval_new_r(R, prec);
                
                eval_free(eval);
            }
            
            eval = eval_new_r(R, prec);
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
            eval = eval_new(P, prec);
            if(eval == NULL) {
                printf("Could not create an evaluator for this polynomial.\n");
                
                return false;
            }
            
            timePP = lap(&ts, time);
            printf("Evaluator created in %s.\n", time);
        } else {
            for (ulong i = 0; i < count - 1; i++) {
                eval = eval_new(P, prec);
                
                eval_free(eval);
            }
            
            eval = eval_new(P, prec);
            if(eval == NULL) {
                printf("Could not create an evaluator for this polynomial.\n");
                
                return false;
            }
            
            timePP = lap(&ts, time);
            printf("Evaluator created %ld times in %s.\n", count, time);
        }
    }
    
    pts = z->len;
    
    array pz = array_new(z->len);
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
    
    comp v;
    comp_init(v, prec);
    
    comp_ptr p;
    // order of iterations avoids eventual optimizations by caches
    for (ulong j = 0; j < count; j++) {
        for (ulong i = 0; i < z->len; i++) {
            p = z->zi + i;
            
            eval_der(v, eval, p);
            
            if(j == 0) {
                array_add(pz, v, prec);
                
                if(err != NULL) {
                    err[i] = eval->derErr;
                    terms[i] = eval->terms;
                }
            }
        }
    }

    timeEV = lap(&ts, time);
    printf("done in %s\n", time);
    
    bool ok = array_write(pz, outFile, 1 + prec / 3.2, true);
    if(ok) {
        printf("The values are written to %s.\n", outFile);
    }
    
    if(ok && err != NULL) {
        if(write_err(err, pz, terms, outError)) {
            printf("The error bounds in bits are written to %s.\n", outError);
        } else {
            printf("Could not write the error bounds to %s.\n", outError);
        }
        
        free(err);
        free(terms);
    }
    
    if(ok && outHorner != NULL) {
        array cz = array_new(z->len);
        if(cz == NULL) {
            printf("Out of memory !\n");
            
            return false;
        }
        
        if(countHorner == 1) {
            printf("Evaluating the derivative by Horner's method in %" PRIu64 " points ... ", z->len);
        } else {
            printf("Evaluating %ld times the derivative by Horner's method in %"
                   FMT_DEG " points ... ", countHorner, z->len);
        }
        lap(&ts, NULL);
        
        poly D = poly_derivative(P);
        
        for (ulong j = 0; j < countHorner; j++) {
            for (ulong i = 0; i < z->len; i++) {
                poly_eval(v, D, z->zi + i);
                
                if(j == 0) {
                    array_add(cz, v, prec);
                }
            }
        }
        
        timeHo = lap(&ts, time);
        printf("done in %s\n", time);
        
        ok = array_write(cz, outHorner, 1 + prec / 3.2, true);
        if(ok) {
            printf("The values by Horner's method are written to %s.\n", outHorner);
        }
        
        poly_free(D);
        array_free(cz);
    }
    
    array_free(z);
    array_free(pz);
    eval_free(eval);
    polyr_free(R);
    poly_free(P);
    
    comp_clear(v);
    
    if (ok) {
        statsmp("StatsD", prec, deg, real, (char *) inFiles[0], (char *) inFiles[1], outFile, outError,
               timePP, count, timeEV, pts * count, outHorner, timeHo, pts * countHorner, 0, 0);
    }
    
    return ok;
}

bool app_eval_n(long prec, poly P, array z, char *outFile, char *outError,
                long count, char *outHorner, long countHorner, const char **inFiles) {
    print_time_res();
    ptime ts;
    char time[100];
    
    deg_t deg, pts;
    double timePP, timeEV, timeHo = 0;
    
    deg = P->degree;
    
    polyr R = is_real(P);
    eval eval;
    bool real = R != NULL;
    
    if(real) {
        printf("Read the real polynomial of degree %" FMT_DEG ".\n", R->degree);
        
        lap(&ts, NULL);
        if(count == 1) {
            eval = eval_new_r(R, prec);
            if(eval == NULL) {
                printf("Could not create an evaluator for this polynomial.\n");
                
                return false;
            }
            
            timePP = lap(&ts, time);
            printf("Evaluator created in %s.\n", time);
        } else {
            for (ulong i = 0; i < count - 1; i++) {
                eval = eval_new_r(R, prec);
                
                eval_free(eval);
            }
            
            eval = eval_new_r(R, prec);
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
            eval = eval_new(P, prec);
            if(eval == NULL) {
                printf("Could not create an evaluator for this polynomial.\n");
                
                return false;
            }
            
            timePP = lap(&ts, time);
            printf("Evaluator created in %s.\n", time);
        } else {
            for (ulong i = 0; i < count - 1; i++) {
                eval = eval_new(P, prec);
                
                eval_free(eval);
            }
            
            eval = eval_new(P, prec);
            if(eval == NULL) {
                printf("Could not create an evaluator for this polynomial.\n");
                
                return false;
            }
            
            timePP = lap(&ts, time);
            printf("Evaluator created %lu times in %s.\n", count, time);
        }
    }
    
    pts = z->len;
    
    array pz = array_new(z->len);
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
    
    comp v;
    comp_init(v, prec);
    
    comp_ptr p;
    // order of iterations avoids eventual optimizations by caches
    for (ulong j = 0; j < count; j++) {
        for (ulong i = 0; i < z->len; i++) {
            p = z->zi + i;
            
            eval_newton(v, eval, p);
            
            if(j == 0) {
                array_add(pz, v, prec);
                
                if(err != NULL) {
                    err[i] = eval->ntErr;
                    terms[i] = eval->terms;
                }
            }
        }
    }

    timeEV = lap(&ts, time);
    printf("done in %s\n", time);
    
    bool ok = array_write(pz, outFile, 1 + prec / 3.2, true);
    if(ok) {
        printf("The values are written to %s.\n", outFile);
    }
    
    if(ok && err != NULL) {
        if(write_err(err, pz, terms, outError)) {
            printf("The error bounds in bits are written to %s.\n", outError);
        } else {
            printf("Could not write the error bounds to %s.\n", outError);
        }
        
        free(err);
        free(terms);
    }
    
    if(ok && outHorner != NULL) {
        array cz = array_new(z->len);
        if(cz == NULL) {
            printf("Out of memory !\n");
            
            return false;
        }
        
        if(countHorner == 1) {
            printf("Evaluating the derivative by Horner's method in %" PRIu64 " points ... ", z->len);
        } else {
            printf("Evaluating %ld times the derivative by Horner's method in %"
                   PRIu64 " points ... ", countHorner, z->len);
        }
        lap(&ts, NULL);
        
        poly D = poly_derivative(P);
        comp d, nt;
        comp_init(d, prec);
        comp_init(nt, prec);
        mpfr_t a, b, c;
        mpfr_init2(a, prec);
        mpfr_init2(b, prec);
        mpfr_init2(c, prec);
        
        for (ulong j = 0; j < countHorner; j++) {
            for (ulong i = 0; i < z->len; i++) {
                poly_eval(v, P, z->zi + i);
                poly_eval(d, D, z->zi + i);
                comp_div(nt, v, d, a, b, c);
                
                if(j == 0) {
                    array_add(cz, nt, prec);
                }
            }
        }
        
        timeHo = lap(&ts, time);
        printf("done in %s\n", time);
        
        ok = array_write(cz, outHorner, 1 + prec / 3.2, true);
        if(ok) {
            printf("The values by Horner's method are written to %s.\n", outHorner);
        }
        
        poly_free(D);
        array_free(cz);
        comp_clear(d);
        comp_clear(nt);
        mpfr_clear(a);
        mpfr_clear(b);
        mpfr_clear(c);
    }
    
    array_free(z);
    array_free(pz);
    eval_free(eval);
    polyr_free(R);
    poly_free(P);
    
    comp_clear(v);
    
    if (ok) {
        statsmp("StatsN", prec, deg, real, (char *) inFiles[0], (char *) inFiles[1], outFile, outError,
               timePP, count, timeEV, pts * count, outHorner, timeHo, pts * countHorner, 0, 0);
    }
    
    return ok;
}

static bool write_errn(long prec, real_t v[], array pz, real_t terms[], long iters[],
                       const char *fileName) {
    if(fileName == NULL || v == NULL || pz == NULL || pz->len == 0 || terms == NULL) {
        return false;
    }
    
    FILE *f = fopen(fileName, "w");
    if(f == NULL) {
        return false;
    }
    
    bool ok = fprintf(f, "Correct bits, log_2(abs error), avg # terms, iterates\n") > 1;
    for (long i = 0; i < pz->len && ok; i++) {
        real_t cb = v[i] == INF_M ? prec : comp_log2(pz->zi + i) - v[i]; // correct bits
        cb = cb < 1 ? 0 : cb;
        ok = fprintf(f, "%.5" FMT_REAL "g, %.5" FMT_REAL "g, %.2" FMT_REAL "f, %ld\n",
                     cb, v[i], terms[i], iters[i]) > 1;
    }
    
    fclose(f);
    
    return ok;
}

bool app_eval_r(long prec, poly P, array z, ulong maxIter, char *outFile, char *outError,
                 long count, char *outHorner, long countHorner, const char **inFiles) {
    print_time_res();
    ptime ts;
    char time[100];
    
    deg_t deg, pts;
    double timePP, timeEV, timeHo = 0;
        
    deg = P->degree;
    
    polyr R = is_real(P);
    eval eval;
    bool real = R != NULL;
    
    if(real) {
        printf("Read the real polynomial of degree %" FMT_DEG ".\n", R->degree);
        
        lap(&ts, NULL);
        if(count == 1) {
            eval = eval_new_r(R, prec);
            if(eval == NULL) {
                printf("Could not create an evaluator for this polynomial.\n");
                
                return false;
            }
            
            timePP = lap(&ts, time);
            printf("Evaluator created in %s.\n", time);
        } else {
            for (ulong i = 0; i < count - 1; i++) {
                eval = eval_new_r(R, prec);
                
                eval_free(eval);
            }
            
            eval = eval_new_r(R, prec);
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
            eval = eval_new(P, prec);
            if(eval == NULL) {
                printf("Could not create an evaluator for this polynomial.\n");
                
                return false;
            }
            
            timePP = lap(&ts, time);
            printf("Evaluator created in %s.\n", time);
        } else {
            for (ulong i = 0; i < count - 1; i++) {
                eval = eval_new(P, prec);
                
                eval_free(eval);
            }
            
            eval = eval_new(P, prec);
            if(eval == NULL) {
                printf("Could not create an evaluator for this polynomial.\n");
                
                return false;
            }
            
            timePP = lap(&ts, time);
            printf("Evaluator created %lu times in %s.\n", count, time);
        }
    }
    
    pts = z->len;
    
    array pz = array_new(z->len);
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
    
    comp v, s;
    comp_init(v, prec);
    comp_init(s, prec);
    
    ulong tit = 0;
    // order of iterations avoids eventual optimizations by caches
    for (ulong j = 0; j < count; j++) {
        for (ulong i = 0; i < z->len; i++) {
            comp_set(s, z->zi + i);
            bool conv = false, out = false;
            real_t vb = 0, sb = comp_log2(s), pb = sb, nb;
            
            long k = 0;
            real_t tr = 0;
            for (; k < maxIter && ! conv && ! out; k++) {
                eval_newton(v, eval, s);
                vb = comp_log2(v);
                nb = eval->ntErr;
                // convergence or loss of precision
                conv = vb == INF_M || nb == INF_M || vb - nb < NEWTON_CONV_BITS;
                out = pb - sb >= NEWTON_ESCAPE_BITS; // divergence
                
                comp_sub(s, s, v);
                pb = comp_log2(s);
                
                tr += eval->terms;
            }
            
            if(j == 0) {
                array_add(pz, s, prec);
                
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
    
    bool ok = array_write(pz, outFile, 1 + prec / 3.2, true);
    if(ok) {
        printf("The values are written to %s.\n", outFile);
    }
    
    if(ok && err != NULL) {
        if(write_errn(prec, err, pz, terms, iters, outError)) {
            printf("The error bounds in bits are written to %s.\n", outError);
        } else {
            printf("Could not write the error bounds to %s.\n", outError);
        }
        
        free(err);
        free(terms);
        free(iters);
    }
    
    ulong tith = 0;
    if(ok && outHorner != NULL) {
        array cz = array_new(z->len);
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
        
        poly D = poly_derivative(P);
        comp d, nt;
        comp_init(d, prec);
        comp_init(nt, prec);
        mpfr_t a, b, c;
        mpfr_init2(a, prec);
        mpfr_init2(b, prec);
        mpfr_init2(c, prec);
        
        real_t lb = plog2(z->len);
        real_t convb = 3 * prec / 4 - lb - NEWTON_CONV_BITS;
        for (ulong j = 0; j < countHorner; j++) {
            for (ulong i = 0; i < z->len; i++) {
                comp_set(s, z->zi + i);
                bool conv = false, out = false;
                real_t vb = 0, sb = comp_log2(s), pb = sb;
                
                ulong k = 0;
                for (; k < maxIter && ! conv && ! out; k++) {
                    poly_eval(v, P, s);
                    poly_eval(d, D, s);
                    comp_div(nt, v, d, a, b, c);
                    
                    vb = comp_log2(nt);
                    conv = pb - vb >= convb || vb < -convb;
                    out = vb - sb >= NEWTON_ESCAPE_BITS;
                    comp_sub(s, s, nt);
                    
                    pb = comp_log2(s);
                }
                
                if(j == 0) {
                    array_add(cz, s, prec);
                    tith += k;
                }
            }
        }
        
        timeHo = lap(&ts, time);
        printf("done in %s\n", time);
        
        ok = array_write(cz, outHorner, 1 + prec / 3.2, true);
        if(ok) {
            printf("The roots by Horner's method are written to %s.\n", outHorner);
        }
        
        poly_free(D);
        array_free(cz);
        comp_clear(d);
        comp_clear(nt);
        mpfr_clear(a);
        mpfr_clear(b);
        mpfr_clear(c);
    }
    
    array_free(z);
    array_free(pz);
    eval_free(eval);
    polyr_free(R);
    poly_free(P);
    
    comp_clear(v);
    comp_clear(s);
    
    if (ok) {
        statsmp("StatsR", prec, deg, real, (char *) inFiles[0], (char *) inFiles[1], outFile, outError,
               timePP, count, timeEV, pts * count, outHorner, timeHo, pts * countHorner, tit, tith);
    }
    
    return ok;
}

static bool write_concave(poly P, concave cc, const char *fileName) {
    if(fileName == NULL || P == NULL || cc == NULL) {
        return false;
    }
    
    FILE *f = fopen(fileName, "w");
    if(f == NULL) {
        return false;
    }
    
    bool ok = fprintf(f, "Degree k, log_2 |a_k|, s(a_k), concave, always ignored\n") > 1;
    char *is_conc, *is_ignored;
    long kc = 0, ki = 0;
    long dc = cc->def->k[kc];
    long di = cc->all->k[kc];
    
    for (long i = 0; i <= P->degree && ok; i++) {
        real_t cb = comp_log2(P->a + i);
        long sk = comp_s(P->a + i);
        
        is_conc = "false";
        if(i == dc) {
            is_conc = "true";
            kc ++;
            dc = kc < cc->def->count ? cc->def->k[kc] : -1;
        }
        
        is_ignored = "true";
        if(i == di) {
            is_ignored = "false";
            ki ++;
            di = ki < cc->all->count ? cc->all->k[ki] : -1;
        }
        
        ok = fprintf(f, "%ld, %.5" FMT_REAL "g, %ld, %s, %s\n", i, cb, sk, is_conc, is_ignored) > 1;
    }
    
    fclose(f);
    
    return ok;
}

static bool write_intervals(list il, concave cc, const char *fileName) {
    if(fileName == NULL || il == NULL || cc == NULL) {
        return false;
    }
    
    FILE *f = fopen(fileName, "w");
    if(f == NULL) {
        return false;
    }
    
    bool ok = fprintf(f, "min(log_2 |z|), max(log_2 |z|), Lowest deg, Highest deg, Deg k of max |a_k z^k|, Terms\n") > 1;
    
    real_t l, r, m;
    deg_t ld, hd, md, t;
    for (long i = il->count; i >= 0 && ok; i--) {
        l = i == il->count ? -HUGE_VALL : il->s[i];
        r = i == 0 ? HUGE_VALL : il->s[i - 1];
        m = i == il->count ? r - 1 : i == 0 ? l + 1 : (l + r) / 2;
        
        conc_range(cc, m);
        ld = cc->all->k[cc->end];
        hd = cc->all->k[cc->start];
        md = cc->def->k[cc->mid];
        t = 1 + cc->start - cc->end;
        
        ok = fprintf(f, "%.20" FMT_REAL "g, %.20" FMT_REAL "g, %" FMT_DEG ", %" FMT_DEG
                     ", %" FMT_DEG ", %" FMT_DEG "\n", l, r, ld, hd, md, t) > 1;
    }
    
    fclose(f);
    
    return ok;
}

bool add_change_points(list chg, concave cc, real_t min, real_t max) {
    if(min >= max) {
        return false;
    }
    
    conc_range(cc, min);
    real_t lmin = cc->end;
    real_t rmin = cc->start;
    
    conc_range(cc, max);
    real_t lmax = cc->end;
    real_t rmax = cc->start;
    
    if(lmin == lmax && rmin == rmax) {
        return true;
    }
    
    real_t mid = (min + max) / 2;
    bool ok = true;
    if(mid == min || mid == max) {
        ok = list_add(chg, 0, mid);
        
        return ok;
    }
    
    ok = ok && add_change_points(chg, cc, min, mid);
    ok = ok && add_change_points(chg, cc, mid, max);
    
    return ok;
}

list find_change_points(concave cc) {
    if(cc == NULL || cc->all->count < 2) {
        printf("There is only one non-zero coefficient, cannot find intervals !\n");
        
        return NULL;
    }
    
    deg_t last = cc->all->count - 1;
    
    list chg = list_new(2 * last + 2);
    real_t min = 0, max = 0, inc = 1;
    
    conc_range(cc, min);
    while (min != -HUGE_VALL && cc->start > 0) {
        min -= inc;
        inc *= 2;
        
        conc_range(cc, min);
    }
    
    if(min == -HUGE_VALL) {
        printf("The concave cover is not consistent for low degrees !\n");
        
        return NULL;
    }
    
    inc = 1;
    conc_range(cc, max);
    while (max != HUGE_VALL && cc->end < last) {
        max += inc;
        inc *= 2;
        
        conc_range(cc, max);
    }
    
    if(max == HUGE_VALL) {
        printf("The concave cover is not consistent for high degrees !\n");
        
        return NULL;
    }
    
    if(! add_change_points(chg, cc, min, max)) {
        list_free(chg);
        
        return NULL;
    }
    
    return chg;
}

bool app_analyse(long prec, poly P, char *outFile, char *outChange) {
    eval eval;
    eval = eval_new(P, prec);
    if(eval == NULL) {
        printf("Could not create an evaluator for the polynomial with %ld bits of precision.\n",
               prec);
        
        return false;
    }
    
    if(! write_concave(P, eval->f, outFile)) {
        printf("Could not write the concave cover definition to %s !\n", outFile);
        
        return false;
    }
    
    printf("Written the concave cover definition to %s\n", outFile);
    
    if(outChange != NULL) {
        list chg_pts = find_change_points(eval->f);
        if(chg_pts == NULL) {
            printf("Could not find the meaningful intervals for this polynomial !\n");
            
            return false;
        }
        
        list_sort(chg_pts);
        
        if(! write_intervals(chg_pts, eval->f, outChange)) {
            printf("Could not write the intervals list to %s !\n", outChange);
            
            return false;
        }
        
        printf("Written the the intervals list to %s\n", outChange);
        
        list_free(chg_pts);
    }
    
    eval_free(eval);
    
    return true;
}
