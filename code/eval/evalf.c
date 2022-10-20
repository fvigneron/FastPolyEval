//
//  evalf.c
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

#include <stdlib.h>

#include "evalf.h"

// MARK: Implementation of public functions

evalf evalf_new(polyf P) {
    if(P == NULL) {
        return NULL;
    }
    
    evalf ev = malloc(sizeof(evalf_struct));
    ev->P = P;
    ev->real = false;
    ev->Q = NULL;
    ev->xn = powsfr_new(sqrt(P->degree) + 3);
    ev->zn = powsf_new(sqrt(P->degree) + 3);
    ev->valErr = NAN;
    ev->derErr = NAN;
    ev->ntErr = NAN;
    ev->terms = 0;
    ev->f = NULL;
    
    if(evalf_analyse(ev)) {
        return ev;
    }
    
    evalf_free(ev);
    
    return NULL;
}

evalf evalf_new_r(polyfr Q) {
    if(Q == NULL) {
        return NULL;
    }
    
    evalf ev = malloc(sizeof(evalf_struct));
    ev->Q = Q;
    ev->real = true;
    ev->P = NULL;
    ev->xn = powsfr_new(sqrt(Q->degree) + 3);
    ev->zn = powsf_new(sqrt(Q->degree) + 3);
    ev->valErr = NAN;
    ev->derErr = NAN;
    ev->ntErr = NAN;
    ev->terms = 0;
    ev->f = NULL;
    
    if(evalf_analyse_r(ev)) {
        return ev;
    }
    
    evalf_free(ev);
    
    return NULL;
}

bool evalf_free(evalf ev) {
    if(ev == NULL) {
        return false;
    }
    
    if(ev->f != NULL) {
        conc_free(ev->f);
    }
    
    if(ev->zn != NULL) {
        powsf_free(ev->zn);
    }
    
    if(ev->xn != NULL) {
        powsfr_free(ev->xn);
    }
    
    free(ev);
    
    return true;
}

bool evalf_val(compf v, evalf ev, compf z) {
    if(v == NULL || ev == NULL || z == NULL) {
        return false;
    }
    
    if(ev->real) { // real poly
        if(ev->Q->degree == 0) {
            v->x = ev->Q->a[0];
            v->y = 0;
            ev->valErr = INF_M;
            
            return true;
        }
        
        if(! evalf_analyse_r(ev)) {
            return false;
        }
        
        if(z->y == 0) {
            v->y = 0;
            
            if(z->x == 0) {
                v->x = ev->Q->a[0];
                ev->valErr = INF_M;
                
                return true;
            }
            
            v->x = evalf_val_rr(ev, z->x);
            
            return ! isnan(v->x);
        }
        
        return evalf_val_rc(v, ev, z);
    } else { // complex poly
        if(ev->P->degree == 0) {
            v->x = ev->P->a[0].x;
            v->y = ev->P->a[0].y;
            ev->valErr = INF_M;
            
            return true;
        }
        
        if(! evalf_analyse(ev)) {
            return false;
        }
        
        if(z->y == 0) {
            if(z->x == 0) {
                v->x = ev->P->a[0].x;
                v->y = ev->P->a[0].y;
                ev->valErr = INF_M;
                
                return true;
            }
            
            return evalf_val_cr(v, ev, z->x);
        }
        
        return evalf_val_cc(v, ev, z);
    }
}

bool evalf_der(compf d, evalf ev, compf z) {
    if(d == NULL || ev == NULL || z == NULL) {
        return false;
    }
       
    if(ev->real) { // real poly
        if(ev->Q->degree == 0) {
            d->x = 0;
            d->y = 0;
            
            return true;
        }
        
        if(! evalf_analyse_r(ev)) {
            return false;
        }
        
        if(z->y == 0) {
            d->x = evalf_der_rr(ev, z->x);
            d->y = 0;
            
            return ! isnan(d->x);
        }
        
        return evalf_der_rc(d, ev, z);
    } else { // complex poly
        if(ev->P->degree == 0) {
            d->x = 0;
            d->y = 0;
            
            return true;
        }
        
        if(! evalf_analyse(ev)) {
            return false;
        }
        
        if(z->y == 0) {
            return evalf_der_cr(d, ev, z->x);
        }
        
        return evalf_der_cc(d, ev, z);
    }
}

bool evalf_val_der(compf v, compf d, evalf ev, compf z) {
    if(v == NULL || d == NULL || ev == NULL || z == NULL) {
        return false;
    }

    if(ev->real) { // real poly
        if(ev->Q->degree == 0) {
            v->x = ev->Q->a[0];
            v->y = 0;
            d->x = 0;
            d->y = 0;
            
            return true;
        }
        
        if(! evalf_analyse_r(ev)) {
            return false;
        }
        
        if(z->y == 0) {
            v->y = 0;
            d->y = 0;
            
            return evalf_val_der_rr(&v->x, &d->x, ev, z->x);
        }
        
        return evalf_val_der_rc(v, d, ev, z);
    } else { // complex poly
        if(ev->P->degree == 0) {
            v->x = ev->P->a[0].x;
            v->y = ev->P->a[0].y;
            d->x = 0;
            d->y = 0;
            
            return true;
        }
        
        if(! evalf_analyse(ev)) {
            return false;
        }
        
        if(z->y == 0) {
            return evalf_val_der_cr(v, d, ev, z->x);
        }
        
        return evalf_val_der_cc(v, d, ev, z);
    }
}

bool evalf_newton(compf nt, evalf ev, compf z) {
    if(nt == NULL || ev == NULL || z == NULL) {
        return false;
    }

    if(ev->real) { // real poly
        if(ev->Q->degree == 0) {
            nt->x = INFINITY;
            nt->y = INFINITY;
            
            return false;
        }
        
        if(! evalf_analyse_r(ev)) {
            return false;
        }
        
        if(z->y == 0) {
            nt->x = evalf_newton_rr(ev, z->x);
            nt->y = 0;
            
            return ! isnan(nt->x);
        }
        
        return evalf_newton_rc(nt, ev, z);
    } else { // complex poly
        if(ev->P->degree == 0) {
            nt->x = INFINITY;
            nt->y = INFINITY;
            
            return false;
        }
        
        if(! evalf_analyse(ev)) {
            return false;
        }
        
        if(z->y == 0) {
            return evalf_newton_cr(nt, ev, z->x);
        }
        
        return evalf_newton_cc(nt, ev, z);
    }
}

// MARK: Bare-metal variants

bool evalf_analyse(evalf ev) {
    if(! (ev->P->modified || ev->f == NULL)) {
        return true;
    }
    
    if(ev->f != NULL) {
        conc_free(ev->f);
        
        ev->f = NULL;
    }
    
    // create the list of log2|a_k|, where a_k are the coefficients of ev->P
    real_t s;
    deg_t deg = ev->P->degree;
    compf_ptr a = ev->P->a;
    
    list l = list_new(deg + 1);
    
    bool ok = true;
    for (deg_t i = 0; i <= deg && ok; i++) {
        s = compf_s(a + i);
        ok = ok && list_add(l, i, s);
    }
    
    if(! ok) {
        list_free(l);
        
        return false;
    }
    
    // create the concave cover of the list and store it into ev
    concave f = conc_new(l, precf);
    list_free(l);
    
    ev->f = f;
    if (f != NULL) {
        ev->P->modified = false;
        
        return true;
    }
    
    return false;
}

bool evalf_analyse_r(evalf ev) {
    if(! (ev->Q->modified || ev->f == NULL)) {
        return true;
    }
    
    if(ev->f != NULL) {
        conc_free(ev->f);
        
        ev->f = NULL;
    }
    
    // create the list of log2|a_k|, where a_k are the coefficients of ev->P
    real_t s;
    deg_t deg = ev->Q->degree;
    coeff_t *a = ev->Q->a;
    
    list l = list_new(deg + 1);
    
    bool ok = true;
    for (deg_t i = 0; i <= deg && ok; i++) {
        s = coeff_s(a[i]);
        ok = ok && list_add(l, i, s);
    }
    
    if(! ok) {
        list_free(l);
        
        return false;
    }
    
    // create the concave cover of the list and store it into ev
    concave f = conc_new(l, precf);
    list_free(l);
    
    ev->f = f;
    
    if(f != NULL) {
        ev->Q->modified = false;
        
        return true;
    }
    
    return false;
}

bool evalf_val_cc(compf v, evalf ev, compf z) {
    if(z->x == 0 && z->y == 0) {
        compf_ptr a = ev->P->a;
        
        *v = a[0];
        ev->valErr = compf_log2(v) - precf;
        
        return true;
    }
    
    real_t la = compf_log2(z), err;
    deg_t pp;
    
    powsf zn = ev->zn;
    powsf_set(zn, z);
    compf_ptr zk;
    
    concave f = ev->f;
    if(! conc_range(f, la)) {
        return false;
    }
    deg_t st = f->start;
    deg_t en = f->end; // en <= st
    ev->terms = st - en + 1;
    
    // sparse Horner's scheme with the list from @c ev
    compf_ptr a = ev->P->a;
    deg_t *pows = f->all->k;
    real_t *s = f->all->s;
    
    pp = pows[st];
    deg_t p;
    compf_set(v, a + pp);
    
    err = s[st] + pp * la - precf;
    
    for (deg_t i = st - 1; i >= en && i <= MAX_DEG; i--) {
        p = pows[i];
        zk = powsf_pow(zn, pp - p);
        
        compf_mul(v, v, zk);
        compf_add(v, v, a + p);
        
        err = bits_sum(err, s[i] + p * la - precf);
        
        pp = p;
    }
    
    if(pp > 0) {
        zk = powsf_pow_once(zn, pp);
        
        compf_mul(v, v, zk);
    }
    
    ev->valErr = err;
    
    return true;
}

bool evalf_der_cc(compf d, evalf ev, compf z) {
    if(z->x == 0 && z->y == 0) {
        compf_ptr a = ev->P->a;
        
        *d = a[1];
        ev->derErr = compf_log2(d) - precf;
        
        return true;
    }
    
    real_t la = compf_log2(z), err;
    deg_t pp;
    
    powsf zn = ev->zn;
    powsf_set(zn, z);
    compf_ptr zk;
    
    concave f = ev->f;
    if(! conc_range_der(f, la)) {
        return false;
    }
    deg_t st = f->start;
    deg_t en = f->end; // en <= st
    ev->terms = st - en + 1;
    
    // sparse Horner's scheme with the list from @c ev
    compf_ptr a = ev->P->a;
    deg_t *pows = f->all->k;
    real_t *s = f->all->s;
    
    pp = pows[st];
    deg_t p;
    compf_set(d, a + pp);
    compf_mulr(d, d, pp);
    
    err = s[st] + (pp - 1) * la + plog2(pp) - precf;
    
    for (deg_t i = st - 1; i >= en && i <= MAX_DEG; i--) {
        p = pows[i];
        zk = powsf_pow(zn, pp - p);
        
        compf_mul(d, d, zk);
        compf_amr(d, a + p, p);
        
        err = bits_sum(err, s[i] + (p - 1) * la + plog2(p) - precf);
        
        pp = p;
    }
    
    if(pp > 1) {
        zk = powsf_pow_once(zn, pp - 1);
        
        compf_mul(d, d, zk);
    }
    
    ev->derErr = err;
    
    return true;
}

bool evalf_val_der_cc(compf v, compf d, evalf ev, compf z) {
    if(z->x == 0 && z->y == 0) {
        compf_ptr a = ev->P->a;
        
        *v = a[0];
        ev->valErr = compf_log2(v) - precf;
        
        *d = a[1];
        ev->derErr = compf_log2(d) - precf;
        
        return true;
    }
    
    real_t la = compf_log2(z), erv, erd;
    deg_t pp, pd;
    
    powsf zn = ev->zn;
    powsf_set(zn, z);
    compf_ptr zk;
    
    // compute d = P'(z)
    concave f = ev->f;
    if(! conc_range_der(f, la)) {
        return false;
    }
    deg_t st = f->start;
    deg_t en = f->end; // en <= st
    ev->terms = st - en + 1;
    
    // sparse Horner's scheme for P' with the list from @c ev
    compf_ptr a = ev->P->a;
    deg_t *pows = f->all->k;
    real_t *s = f->all->s;
    
    pp = pows[st];
    deg_t p;
    compf_set(d, a + pp);
    compf_mulr(d, d, pp);
    
    erd = s[st] + (pp - 1) * la + plog2(pp) - precf;
    
    for (deg_t i = st - 1; i >= en && i <= MAX_DEG; i--) {
        p = pows[i];
        zk = powsf_pow(zn, pp - p);
        
        compf_mul(d, d, zk);
        compf_amr(d, a + p, p);
        
        erd = bits_sum(erd, s[i] + (p - 1) * la + plog2(p) - precf);
        
        pp = p;
    }
    
    pd = pp;
    
    // compute v = P(z)
    if(! conc_range(f, la)) {
        return false;
    }
    st = f->start;
    en = f->end; // en <= st
    ev->terms += st - en + 1;
    
    // sparse Horner's scheme for P with the list from @c ev
    pp = pows[st];
    compf_set(v, a + pp);
    
    erv = s[st] + pp * la - precf;
    
    for (deg_t i = st - 1; i >= en && i <= MAX_DEG; i--) {
        p = pows[i];
        zk = powsf_pow(zn, pp - p);
        
        compf_mul(v, v, zk);
        compf_add(v, v, a + p);
        
        erv = bits_sum(erv, s[i] + p * la - precf);
        
        pp = p;
    }
    
    if(pd > 1) {
        zk = powsf_pow_once(zn, pd - 1);
        
        compf_mul(d, d, zk);
    }
    
    if(pp > 0) {
        zk = powsf_pow_once(zn, pp);
        
        compf_mul(v, v, zk);
    }
    
    ev->valErr = erv;
    ev->derErr = erd;
    
    return true;
}

bool evalf_newton_cc(compf nt, evalf ev, compf z) {
    if(z->x == 0 && z->y == 0) {
        compf_ptr a = ev->P->a;
        compf_div(nt, a, a + 1);
        ev->ntErr = compf_log2(nt) - precf;
        
        return true;
    }
    
    compf v, d;
    deg_t pv, pd, mp = 0;
    
    real_t la = compf_log2(z), erv, erd;
    
    powsf zn = ev->zn;
    powsf_set(zn, z);
    compf_ptr zk;
    
    concave f = ev->f;
    if(! conc_range_der(f, la)) {
        return false;
    }
    deg_t st = f->start;
    deg_t en = f->end; // en <= st
    ev->terms = st - en + 1;
    
    // sparse Horner's scheme for P' with the list from @c ev
    compf_ptr a = ev->P->a;
    deg_t *pows = f->all->k;
    real_t *s = f->all->s;
    
    deg_t pp = pows[st];
    deg_t p;
    compf_set(d, a + pp);
    compf_mulr(d, d, pp);
    
    erd = s[st] + (pp - 1) * la + plog2(pp) - precf;
    
    for (deg_t i = st - 1; i >= en && i <= MAX_DEG; i--) {
        p = pows[i];
        zk = powsf_pow(zn, pp - p);
        
        compf_mul(d, d, zk);
        compf_amr(d, a + p, p);
        
        erd = bits_sum(erd, s[i] + (p - 1) * la + plog2(p) - precf);
        
        pp = p;
    }
    
    pd = pp -1;
    
    if(! conc_range(f, la)) {
        return false;
    }
    st = f->start;
    en = f->end; // en <= st
    ev->terms += st - en + 1;
    
    // sparse Horner's scheme for P with the list from @c ev
    pp = pows[st];
    compf_set(v, a + pp);
    
    erv = s[st] + pp * la - precf;
    
    for (deg_t i = st - 1; i >= en && i <= MAX_DEG; i--) {
        p = pows[i];
        zk = powsf_pow(zn, pp - p);
        
        compf_mul(v, v, zk);
        compf_add(v, v, a + p);
        
        erv = bits_sum(erv, s[i] + p * la - precf);
        
        pp = p;
    }
    
    pv = pp;
    
    // Newton term with power simplification
    if(pv > pd) {
        zk = powsf_pow_once(zn, pv - pd);
        
        compf_mul(v, v, zk);
        mp = pd;
    } else if(pd > pv) {
        zk = powsf_pow_once(zn, pd - pv);
        
        compf_mul(d, d, zk);
        mp = pv;
    }
    
    compf_div(nt, v, d);
    
    real_t sc = mp * la;
    ev->ntErr = nt_err(compf_log2(v) + sc, erv, compf_log2(d) + sc, erd);
    
    return true;
}

bool evalf_val_rc(compf v, evalf ev, compf z) {
    if(z->x == 0 && z->y == 0) {
        coeff_t *a = ev->Q->a;
        
        v->x = a[0];
        v->y = 0;
        ev->valErr = flog2(v->x) - precf;
        
        return true;
    }
    
    real_t la = compf_log2(z), err;
    deg_t pp;
    
    powsf zn = ev->zn;
    powsf_set(zn, z);
    compf_ptr zk;
    
    concave f = ev->f;
    if(! conc_range(f, la)) {
        return false;
    }
    deg_t st = f->start;
    deg_t en = f->end; // en <= st
    ev->terms = st - en + 1;
    
    // sparse Horner's scheme with the list from @c ev
    coeff_t *a = ev->Q->a;
    deg_t *pows = f->all->k;
    real_t *s = f->all->s;
    
    pp = pows[st];
    deg_t p;
    compf_setr(v, a[pp]);
    
    err = s[st] + pp * la - precf;
    
    for (deg_t i = st - 1; i >= en && i <= MAX_DEG; i--) {
        p = pows[i];
        zk = powsf_pow(zn, pp - p);
        
        compf_mul(v, v, zk);
        v->x += a[p];
        
        err = bits_sum(err, s[i] + p * la - precf);
        
        pp = p;
    }
    
    if(pp > 0) {
        zk = powsf_pow_once(zn, pp);
        
        compf_mul(v, v, zk);
    }
    
    ev->valErr = err;
    
    return true;
}

bool evalf_der_rc(compf d, evalf ev, compf z) {
    if(z->x == 0 && z->y == 0) {
        coeff_t *a = ev->Q->a;
        
        d->x = a[1];
        d->y = 0;
        ev->derErr = flog2(d->x) - precf;
        
        return true;
    }
    
    real_t la = compf_log2(z), err;
    deg_t pp;
    
    powsf zn = ev->zn;
    powsf_set(zn, z);
    compf_ptr zk;
    
    concave f = ev->f;
    if(! conc_range_der(f, la)) {
        return false;
    }
    deg_t st = f->start;
    deg_t en = f->end; // en <= st
    ev->terms = st - en + 1;
    
    // sparse Horner's scheme with the list from @c ev
    coeff_t *a = ev->Q->a;
    deg_t *pows = f->all->k;
    real_t *s = f->all->s;
    
    pp = pows[st];
    deg_t p;
    compf_setr(d, a[pp] * pp);
    
    err = s[st] + (pp - 1) * la + plog2(pp) - precf;
    
    for (deg_t i = st - 1; i >= en && i <= MAX_DEG; i--) {
        p = pows[i];
        zk = powsf_pow(zn, pp - p);
        
        compf_mul(d, d, zk);
        d->x += a[p] * p;
        
        err = bits_sum(err, s[i] + (p - 1) * la + plog2(p) - precf);
        
        pp = p;
    }
    
    if(pp > 1) {
        zk = powsf_pow_once(zn, pp - 1);
        
        compf_mul(d, d, zk);
    }
    
    ev->derErr = err;
    
    return true;
}

bool evalf_val_der_rc(compf v, compf d, evalf ev, compf z) {
    if(z->x == 0 && z->y == 0) {
        coeff_t *a = ev->Q->a;
        
        v->x = a[0];
        v->y = 0;
        ev->valErr = flog2(v->x) - precf;
        
        d->x = a[1];
        d->y = 0;
        ev->derErr = flog2(d->x) - precf;
        
        return true;
    }
    
    real_t la = compf_log2(z), erv, erd;
    deg_t pp, pd;
    
    powsf zn = ev->zn;
    powsf_set(zn, z);
    compf_ptr zk;
    
    // compute d=Q'(z)
    concave f = ev->f;
    if(! conc_range_der(f, la)) {
        return false;
    }
    deg_t st = f->start;
    deg_t en = f->end; // en <= st
    ev->terms = st - en + 1;
    
    // sparse Horner's scheme for Q' with the list from @c ev
    coeff_t *a = ev->Q->a;
    deg_t *pows = f->all->k;
    real_t *s = f->all->s;
    
    pp = pows[st];
    deg_t p;
    compf_setr(d, a[pp] * pp);
    
    erd = s[st] + (pp - 1) * la + plog2(pp) - precf;
    
    for (deg_t i = st - 1; i >= en && i <= MAX_DEG; i--) {
        p = pows[i];
        zk = powsf_pow(zn, pp - p);
        
        compf_mul(d, d, zk);
        d->x += a[p] * p;
        
        erd = bits_sum(erd, s[i] + (p - 1) * la + plog2(p) - precf);
        
        pp = p;
    }
    
    pd = pp;
    
    // compute v=Q(z)
    if(! conc_range(f, la)) {
        return false;
    }
    st = f->start;
    en = f->end; // en <= st
    ev->terms += st - en + 1;
    
    // sparse Horner's scheme with the list from @c ev
    pows = f->all->k;
    s = f->all->s;
    
    pp = pows[st];
    compf_setr(v, a[pp]);
    
    erv = s[st] + pp * la - precf;
    
    for (deg_t i = st - 1; i >= en && i <= MAX_DEG; i--) {
        p = pows[i];
        zk = powsf_pow(zn, pp - p);
        
        compf_mul(v, v, zk);
        v->x += a[p];
        
        erv = bits_sum(erv, s[i] + p * la - precf);
        
        pp = p;
    }
    
    if(pd > 1) {
        zk = powsf_pow_once(zn, pd - 1);
        
        compf_mul(d, d, zk);
    }
    
    if(pp > 0) {
        zk = powsf_pow_once(zn, pp);
        
        compf_mul(v, v, zk);
    }
    
    ev->valErr = erv;
    ev->derErr = erd;
    
    return true;
}

bool evalf_newton_rc(compf nt, evalf ev, compf z) {
    if(z->x == 0 && z->y == 0) {
        coeff_t *a = ev->Q->a;
        nt->x = a[0] / a[1];
        nt->y = 0;
        
        ev->ntErr = flog2(nt->x) - precf;
        
        return true;
    }
    
    compf v, d;
    deg_t pv, pd, mp = 0;
    
    real_t la = compf_log2(z), erv, erd;
    
    powsf zn = ev->zn;
    powsf_set(zn, z);
    compf_ptr zk;
    
    // compute d=Q'(z)
    concave f = ev->f;
    if(! conc_range_der(f, la)) {
        return false;
    }
    deg_t st = f->start;
    deg_t en = f->end; // en <= st
    ev->terms = st - en + 1;
    
    // sparse Horner's scheme for Q' with the list from @c ev
    coeff_t *a = ev->Q->a;
    deg_t *pows = f->all->k;
    real_t *s = f->all->s;
    
    deg_t pp = pows[st];
    deg_t p;
    compf_setr(d, a[pp] * pp);
    
    erd = s[st] + (pp - 1) * la + plog2(pp) - precf;
    
    for (deg_t i = st - 1; i >= en && i <= MAX_DEG; i--) {
        p = pows[i];
        zk = powsf_pow(zn, pp - p);
        
        compf_mul(d, d, zk);
        d->x += a[p] * p;
        
        erd = bits_sum(erd, s[i] + (p - 1) * la + plog2(p) - precf);
        
        pp = p;
    }
    
    pd = pp - 1;
    
    // compute v=Q(z)
    if(! conc_range(f, la)) {
        return false;
    }
    st = f->start;
    en = f->end; // en <= st
    ev->terms += st - en + 1;
    
    // sparse Horner's scheme with the list from @c ev
    pp = pows[st];
    compf_setr(v, a[pp]);
    
    erv = s[st] + pp * la - precf;
    
    for (deg_t i = st - 1; i >= en && i <= MAX_DEG; i--) {
        p = pows[i];
        zk = powsf_pow(zn, pp - p);
        
        compf_mul(v, v, zk);
        v->x += a[p];
        
        erv = bits_sum(erv, s[i] + p * la - precf);
        
        pp = p;
    }
    
    pv = pp;
    
    // Newton term with power simplification
    if(pv > pd) {
        zk = powsf_pow_once(zn, pv - pd);
        
        compf_mul(v, v, zk);
        mp = pd;
    } else if(pd > pv) {
        zk = powsf_pow_once(zn, pd - pv);
        
        compf_mul(d, d, zk);
        mp = pv;
    }
    
    compf_div(nt, v, d);
    
    real_t sc = mp * la;
    ev->ntErr = nt_err(compf_log2(v) + sc, erv, compf_log2(d) + sc, erd);
    
    return true;
}

bool evalf_val_cr(compf v, evalf ev, coeff_t x) {
    if(x == 0) {
        compf_ptr a = ev->P->a;
        
        *v = a[0];
        ev->valErr = compf_log2(v) - precf;
        
        return true;
    }
    
    real_t la = coeff_log2(x), err;
    deg_t pp;
    
    powsfr xn = ev->xn;
    powsfr_set(xn, x);
    coeff_t xk;
    
    concave f = ev->f;
    if(! conc_range(f, la)) {
        return false;
    }
    deg_t st = f->start;
    deg_t en = f->end; // en <= st
    ev->terms = st - en + 1;
    
    // sparse Horner's scheme with the list from @c ev
    compf_ptr a = ev->P->a;
    deg_t *pows = f->all->k;
    real_t *s = f->all->s;
    
    pp = pows[st];
    deg_t p;
    compf_set(v, a + pp);
    
    err = s[st] + pp * la - precf;
    
    for (deg_t i = st - 1; i >= en && i <= MAX_DEG; i--) {
        p = pows[i];
        xk = powsfr_pow(xn, pp - p);
        
        compf_mulr(v, v, xk);
        compf_add(v, v, a + p);
        
        err = bits_sum(err, s[i] + p * la - precf);
        
        pp = p;
    }
    
    if(pp > 0) {
        xk = powsfr_pow_once(xn, pp);
        
        compf_mulr(v, v, xk);
    }
    
    ev->valErr = err;
    
    return true;
}

bool evalf_der_cr(compf d, evalf ev, coeff_t x) {
    if(x == 0) {
        compf_ptr a = ev->P->a;
        
        *d = a[1];
        ev->derErr = compf_log2(d) - precf;
        
        return true;
    }
    
    real_t la = coeff_log2(x), err;
    deg_t pp;
    
    powsfr xn = ev->xn;
    powsfr_set(xn, x);
    coeff_t xk;
    
    concave f = ev->f;
    if(! conc_range_der(f, la)) {
        return false;
    }
    deg_t st = f->start;
    deg_t en = f->end; // en <= st
    ev->terms = st - en + 1;
    
    // sparse Horner's scheme with the list from @c ev
    compf_ptr a = ev->P->a;
    deg_t *pows = f->all->k;
    real_t *s = f->all->s;
    
    pp = pows[st];
    deg_t p;
    compf_set(d, a + pp);
    compf_mulr(d, d, pp);
    
    err = s[st] + (pp - 1) * la + plog2(pp) - precf;
    
    for (deg_t i = st - 1; i >= en && i <= MAX_DEG; i--) {
        p = pows[i];
        xk = powsfr_pow(xn, pp - p);
        
        compf_mulr(d, d, xk);
        compf_amr(d, a + p, p);
        
        err = bits_sum(err, s[i] + (p - 1) * la + plog2(p) - precf);
        
        pp = p;
    }
    
    if(pp > 1) {
        xk = powsfr_pow_once(xn, pp - 1);
        
        compf_mulr(d, d, xk);
    }
    
    ev->derErr = err;
    
    return true;
}

bool evalf_val_der_cr(compf v, compf d, evalf ev, coeff_t x) {
    if(x == 0) {
        compf_ptr a = ev->P->a;
        
        *v = a[0];
        ev->valErr = compf_log2(v) - precf;
        
        *d = a[1];
        ev->derErr = compf_log2(d) - precf;
        
        return true;
    }
    
    real_t la = coeff_log2(x), erv, erd;
    deg_t pp, pd;
    
    powsfr xn = ev->xn;
    powsfr_set(xn, x);
    coeff_t xk;
    
    concave f = ev->f;
    if(! conc_range_der(f, la)) {
        return false;
    }
    deg_t st = f->start;
    deg_t en = f->end; // en <= st
    ev->terms = st - en + 1;
    
    // sparse Horner's scheme for P' with the list from @c ev
    compf_ptr a = ev->P->a;
    deg_t *pows = f->all->k;
    real_t *s = f->all->s;
    
    pp = pows[st];
    deg_t p;
    compf_set(d, a + pp);
    compf_mulr(d, d, pp);
    
    erd = s[st] + (pp - 1) * la + plog2(pp) - precf;
    
    for (deg_t i = st - 1; i >= en && i <= MAX_DEG; i--) {
        p = pows[i];
        xk = powsfr_pow(xn, pp - p);
        
        compf_mulr(d, d, xk);
        compf_amr(d, a + p, p);
        
        erd = bits_sum(erd, s[i] + (p - 1) * la + plog2(p) - precf);
        
        pp = p;
    }
    
    pd = pp;
    
    if(! conc_range(f, la)) {
        return false;
    }
    st = f->start;
    en = f->end; // en <= st
    ev->terms += st - en + 1;
    
    // sparse Horner's scheme for P with the list from @c ev
    pp = pows[st];
    compf_set(v, a + pp);
    
    erv = s[st] + pp * la - precf;
    
    for (deg_t i = st - 1; i >= en && i <= MAX_DEG; i--) {
        p = pows[i];
        xk = powsfr_pow(xn, pp - p);
        
        compf_mulr(v, v, xk);
        compf_add(v, v, a + p);
        
        erv = bits_sum(erv, s[i] + p * la - precf);
        
        pp = p;
    }
    
    if(pd > 1) {
        xk = powsfr_pow_once(xn, pd - 1);
        
        compf_mulr(d, d, xk);
    }
    
    if(pp > 0) {
        xk = powsfr_pow_once(xn, pp);
        
        compf_mulr(v, v, xk);
    }
    
    ev->valErr = erv;
    ev->derErr = erd;
    
    return true;
}

bool evalf_newton_cr(compf nt, evalf ev, coeff_t x) {
    if(x == 0) {
        compf_ptr a = ev->P->a;
        compf_div(nt, a, a + 1);
        ev->ntErr = compf_log2(nt) - precf;
        
        return true;
    }
    
    compf v, d;
    deg_t pv, pd, mp = 0;
    
    real_t la = coeff_log2(x), erv, erd;
    
    powsfr xn = ev->xn;
    powsfr_set(xn, x);
    coeff_t xk;
    
    concave f = ev->f;
    if(! conc_range_der(f, la)) {
        return false;
    }
    deg_t st = f->start;
    deg_t en = f->end; // en <= st
    ev->terms = st - en + 1;
    
    // sparse Horner's scheme for P' with the list from @c ev
    compf_ptr a = ev->P->a;
    deg_t *pows = f->all->k;
    real_t *s = f->all->s;
    
    deg_t pp = pows[st];
    deg_t p;
    compf_set(d, a + pp);
    compf_mulr(d, d, pp);
    
    erd = s[st] + (pp - 1) * la + plog2(pp) - precf;
    
    for (deg_t i = st - 1; i >= en && i <= MAX_DEG; i--) {
        p = pows[i];
        xk = powsfr_pow(xn, pp - p);
        
        compf_mulr(d, d, xk);
        compf_amr(d, a + p, p);
        
        erd = bits_sum(erd, s[i] + (p - 1) * la + plog2(p) - precf);
        
        pp = p;
    }
    
    pd = pp -1;
    
    if(! conc_range(f, la)) {
        return false;
    }
    st = f->start;
    en = f->end; // en <= st
    ev->terms += st - en + 1;
    
    // sparse Horner's scheme for P with the list from @c ev
    pp = pows[st];
    compf_set(v, a + pp);
    
    erv = s[st] + pp * la - precf;
    
    for (deg_t i = st - 1; i >= en && i <= MAX_DEG; i--) {
        p = pows[i];
        xk = powsfr_pow(xn, pp - p);
        
        compf_mulr(v, v, xk);
        compf_add(v, v, a + p);
        
        erv = bits_sum(erv, s[i] + p * la - precf);
        
        pp = p;
    }
    
    pv = pp;
    
    // Newton term with power simplification
    if(pv > pd) {
        xk = powsfr_pow_once(xn, pv - pd);
        
        compf_mulr(v, v, xk);
        mp = pd;
    } else if(pd > pv) {
        xk = powsfr_pow_once(xn, pd - pv);
        
        compf_mulr(d, d, xk);
        mp = pv;
    }
    
    compf_div(nt, v, d);
    
    real_t sc = mp * la;
    ev->ntErr = nt_err(compf_log2(v) + sc, erv, compf_log2(d) + sc, erd);
    
    return true;
}

coeff_t evalf_val_rr(evalf ev, coeff_t x) {
    if(x == 0) {
        coeff_t *a = ev->Q->a;
        
        ev->valErr = flog2(a[0]) - precf;
        
        return a[0];
    }
    
    real_t la = coeff_log2(x), err;
    deg_t pp;
    
    powsfr xn = ev->xn;
    powsfr_set(xn, x);
    coeff_t v;
    
    concave f = ev->f;
    if(! conc_range(f, la)) {
        return false;
    }
    deg_t st = f->start;
    deg_t en = f->end; // en <= st
    ev->terms = st - en + 1;
    
    // sparse Horner's scheme with the list from @c ev
    coeff_t *a = ev->Q->a;
    deg_t *pows = f->all->k;
    real_t *s = f->all->s;
    
    pp = pows[st];
    deg_t p;
    v = a[pp];
    
    err = s[st] + pp * la - precf;
    
    for (deg_t i = st - 1; i >= en && i <= MAX_DEG; i--) {
        p = pows[i];
        v *= powsfr_pow(xn, pp - p);
        v += a[p];
        
        err = bits_sum(err, s[i] + p * la - precf);
        
        pp = p;
    }
        
    if(pp > 0) {
        v *= powsfr_pow_once(xn, pp);
    }
    
    ev->valErr = err;
    
    return v;
}

coeff_t evalf_der_rr(evalf ev, coeff_t x) {
    if(x == 0) {
        coeff_t *a = ev->Q->a;
        
        ev->derErr = flog2(a[1]) - precf;
        
        return a[1];
    }
    
    real_t la = coeff_log2(x), err;
    deg_t pp;
    
    powsfr xn = ev->xn;
    powsfr_set(xn, x);
    coeff_t d;
    
    concave f = ev->f;
    if(! conc_range_der(f, la)) {
        return false;
    }
    deg_t st = f->start;
    deg_t en = f->end; // en <= st
    ev->terms = st - en + 1;
    
    // sparse Horner's scheme with the list from @c ev
    coeff_t *a = ev->Q->a;
    deg_t *pows = f->all->k;
    real_t *s = f->all->s;
    
    pp = pows[st];
    deg_t p;
    d = a[pp] * pp;
    
    err = s[st] + (pp - 1) * la + plog2(pp) - precf;
    
    for (deg_t i = st - 1; i >= en && i <= MAX_DEG; i--) {
        p = pows[i];
        d *= powsfr_pow(xn, pp - p);
        d += a[p] * p;
        
        err = bits_sum(err, s[i] + (p - 1) * la + plog2(p) - precf);
        
        pp = p;
    }

    if(pp > 1) {
        d *= powsfr_pow_once(xn, pp - 1);
    }
    
    ev->derErr = err;
    
    return d;
}

bool evalf_val_der_rr(coeff_t *v, coeff_t *d, evalf ev, coeff_t x) {
    if(x == 0) {
        coeff_t *a = ev->Q->a;
        
        *v = a[0];
        ev->valErr = flog2(*v) - precf;
        
        *d = a[1];
        ev->derErr = flog2(*d) - precf;
        
        return true;
    }
    
    real_t la = coeff_log2(x), erv, erd;
    deg_t pp, pd;
    
    powsfr xn = ev->xn;
    powsfr_set(xn, x);
    coeff_t vx, dx;
    
    // compute d=Q'(z)
    concave f = ev->f;
    if(! conc_range_der(f, la)) {
        return false;
    }
    deg_t st = f->start;
    deg_t en = f->end; // en <= st
    ev->terms = st - en + 1;
    
    // sparse Horner's scheme for Q' with the list from @c ev
    coeff_t *a = ev->Q->a;
    deg_t *pows = f->all->k;
    real_t *s = f->all->s;
    
    pp = pows[st];
    deg_t p;
    dx = a[pp] * pp;
    
    erd = s[st] + (pp - 1) * la + plog2(pp) - precf;
    
    for (deg_t i = st - 1; i >= en && i <= MAX_DEG; i--) {
        p = pows[i];
        dx *= powsfr_pow(xn, pp - p);
        dx += a[p] * p;
        
        erd = bits_sum(erd, s[i] + (p - 1) * la + plog2(p) - precf);
        
        pp = p;
    }
    
    pd = pp;
    
    // compute v=Q(z)
    if(! conc_range(f, la)) {
        return false;
    }
    st = f->start;
    en = f->end; // en <= st
    ev->terms += st - en + 1;
    
    // sparse Horner's scheme with the list from @c ev
    pp = pows[st];
    vx = a[pp];
    
    erv = s[st] + pp * la - precf;
    
    for (deg_t i = st - 1; i >= en && i <= MAX_DEG; i--) {
        p = pows[i];
        vx *= powsfr_pow(xn, pp - p);
        vx += a[p];
        
        erv = bits_sum(erv, s[i] + p * la - precf);
        
        pp = p;
    }
    
    if(pd > 1) {
        dx *= powsfr_pow_once(xn, pd - 1);
    }
    
    if(pp > 0) {
        vx *= powsfr_pow_once(xn, pp);
    }
    
    *v = vx;
    *d = dx;
    
    ev->valErr = erv;
    ev->derErr = erd;
    
    return true;
}

coeff_t evalf_newton_rr(evalf ev, coeff_t x) {
    if(x == 0) {
        coeff_t *a = ev->Q->a;
        coeff_t nt = a[0] / a[1];
        ev->ntErr = flog2(nt) - precf;
        
        return nt;
    }
    
    coeff_t v, d;
    deg_t pv, pd, mp = 0;
    
    real_t la = coeff_log2(x), erv, erd;
    
    powsfr xn = ev->xn;
    powsfr_set(xn, x);
    
    // compute d=Q'(z)
    concave f = ev->f;
    if(! conc_range_der(f, la)) {
        return false;
    }
    deg_t st = f->start;
    deg_t en = f->end; // en <= st
    ev->terms = st - en + 1;
    
    // sparse Horner's scheme for Q' with the list from @c ev
    coeff_t *a = ev->Q->a;
    deg_t *pows = f->all->k;
    real_t *s = f->all->s;
    
    deg_t pp = pows[st];
    deg_t p;
    d = a[pp] * pp;
    
    erd = s[st] + (pp - 1) * la + plog2(pp) - precf;
    
    for (deg_t i = st - 1; i >= en && i <= MAX_DEG; i--) {
        p = pows[i];
        d *= powsfr_pow(xn, pp - p);
        d += a[p] * p;
        
        erd = bits_sum(erd, s[i] + (p - 1) * la + plog2(p) - precf);
        
        pp = p;
    }
    
    pd = pp - 1;
    
    // compute v=Q(z)
    if(! conc_range(f, la)) {
        return false;
    }
    st = f->start;
    en = f->end; // en <= st
    ev->terms += st - en + 1;
    
    // sparse Horner's scheme with the list from @c ev
    pp = pows[st];
    v = a[pp];
    
    erv = s[st] + pp * la - precf;
    
    for (deg_t i = st - 1; i >= en && i <= MAX_DEG; i--) {
        p = pows[i];
        v *= powsfr_pow(xn, pp - p);
        v += a[p];
        
        erv = bits_sum(erv, s[i] + p * la - precf);
        
        pp = p;
    }
    
    pv = pp;
    
    // Newton term with power simplification
    if(pv > pd) {
        v *= powsfr_pow_once(xn, pv - pd);
        mp = pd;
    } else if(pd > pv) {
        d *= powsfr_pow_once(xn, pd - pv);
        mp = pv;
    }
    
    real_t sc = mp * la;
    ev->ntErr = nt_err(coeff_log2(v) + sc, erv, coeff_log2(d) + sc, erd);
    
    return v / d;
}
