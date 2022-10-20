//
//  eval.c
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

#include "eval.h"

// MARK: Implementation of public functions

eval eval_new(poly P, prec_t prec) {
    if(P == NULL || prec < MIN_EVAL_PREC) {
        return NULL;
    }
    
    eval ev = malloc(sizeof(eval_struct));
    ev->P = P;
    ev->real = false;
    ev->Q = NULL;
    ev->xn = powsr_new(prec, sqrt(P->degree) + 3);
    ev->zn = pows_new(prec, sqrt(P->degree) + 3);
    ev->valErr = NAN;
    ev->derErr = NAN;
    ev->ntErr = NAN;
    ev->terms = 0;
    ev->f = NULL;
    ev->prec = prec;
    
    comp_init(ev->buf, prec);
    mpfr_init2(ev->br, prec);
    
    if(eval_analyse(ev)) {
        return ev;
    }
    
    eval_free(ev);
    
    return NULL;
}

eval eval_new_r(polyr Q, prec_t prec) {
    if(Q == NULL || prec < MIN_EVAL_PREC) {
        return NULL;
    }
    
    eval ev = malloc(sizeof(eval_struct));
    ev->Q = Q;
    ev->real = true;
    ev->P = NULL;
    ev->xn = powsr_new(prec, sqrt(Q->degree) + 3);
    ev->zn = pows_new(prec, sqrt(Q->degree) + 3);
    ev->valErr = NAN;
    ev->derErr = NAN;
    ev->ntErr = NAN;
    ev->terms = 0;
    ev->f = NULL;
    ev->prec = prec;
    
    comp_init(ev->buf, prec);
    mpfr_init2(ev->br, prec);
    
    if(eval_analyse_r(ev)) {
        return ev;
    }
    
    eval_free(ev);
    
    return NULL;
}

bool eval_free(eval ev) {
    if(ev == NULL) {
        return false;
    }
    
    if(ev->f != NULL) {
        conc_free(ev->f);
    }
    
    if(ev->zn != NULL) {
        pows_free(ev->zn);
    }
    
    if(ev->xn != NULL) {
        powsr_free(ev->xn);
    }
    
    comp_clear(ev->buf);
    mpfr_clear(ev->br);
    
    free(ev);
    
    return true;
}

bool eval_val(comp v, eval ev, comp z) {
    if(v == NULL || ev == NULL || z == NULL) {
        return false;
    }
    
    if(ev->real) { // real poly
        if(ev->Q->degree == 0) {
            mpfr_set(v->x, ev->Q->a, MPFR_RNDN);
            mpfr_set_zero(v->y, 1);
            ev->valErr = INF_M;
            
            return true;
        }
        
        if(! eval_analyse_r(ev)) {
            return false;
        }
        
        if(mpfr_zero_p(z->y)) {
            mpfr_set_zero(v->y, 1);
            
            if(mpfr_zero_p(z->x)) {
                mpfr_set(v->x, ev->Q->a, MPFR_RNDN);
                ev->valErr = INF_M;
                
                return true;
            }
            
            return eval_val_rr(v->x, ev, z->x);
        }
        
        return eval_val_rc(v, ev, z);
    } else { // complex poly
        if(ev->P->degree == 0) {
            mpfr_set(v->x, ev->P->a[0].x, MPFR_RNDN);
            mpfr_set(v->y, ev->P->a[0].y, MPFR_RNDN);
            ev->valErr = INF_M;
            
            return true;
        }
        
        if(! eval_analyse(ev)) {
            return false;
        }
        
        if(mpfr_zero_p(z->y)) {
            if(mpfr_zero_p(z->x)) {
                mpfr_set(v->x, ev->P->a[0].x, MPFR_RNDN);
                mpfr_set(v->y, ev->P->a[0].y, MPFR_RNDN);
                ev->valErr = INF_M;
                
                return true;
            }
            
            return eval_val_cr(v, ev, z->x);
        }
        
        return eval_val_cc(v, ev, z);
    }
}

bool eval_der(comp d, eval ev, comp z) {
    if(d == NULL || ev == NULL || z == NULL) {
        return false;
    }
       
    if(ev->real) { // real poly
        if(ev->Q->degree == 0) {
            mpfr_set_zero(d->x, 1);
            mpfr_set_zero(d->y, 1);
            
            return true;
        }
        
        if(! eval_analyse_r(ev)) {
            return false;
        }
        
        if(mpfr_zero_p(z->y)) {
            mpfr_set_zero(d->y, 1);
            
            return eval_der_rr(d->x, ev, z->x);
        }
        
        return eval_der_rc(d, ev, z);
    } else { // complex poly
        if(ev->P->degree == 0) {
            mpfr_set_zero(d->x, 1);
            mpfr_set_zero(d->y, 1);
            
            return true;
        }
        
        if(! eval_analyse(ev)) {
            return false;
        }
        
        if(mpfr_zero_p(z->y)) {
            return eval_der_cr(d, ev, z->x);
        }
        
        return eval_der_cc(d, ev, z);
    }
}

bool eval_val_der(comp v, comp d, eval ev, comp z) {
    if(v == NULL || d == NULL || ev == NULL || z == NULL) {
        return false;
    }

    if(ev->real) { // real poly
        if(ev->Q->degree == 0) {
            mpfr_set(v->x, ev->Q->a, MPFR_RNDN);
            mpfr_set_zero(v->y, 1);
            mpfr_set_zero(d->x, 1);
            mpfr_set_zero(d->y, 1);
            
            return true;
        }
        
        if(! eval_analyse_r(ev)) {
            return false;
        }
        
        if(mpfr_zero_p(z->y)) {
            mpfr_set_zero(v->y, 1);
            mpfr_set_zero(d->y, 1);
            
            return eval_val_der_rr(v->x, d->x, ev, z->x);
        }
        
        return eval_val_der_rc(v, d, ev, z);
    } else { // complex poly
        if(ev->P->degree == 0) {
            mpfr_set(v->x, ev->P->a[0].x, MPFR_RNDN);
            mpfr_set(v->y, ev->P->a[0].y, MPFR_RNDN);
            mpfr_set_zero(d->x, 1);
            mpfr_set_zero(d->y, 1);
            
            return true;
        }
        
        if(! eval_analyse(ev)) {
            return false;
        }
        
        if(mpfr_zero_p(z->y)) {
            return eval_val_der_cr(v, d, ev, z->x);
        }
        
        return eval_val_der_cc(v, d, ev, z);
    }
}

bool eval_newton(comp nt, eval ev, comp z) {
    if(nt == NULL || ev == NULL || z == NULL) {
        return false;
    }

    if(ev->real) { // real poly
        if(ev->Q->degree == 0) {
            mpfr_set_inf(nt->x, 1);
            mpfr_set_inf(nt->y, 1);
            
            return false;
        }
        
        if(! eval_analyse_r(ev)) {
            return false;
        }
        
        if(mpfr_zero_p(z->y)) {
            mpfr_set_zero(nt->y, 1);
            
            return eval_newton_rr(nt->x, ev, z->x);
        }
        
        return eval_newton_rc(nt, ev, z);
    } else { // complex poly
        if(ev->P->degree == 0) {
            mpfr_set_inf(nt->x, 1);
            mpfr_set_inf(nt->y, 1);
            
            return false;
        }
        
        if(! eval_analyse(ev)) {
            return false;
        }
        
        if(mpfr_zero_p(z->y)) {
            return eval_newton_cr(nt, ev, z->x);
        }
        
        return eval_newton_cc(nt, ev, z);
    }
}

// MARK: Bare-metal variants

bool eval_analyse(eval ev) {
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
    comp_ptr a = ev->P->a;
    
    list l = list_new(deg + 1);
    
    bool ok = true;
    for (deg_t i = 0; i <= deg && ok; i++) {
        s = comp_s(a + i);
        ok = ok && list_add(l, i, s);
    }
    
    if(! ok) {
        list_free(l);
        
        return false;
    }
    
    // create the concave cover of the list and store it into ev
    concave f = conc_new(l, ev->prec);
    list_free(l);
    
    ev->f = f;
    
    if (f != NULL) {
        ev->P->modified = false;
        
        return true;
    }
    
    return false;
}

bool eval_analyse_r(eval ev) {
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
    mpfr_ptr a = ev->Q->a;
    
    list l = list_new(deg + 1);
    
    bool ok = true;
    for (deg_t i = 0; i <= deg && ok; i++) {
        s = mpfr_s(a + i);
        ok = ok && list_add(l, i, s);
    }
    
    if(! ok) {
        list_free(l);
        
        return false;
    }
    
    // create the concave cover of the list and store it into ev
    concave f = conc_new(l, ev->prec);
    list_free(l);
    
    ev->f = f;
    
    if (f != NULL) {
        ev->Q->modified = false;
        
        return true;
    }
    
    return false;
}

bool eval_val_cc(comp v, eval ev, comp z) {
    if(mpfr_zero_p(z->x) && mpfr_zero_p(z->y)) {
        comp_ptr a = ev->P->a;
        
        comp_set(v, a);
        ev->valErr = comp_log2(v) - ev->prec;
        
        return true;
    }
    
    real_t la = comp_log2(z), err;
    deg_t pp;
    
    pows zn = ev->zn;
    pows_set(zn, z);
    comp_ptr zk;
    
    mpfr_ptr b1 = ev->buf->x;
    mpfr_ptr b2 = ev->buf->y;
    
    concave f = ev->f;
    if(! conc_range(f, la)) {
        return false;
    }
    deg_t st = f->start;
    deg_t en = f->end; // en <= st
    ev->terms = st - en + 1;
    
    // sparse Horner's scheme with the list from @c ev
    comp_ptr a = ev->P->a;
    deg_t *pows = f->all->k;
    real_t *s = f->all->s;
    
    pp = pows[st];
    deg_t p;
    comp_set(v, a + pp);
    
    long prec = f->prec;
    err = s[st] + pp * la - prec;
    
    for (deg_t i = st - 1; i >= en && i <= MAX_DEG; i--) {
        p = pows[i];
        zk = pows_pow(zn, pp - p);
        
        comp_mul(v, v, zk, b1, b2);
        comp_add(v, v, a + p);
        
        err = bits_sum(err, s[i] + p * la - prec);
        
        pp = p;
    }
    
    if(pp > 0) {
        zk = pows_pow_once(zn, pp);
        
        comp_mul(v, v, zk, b1, b2);
    }
    
    ev->valErr = err;
    
    return true;
}

bool eval_der_cc(comp d, eval ev, comp z) {
    if(mpfr_zero_p(z->x) && mpfr_zero_p(z->y)) {
        comp_ptr a = ev->P->a;
        
        comp_set(d, a + 1);
        ev->derErr = comp_log2(d) - ev->prec;
        
        return true;
    }
    
    real_t la = comp_log2(z), err;
    deg_t pp;
    
    pows zn = ev->zn;
    pows_set(zn, z);
    comp_ptr zk;
    
    mpfr_ptr b1 = ev->buf->x;
    mpfr_ptr b2 = ev->buf->y;
    
    concave f = ev->f;
    if(! conc_range_der(f, la)) {
        return false;
    }
    deg_t st = f->start;
    deg_t en = f->end; // en <= st
    ev->terms = st - en + 1;
    
    // sparse Horner's scheme with the list from @c ev
    comp_ptr a = ev->P->a;
    deg_t *pows = f->all->k;
    real_t *s = f->all->s;
    
    pp = pows[st];
    deg_t p;
    comp_set(d, a + pp);
    comp_mulu(d, d, pp);
    
    long prec = f->prec;
    err = s[st] + (pp - 1) * la + plog2(pp) - prec;
    
    for (deg_t i = st - 1; i >= en && i <= MAX_DEG; i--) {
        p = pows[i];
        zk = pows_pow(zn, pp - p);
        
        comp_mul(d, d, zk, b1, b2);
        comp_amu(d, a + p, p, b1);
        
        err = bits_sum(err, s[i] + (p - 1) * la + plog2(p) - prec);
        
        pp = p;
    }
    
    if(pp > 1) {
        zk = pows_pow_once(zn, pp - 1);
        
        comp_mul(d, d, zk, b1, b2);
    }
    
    ev->derErr = err;
    
    return true;
}

bool eval_val_der_cc(comp v, comp d, eval ev, comp z) {
    if(mpfr_zero_p(z->x) && mpfr_zero_p(z->y)) {
        comp_ptr a = ev->P->a;
        
        comp_set(v, a);
        ev->valErr = comp_log2(v) - ev->prec;
        
        comp_set(d, a + 1);
        ev->derErr = comp_log2(d) - ev->prec;
        
        return true;
    }
    
    real_t la = comp_log2(z), erv, erd;
    deg_t pp, pd;
    
    pows zn = ev->zn;
    pows_set(zn, z);
    comp_ptr zk;
    
    mpfr_ptr b1 = ev->buf->x;
    mpfr_ptr b2 = ev->buf->y;
    
    // compute d = P'(z)
    concave f = ev->f;
    if(! conc_range_der(f, la)) {
        return false;
    }
    deg_t st = f->start;
    deg_t en = f->end; // en <= st
    ev->terms = st - en + 1;
    
    // sparse Horner's scheme for P' with the list from @c ev
    comp_ptr a = ev->P->a;
    deg_t *pows = f->all->k;
    real_t *s = f->all->s;
    
    pp = pows[st];
    deg_t p;
    comp_set(d, a + pp);
    comp_mulu(d, d, pp);
    
    long prec = f->prec;
    erd = s[st] + (pp - 1) * la + plog2(pp) - prec;
    
    for (deg_t i = st - 1; i >= en && i <= MAX_DEG; i--) {
        p = pows[i];
        zk = pows_pow(zn, pp - p);
        
        comp_mul(d, d, zk, b1, b2);
        comp_amu(d, a + p, p, b1);
        
        erd = bits_sum(erd, s[i] + (p - 1) * la + plog2(p) - prec);
        
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
    comp_set(v, a + pp);
    
    erv = s[st] + pp * la - prec;
    
    for (deg_t i = st - 1; i >= en && i <= MAX_DEG; i--) {
        p = pows[i];
        zk = pows_pow(zn, pp - p);
        
        comp_mul(v, v, zk, b1, b2);
        comp_add(v, v, a + p);
        
        erv = bits_sum(erv, s[i] + p * la - prec);
        
        pp = p;
    }
    
    if(pd > 1) {
        zk = pows_pow_once(zn, pd - 1);
        
        comp_mul(d, d, zk, b1, b2);
    }
    
    if(pp > 0) {
        zk = pows_pow_once(zn, pp);
        
        comp_mul(v, v, zk, b1, b2);
    }
    
    ev->valErr = erv;
    ev->derErr = erd;
    
    return true;
}

bool eval_newton_cc(comp nt, eval ev, comp z) {
    mpfr_ptr b1 = ev->buf->x;
    mpfr_ptr b2 = ev->buf->y;
    
    if(mpfr_zero_p(z->x) && mpfr_zero_p(z->y)) {
        comp_ptr a = ev->P->a;
        comp_div(nt, a, a + 1, b1, b2, ev->br);
        ev->ntErr = comp_log2(nt) - ev->prec;
        
        return true;
    }
    
    comp v, d;
    comp_init(v, ev->prec);
    comp_init(d, ev->prec);
    deg_t pv, pd, mp = 0;
    
    real_t la = comp_log2(z), erv, erd;
    
    pows zn = ev->zn;
    pows_set(zn, z);
    comp_ptr zk;
    
    concave f = ev->f;
    if(! conc_range_der(f, la)) {
        return false;
    }
    deg_t st = f->start;
    deg_t en = f->end; // en <= st
    ev->terms = st - en + 1;
    
    // sparse Horner's scheme with the list from @c ev
    comp_ptr a = ev->P->a;
    deg_t *pows = f->all->k;
    real_t *s = f->all->s;
    
    deg_t pp = pows[st];
    deg_t p;
    comp_set(d, a + pp);
    comp_mulu(d, d, pp);
    
    long prec = f->prec;
    erd = s[st] + (pp - 1) * la + plog2(pp) - prec;
    
    for (deg_t i = st - 1; i >= en && i <= MAX_DEG; i--) {
        p = pows[i];
        zk = pows_pow(zn, pp - p);
        
        comp_mul(d, d, zk, b1, b2);
        comp_amu(d, a + p, p, b1);
        
        erd = bits_sum(erd, s[i] + (p - 1) * la + plog2(p) - prec);
        
        pp = p;
    }
    
    pd = pp - 1;
    
    if(! conc_range(f, la)) {
        return false;
    }
    st = f->start;
    en = f->end; // en <= st
    ev->terms += st - en + 1;
    
    // sparse Horner's scheme for P with the list from @c ev
    pp = pows[st];
    comp_set(v, a + pp);
    
    erv = s[st] + pp * la - prec;
    
    for (deg_t i = st - 1; i >= en && i <= MAX_DEG; i--) {
        p = pows[i];
        zk = pows_pow(zn, pp - p);
        
        comp_mul(v, v, zk, b1, b2);
        comp_add(v, v, a + p);
        
        erv = bits_sum(erv, s[i] + p * la - prec);
        
        pp = p;
    }
    
    pv = pp;
    
    // Newton term with power simplification
    if(pv > pd) {
        zk = pows_pow_once(zn, pv - pd);
        
        comp_mul(v, v, zk, b1, b2);
        mp = pd;
    } else if(pd > pv) {
        zk = pows_pow_once(zn, pd - pv);
        
        comp_mul(d, d, zk, b1, b2);
        mp = pv;
    } else {
        mp = pd;
    }
    
    comp_div(nt, v, d, b1, b2, ev->br);
    
    real_t sc = mp * la;
    ev->ntErr = nt_err(comp_log2(v) + sc, erv, comp_log2(d) + sc, erd);
    
    comp_clear(v);
    comp_clear(d);
    
    return true;
}

bool eval_val_rc(comp v, eval ev, comp z) {
    if(mpfr_zero_p(z->x) && mpfr_zero_p(z->y)) {
        mpfr_ptr a = ev->Q->a;
        
        mpfr_set(v->x, a, MPFR_RNDN);
        mpfr_set_zero(v->y, 1);
        ev->valErr = real_log2(v->x) - ev->prec;
        
        return true;
    }
    
    real_t la = comp_log2(z), err;
    deg_t pp;
    
    pows zn = ev->zn;
    pows_set(zn, z);
    comp_ptr zk;
    
    mpfr_ptr b1 = ev->buf->x;
    mpfr_ptr b2 = ev->buf->y;
    
    concave f = ev->f;
    if(! conc_range(f, la)) {
        return false;
    }
    deg_t st = f->start;
    deg_t en = f->end; // en <= st
    ev->terms = st - en + 1;
    
    // sparse Horner's scheme with the list from @c ev
    mpfr_ptr a = ev->Q->a;
    deg_t *pows = f->all->k;
    real_t *s = f->all->s;
    
    pp = pows[st];
    deg_t p;
    comp_setr(v, a + pp);
    
    long prec = f->prec;
    err = s[st] + pp * la - prec;
    
    for (deg_t i = st - 1; i >= en && i <= MAX_DEG; i--) {
        p = pows[i];
        zk = pows_pow(zn, pp - p);
        
        comp_mul(v, v, zk, b1, b2);
        mpfr_add(v->x, v->x, a + p, MPFR_RNDN);
        
        err = bits_sum(err, s[i] + p * la - prec);
        
        pp = p;
    }
    
    if(pp > 0) {
        zk = pows_pow_once(zn, pp);
        
        comp_mul(v, v, zk, b1, b2);
    }
    
    ev->valErr = err;
    
    return true;
}

bool eval_der_rc(comp d, eval ev, comp z) {
    if(mpfr_zero_p(z->x) && mpfr_zero_p(z->y)) {
        mpfr_ptr a = ev->Q->a;
        
        mpfr_set(d->x, a + 1, MPFR_RNDN);
        mpfr_set_zero(d->y, 1);
        ev->derErr = real_log2(d->x) - ev->prec;
        
        return true;
    }
    
    real_t la = comp_log2(z), err;
    deg_t pp;
    
    pows zn = ev->zn;
    pows_set(zn, z);
    comp_ptr zk;
    
    mpfr_ptr b1 = ev->buf->x;
    mpfr_ptr b2 = ev->buf->y;
    
    concave f = ev->f;
    if(! conc_range_der(f, la)) {
        return false;
    }
    deg_t st = f->start;
    deg_t en = f->end; // en <= st
    ev->terms = st - en + 1;
    
    // sparse Horner's scheme with the list from @c ev
    mpfr_ptr a = ev->Q->a;
    deg_t *pows = f->all->k;
    real_t *s = f->all->s;
    
    pp = pows[st];
    deg_t p;
    comp_setr(d, a + pp);
    mpfr_mul_ui(d->x, d->x, pp, MPFR_RNDN);
    
    long prec = f->prec;
    err = s[st] + (pp - 1) * la + plog2(pp) - prec;
    
    for (deg_t i = st - 1; i >= en && i <= MAX_DEG; i--) {
        p = pows[i];
        zk = pows_pow(zn, pp - p);
        
        comp_mul(d, d, zk, b1, b2);
        mpfr_mul_ui(b1, a + p, p, MPFR_RNDN);
        mpfr_add(d->x, d->x, b1, MPFR_RNDN);
        
        err = bits_sum(err, s[i] + (p - 1) * la + plog2(p) - prec);
        
        pp = p;
    }
    
    if(pp > 1) {
        zk = pows_pow_once(zn, pp - 1);
        
        comp_mul(d, d, zk, b1, b2);
    }
    
    ev->derErr = err;
    
    return true;
}

bool eval_val_der_rc(comp v, comp d, eval ev, comp z) {
    if(mpfr_zero_p(z->x) && mpfr_zero_p(z->y)) {
        mpfr_ptr a = ev->Q->a;
        
        mpfr_set(v->x, a, MPFR_RNDN);
        mpfr_set_zero(v->y, 1);
        ev->valErr = real_log2(v->x) - ev->prec;
        
        mpfr_set(d->x, a + 1, MPFR_RNDN);
        mpfr_set_zero(d->y, 1);
        ev->derErr = real_log2(d->x) - ev->prec;
        
        return true;
    }
    
    real_t la = comp_log2(z), erv, erd;
    deg_t pp, pd;
    
    pows zn = ev->zn;
    pows_set(zn, z);
    comp_ptr zk;
    
    mpfr_ptr b1 = ev->buf->x;
    mpfr_ptr b2 = ev->buf->y;
    
    // compute d=Q'(z)
    concave f = ev->f;
    if(! conc_range_der(f, la)) {
        return false;
    }
    deg_t st = f->start;
    deg_t en = f->end; // en <= st
    ev->terms = st - en + 1;
    
    // sparse Horner's scheme for Q' with the list from @c ev
    mpfr_ptr a = ev->Q->a;
    deg_t *pows = f->all->k;
    real_t *s = f->all->s;
    
    pp = pows[st];
    deg_t p;
    comp_setr(d, a + pp);
    mpfr_mul_ui(d->x, d->x, pp, MPFR_RNDN);
    
    long prec = f->prec;
    erd = s[st] + (pp - 1) * la + plog2(pp) - prec;
    
    for (deg_t i = st - 1; i >= en && i <= MAX_DEG; i--) {
        p = pows[i];
        zk = pows_pow(zn, pp - p);
        
        comp_mul(d, d, zk, b1, b2);
        mpfr_mul_ui(b1, a + p, p, MPFR_RNDN);
        mpfr_add(d->x, d->x, b1, MPFR_RNDN);
        
        erd = bits_sum(erd, s[i] + (p - 1) * la + plog2(p) - prec);
        
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
    comp_setr(v, a + pp);
    
    erv = s[st] + pp * la - prec;
    
    for (deg_t i = st - 1; i >= en && i <= MAX_DEG; i--) {
        p = pows[i];
        zk = pows_pow(zn, pp - p);
        
        comp_mul(v, v, zk, b1, b2);
        mpfr_add(v->x, v->x, a + p, MPFR_RNDN);
        
        erv = bits_sum(erv, s[i] + p * la - prec);
        
        pp = p;
    }
    
    if(pd > 1) {
        zk = pows_pow_once(zn, pd - 1);
        
        comp_mul(d, d, zk, b1, b2);
    }
    
    if(pp > 0) {
        zk = pows_pow_once(zn, pp);
        
        comp_mul(v, v, zk, b1, b2);
    }
    
    ev->valErr = erv;
    ev->derErr = erd;
    
    return true;
}

bool eval_newton_rc(comp nt, eval ev, comp z) {
    if(mpfr_zero_p(z->x) && mpfr_zero_p(z->y)) {
        mpfr_ptr a = ev->Q->a;
        mpfr_div(nt->x, a, a + 1, MPFR_RNDN);
        ev->ntErr = real_log2(nt->x) - ev->prec;
        mpfr_set_zero(nt->y, 1);
        
        return true;
    }
    
    comp v, d;
    comp_init(v, ev->prec);
    comp_init(d, ev->prec);
    deg_t pv, pd, mp = 0;
    
    real_t la = comp_log2(z), erv, erd;
    
    pows zn = ev->zn;
    pows_set(zn, z);
    comp_ptr zk;
    
    mpfr_ptr b1 = ev->buf->x;
    mpfr_ptr b2 = ev->buf->y;
    
    // compute d=Q'(z)
    concave f = ev->f;
    if(! conc_range_der(f, la)) {
        return false;
    }
    deg_t st = f->start;
    deg_t en = f->end; // en <= st
    ev->terms = st - en + 1;
    
    // sparse Horner's scheme for Q' with the list from @c ev
    mpfr_ptr a = ev->Q->a;
    deg_t *pows = f->all->k;
    real_t *s = f->all->s;
    
    deg_t pp = pows[st];
    deg_t p;
    comp_setr(d, a + pp);
    mpfr_mul_ui(d->x, d->x, pp, MPFR_RNDN);
    
    long prec = f->prec;
    erd = s[st] + (pp - 1) * la + plog2(pp) - prec;
    
    for (deg_t i = st - 1; i >= en && i <= MAX_DEG; i--) {
        p = pows[i];
        zk = pows_pow(zn, pp - p);
        
        comp_mul(d, d, zk, b1, b2);
        mpfr_mul_ui(b1, a + p, p, MPFR_RNDN);
        mpfr_add(d->x, d->x, b1, MPFR_RNDN);
        
        erd = bits_sum(erd, s[i] + (p - 1) * la + plog2(p) - prec);
        
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
    comp_setr(v, a + pp);
    
    erv = s[st] + pp * la - prec;
    
    for (deg_t i = st - 1; i >= en && i <= MAX_DEG; i--) {
        p = pows[i];
        zk = pows_pow(zn, pp - p);
        
        comp_mul(v, v, zk, b1, b2);
        mpfr_add(v->x, v->x, a + p, MPFR_RNDN);
        
        erv = bits_sum(erv, s[i] + p * la - prec);
        
        pp = p;
    }
    
    pv = pp;
    
    // Newton term with power simplification
    if(pv > pd) {
        zk = pows_pow_once(zn, pv - pd);
        
        comp_mul(v, v, zk, b1, b2);
        mp = pd;
    } else if(pd > pv) {
        zk = pows_pow_once(zn, pd - pv);
        
        comp_mul(d, d, zk, b1, b2);
        mp = pv;
    }
    
    comp_div(nt, v, d, b1, b2, ev->br);
    
    real_t sc = mp * la;
    ev->ntErr = nt_err(comp_log2(v) + sc, erv, comp_log2(d) + sc, erd);
    
    comp_clear(v);
    comp_clear(d);
    
    return true;
}

bool eval_val_cr(comp v, eval ev, mpfr_t x) {
    if(mpfr_zero_p(x)) {
        comp_ptr a = ev->P->a;
        
        comp_set(v, a);
        ev->valErr = comp_log2(v) - ev->prec;
        
        return true;
    }
    
    real_t la = real_log2(x), err;
    deg_t pp;
    
    powsr xn = ev->xn;
    powsr_set(xn, x);
    mpfr_ptr xk;
    
    concave f = ev->f;
    if(! conc_range(f, la)) {
        return false;
    }
    deg_t st = f->start;
    deg_t en = f->end; // en <= st
    ev->terms = st - en + 1;
    
    // sparse Horner's scheme with the list from @c ev
    comp_ptr a = ev->P->a;
    deg_t *pows = f->all->k;
    real_t *s = f->all->s;
    
    pp = pows[st];
    deg_t p;
    comp_set(v, a + pp);
    
    long prec = f->prec;
    err = s[st] + pp * la - prec;
    
    for (deg_t i = st - 1; i >= en && i <= MAX_DEG; i--) {
        p = pows[i];
        xk = powsr_pow(xn, pp - p);
        
        comp_mulr(v, v, xk);
        comp_add(v, v, a + p);
        
        err = bits_sum(err, s[i] + p * la - prec);
        
        pp = p;
    }
    
    if(pp > 0) {
        xk = powsr_pow_once(xn, pp);
        
        comp_mulr(v, v, xk);
    }
    
    ev->valErr = err;
    
    return true;
}

bool eval_der_cr(comp d, eval ev, mpfr_t x) {
    if(mpfr_zero_p(x)) {
        comp_ptr a = ev->P->a;
        
        comp_set(d, a + 1);
        ev->derErr = comp_log2(d) - ev->prec;
        
        return true;
    }
    
    real_t la = real_log2(x), err;
    deg_t pp;
    
    powsr xn = ev->xn;
    powsr_set(xn, x);
    mpfr_ptr xk;
    mpfr_ptr b2 = ev->buf->y;
    
    concave f = ev->f;
    if(! conc_range_der(f, la)) {
        return false;
    }
    deg_t st = f->start;
    deg_t en = f->end; // en <= st
    ev->terms = st - en + 1;
    
    // sparse Horner's scheme with the list from @c ev
    comp_ptr a = ev->P->a;
    deg_t *pows = f->all->k;
    real_t *s = f->all->s;
    
    pp = pows[st];
    deg_t p;
    comp_set(d, a + pp);
    comp_mulu(d, d, pp);
    
    long prec = f->prec;
    err = s[st] + (pp - 1) * la + plog2(pp) - prec;
    
    for (deg_t i = st - 1; i >= en && i <= MAX_DEG; i--) {
        p = pows[i];
        xk = powsr_pow(xn, pp - p);
        
        comp_mulr(d, d, xk);
        comp_amu(d, a + p, p, b2);
        
        err = bits_sum(err, s[i] + (p - 1) * la + plog2(p) - prec);
        
        pp = p;
    }
    
    if(pp > 1) {
        xk = powsr_pow_once(xn, pp - 1);
        
        comp_mulr(d, d, xk);
    }
    
    ev->derErr = err;
    
    return true;
}

bool eval_val_der_cr(comp v, comp d, eval ev, mpfr_t x) {
    if(mpfr_zero_p(x)) {
        comp_ptr a = ev->P->a;
        
        comp_set(v, a);
        ev->valErr = comp_log2(v) - ev->prec;
        
        comp_set(d, a + 1);
        ev->derErr = comp_log2(d) - ev->prec;
        
        return true;
    }
    
    real_t la = real_log2(x), erv, erd;
    deg_t pp, pd;
    
    powsr xn = ev->xn;
    powsr_set(xn, x);
    mpfr_ptr xk;
    mpfr_ptr b2 = ev->buf->y;
    
    concave f = ev->f;
    if(! conc_range_der(f, la)) {
        return false;
    }
    deg_t st = f->start;
    deg_t en = f->end; // en <= st
    ev->terms = st - en + 1;
    
    // sparse Horner's scheme for P' with the list from @c ev
    comp_ptr a = ev->P->a;
    deg_t *pows = f->all->k;
    real_t *s = f->all->s;
    
    pp = pows[st];
    deg_t p;
    comp_set(d, a + pp);
    comp_mulu(d, d, pp);
    
    long prec = f->prec;
    erd = s[st] + (pp - 1) * la + plog2(pp) - prec;
    
    for (deg_t i = st - 1; i >= en && i <= MAX_DEG; i--) {
        p = pows[i];
        xk = powsr_pow(xn, pp - p);
        
        comp_mulr(d, d, xk);
        comp_amu(d, a + p, p, b2);
        
        erd = bits_sum(erd, s[i] + (p - 1) * la + plog2(p) - prec);
        
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
    comp_set(v, a + pp);
    
    erv = s[st] + pp * la - prec;
    
    for (deg_t i = st - 1; i >= en && i <= MAX_DEG; i--) {
        p = pows[i];
        xk = powsr_pow(xn, pp - p);
        
        comp_mulr(v, v, xk);
        comp_add(v, v, a + p);
        
        erv = bits_sum(erv, s[i] + p * la - prec);
        
        pp = p;
    }
    
    if(pd > 1) {
        xk = powsr_pow(xn, pd - 1);
        
        comp_mulr(d, d, xk);
    }
    
    if(pp > 0) {
        xk = powsr_pow_once(xn, pp);
        
        comp_mulr(v, v, xk);
    }
    
    ev->valErr = erv;
    ev->derErr = erd;
    
    return true;
}

bool eval_newton_cr(comp nt, eval ev, mpfr_t x) {
    mpfr_ptr b1 = ev->buf->x;
    mpfr_ptr b2 = ev->buf->y;
    
    if(mpfr_zero_p(x)) {
        comp_ptr a = ev->P->a;
        comp_div(nt, a, a + 1, b1, b2, ev->br);
        ev->ntErr = comp_log2(nt) - ev->prec;
        
        return true;
    }
    
    comp v, d;
    comp_init(v, ev->prec);
    comp_init(d, ev->prec);
    deg_t pv, pd, mp = 0;
    
    real_t la = real_log2(x), erv, erd;
    
    powsr xn = ev->xn;
    powsr_set(xn, x);
    mpfr_ptr xk;
    
    concave f = ev->f;
    if(! conc_range_der(f, la)) {
        return false;
    }
    deg_t st = f->start;
    deg_t en = f->end; // en <= st
    ev->terms = st - en + 1;
    
    // sparse Horner's scheme for P' with the list from @c ev
    comp_ptr a = ev->P->a;
    deg_t *pows = f->all->k;
    real_t *s = f->all->s;
    
    deg_t pp = pows[st];
    deg_t p;
    comp_set(d, a + pp);
    comp_mulu(d, d, pp);
    
    long prec = f->prec;
    erd = s[st] + (pp - 1) * la + plog2(pp) - prec;
    
    for (deg_t i = st - 1; i >= en && i <= MAX_DEG; i--) {
        p = pows[i];
        xk = powsr_pow(xn, pp - p);
        
        comp_mulr(d, d, xk);
        comp_amu(d, a + p, p, b2);
        
        erd = bits_sum(erd, s[i] + (p - 1) * la + plog2(p) - prec);
        
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
    comp_set(v, a + pp);
    
    erv = s[st] + pp * la - prec;
    
    for (deg_t i = st - 1; i >= en && i <= MAX_DEG; i--) {
        p = pows[i];
        xk = powsr_pow(xn, pp - p);
        
        comp_mulr(v, v, xk);
        comp_add(v, v, a + p);
        
        erv = bits_sum(erv, s[i] + p * la - prec);
        
        pp = p;
    }
    
    pv = pp;
    
    // Newton term with power simplification
    if(pv > pd) {
        xk = powsr_pow_once(xn, pv - pd);
        
        comp_mulr(v, v, xk);
        mp = pd;
    } else if(pd > pv) {
        xk = powsr_pow_once(xn, pd - pv);
        
        comp_mulr(d, d, xk);
        mp = pv;
    }
    
    comp_div(nt, v, d, b1, b2, ev->br);
    
    real_t sc = mp * la;
    ev->ntErr = nt_err(comp_log2(v) + sc, erv, comp_log2(d) + sc, erd);
    
    comp_clear(v);
    comp_clear(d);
    
    return true;
}

bool eval_val_rr(mpfr_t v, eval ev, mpfr_t x) {
    if(mpfr_zero_p(x)) {
        mpfr_ptr a = ev->Q->a;
        
        mpfr_set(v, a, MPFR_RNDN);
        ev->valErr = real_log2(v) - ev->prec;
        
        return true;
    }
    
    real_t la = real_log2(x), err;
    deg_t pp;
    
    powsr xn = ev->xn;
    powsr_set(xn, x);
    mpfr_ptr xk;
    
    concave f = ev->f;
    if(! conc_range(f, la)) {
        return false;
    }
    deg_t st = f->start;
    deg_t en = f->end; // en <= st
    ev->terms = st - en + 1;
    
    // sparse Horner's scheme with the list from @c ev
    mpfr_ptr a = ev->Q->a;
    deg_t *pows = f->all->k;
    real_t *s = f->all->s;
    
    pp = pows[st];
    deg_t p;
    mpfr_set(v, a + pp, MPFR_RNDN);
    
    long prec = f->prec;
    err = s[st] + pp * la - prec;
    
    for (deg_t i = st - 1; i >= en && i <= MAX_DEG; i--) {
        p = pows[i];
        xk = powsr_pow(xn, pp - p);
        
        mpfr_mul(v, v, xk, MPFR_RNDN);
        mpfr_add(v, v, a + p, MPFR_RNDN);
        
        err = bits_sum(err, s[i] + p * la - prec);
        
        pp = p;
    }
        
    if(pp > 0) {
        xk = powsr_pow_once(xn, pp);
        
        mpfr_mul(v, v, xk, MPFR_RNDN);
    }
    
    ev->valErr = err;
    
    return true;
}

bool eval_der_rr(mpfr_t d, eval ev, mpfr_t x) {
    if(mpfr_zero_p(x)) {
        mpfr_ptr a = ev->Q->a;
        
        mpfr_set(d, a + 1, MPFR_RNDN);
        ev->derErr = real_log2(d) - ev->prec;
        
        return true;
    }
    
    real_t la = real_log2(x), err;
    deg_t pp;
    
    powsr xn = ev->xn;
    powsr_set(xn, x);
    mpfr_ptr xk;
    
    mpfr_ptr b1 = ev->buf->x;
    
    concave f = ev->f;
    if(! conc_range_der(f, la)) {
        return false;
    }
    deg_t st = f->start;
    deg_t en = f->end; // en <= st
    ev->terms = st - en + 1;
    
    // sparse Horner's scheme with the list from @c ev
    mpfr_ptr a = ev->Q->a;
    deg_t *pows = f->all->k;
    real_t *s = f->all->s;
    
    pp = pows[st];
    deg_t p;
    mpfr_set(d, a + pp, MPFR_RNDN);
    mpfr_mul_ui(d, d, pp, MPFR_RNDN);
    
    long prec = f->prec;
    err = s[st] + (pp - 1) * la + plog2(pp) - prec;
    
    for (deg_t i = st - 1; i >= en && i <= MAX_DEG; i--) {
        p = pows[i];
        xk = powsr_pow(xn, pp - p);
        
        mpfr_mul(d, d, xk, MPFR_RNDN);
        mpfr_mul_ui(b1, a + p, p, MPFR_RNDN);
        mpfr_add(d, d, b1, MPFR_RNDN);
        
        err = bits_sum(err, s[i] + (p - 1) * la + plog2(p) - prec);
        
        pp = p;
    }

    if(pp > 1) {
        xk = powsr_pow_once(xn, pp - 1);
        
        mpfr_mul(d, d, xk, MPFR_RNDN);
    }
    
    ev->derErr = err;
    
    return true;
}

bool eval_val_der_rr(mpfr_t v, mpfr_t d, eval ev, mpfr_t x) {
    if(mpfr_zero_p(x)) {
        mpfr_ptr a = ev->Q->a;
        
        mpfr_set(v, a, MPFR_RNDN);
        ev->valErr = real_log2(v) - ev->prec;
        
        mpfr_set(d, a + 1, MPFR_RNDN);
        ev->derErr = real_log2(d) - ev->prec;
        
        return true;
    }
    
    real_t la = real_log2(x), erv, erd;
    deg_t pp, pd;
    
    powsr xn = ev->xn;
    powsr_set(xn, x);
    mpfr_ptr xk;
    
    mpfr_ptr b1 = ev->buf->x;
    
    // compute d=Q'(z)
    concave f = ev->f;
    if(! conc_range_der(f, la)) {
        return false;
    }
    deg_t st = f->start;
    deg_t en = f->end; // en <= st
    ev->terms = st - en + 1;
    
    // sparse Horner's scheme for Q' with the list from @c ev
    mpfr_ptr a = ev->Q->a;
    deg_t *pows = f->all->k;
    real_t *s = f->all->s;
    
    pp = pows[st];
    deg_t p;
    mpfr_set(d, a + pp, MPFR_RNDN);
    mpfr_mul_ui(d, d, pp, MPFR_RNDN);
    
    long prec = f->prec;
    erd = s[st] + (pp - 1) * la + plog2(pp) - prec;
    
    for (deg_t i = st - 1; i >= en && i <= MAX_DEG; i--) {
        p = pows[i];
        xk = powsr_pow(xn, pp - p);
        
        mpfr_mul(d, d, xk, MPFR_RNDN);
        mpfr_mul_ui(b1, a + p, p, MPFR_RNDN);
        mpfr_add(d, d, b1, MPFR_RNDN);
        
        erd = bits_sum(erd, s[i] + (p - 1) * la + plog2(p) - prec);
        
        pp = p;
    }
    
    pd = pp;
    
    // compute v=Q(x)
    if(! conc_range(f, la)) {
        return false;
    }
    st = f->start;
    en = f->end; // en <= st
    ev->terms += st - en + 1;
    
    // sparse Horner's scheme with the list from @c ev
    pp = pows[st];
    mpfr_set(v, a + pp, MPFR_RNDN);
    
    erv = s[st] + pp * la - prec;
    
    for (deg_t i = st - 1; i >= en && i <= MAX_DEG; i--) {
        p = pows[i];
        xk = powsr_pow(xn, pp - p);
        
        mpfr_mul(v, v, xk, MPFR_RNDN);
        mpfr_add(v, v, a + p, MPFR_RNDN);
        
        erv = bits_sum(erv, s[i] + p * la - prec);
        
        pp = p;
    }
    
    if(pd > 1) {
        xk = powsr_pow(xn, pp - 1);
        
        mpfr_mul(d, d, xk, MPFR_RNDN);
    }
    
    if(pp > 0) {
        xk = powsr_pow_once(xn, pp);
        
        mpfr_mul(v, v, xk, MPFR_RNDN);
    }
    
    ev->valErr = erv;
    ev->derErr = erd;
    
    return true;
}

bool eval_newton_rr(mpfr_t nt, eval ev, mpfr_t x) {
    if(mpfr_zero_p(x)) {
        mpfr_ptr a = ev->Q->a;
        mpfr_div(nt, a, a + 1, MPFR_RNDN);
        ev->ntErr = real_log2(nt) - ev->prec;
        
        return true;
    }
    
    mpfr_t v, d;
    mpfr_init2(v, ev->prec);
    mpfr_init2(d, ev->prec);
    
    deg_t pv, pd, mp = 0;
    
    real_t la = real_log2(x), erv, erd;
    
    powsr xn = ev->xn;
    powsr_set(xn, x);
    mpfr_ptr xk;
    
    mpfr_ptr b1 = ev->buf->x;
    
    // compute d=Q'(z)
    concave f = ev->f;
    if(! conc_range_der(f, la)) {
        return false;
    }
    deg_t st = f->start;
    deg_t en = f->end; // en <= st
    ev->terms = st - en + 1;
    
    // sparse Horner's scheme for Q' with the list from @c ev
    mpfr_ptr a = ev->Q->a;
    deg_t *pows = f->all->k;
    real_t *s = f->all->s;
    
    deg_t pp = pows[st];
    deg_t p;
    mpfr_set(d, a + pp, MPFR_RNDN);
    mpfr_mul_ui(d, d, pp, MPFR_RNDN);
    
    long prec = f->prec;
    erd = s[st] + (pp - 1) * la + plog2(pp) - prec;
    
    for (deg_t i = st - 1; i >= en && i <= MAX_DEG; i--) {
        p = pows[i];
        xk = powsr_pow(xn, pp - p);
        
        mpfr_mul(d, d, xk, MPFR_RNDN);
        mpfr_mul_ui(b1, a + p, p, MPFR_RNDN);
        mpfr_add(d, d, b1, MPFR_RNDN);
        
        erd = bits_sum(erd, s[i] + (p - 1) * la + plog2(p) - prec);
        
        pp = p;
    }
    
    pd = pp - 1;
    
    // compute v=Q(x)
    if(! conc_range(f, la)) {
        return false;
    }
    st = f->start;
    en = f->end; // en <= st
    ev->terms += st - en + 1;
    
    // sparse Horner's scheme with the list from @c ev
    pp = pows[st];
    mpfr_set(v, a + pp, MPFR_RNDN);
    
    erv = s[st] + pp * la - prec;
    
    for (deg_t i = st - 1; i >= en && i <= MAX_DEG; i--) {
        p = pows[i];
        xk = powsr_pow(xn, pp - p);
        
        mpfr_mul(v, v, xk, MPFR_RNDN);
        mpfr_add(v, v, a + p, MPFR_RNDN);
        
        erv = bits_sum(erv, s[i] + p * la - prec);
        
        pp = p;
    }
    
    pv = pp;
    
    // Newton term with power simplification
    if(pv > pd) {
        xk = powsr_pow_once(xn, pv - pd);
        
        mpfr_mul(v, v, xk, MPFR_RNDN);
        mp = pd;
    } else if(pd > pv) {
        xk = powsr_pow_once(xn, pd - pv);
        
        mpfr_mul(d, d, xk, MPFR_RNDN);
        mp = pv;
    }
    
    mpfr_div(nt, v, d, MPFR_RNDN);
    
    real_t sc = mp * la;
    ev->ntErr = nt_err(real_log2(v) + sc, erv, real_log2(d) + sc, erd);
    
    mpfr_clear(d);
    mpfr_clear(v);
    
    return true;
}
