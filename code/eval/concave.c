//
//  concave.c
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
#include <stdio.h>

#include "concave.h"

// MARK: Static functions

/// @brief Computes the concave cover @c f of the @b sorted list of scales @c l of the coefficients of some complex (or real) polynomial.
///
/// Only initializes @c f->def, not @c f->all, use @c conc_filter() for this second list.
///
/// @param f the result
/// @param l the list of scales
///
/// @return @ref true is successfull, @ref false otherwise
static inline bool conc_build(concave f, list l) {
    list_t ll, lr; // list left, list right
    list_init(ll, f->def->size);
    list_init(lr, f->def->size);
    
    deg_t *kl = l->k, *ks;
    real_t *sl = l->s, *ss, s = sl[0]; // s is the largest value in l
    // init n = total length, k - the current index; k0 = the index of the largest value in l
    // lkl - last k in ll; lkr - last k in lr
    deg_t n = l->count, k = kl[0], k0 = k, lkl = k, lkr = k;
    
    list_add(ll, k, s);
    list_add(lr, k, s);
    // def. of li/ ri : left / right list has at least 2 points => (a, b) affine coeff.
    bool li = false, ri = false;
    
    // (a, b) for left / right list: the last / first(0) segment is on the line x -> a * x + b
    real_t ar = 0, br = 0, al = 0, bl = 0, ar0 = 0, br0 = 0, al0 = 0, bl0 = 0;
    real_t ymin, ymax, y, y1, y2;
    deg_t imin, imax, imid, k1, k2;
    
    for (deg_t i = 1; i < n; i++) { // index in the list l
        k = kl[i];
        s = sl[i];
        
        if(s <= INF_M) { // done, no more useful coefficients
            break;
        }
        
        if(k > k0) { // right list
            if(! ri) { // initialize the coefficients (a, b) for the right list
                list_add(lr, k, s);
                ri = true;
                lkr = k;
                
                ar = (s - sl[0]) / (k - k0); // k > k0, unsigned !
                br = s - ar * k;
                
                ar0 = ar;
                br0 = br;
                
                continue;
            }
            if(k < lkr) {
                continue;
            }
            
            ymin = ar * k + br;
            if(s < ymin) {
                ar = (s - lr->s[lr->count - 1]) / (k - lkr); // k > lkr, unsigned !
                br = s - ar * k;
                
                list_add(lr, k, s);
                lkr = k;
                
                continue;
            }
            
            ymax = ar0 * k + br0;
            if(s >= ymax) { // larger than the first slope
                lr->count = 1;
                ar = (s - lr->s[0]) / (k - k0);  // k > k0, unsigned !
                br = s - ar * k;
                
                ar0 = ar;
                br0 = br;
                
                list_add(lr, k, s);
                lkr = k;
                
                continue;
            }
            
            // binary search of the last line above @c (k,s)
            imin = 0;
            imax = lr->count - 2;
            ks = lr->k;
            ss = lr->s;
            
            while(imax - imin > 1) {
                imid = (imin + imax) >> 1;
                k1 = ks[imid];
                k2 = ks[imid + 1];
                y1 = ss[imid];
                y2 = ss[imid + 1];
                
                y = (y2 * (k - k1) - y1 * (k - k2)) / (k2 - k1); // k1 < k2 < k, unsigned !
                if(s >= y) {
                    imax = imid;
                } else {
                    imin = imid;
                }
            }
            
            lr->count = imax + 1;
            
            ar = (s - lr->s[imax]) / (k - lr->k[imax]);  // k > lr->k[imax], unsigned !
            br = s - ar * k;
            
            if(imax == 0) {
                ar0 = ar;
                br0 = br;
            }
            
            list_add(lr, k, s);
            lkr = k;
        } else { // left list
            if(! li) { // initialize the coefficients (a, b) for the left list
                list_add(ll, k, s);
                li = true;
                lkl = k;
                
                al = (sl[0] - s) / (k0 - k); // k0 > k, unsigned !
                bl = s - al * k;
                
                al0 = al;
                bl0 = bl;
                
                continue;
            }
            
            if(k > lkl) {
                continue;
            }
            
            ymin = al * k + bl;
            if(s < ymin) {
                al = (ll->s[ll->count - 1] - s) / (lkl - k); // lkl > k, unsigned !
                bl = s - al * k;
                
                list_add(ll, k, s);
                lkl = k;
                
                continue;
            }
            
            ymax = al0 * k + bl0;
            if(s >= ymax) { // larger than the first slope
                ll->count = 1;
                
                al = (ll->s[0] - s) / (k0 - k); // k0 > k, unsigned !
                bl = s - al * k;
                
                al0 = al;
                bl0 = bl;
                
                list_add(ll, k, s);
                lkl = k;
                
                continue;
            }
            
            // binary search of the last line above @c (k,s)
            imin = 0;
            imax = ll->count - 2;
            ks = ll->k;
            ss = ll->s;
            
            while(imax - imin > 1) {
                imid = (imin + imax) >> 1;
                k1 = ks[imid];
                k2 = ks[imid + 1];
                y1 = ss[imid];
                y2 = ss[imid + 1];
                
                y = (y2 * (k1 - k) - y1 * (k2 - k)) / (k1 - k2); // k1 > k2 > k, unsigned !
                if(s >= y) {
                    imax = imid;
                } else {
                    imin = imid;
                }
            }
            
            ll->count = imax + 1;
            
            al = (ll->s[imax] - s) / (ll->k[imax] - k);  // ll->k[imax] > k, unsigned !
            bl = s - al * k;
            
            if(imax == 0) {
                al0 = al;
                bl0 = bl;
            }
            
            list_add(ll, k, s);
            lkl = k;
        }
    }
    
    list def = f->def;
    kl = ll->k;
    sl = ll->s;
    
    for (deg_t i = ll->count - 1; i <= MAX_DEG; i--) {
        list_add(def, kl[i], sl[i]);
    }
    
    kl = lr->k;
    sl = lr->s;
    n = lr->count;
    for (deg_t i = 1; i < n; i++) {
        list_add(def, kl[i], sl[i]);
    }
    
    list_trim(f->def);
    
    list_clear(ll);
    list_clear(lr);
    
    return true;
}

/// @brief Filters the list @c f->all from the @b sorted list @c l and using the definition @c f->def of the concave cover of @c l.
///
/// @note Use only after @c conc_build().
///
/// @param f the concave map
/// @param l the list of scales of coeffcients
/// @param p the inverse permutation of elements in @c l
static inline void conc_filter(concave f, list l, deg_t *p) {
    deg_t count = f->def->count;
    deg_t *k = f->def->k;
    real_t *s = f->def->s;
    
    deg_t st = 0;
    real_t sl = s[st], sr;
    while(sl <= INF_M && sl < count) {
        sl = s[++ st];
    }
    
    deg_t lf = k[st], rt;
    list_add(f->all, lf, sl);
    prec_t del = f->prec + f->extraBits;
    real_t l2l = sl + plog2(lf), l2r;
    
    real_t *sa = l->s;
    
    for (deg_t i = 1; i < count; i++) {
        rt = k[i];
        sr = s[i];
        l2r = sr + plog2(rt);
        
        // check for values of the polynomial
        deg_t d = rt - lf;
        real_t a = (sr - sl) / d;
        real_t b = sr - a * rt;
        
        // check for values of the derivative
        real_t ad = lf == 0 ? 1 : (l2r - l2l) / d;
        real_t bd = lf == 0 ? INF_M : l2r - ad * rt;
        
        real_t y, yd, sj;
        for (deg_t j = lf + 1; j < rt; j++) {
            sj = sa[p[j]];
            y = a * j + b;
            yd = ad * j + bd;
            
            if(sj >= y - del || sj + plog2(j) + 1 >= yd - del) {
                list_add(f->all, j, sj);
            }
        }
        
        lf = rt;
        sl = sr;
        l2l = l2r;
        list_add(f->all, lf, sl);
    }
    
    list_trim(f->all);
}

/// @brief Computes the index in @c f->def that will give the largest modulus of @c a_k*z^k, where @c la=log_2|z|.
///
/// @param f the concave cover of the scales of the coefficients of the polynomial
/// @param la the modulus of @c z in @c log_2 scale
///
/// @return the index in @c f->def.
static inline deg_t conc_argmax_val(concave f, real_t la) {
    // find j = argmax |a_k[j] z^k[j]| ~ argmax s_k[j] + k[j] * la
    list def = f->def;
    if(def->count == 1) {
        return 0;
    }
    
    deg_t *kl = def->k;
    real_t *sl = def->s;
    
    deg_t l = 0, r = def->count - 2; // left, right
    deg_t w = kl[l + 1] - kl[l];
    real_t t = sl[l + 1] - sl[l] + w * la;
    
    if(t <= 0) {
        return l;
    }
    
    w = kl[r + 1] - kl[r];
    t = sl[r + 1] - sl[r] + w * la;
    
    if(t > 0) {
        return r + 1;
    }
    
    deg_t m;
    while(r - l > 1) {
        m = (l + r) >> 1;
        w = kl[m + 1] - kl[m];
        t = sl[m + 1] - sl[m] + w * la;
        
        if(t <= 0) {
            r = m;
        } else {
            l = m;
        }
    }
    
    return r;
}

/// @brief Computes the solution of @c y=affine(f)(t)+la*t, where @c affine(f) is the piecewise affine map
/// of the list @c f(st),...,f(en).
///
/// The list @c f(st),...,f(en) has increasing y-coordinates, both cases @c st<en and @c st>en are used.
///
/// @param f the list of planar points
/// @param st the first index to use in @c f
/// @param en the last index to use in @c f
/// @param la the linear slope correction (equal to @c log_2|z|)
/// @param y the value of the map to search for
///
/// @return the real value @c t, the solution of the above equation.
static inline real_t conc_inverse(list f, deg_t st, deg_t en, real_t la, real_t y) {
    if(st == en) {
        return st;
    }
    
    deg_t *kl = f->k;
    real_t *sl = f->s;
    
    if(y <= sl[st] + kl[st] * la) { // it is "guaranteed" that end value is larger than y
        return kl[st];
    }
    
    if(st < en) { // increasing
        // binary search for the segment of the graph of f(x) + la * x that crosses height y
        deg_t l = st, r = en - 1, m;
        while(r - l > 1) {
            m = (l + r) >> 1;
            if(y > sl[m] + kl[m] * la) {
                l = m;
            } else {
                r = m;
            }
        }
        
        // here y is between the values in l and l + 1
        real_t y1 = sl[l], y2 = sl[l + 1];
        deg_t k1 = kl[l], k2 = kl[l + 1], w = k2 - k1;
        real_t den = y2 - y1 + la * w;
        
        if(den <= PDEL) { // slope very close to 0 or some error, return a safe value
            return k1;
        }
        
        return (y * w + k1 * y2 - k2 * y1) / den;
    } else { // decreasing
        // binary search for the segment of the graph of f(x) + la * x that crosses height y
        deg_t l = en, r = st - 1, m;
        while(r - l > 1) {
            m = (l + r) >> 1;
            if(y <= sl[m] + kl[m] * la) {
                l = m;
            } else {
                r = m;
            }
        }
        
        // here y is between the values in l and l + 1
        real_t y1 = sl[l], y2 = sl[l + 1];
        deg_t k1 = kl[l], k2 = kl[l + 1], w = k2 - k1;
        real_t den = y2 - y1 + la * w;
        
        if(den >= -PDEL) { // slope very close to 0 or some error, return a safe value
            return k2;
        }
        
        return (y * w + k1 * y2 - k2 * y1) / den;
    }
}

/// @brief Binary search of @c k in the increasing list @c l.
///
/// @param f the list
/// @param k the value to search for
/// @param left @ref true to choose the value at the left of @c k if @c k is not in @c l, @ref false to choose the value at the right
///
/// @retuen the position in the increasing list.
static inline deg_t conc_search(list f, deg_t k, bool left) {
    deg_t l = 0, r = f->count - 1, m;
    deg_t *kl = f->k;
    
    if(k <= kl[l]) {
        return l;
    }
    
    if(k >= kl[r]) {
        return r;
    }
    
    while(r - l > 1) {
        m = (r + l) >> 1;
        if(k >= kl[m]) {
            l = m;
        } else {
            r = m;
        }
    }
    
    return k == kl[l] || left ? l : l + 1;
}

/// @brief Computes the index in @c f->def that will give the largest modulus of @c k*a_k*z^{k-1}.
///
/// @param f the concave cover of the scales of the coefficients of the polynomial
/// @param la the modulus of @c z in @c log_2 scale
///
/// @return the index in @c f->def.
static inline deg_t conc_argmax_der(concave f, real_t la) {
    // find argmax |k a_k z^(k - 1)| ~ argmax s_k + (k - 1) * la + log_2(k)
    list def = f->def;
    if(def->count <= 1) {
        return 0;
    }
    
    deg_t *kl = def->k;
    
    deg_t l = kl[0] == 0 ? 1 : 0, r = def->count - 2;
    if(r < l) {
        return l;
    }
    
    real_t *sl = def->s;
    deg_t k1 = kl[l], k2 = kl[l + 1], w = k2 - k1;
    real_t t = sl[l + 1] - sl[l] + w * la + plog2(((real_t) k2) / k1);
    
    if(t <= 0) {
        return l;
    }
    
    k1 = kl[r];
    k2 = kl[r + 1];
    w = k2 - k1;
    t = sl[r + 1] - sl[r] + w * la + plog2(((real_t) k2) / k1);
    
    if(t > 0) {
        return r + 1;
    }
    
    deg_t m;
    while(r - l > 1) {
        m = (l + r) >> 1;
        k1 = kl[m];
        k2 = kl[m + 1];
        w = k2 - k1;
        t = sl[m + 1] - sl[m] + w * la + plog2(((real_t) k2) / k1);
        
        if(t <= 0) {
            r = m;
        } else {
            l = m;
        }
    }
    
    return r;
}

/// @brief Computes the solution of @c y=affine(f)(t)+la*(t-1)+log_2(t), where @c affine(f) is the
/// piecewise affine map of the list @c f(st),...,f(en).
///
/// The list @c f(st),...,f(en) has increasing y-coordinates, both cases @c st<en and @c st>en are used.
///
/// @param f the list of planar points
/// @param st the first index to use in @c f
/// @param en the last index to use in @c f
/// @param la the linear slope correction (equal to @c log_2|z|)
/// @param y the value of the map to search for
///
/// @return the real value @c t, the solution of the above equation.
static inline real_t conc_inverse_der(list f, deg_t st, deg_t en, real_t la, real_t y) {
    if(st == en) {
        return st;
    }
    
    deg_t *kl = f->k, k = kl[st];
    real_t *sl = f->s;
    
    if(y <= sl[st] + (k - 1) * la + plog2(k)) { // it is "guaranteed" that end value is larger than y
        return k;
    }
    
    if(st < en) { // increasing
        // binary search for the segment of the graph of f(x) + la * x that crosses height y
        deg_t l = st, r = en - 1, m;
        while(r - l > 1) {
            m = (l + r) >> 1;
            k = kl[m];
            
            if(y > sl[m] + (k - 1) * la + plog2(k)) {
                l = m;
            } else {
                r = m;
            }
        }
        
        // here y is between the values in l and l + 1
        deg_t k1 = kl[l], k2 = kl[l + 1], w = k2 - k1;
        if(k1 == 0) { // log(0) creates some issues
            return 0.5;
        }
        
        real_t y1 = sl[l] + plog2(k1), y2 = sl[l + 1] + plog2(k2);
        real_t den = y2 - y1 + la * w;
        
        if(den <= PDEL) { // slope very close to 0 or some error, return a safe value
            return k1;
        }
        
        // the map is locally strictly convex because of log_2, adjust accordingly
        real_t kr = (y * w + k1 * y2 - k2 * y1) / den;
        k = kr < 1 ? 0 : pfloor(kr); // works even if negative, while k is unsigned
        y1 = sl[l];
        y2 = sl[l + 1];
        real_t a = (y2 - y1) / w;
        real_t b = y2 - a * (k2 - 1);
        
        while(k > k1 && y <= a * (k - 1) + b + plog2(k)) {
            k --;
        }
        
        return k;
    } else { // decreasing
        // binary search for the segment of the graph of f(x) + la * x that crosses height y
        deg_t l = en, r = st - 1, m;
        while(r - l > 1) {
            m = (l + r) >> 1;
            k = kl[m];
            
            if(y <= sl[m] + (k - 1) * la + plog2(k)) {
                l = m;
            } else {
                r = m;
            }
        }
        
        // here y is between the values in l and l + 1
        deg_t k1 = kl[l], k2 = kl[l + 1], w = k2 - k1;
        real_t y1 = sl[l] + plog2(k1), y2 = sl[l + 1] + plog2(k2);
        real_t den = y2 - y1 + la * w;
        
        if(den >= -PDEL) { // slope very close to 0 or some error, return a safe value
            return k2;
        }
        
        // the map is locally strictly convex because of log_2, adjust accordingly
        k = pceil((y * w + k1 * y2 - k2 * y1) / den);
        y1 = sl[l];
        y2 = sl[l + 1];
        real_t a = (y2 - y1) / w;
        real_t b = y2 - a * (k2 - 1);
        
        while(k < k2 && y <= a * (k - 1) + b + plog2(k)) {
            k ++;
        }
        
        return k;
    }
}

// MARK: Implementation of public functions

concave conc_new(list l, prec_t prec) {
    if(l == NULL || l->count == 0 || prec < precf) {
        return NULL;
    }
    
    if(! l->sorted) {
        list_sort(l);
    }
    
    if(l->s[0] <= INF_M) {
        return NULL;
    }
    
    deg_t n = l->size;
    deg_t size = n / 7 + 100;
    
    concave f = malloc(sizeof(concave_struct));
    list_init(f->def, size);
    list_init(f->all, size);
    
    f->prec = prec;
    f->extraBits = pceil(plog2(n)) + BITS_GUARD;
    
    if(! conc_build(f, l)) {
        conc_free(f);
        
        return NULL;
    }
    
    // inverse permutation
    deg_t *p = malloc(sizeof(deg_t) * n);
    deg_t *k = l->k;
    for (deg_t i = 0; i < n; i++) {
        p[k[i]] = i;
    }
    
    conc_filter(f, l, p);
    
    free(p);
    
    return f;
}

bool conc_free(concave f) {
    if(f == NULL) {
        return false;
    }
    
    f->prec = 0;
    bool ok = list_clear(f->def);
    ok = ok && list_clear(f->all);
    
    free(f);
    
    return ok;
}

bool conc_range(concave f, real_t la) {
    if(isnan(la) || la == INF_M || la == INF_P) {
        return false;
    }
    
    if(f->def->count <= 0) {
        f->mid = 1;
        f->start = 0;
        f->end = 1;
        
        return true;
    }
    
    deg_t m = conc_argmax_val(f, la);
    
    list def = f->def;
    deg_t *kl = def->k;
    real_t *sl = def->s;
    deg_t k0 = kl[m];
    
    real_t ml2 = sl[m] + k0 * la - f->prec - f->extraBits;
    
    deg_t d = def->count - 1;
    deg_t lk = k0 == 0 ? 0 : pceil(conc_inverse(def, 0, m, la, ml2));
    deg_t rk = m == d ? kl[d] : pfloor(conc_inverse(def, d, m, la, ml2));
    
    list all = f->all;
    deg_t li = conc_search(all, lk, false);
    deg_t ri = conc_search(all, rk, true);
    
    kl = all->k;
    sl = all->s;
    deg_t ki = kl[li];
    while(ki < k0 && sl[li] + ki * la < ml2) {
        ki = kl[++ li];
    }
    
    ki = kl[ri];
    while(ki > k0 && sl[ri] + ki * la < ml2) {
        ki = kl[-- ri];
    }
    
    f->mid = m;
    f->start = ri;
    f->end = li;
    
    return kl[li] <= k0 && kl[ri] >= k0;
}

bool conc_range_der(concave f, real_t la) {
    if(isnan(la) || la == INF_M || la == INF_P) {
        return false;
    }
    
    if(f->def->count <= 0 || (f->def->count == 1 && f->def->k[0] == 0)) {
        f->mid = 1;
        f->start = 0;
        f->end = 1;
        
        return true;
    }
    
    deg_t m = conc_argmax_der(f, la);
    
    list def = f->def;
    deg_t *kl = def->k;
    real_t *sl = def->s;
    deg_t k0 = kl[m];
    
    real_t ml2 = sl[m] + (k0 - 1) * la + plog2(k0) - f->prec - f->extraBits;
    
    deg_t d = def->count - 1;
    deg_t lk = k0 <= 1 ? 1 : pceil(conc_inverse_der(def, 0, m, la, ml2));
    deg_t rk = m >= d ? kl[d] : pfloor(conc_inverse_der(def, d, m, la, ml2));
    
    list all = f->all;
    deg_t li = conc_search(all, lk, false);
    deg_t ri = conc_search(all, rk, true);
    
    kl = all->k;
    sl = all->s;
    deg_t ki = kl[li];
    while(ki < k0 && sl[li] + (ki - 1) * la + plog2(ki) < ml2) {
        ki = kl[++ li];
    }
    
    ki = kl[ri];
    while(ki > k0 && sl[ri] + (ki - 1) * la + plog2(ki) < ml2) {
        ki = kl[-- ri];
    }
    
    f->mid = m;
    f->start = ri;
    f->end = li;
    
    return kl[li] <= k0 && kl[ri] >= k0;
}
