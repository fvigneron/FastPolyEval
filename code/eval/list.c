//
//  list.c
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
#include <string.h>

#include "list.h"

// MARK: Static functions

/// @brief Sorts in descending order the vector @c s with and permutes accordingly @c k of length @c len.
///
/// It uses @c bs and @c bk as buffers.
///
/// @warning @c bs and resp. @c bk must be copies of @c b and resp. @c k at the top call.
///
/// @param s the list of numbers to sort
/// @param k some indexes of numbers in @c s
/// @param len the length of this list
/// @param bs a buffer for real numbers
/// @param bk a buffer for indexes
static void list_psortd(real_t *s, deg_t *k, deg_t len, real_t *bs, deg_t *bk) {
    deg_t l = len >> 1;
    deg_t h = len - l;
    
    // here s == bs, k == bk
    switch(len) {
        case 0:
        case 1:
            return;
        case 2: // here s == bs, k == bk, doing quick merge in s, k
            if(s[0] < s[1]) {
                real_t ts = s[0];
                s[0] = s[1];
                s[1] = ts;
                
                deg_t tk = k[0];
                k[0] = k[1];
                k[1] = tk;
            }
            
            return;
            
        case 3: // here l == 1 and h == 2
            if(bs[1] < bs[2]) {
                real_t ts = bs[1];
                bs[1] = bs[2];
                bs[2] = ts;
                
                deg_t tk = bk[1];
                bk[1] = bk[2];
                bk[2] = tk;
            }
            
            break;
            
        case 4: // here l == 2 and h == 2
            if(bs[0] < bs[1]) {
                real_t ts = bs[0];
                bs[0] = bs[1];
                bs[1] = ts;
                
                deg_t tk = bk[0];
                bk[0] = bk[1];
                bk[1] = tk;
            }
            
            if(bs[2] < bs[3]) {
                real_t ts = bs[2];
                bs[2] = bs[3];
                bs[3] = ts;
                
                deg_t tk = bk[2];
                bk[2] = bk[3];
                bk[3] = tk;
            }
            
            break;
            
        case 5: // here l == 2 and h == 3
            if(bs[0] < bs[1]) {
                real_t ts = bs[0];
                bs[0] = bs[1];
                bs[1] = ts;
                
                deg_t tk = bk[0];
                bk[0] = bk[1];
                bk[1] = tk;
            }
            
            list_psortd(bs + 2, bk + 2, h, s + 2, k + 2);
            
            break;
            
        default:
            list_psortd(bs, bk, l, s, k);
            list_psortd(bs + l, bk + l, h, s + l, k + l);
    }
    
    // merge the two lists from bs to s (resp. from bk to k)
    deg_t i = 0, j = l, d = 0;
    real_t vl = bs[i];
    real_t vr = bs[j];
    while(i < l && j < len) {
        if(vl >= vr) {
            s[d] = vl;
            k[d ++] = bk[i ++];
            vl = bs[i];
        } else {
            s[d] = vr;
            k[d ++] = bk[j ++];
            vr = bs[j];
        }
    }
    
    // either this, or the next for, but not both
    for (; i < l; i++) {
        s[d] = bs[i];
        k[d ++] = bk[i];
    }
    
    for (; j < len; j++) {
        s[d] = bs[j];
        k[d ++] = bk[j];
    }
}

// MARK: Implementation of public functions

list list_new(deg_t size) {
    if(size == 0) {
        return NULL;
    }
    
    list l = malloc(sizeof(list_struct));
    list_init(l, size);
    
    return l;
}

bool list_init(list l, deg_t size) {
    if(l == NULL || size == 0) {
        return false;
    }
    
    l->size = size;
    l->count = 0;
    l->sorted = false;
    
    l->s = malloc(sizeof(real_t) * size);
    l->k = malloc(sizeof(deg_t) * size);
    
    return true;
}

list list_clone(list l) {
    if(l == NULL) {
        return NULL;
    }
    
    deg_t size = l->size;
    list c = list_new(size);
    memcpy(c->s, l->s, sizeof(real_t) * size);
    memcpy(c->k, l->k, sizeof(deg_t) * size);
    
    c->count = l->count;
    c->sorted = l->sorted;
    
    return c;
}

bool list_free(list l) {
    if(l == NULL) {
        return false;
    }
    
    bool ok = list_clear(l);
    free(l);
    
    return ok;
}

bool list_clear(list l) {
    if(l == NULL) {
        return false;
    }
    
    int count = 0;
    if(l->s != NULL) {
        free(l->s);
        
        count ++;
    }
    
    if(l->k != NULL) {
        free(l->k);
        
        count ++;
    }
    
    return count == 2;
}

bool list_add(list l, deg_t k, real_t s) {
    if(l == NULL || isnan(s) || l->count > MAX_DEG) {
        return false;
    }

    if(l->count == l->size) { // slow operation, increase memory exponentially fast
        l->size *= LIST_FACTOR;
        l->size += LIST_TERM;
        
        l->k = realloc(l->k, sizeof(deg_t) * l->size);
        l->s = realloc(l->s, sizeof(real_t) * l->size);
    }
    
    l->k[l->count] = k;
    l->s[l->count ++] = s;
    l->sorted = false;
    
    return true;
}

bool list_trim(list l) {
    if(l == NULL || l->s == NULL || l->k == NULL) {
        return false;
    }

    if(l->size > l->count) {
        l->size = l->count;
        
        l->k = realloc(l->k, sizeof(deg_t) * l->size);
        l->s = realloc(l->s, sizeof(real_t) * l->size);
    }
    
    return true;
}

bool list_sort(list l) {
    if(l == NULL) {
        return false;
    }
    
    l->sorted = true;
    
    if(l->size <= 1) {        
        return true;
    }
    
    list b = list_clone(l);
    list_psortd(l->s, l->k, l->count, b->s, b->k);
    list_free(b);
    
    return true;
}
