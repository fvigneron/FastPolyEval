//
//  chrono.c
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


#include "chrono.h"


void millis(long df, char *str) {
    if(df > 100) {    // pretty print durations that exceed 100 milliseconds
        double s = df - (df / MINUTE) * MINUTE;
        s = s / 1000;
        int ss = 2 * s + 1;
        ss = ss / 2;
        
        long dd = df / DAY;
        long hh = (df - dd * DAY) / HOUR;
        long r = df - dd * DAY - hh * HOUR;
        long mm = r / MINUTE;
        long mmm = r / MINUTE + (r % MINUTE > MINUTE / 2 ? 1 : 0);
        
        if(df < 1000) {
            snprintf(str, 90, "%.3lf s", s);
        } else if(df < MINUTE) {
            snprintf(str, 90, "%.2lf s", s); // shorter than a minute : seconds, rounded up to 10 milliseconds
        } else if(df < HOUR) {
                snprintf(str, 90, "%ld m %.1lf s", mm, s); // minute long : rounded up to 0.1 seconds
        } else if(df < DAY) {
                snprintf(str, 90, "%ld h %ldm %d s", hh, mm, ss); // hour long : rounded up to 1 second
        } else {
            snprintf(str, 90, "%ld d %ld h %ld m", dd, hh, mmm); // day long : rounded up to 1 minute
        }
    } else {
        snprintf(str, 90, "%ld ms", df);  // raw if shorter than 100 milliseconds
    }
}

void nanos(long ns, char *str) {
    if(ns > 100 * MILLION) {
        millis((ns + MILLION / 2) / MILLION, str);
        
        return;
    } else if(ns > 10 * MILLION) {
        double dns = ns;
        dns /= MILLION;
        
        snprintf(str, 90, "%.1lf ms", dns);
    } else if(ns > MILLION) {
        double dns = ns;
        dns /= MILLION;
        
        snprintf(str, 90, "%.2lf ms", dns);
    } else if(ns > MILLION / 100) {
        double dns = ns;
        dns /= 1000;
        
        snprintf(str, 90, "%.1lf us", dns);
    } else if(ns > 1000) {
        double dns = ns;
        dns /= 1000;
        
        snprintf(str, 90, "%.2lf us", dns);
    } else if(ns == 1000) {
        snprintf(str, 90, "1 us");
    } else {
        snprintf(str, 90, "%ld ns", ns);
    }
}

ulong lap(ptime *ts, char *str) {
    ptime now;
    clock_gettime(USED_CLOCK, &now);
    
    long ns = (now.tv_sec - ts->tv_sec) * BILLION + now.tv_nsec - ts->tv_nsec;
    
    *ts = now;
    
    if(str != NULL) {
        nanos(ns, str);
    }
    
    return ns;
}

void print_time_res() {
    ptime pt;
    clock_getres(USED_CLOCK, &pt);
    
    char nns[100];
    nanos(pt.tv_nsec, nns);
    printf("Timer resolution: %s\n", nns);
}
