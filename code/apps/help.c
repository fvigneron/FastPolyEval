//
//  help.c
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

#include "help.h"

static char *main_comm[APP_COUNT] = {
    "-hyperbolic",
    "-Chebyshev",
    "-Legendre",
    "-Hermite",
    "-Laguerre",
    "-sum",
    "-diff",
    "-prod",
    "-cat",
    "-join",
    "-grid",
    "-exp",
    "-rot",
    "-polar",
    "-roots",
    "-der",
    "-unif",
    "-rand",
    "-normal",
    "-sphere",
    "-eval",
    "-evalD",
    "-evalN",
    "-iterN",
    "-comp",
    "-re",
    "-im",
    "-conj",
    "-tensor",
    "-analyse"
};

static char *main_param[APP_COUNT] = {
    "prec period outFile",
    "prec degree outFile",
    "prec degree outFile",
    "prec degree outFile",
    "prec degree outFile",
    "prec inFile1 inFile2 outFile",
    "prec inFile1 inFile2 outFile",
    "prec inFile1 inFile2 outFile",
    "prec inFile1 inFile2 outFile",
    "prec inFile1 inFile2 outFile",
    "prec inFile1 inFile2 outFile",
    "prec inFile outFile",
    "prec inFile outFile",
    "prec inFile outFile",
    "prec inFile outFile",
    "prec inFile outFile",
    "prec count outFile [start] [end]",
    "prec count outFile [start] [end]",
    "prec count outFile [center] [var]",
    "prec count outFile",
    "prec polyFile pointsFile outFile [checkFile] [count]",
    "prec polyFile pointsFile outFile [checkFile] [count]",
    "prec polyFile pointsFile outFile [checkFile] [count]",
    "prec polyFile pointsFile outFile [checkFile] [count]",
    "prec inFile1 inFile2 [outFile]",
    "prec inFile outFile",
    "prec inFile outFile",
    "prec inFile outFile",
    "prec inFile1 inFile2 outFile",
    "prec polyFile outFile1 [outFile2]"
};

static char *main_desc[APP_COUNT] = {
    "writes the coefficients of the hyperbolic polynomial to a CSV file",
    "writes the coefficients of the Chebyshev polynomial to a CSV file",
    "writes the coefficients of the Legendre polynomial to a CSV file",
    "writes the coefficients of the Hermite polynomial to a CSV file",
    "writes the coefficients of the Laguerre polynomial to a CSV file",
    "computes the sum of two polynomials and writes the result to a CSV file",
    "computes the difference of two polynomials and writes the result to a CSV file",
    "computes the product of two polynomials and writes the result to a CSV file",
    "concatenates two CSV files containing complex numbers",
    "joins the real part of two sequnces into one sequence of complex numbers",
    "computes the set product of the real parts of two sequences",
    "computes the complex exponential of a list of points and writes the results to a CSV file",
    "maps complex numbers (a, b) to a*exp(ib)",
    "computes the points given by polar coordinates on the sphere and writes the result to a CSV file",
    "computes the polynomial with a given list of roots and writes the result to a CSV file",
    "computes the derivative of a polynomials and writes the result to a CSV file",
    "writes count real numbers in arithmetic progression to a CSV file",
    "writes count real random numbers uniformly distributed in an interval to a CSV file",
    "writes count real random numbers with Gaussian distibution to a CSV file",
    "writes polar coordinates approximating an uniform distribution on the sphere to a CSV file",
    "quickly evaluates a polynomial on a set of points",
    "quickly evaluates the derivative of a polynomial on a set of points",
    "quickly evaluates one Newton step of a polynomial on a set of points",
    "quickly iterates and searches for limits of the Newton method (roots of the polynomial)",
    "compares two lists of points",
    "writes the real part of the list of complex numbers",
    "writes the imaginary part of the list of complex numbers",
    "writes the conjugates of the list of complex numbers",
    "computes the tensorial product of the two lists of numbers (c_i = a_i * b_i)",
    "computes the concave cover and the intervals of |z| for which the evaluation strategy changes"
};

static char *main_head[] = {"Task", "Parameters", "Short description"};
static int main_widths[] = {15, 55};
static char **main_lines[] = {main_comm, main_param, main_desc};

static help_struct help_main = {
    .linesCount = APP_COUNT,
    .columnCount = 3,
    .before = "\nThe following tasks can be performed with the corresponding command line arguments below:\n\n",
    .columnWidths = main_widths,
    .headers = main_head,
    .lines = main_lines,
    .after = "\nFor each of the tasks above, use -task -help for more detailed information.\n"
             "If the precision is at most " PRECF_STR ", machine floating numbers are used ("
             TYPEFP_STR ", " TYPEF_STR "),\notherwise MPFR floating numbers are used."
};


static char *hyp_param[] = {
    "precision",
    "period",
    "outFile"
};

static char *poly_param[] = {
    "precision",
    "degree",
    "outFile"
};

static char *hyp_type[] = {
    "required",
    "required",
    "required"
};

static char *hyp_desc[] = {
    "the precision in bits of the coefficients (and intermediary variables)",
    "the period of the hyperbolic centers that are the roots of the polynomial",
    "the path of the output CSV file"
};

static char *poly_desc[] = {
    "the precision in bits of the coefficients (and intermediary variables)",
    "the degree of the polynomial",
    "the path of the output CSV file"
};

static int app_widths_no[] = {18, 18};
static char *app_head_no[] = {"Parameter", "Type", "Description"};
static char **hyp_lines[] = {hyp_param, hyp_type, hyp_desc};
static char **poly_lines[] = {poly_param, hyp_type, poly_desc};

static help_struct help_hyp = {
    .linesCount = 3,
    .columnCount = 3,
    .before = "\nGenerates a CSV file with the coefficients of the hyperbolic polynomial:\n\n",
    .columnWidths = app_widths_no,
    .headers = app_head_no,
    .lines = hyp_lines,
    .after = "\nHyperbolic polynomials coefficients and the app -diff can be combined to obtain Misiurewicz polynomials.\n"
};

static help_struct help_cheb = {
    .linesCount = 3,
    .columnCount = 3,
    .before = "\nGenerates a CSV file with the coefficients of the Chebyshev polynomial:\n\n",
    .columnWidths = app_widths_no,
    .headers = app_head_no,
    .lines = poly_lines,
    .after = "\n"
};

static help_struct help_leg = {
    .linesCount = 3,
    .columnCount = 3,
    .before = "\nGenerates a CSV file with the coefficients of the Legendre polynomial:\n\n",
    .columnWidths = app_widths_no,
    .headers = app_head_no,
    .lines = poly_lines,
    .after = "\n"
};

static help_struct help_her = {
    .linesCount = 3,
    .columnCount = 3,
    .before = "\nGenerates a CSV file with the coefficients of the Hermite polynomial:\n\n",
    .columnWidths = app_widths_no,
    .headers = app_head_no,
    .lines = poly_lines,
    .after = "\n"
};

static help_struct help_lag = {
    .linesCount = 3,
    .columnCount = 3,
    .before = "\nGenerates a CSV file with the coefficients of the Laguerre polynomial:\n\n",
    .columnWidths = app_widths_no,
    .headers = app_head_no,
    .lines = poly_lines,
    .after = "\n"
};

static char *sum_param[] = {
    "precision",
    "inFile1",
    "inFile2",
    "outFile"
};

static char *sum_type[] = {
    "required",
    "required",
    "required",
    "required"
};

static char *sum_desc[] = {
    "precision in bits of the coefficients (and intermediary variables)",
    "the path of the input CSV file containing the coefficients of the first polynomial",
    "the path of the input CSV file containing the coefficients of the second polynomial",
    "the path of the output CSV file"
};

static char **sum_lines[] = {sum_param, sum_type, sum_desc};

static help_struct help_sum = {
    .linesCount = 4,
    .columnCount = 3,
    .before = "\nComputes the sum of two polynomials and writes the result to a CSV file:\n\n",
    .columnWidths = app_widths_no,
    .headers = app_head_no,
    .lines = sum_lines,
    .after = "\n"
};

static help_struct help_diff = {
    .linesCount = 4,
    .columnCount = 3,
    .before = "\nComputes the difference of two polynomials and writes the result to a CSV file:\n\n",
    .columnWidths = app_widths_no,
    .headers = app_head_no,
    .lines = sum_lines,
    .after = "\n"
};

static help_struct help_prod = {
    .linesCount = 4,
    .columnCount = 3,
    .before = "\nComputes the product of two polynomials and writes the result to a CSV file:\n\n",
    .columnWidths = app_widths_no,
    .headers = app_head_no,
    .lines = sum_lines,
    .after = "\n"
};

static char *concat_desc[] = {
    "precision in bits of the complex numbers (and intermediary variables)",
    "the path of the input CSV file containing the first list of points",
    "the path of the input CSV file containing the second list of points",
    "the path of the output CSV file"
};

static char **concat_lines[] = {sum_param, sum_type, concat_desc};

static help_struct help_concat = {
    .linesCount = 4,
    .columnCount = 3,
    .before = "\nConcatenates two CSV files containing complex numbers and writes the result to a CSV file:\n\n",
    .columnWidths = app_widths_no,
    .headers = app_head_no,
    .lines = concat_lines,
    .after = "\n"
};

static help_struct help_join = {
    .linesCount = 4,
    .columnCount = 3,
    .before = "\nJoins the real part of two sequnces into one sequence of complex numbers and writes the result to a CSV file:\n\n",
    .columnWidths = app_widths_no,
    .headers = app_head_no,
    .lines = concat_lines,
    .after = "\n"
};

static help_struct help_grid = {
    .linesCount = 4,
    .columnCount = 3,
    .before = "\nComputes the set product of the real parts of two sequences and writes the result to a CSV file:\n\n",
    .columnWidths = app_widths_no,
    .headers = app_head_no,
    .lines = concat_lines,
    .after = "\n"
};

static char *exp_param[] = {
    "precision",
    "inFile",
    "outFile"
};

static char *exp_type[] = {
    "required",
    "required",
    "required"
};

static char *exp_desc[] = {
    "precision in bits of the coefficients (and intermediary variables)",
    "the path of the input CSV file containing the coordinates",
    "the path of the output CSV file"
};

static char **exp_lines[] = {exp_param, exp_type, exp_desc};

static help_struct help_exp = {
    .linesCount = 3,
    .columnCount = 3,
    .before = "\nComputes the complex exponential of the points in the input file and writes the result to a CSV file:\n\n",
    .columnWidths = app_widths_no,
    .headers = app_head_no,
    .lines = exp_lines,
    .after = "\n"
};

static help_struct help_rot = {
    .linesCount = 3,
    .columnCount = 3,
    .before = "\nMaps complex numbers (a, b) to a exp(ib) and writes the result to a CSV file:\n\n",
    .columnWidths = app_widths_no,
    .headers = app_head_no,
    .lines = exp_lines,
    .after = "\n"
};

static help_struct help_polar = {
    .linesCount = 3,
    .columnCount = 3,
    .before = "\nComputes the stereographical projection of the angles (vertical, horizontal) in the input file\nand writes the result to a CSV file:\n\n",
    .columnWidths = app_widths_no,
    .headers = app_head_no,
    .lines = exp_lines,
    .after = "\n"
};

static char *roots_desc[] = {
    "precision in bits of the coefficients (and intermediary variables)",
    "the path of the input CSV file containing the roots",
    "the path of the output CSV file"
};

static char **roots_lines[] = {exp_param, exp_type, roots_desc};

static help_struct help_roots = {
    .linesCount = 3,
    .columnCount = 3,
    .before = "\nComputes the coefficients of a polynomial with given roots and writes the result to a CSV file:\n\n",
    .columnWidths = app_widths_no,
    .headers = app_head_no,
    .lines = roots_lines,
    .after = "\n"
};

static char *der_desc[] = {
    "precision in bits of the coefficients (and intermediary variables)",
    "the path of the input CSV file containing the coefficients of the polynomial",
    "the path of the output CSV file"
};

static char **der_lines[] = {exp_param, exp_type, der_desc};

static help_struct help_der = {
    .linesCount = 3,
    .columnCount = 3,
    .before = "\nComputes the derivative of a polynomial and writes the result to a CSV file:\n\n",
    .columnWidths = app_widths_no,
    .headers = app_head_no,
    .lines = der_lines,
    .after = "\n"
};

static int app_widths[] = {15, 15, 15};
static char *app_head[] = {"Parameter", "Type", "Default", "Description"};

static char *unif_param[] = {
    "precision",
    "count",
    "outFile",
    "start",
    "end"
};

static char *unif_type[] = {
    "required",
    "required",
    "required",
    "optional",
    "optional"
};

static char *unif_desc[] = {
    "precision in bits of the coefficients (and intermediary variables)",
    "the number of points to generate (at least 1 for -rand, -normal and 2 for -unif)",
    "the path of the output CSV file",
    "the left side of the interval",
    "the right side of the interval"
};

static char *unif_def[] = {
    "",
    "",
    "",
    "0",
    "1"
};

static char **unif_lines[] = {unif_param, unif_type, unif_def, unif_desc};

static help_struct help_unif = {
    .linesCount = 5,
    .columnCount = 4,
    .before = "\nWrites n real numbers in arithmetic progression in the interval [start, end] to a CSV file:\n\n",
    .columnWidths = app_widths,
    .headers = app_head,
    .lines = unif_lines,
    .after = "\n"
};

static help_struct help_rand = {
    .linesCount = 5,
    .columnCount = 4,
    .before = "\nWrites n real random numbers uniformly distributed in the interval [start, end) to a CSV file:\n\n",
    .columnWidths = app_widths,
    .headers = app_head,
    .lines = unif_lines,
    .after = "\n"
};

static char *norm_param[] = {
    "precision",
    "count",
    "outFile",
    "center",
    "deviation"
};

static char *norm_desc[] = {
    "precision in bits of the coefficients (and intermediary variables)",
    "the number of points to generate (at least 1 for -rand, -normal and 2 for -unif)",
    "the path of the output CSV file",
    "the center of the normal distribution",
    "the standard deviation of the normal distribution"
};

static char **norm_lines[] = {norm_param, unif_type, unif_def, norm_desc};

static help_struct help_norm = {
    .linesCount = 5,
    .columnCount = 4,
    .before = "\nWrites n real random numbers with Gaussian distibution N(center, deviation^2) to a CSV file:\n\n",
    .columnWidths = app_widths,
    .headers = app_head,
    .lines = norm_lines,
    .after = "\n"
};

static char *sphere_param[] = {
    "precision",
    "count",
    "outFile"
};

static char *sphere_desc[] = {
    "precision in bits of the coefficients (and intermediary variables)",
    "the number of points to generate (at the equator)",
    "the path of the output CSV file"
};

static char **sphere_lines[] = {sphere_param, exp_type, sphere_desc};

static help_struct help_sphere = {
    .linesCount = 3,
    .columnCount = 3,
    .before = "\nWrites polar coordinates approximating an uniform distribution on the sphere,\nsuch that on the equator there are exactly count points, to a CSV file:\n\n",
    .columnWidths = app_widths_no,
    .headers = app_head_no,
    .lines = sphere_lines,
    .after = "The North and South poles will be omitted.\nA point (a, b) represents (cos a cos b, cos a sin b, sin a) in S^2 or R^3.\n"
};

static char *eval_param[] = {
    "precision",
    "polyFile",
    "pointsFile",
    "outFile",
    "errorsFile",
    "count",
    "checkFile",
    "Horner count"
};

static char *limn_param[] = {
    "precision",
    "polyFile",
    "pointsFile",
    "maxIter",
    "outFile",
    "errorsFile",
    "count",
    "checkFile",
    "Horner count"
};

static char *eval_type[] = {
    "required",
    "required",
    "required",
    "required",
    "optional",
    "optional",
    "optional",
    "optional"
};

static char *limn_type[] = {
    "required",
    "required",
    "required",
    "required",
    "required",
    "optional",
    "optional",
    "optional",
    "optional"
};

static char *eval_desc[] = {
    "precision in bits of the coefficients and points (and intermediary variables)",
    "the path of the input CSV file containing the coefficients of the polynomial",
    "the path of the input CSV file containing the coordinates of the points",
    "the path of the output CSV file",
    "the path of the output file containing the bounds for evaluation errors in bits (log2 scale)",
    "the count of the pre-conditionning and of the quick evaluations (for timing measurements)",
    "the path of the output CSV file containing the results by the Horner's method",
    "the count of evaluations of the polynomial by Horner's method"
};

static char *limn_desc[] = {
    "precision in bits of the coefficients and points (and intermediary variables)",
    "the path of the input CSV file containing the coefficients of the polynomial",
    "the path of the input CSV file containing the coordinates of the starting points",
    "the maximal number of iterates of the Newton method",
    "the path of the output CSV file",
    "the path of the output file containing the bounds for evaluation errors in bits (log2 scale)",
    "the count of the pre-conditionning and of the quick evaluations (for timing measurements)",
    "the path of the output CSV file containing the results by the Horner's method",
    "the count of evaluations of the polynomial by Horner's method"
};

static char *eval_def[] = {
    "",
    "",
    "",
    "",
    "",
    "1",
    "",
    "1"
};

static char *limn_def[] = {
    "",
    "",
    "",
    "",
    "",
    "",
    "1",
    "",
    "1"
};

static char **eval_lines[] = {eval_param, eval_type, eval_def, eval_desc};

static char **limn_lines[] = {limn_param, limn_type, limn_def, limn_desc};

static help_struct help_eval = {
    .linesCount = 8,
    .columnCount = 4,
    .before = "\nQuickly evaluates a polynomial on a set of points and writes the result to a CSV file:\n\n",
    .columnWidths = app_widths,
    .headers = app_head,
    .lines = eval_lines,
    .after = "\n"
};

static help_struct help_evald = {
    .linesCount = 8,
    .columnCount = 4,
    .before = "\nQuickly evaluates the derivative of a polynomial on a set of points and writes the result to a CSV file:\n\n",
    .columnWidths = app_widths,
    .headers = app_head,
    .lines = eval_lines,
    .after = "\n"
};

static help_struct help_evaln = {
    .linesCount = 8,
    .columnCount = 4,
    .before = "\nQuickly evaluates one Newton step of a polynomial on a set of points and writes the result to a CSV file:\n\n",
    .columnWidths = app_widths,
    .headers = app_head,
    .lines = eval_lines,
    .after = "\n"
};

static help_struct help_evalr = {
    .linesCount = 9,
    .columnCount = 4,
    .before = "\nQuickly iterates the Newton method (and searches for eventual limits) of the polynomial and writes the result to a CSV file:\n\n",
    .columnWidths = app_widths,
    .headers = app_head,
    .lines = limn_lines,
    .after = "\nFor best results, we recommend successive calls where the precision is gradually increased and maxIter is reduced.\nFor guidance in the choice of the starting points, please refer to J.H. Hubbard, D. Schleicher, S. Sutherland:\nHow to find all roots of complex polynomials by Newton's method. Invent. math., 146:1-33, 2001.\n"
};

static char *comp_desc[] = {
    "precision in bits of the complex numbers (and intermediary variables)",
    "the path of the input CSV file containing the first list of points",
    "the path of the input CSV file containing the second list of points",
    "the path of the output CSV file containg the relative errors"
};

static char **comp_lines[] = {sum_param, unif_type, comp_desc};

static help_struct help_comp = {
    .linesCount = 4,
    .columnCount = 3,
    .before = "\nCompares two sequences of complex numbers:\n\n",
    .columnWidths = app_widths_no,
    .headers = app_head_no,
    .lines = comp_lines,
    .after = "\n"
};


static char *analyse_param[] = {
    "precision",
    "polyFile",
    "concaveFile",
    "coeffFile"
};

static char *analyse_type[] = {
    "required",
    "required",
    "required",
    "optional"
};

static char *analyse_desc[] = {
    "precision in bits of the complex numbers (and intermediary variables)",
    "the path of the input CSV file containing the coefficients of the polynomial",
    "the path of the output CSV file containing the concave cover",
    "the path of the output CSV file containg the intervals for |z| ans associated coefficients ranges"
};

static char **analyse_lines[] = {analyse_param, analyse_type, analyse_desc};

static help_struct help_analyse = {
    .linesCount = 4,
    .columnCount = 3,
    .before = "\nComputes the concave cover used for quick evaluations. Optionally,\nit outputs the intervals for which the ignored terms change:\n\n",
    .columnWidths = app_widths_no,
    .headers = app_head_no,
    .lines = analyse_lines,
    .after = "The precision is used only for determining which terms to be ignored, not for computing the concave cover.\n"
};

static help_struct help_re = {
    .linesCount = 3,
    .columnCount = 3,
    .before = "\nComputes the real part of the points in the input file and writes the result to a CSV file:\n\n",
    .columnWidths = app_widths_no,
    .headers = app_head_no,
    .lines = exp_lines,
    .after = "The values are written as complex numbers, but their imaginary part is zero.\n"
};

static help_struct help_im = {
    .linesCount = 3,
    .columnCount = 3,
    .before = "\nComputes the imaginary part of the points in the input file and writes the result to a CSV file:\n\n",
    .columnWidths = app_widths_no,
    .headers = app_head_no,
    .lines = exp_lines,
    .after = "The values are written as complex numbers, but their imaginary part is zero.\n"
};

static help_struct help_conj = {
    .linesCount = 3,
    .columnCount = 3,
    .before = "\nComputes the conjugates of the points in the input file and writes the result to a CSV file:\n\n",
    .columnWidths = app_widths_no,
    .headers = app_head_no,
    .lines = exp_lines,
    .after = "\n"
};

static help_struct help_tensor = {
    .linesCount = 4,
    .columnCount = 3,
    .before = "\nComputes pointwise product of two lists of points and writes the result to a CSV file:\n\n",
    .columnWidths = app_widths_no,
    .headers = app_head_no,
    .lines = concat_lines,
    .after = "\n"
};

static help help_list[HELP_COUNT] = {
    &help_main,
    &help_hyp,
    &help_cheb,
    &help_leg,
    &help_her,
    &help_lag,
    &help_sum,
    &help_diff,
    &help_prod,
    &help_concat,
    &help_join,
    &help_grid,
    &help_exp,
    &help_rot,
    &help_polar,
    &help_roots,
    &help_der,
    &help_unif,
    &help_rand,
    &help_norm,
    &help_sphere,
    &help_eval,
    &help_evald,
    &help_evaln,
    &help_evalr,
    &help_comp,
    &help_re,
    &help_im,
    &help_conj,
    &help_tensor,
    &help_analyse
};

bool fpe_help_print(help h) {
    if(h == NULL) {
        printf("Could not print NULL help system !\n");
        
        return false;
    }
    
    printf("%s", h->before);
    
    char format[50];
    printf("   ");
    for (int i = 0; i < h->columnCount - 1; i++) {
        snprintf(format, 45, " %%-%ds", h->columnWidths[i]);
        printf(format, h->headers[i]);
    }
        
    printf(" %s\n", h->headers[h->columnCount - 1]);
    printf("\n");
    
    for(int i = 0; i < h->linesCount; i++) {
        printf("   ");
        for (int j = 0; j < h->columnCount - 1; j++) {
            snprintf(format, 45, " %%-%ds", h->columnWidths[j]);
            printf(format, h->lines[j][i]);
        }
            
        printf(" %s\n", h->lines[h->columnCount - 1][i]);
    }
    
    printf("\n%s\n", h->after);
    
    return true;
}

help fpe_help_get(int ind) {
    if(ind < 0 || ind >= HELP_COUNT) {
        return NULL;
    }
    
    return help_list[ind];
}

void stats_print(char *mes, long prec, deg_t deg, char *numType, bool real, char *polyFile, char *ptsFile,
                 char *outFile, char *outError, double timePP, long count, double timeEV,
                 ulong tpts, char *outHorner, double timeHo, ulong tptsho, ulong tit, ulong titho) {
    int b = 1000000000;
    
    printf("\n%s: %ld, %s, %" FMT_DEG ", %s, %s, %s",
           mes, prec, numType, deg, polyFile, ptsFile, outFile);
    
    if(outError != NULL) {
        printf(", %s", outError);
    }
    
    printf(", %.6f, %lu, %.6f, %" PRIu64, timePP / b, count, timeEV / b, tpts);
    
    if(outHorner != NULL) {
        printf(", %s, %.6f, %" FMT_DEG "\n", outHorner, timeHo / b, tptsho);
    }
    
    if(tit > 0) {
        printf(", %" PRIu64, tit);
    }
    
    if(titho > 0) {
        printf(", %" PRIu64, titho);
    }
        
    printf("\n");
}
