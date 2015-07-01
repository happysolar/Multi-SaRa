#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <Rdefines.h>
double quick_select(double arr[], int n);

int which_min(double values[], int start, int end){
	int i, index;
	double min;
	if (start > end) {
		return 0;
	}
	index = start;
	min = values[start - 1];
	for (i = start - 1; i <= end - 1; i++){
		if (values[i] < min){
			min = values[i];
			index = i + 1;
		}
	}
	return index;
}

SEXP backward_del_cutoff(SEXP means, SEXP lengths, SEXP candidates, SEXP num_segments, SEXP cutoff){
	int J = INTEGER_VALUE(num_segments);
	int i, index_min;
	int current_length = J;
	double min;
	double m[J], l[J], d[J - 1];
	int c[J], candidate_index[J];
	for (i = 0; i <= J - 1; i++){
		m[i] = REAL(means)[i];
		l[i] = REAL(lengths)[i];
		c[i] = INTEGER(candidates)[i];
		candidate_index[i] = i;
		if (i <= (J - 2)){
			d[i] = fabs(REAL(means)[i + 1] - REAL(means)[i]);
		}
	}
	index_min = which_min(d, 1, J - 1);
	min = d[index_min - 1];
	while (min < NUMERIC_VALUE(cutoff)){
		current_length = current_length - 1;
		if (current_length == 1) break;
		/*for (i = 0; i <= index_min - 2; i++){
			m stays;
			l stays;
		}*/
		m[index_min - 1] = (m[index_min - 1] * l[index_min - 1] + m[index_min] * l[index_min])/(l[index_min - 1] + l[index_min]);
		l[index_min - 1] = l[index_min - 1] + l[index_min];
		for (i = index_min; i <= current_length - 1; i++){
			m[i] = m[i + 1];
			l[i] = l[i + 1];
			c[i] = c[i + 1];
			candidate_index[i] = candidate_index[i + 1];
		}
		/*for (i = 0; i <= index_min - 3; i++){
			d stays;
		}*/
		d[index_min - 2] = fabs(m[index_min - 1] - m[index_min - 2]);
		if(index_min <= current_length - 1){
			d[index_min - 1] = fabs(m[index_min] - m[index_min - 1]);
		}
		for (i = index_min; i <= current_length - 2; i++){
			d[i] = d[i + 1];
		}
		index_min = which_min(d, 1, current_length - 1);
		min = d[index_min - 1];
	}
	SEXP output_candidate = PROTECT(allocVector(INTSXP, current_length));
	SEXP output_index = PROTECT(allocVector(INTSXP, current_length));
	SEXP output_diff = PROTECT(allocVector(REALSXP, current_length));
	REAL(output_diff)[0] = 0;
	for (i = 0; i <= current_length - 1; i++){
		INTEGER(output_candidate)[i] = c[i];
		INTEGER(output_index)[i] = candidate_index[i];
		if (i >= 1){
			REAL(output_diff)[i] = m[i] - m[i - 1];
		}
	}
	SEXP output = PROTECT(allocVector(VECSXP, 3));
	SET_VECTOR_ELT(output, 0, output_candidate);
	SET_VECTOR_ELT(output, 1, output_index);
	SET_VECTOR_ELT(output, 2, output_diff);
	UNPROTECT(4);
	return output;
}

SEXP median_1(SEXP original, SEXP vec_length, SEXP flank_size)
{
	int i, j;
	int len = INTEGER_VALUE(vec_length);
	int n = INTEGER_VALUE(flank_size) * 2 + 1;
	int flank = INTEGER_VALUE(flank_size);
	double median;
	double temp[n];
        double mean, sd, lower, upper, sum, sum_sq;
	SEXP vec_sum = PROTECT(allocVector(REALSXP, (len - n + 1)));
	SEXP vec_ss = PROTECT(allocVector(REALSXP, (len - n + 1)));
	SEXP output = PROTECT(allocVector(REALSXP, (len - n + 1)));
	sum = 0;
	sum_sq = 0;
	for (i = 0; i < (len - n + 1); i++){
		if (i == 0){
			for (j = 0; j < n; j++){
				sum += REAL(original)[i + j];
			}
			for (j = 0; j < n; j++){
				sum_sq += REAL(original)[i + j] * REAL(original)[i + j];
			}
		}
		else {
			sum = sum - REAL(original)[i-1] + REAL(original)[i + n - 1];
			sum_sq = sum_sq - REAL(original)[i-1] * REAL(original)[i-1] + REAL(original)[i + n - 1] * REAL(original)[i + n - 1];
		}
		REAL(vec_sum)[i] = sum;
		REAL(vec_ss)[i] = sum_sq;
		mean = sum/n;
		sd = sqrt((sum_sq - sum * sum/n)/(n-1));
		upper = mean + 2 * sd;
		lower = mean - 2 * sd;
		if (REAL(original)[flank + i] < upper && REAL(original)[flank + i] > lower){
			REAL(output)[i] = REAL(original)[flank + i];
		}
		else{
			for (j = 0; j < n; j++){
				temp[j] = REAL(original)[i + j];
			}
			median = quick_select(temp, n);
			REAL(output)[i] = median;
		}
	}
	SEXP vec = PROTECT(allocVector(VECSXP, 3));
	SET_VECTOR_ELT(vec, 0, output);
	SET_VECTOR_ELT(vec, 1, vec_sum);
	SET_VECTOR_ELT(vec, 2, vec_ss);
	UNPROTECT(4);
	return vec;
}

SEXP median_2(SEXP original, SEXP vec_length, SEXP flank_size, SEXP multiplier)
{
	int i, j;
	int len = INTEGER_VALUE(vec_length);
	int n = INTEGER_VALUE(flank_size) * 2 + 1;
	int flank = INTEGER_VALUE(flank_size);
	double k = NUMERIC_VALUE(multiplier);
	double temp[n], temp2[n];
        double median, mad, lower, upper;
/*	SEXP vec_median = PROTECT(allocVector(REALSXP, len));
	SEXP vec_mad = PROTECT(allocVector(REALSXP, len));
	SEXP vec_lower = PROTECT(allocVector(REALSXP, len));
	SEXP vec_upper = PROTECT(allocVector(REALSXP, len));
*/
	SEXP output = PROTECT(allocVector(REALSXP, len));
	for (i = 0; i < len; i++){
		if (i < flank || i > (len - flank - 1)){
			REAL(output)[i] = REAL(original)[i];
/*			REAL(vec_median)[i] = 0;
			REAL(vec_mad)[i] = 0;
			REAL(vec_lower)[i] = 0;
			REAL(vec_upper)[i] = 0;
*/
		}
		else{
			for (j = 0; j < n; j++){
				temp[j] = REAL(original)[i + j - flank];
			}
			median = quick_select(temp, n);
			for (j = 0; j < n; j++){
				temp2[j] = fabs(REAL(original)[i + j - flank] - median);
			}
			mad = quick_select(temp2, n);
			upper = median + k * mad;
			lower = median - k * mad;
/*			REAL(vec_median)[i] = median;
			REAL(vec_mad)[i] = mad;
			REAL(vec_upper)[i] = upper;
			REAL(vec_lower)[i] = lower;
*/
			if (REAL(original)[i] < upper && REAL(original)[i] > lower){
				REAL(output)[i] = REAL(original)[i];
			}
			else{
				REAL(output)[i] = median;
			}
		}
	}
/*	SEXP vec = PROTECT(allocVector(VECSXP, 5));
	SET_VECTOR_ELT(vec, 0, output);
	SET_VECTOR_ELT(vec, 1, vec_median);
	SET_VECTOR_ELT(vec, 2, vec_mad);
	SET_VECTOR_ELT(vec, 3, vec_upper);
	SET_VECTOR_ELT(vec, 4, vec_lower);
*/
	UNPROTECT(1);
	return output;
}


/*
 *  This Quickselect routine is based on the algorithm described in
 *  "Numerical recipes in C", Second Edition,
 *  Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
 *  This code by Nicolas Devillard - 1998. Public domain.
 */


#define ELEM_SWAP(a,b) { register double t=(a);(a)=(b);(b)=t; }

double quick_select(double arr[], int n) 
{
    int low, high ;
    int median;
    int middle, ll, hh;

    low = 0 ; high = n-1 ; median = (low + high) / 2;
    for (;;) {
        if (high <= low) /* One element only */
            return arr[median] ;

        if (high == low + 1) {  /* Two elements only */
            if (arr[low] > arr[high])
                ELEM_SWAP(arr[low], arr[high]) ;
            return arr[median] ;
        }

    /* Find median of low, middle and high items; swap into position low */
    middle = (low + high) / 2;
    if (arr[middle] > arr[high])    ELEM_SWAP(arr[middle], arr[high]) ;
    if (arr[low] > arr[high])       ELEM_SWAP(arr[low], arr[high]) ;
    if (arr[middle] > arr[low])     ELEM_SWAP(arr[middle], arr[low]) ;

    /* Swap low item (now in position middle) into position (low+1) */
    ELEM_SWAP(arr[middle], arr[low+1]) ;

    /* Nibble from each end towards middle, swapping items when stuck */
    ll = low + 1;
    hh = high;
    for (;;) {
        do ll++; while (arr[low] > arr[ll]) ;
        do hh--; while (arr[hh]  > arr[low]) ;

        if (hh < ll)
        break;

        ELEM_SWAP(arr[ll], arr[hh]) ;
    }

    /* Swap middle item (in position low) back into correct position */
    ELEM_SWAP(arr[low], arr[hh]) ;

    /* Re-set active partition */
    if (hh <= median)
        low = ll;
        if (hh >= median)
        high = hh - 1;
    }
}

#undef ELEM_SWAP

