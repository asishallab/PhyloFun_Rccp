#ifndef _PhyloFun_RCPP_PHYLO_FUN_H
#define _PhyloFun_RCPP_PHYLO_FUN_H

#include <Rcpp.h>
#include <iostream>

#ifdef _OPENMP
   #include <omp.h>
#else
   #define omp_set_num_thread(x) 1
#endif

using namespace Rcpp ;

double findMatchingRow( NumericMatrix table, double value, int columnIndex ) ;

double mutationProbability( SEXP compositeAnnotation, double branchLength, SEXP
    annotsMutProbTables, int distanceColumnIndx ) ;

SEXP conditionalProbabilityTable( SEXP branchLength, SEXP annos,
    SEXP stringifiedAnnotations, SEXP annotsMutationProbTables, SEXP
    mutTblLengthColIndx );

RcppExport SEXP conditionalProbabilityTables( SEXP uniqueEdgeLengths, SEXP annos, SEXP
    stringifiedAnnotations, SEXP annotsMutationProbTableList, SEXP
    mutTblLengthColIndx, SEXP nThreads );

RcppExport SEXP r_convertListOfMatrices( SEXP sList );

#endif
