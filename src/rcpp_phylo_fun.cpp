#include "rcpp_phylo_fun.h"
using namespace Rcpp ;

NumericVector findMatchingRow( NumericMatrix table, SEXP val, SEXP colInd ){
  BEGIN_RCPP

    NumericVector value = NumericVector( val );
    NumericVector columnIndex = NumericVector( colInd );
    NumericVector clmn = table( _ , columnIndex( 0 ) );
    int indx = clmn.size() - 1;
    for( int i=0; i<clmn.size(); i++ ) {
      if( clmn( i ) >= value( 0 ) ) {
        indx = i;
        break;
      }
    }
    return( table( indx, _ ) );

  END_RCPP
}

NumericVector mutationProbability( SEXP compositeAnnotation, SEXP branchLength, SEXP
    annotsMutProbTables, SEXP distanceColumnIndx ) {
  BEGIN_RCPP

    CharacterVector compAnnos( compositeAnnotation );
    List anMutProbTbls = List( annotsMutProbTables );
    double mutProb = 0.0;
    for ( int i = 0; i < compAnnos.size(); i++ ) {
      std::string singlAnno( compAnnos( i ) );
      NumericMatrix pMutTbl = anMutProbTbls( singlAnno );
      NumericVector pMutRow = findMatchingRow( pMutTbl, branchLength, distanceColumnIndx );
      if ( pMutRow( 0 ) > mutProb ) {
        mutProb = pMutRow( 0 );
      }
    }
    return( NumericVector( 1, mutProb ) );

  END_RCPP
}

SEXP conditionalProbabilityTable( SEXP branchLength, SEXP annos,
    SEXP stringifiedAnnotations, SEXP annotsMutationProbTables, SEXP
    mutTblLengthColIndx ) {
  BEGIN_RCPP

    List annosLst( annos );
    CharacterVector annosAsStrs( stringifiedAnnotations );
    NumericMatrix cpt = NumericMatrix( annosLst.size(), annosLst.size() );
    cpt.attr( "dimnames" ) = List::create( annosAsStrs, annosAsStrs );

    for ( int i = 0; i < annosLst.size(); i++ ) {
      CharacterVector compositeAnnotation = annosLst( i );
      double compAnnoMutProb = 1.0;
      std::string ua = "unknown";
      std::string caFirst = as<std::string>( compositeAnnotation( 0 ) );
      if ( ua != caFirst ) {
        compAnnoMutProb = mutationProbability( compositeAnnotation,
            branchLength, annotsMutationProbTables, mutTblLengthColIndx )( 0 );
      }
      double mutToOtherAnnosProb = compAnnoMutProb / ( annosLst.size() - 1 );
      NumericVector colmn( annosLst.size(), mutToOtherAnnosProb );
      colmn( i ) = 1.0 - compAnnoMutProb;
      cpt( _, i ) = colmn;
    }

    return( wrap( cpt ) );

  END_RCPP
}

SEXP conditionalProbabilityTables( SEXP uniqueEdgeLengths, SEXP annos, SEXP
    stringifiedAnnotations, SEXP annotsMutationProbTableList, SEXP
    mutTblLengthColIndx, SEXP nThreads ) {
  BEGIN_RCPP

    std::cout << "1" << "\n";
    NumericVector numberThreads = NumericVector( nThreads );
    std::cout << "2" << "\n";
    omp_set_num_threads( numberThreads(0) );
    std::cout << "3" << "\n";

    NumericVector edgeLengths = NumericVector( uniqueEdgeLengths );
    std::cout << "4" << "\n";
    CharacterVector edgeLengthsAsStrs = as<CharacterVector>( edgeLengths );
    std::cout << "5" << "\n";
    List cpts = List( 0 );
    List& rCpts = cpts;
    std::cout << "6" << "\n";

    for ( int i = 0; i < edgeLengths.size(); i++ )
    {
      NumericVector currBranchLength( 1, edgeLengths( i ) );
      NumericMatrix cpt = conditionalProbabilityTable( currBranchLength, annos,
          stringifiedAnnotations, annotsMutationProbTableList,
          mutTblLengthColIndx );
      rCpts.push_back( cpt, std::string( edgeLengthsAsStrs( i ) ) );
    }
    return( wrap( cpts ) );

  END_RCPP
}
