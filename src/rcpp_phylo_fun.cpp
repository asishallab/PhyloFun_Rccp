#include "rcpp_phylo_fun.h"
using namespace Rcpp ;

double findMatchingRow( NumericMatrix table, double value, int columnIndex ){

    NumericVector clmn = table( _ , columnIndex );
    int indx = clmn.size() - 1;
    for( int i=0; i < clmn.size(); i++ ) {
      if( clmn( i ) >= value ) {
        indx = i;
        break;
      }
    }
    return( table( indx, _ )( 0 ) );

}

double mutationProbability( SEXP compositeAnnotation, double branchLength, SEXP
    annotsMutProbTables, int distanceColumnIndx ) {

    CharacterVector compAnnos( compositeAnnotation );
    List anMutProbTbls = List( annotsMutProbTables );
    double mutProb = 0.0;
    for ( int i = 0; i < compAnnos.size(); i++ ) {
      std::string singlAnno( compAnnos( i ) );
      NumericMatrix pMutTbl = anMutProbTbls( singlAnno );
      double pMut = findMatchingRow( pMutTbl, branchLength, distanceColumnIndx );
      if ( pMut > mutProb ) {
        mutProb = pMut;
      }
    }
    return( mutProb );

}

std::vector< std::vector< double > > conditionalProbabilityTable( double
    branchLength, List annosLst, List annotsMutationProbTables, int
    mutTblLengthColIndx ) {

  std::vector< std::vector <double> > cpt( annosLst.size() );
  for ( int i = 0; i < annosLst.size(); i++ ) {
    CharacterVector compositeAnnotation = annosLst( i );
    double compAnnoMutProb = 1.0;
    std::string ua = "unknown";
    std::string caFirst = as<std::string>( compositeAnnotation( 0 ) );
    if ( ua != caFirst ) {
      compAnnoMutProb = mutationProbability( compositeAnnotation,
          branchLength, annotsMutationProbTables, mutTblLengthColIndx );
    }
    double mutToOtherAnnosProb = compAnnoMutProb / ( annosLst.size() - 1 );
    std::vector< double > colmn( annosLst.size(), mutToOtherAnnosProb );
    colmn[ i ] = 1.0 - compAnnoMutProb;
    cpt[ i ] = colmn;
  }

  return( cpt );
}

SEXP conditionalProbabilityTables( SEXP sUniqueEdgeLengths, SEXP sAnnos, SEXP
    sStringifiedAnnotations, SEXP sAnnotsMutationProbTableList, SEXP
    sMutTblLengthColIndx, SEXP sNThreads ) {
  BEGIN_RCPP

    NumericVector nThreads = NumericVector( sNThreads );
    //std::cout << "2" << "\n";
    omp_set_num_threads( nThreads(0) );
    //std::cout << "3" << "\n";

    NumericVector uniqueEdgeLengths = NumericVector( sUniqueEdgeLengths );
    CharacterVector edgeLengthsAsStrs = as<CharacterVector>( edgeLengths );

    List annos( sAnnos );
    List annotsMutationProbTables( sAnnotsMutationProbTableList );
    NumericVector rMutTblLengthColIndx( sMutTblLengthColIndx );
    int mutTblLengthColIndx = rMutTblLengthColIndx( 0 );
    //std::cout << "5" << "\n";
    /* List cpts = List( 0 ); */
    /* std::map< std::string, std::vector<double> > cpts; */
    //std::cout << "6" << "\n";

    std::vector< std::vector< std::vector< double > > > cpts( uniqueEdgeLengths.size() );
    #pragma omp parallel for
    for ( int i = 0; i < uniqueEdgeLengths.size(); i++ )
    {
      double currBranchLength = uniqueEdgeLengths( i );
      std::vector< std::vector< double > > cpt = c_conditionalProbabilityTable(
          currBranchLength, annos, 
          annotsMutationProbTables, mutTblLengthColIndx );
      cpts[ i ] = cpt;
    }
    return( wrap( cpts ) );

  END_RCPP
}
