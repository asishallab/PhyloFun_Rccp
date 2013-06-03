#include "rcpp_phylo_fun.h"
using namespace Rcpp ;

double findMatchingRow( std::vector< std::vector< double > > table, double
    value, int columnIndex ){

  std::vector< double > clmn = table[ columnIndex ];
  int indx = clmn.size() - 1;
  for( int i=0; i < clmn.size(); ++i ) {
    if( clmn[ i ] >= value ) {
      indx = i;
      break;
    }
  }
  return( clmn[ indx ] );

}

double mutationProbability( std::vector< std::string > compositeAnnotation,
    double branchLength, std::map< std::string, std::vector< std::vector<
    double > > > annotsMutationProbTables, int distanceColumnIndx ) {

    double mutProb = 0.0;
    for ( int i = 0; i < compositeAnnotation.size(); ++i ) {
      std::string singlAnno = compositeAnnotation[ i ];
      std::vector< std::vector< double > > pMutTbl = annotsMutationProbTables[ singlAnno ];
      double pMut = findMatchingRow( pMutTbl, branchLength, distanceColumnIndx );
      if ( pMut > mutProb ) {
        mutProb = pMut;
      }
    }
    return( mutProb );

}

std::vector< std::vector< double > > conditionalProbabilityTable( double
    branchLength, std::vector< std::vector< std::string > >annos, std::map<
    std::string, std::vector< std::vector< double > > >
    annotsMutationProbTables, int mutTblLengthColIndx ) {

  std::vector< std::vector <double> > cpt( annos.size() );
  for ( int i = 0; i < annos.size(); ++i )
  {
    std::vector< std::string > compositeAnnotation = annos[ i ];
    double compAnnoMutProb = 1.0;
    std::string ua = "unknown";
    std::string caFirst = compositeAnnotation[ 0 ];
    if ( ua != caFirst ) {
      compAnnoMutProb = mutationProbability( compositeAnnotation,
          branchLength, annotsMutationProbTables, mutTblLengthColIndx );
    }
    double mutToOtherAnnosProb = compAnnoMutProb / ( annos.size() - 1 );
    std::vector< double > colmn( annos.size(), mutToOtherAnnosProb );
    colmn[ i ] = 1.0 - compAnnoMutProb;
    cpt[ i ] = colmn;
  }

  return( cpt );
}

std::vector< std::vector< std::string > > convertListOfCharacters( List list ) {
  std::vector< std::vector< std::string > > cnvrt( list.size() );
  for (int i = 0; i < list.size(); ++i)
  {
    CharacterVector cv = list( i );
    std::vector< std::string > v( cv.size() );
    cnvrt[ i ] = as< std::vector< std::string > >( list( i ) );
  }
  return( cnvrt );
}

std::vector< std::vector <double> > convertNumericMatrix( NumericMatrix mtrx ) {
  std::vector< std::vector <double> > clmns( mtrx.ncol() );
  for ( int i = 0; i < mtrx.ncol(); ++i )
  {
    NumericVector rClmn = mtrx( _, i );
    clmns[ i ] = as< std::vector< double > >( rClmn );
  }
  return( clmns );
}

std::map< std::string, std::vector< std::vector< double > > > convertListOfMatrices( List list ) {
  std::map< std::string, std::vector< std::vector< double > > > cnvrt;
  CharacterVector nms = list.names();
  for ( int i = 0; i < list.size(); ++i )
  {
    std::string key( nms( i ) ); 
    NumericMatrix m = list( i );
    cnvrt[ key ] = convertNumericMatrix( m );
  }
  return( cnvrt );
}

SEXP r_convertListOfMatrices( SEXP sList ) {
  return( wrap( convertListOfMatrices( List( sList ) ) ) );
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
    CharacterVector edgeLengthsAsStrs = as<CharacterVector>( uniqueEdgeLengths );

    List rAnnos( sAnnos );
    std::vector< std::vector< std::string > > annos = convertListOfCharacters( rAnnos );
    List rAnnotsMutationProbTables( sAnnotsMutationProbTableList );
    std::map< std::string, std::vector< std::vector< double > > >
      annotsMutationProbTables = convertListOfMatrices(
          rAnnotsMutationProbTables );
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
      std::vector< std::vector< double > > cpt = conditionalProbabilityTable(
          currBranchLength, annos, annotsMutationProbTables,
          mutTblLengthColIndx );
      cpts[ i ] = cpt;
    }
    return( wrap( cpts ) );

  END_RCPP
}
