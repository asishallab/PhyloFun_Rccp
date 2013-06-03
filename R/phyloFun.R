convListOfMatrcs <- function( lst ) {
  ret <- .Call( "r_convertListOfMatrices", lst, PACKAGE="PhyloFun" )
  return( ret )
}

mutProb <- function( comp.anno, brnch.ln, mutProbTbls ) {
  ret <- .Call( "mutationProbability", comp.anno, brnch.ln, mutProbTbls, 4,
               PACKAGE="PhyloFun" )
  return( ret )
}

condProbTbl <- function( branchLength, annos, stringifiedAnnotations,
                        annotsMutationProbTables, mutTblLengthColIndx ) {
  ret <- .Call( "conditionalProbabilityTable", branchLength, annos,
               stringifiedAnnotations, annotsMutationProbTables,
               mutTblLengthColIndx, PACKAGE="PhyloFun" )
  return( ret )
}

condProbsTbls <- function( uniqueEdgeLengths, annos, stringifiedAnnotations,
                          annotsMutationProbTableList, mutTblLengthColIndx,
                          nThreads ) {
  # Pass only the required mutation probability tables as matrices:
  uniq.annos <- intersect(
    unique( unlist( annos ) ), 
    names( annotsMutationProbTableList )
  )
  annotsMutProbTbls <- lapply( annotsMutationProbTableList[ uniq.annos ], as.matrix )
  ret <- .Call( "conditionalProbabilityTables", uniqueEdgeLengths, annos,
               stringifiedAnnotations, annotsMutProbTbls,
               mutTblLengthColIndx, nThreads, PACKAGE="PhyloFun" )
  return( ret )
}

# testOpenMP <- function( n.threads ) {
#   ret <- .Call( "testOpenMP", n.threads, PACKAGE="PhyloFun" )
#   return( ret )
# }
