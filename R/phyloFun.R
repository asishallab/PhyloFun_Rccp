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
  ret <- .Call( "conditionalProbabilityTables", uniqueEdgeLengths, annos,
               stringifiedAnnotations, annotsMutationProbTableList,
               mutTblLengthColIndx, nThreads, PACKAGE="PhyloFun" )
  return( ret )
}

# testOpenMP <- function( n.threads ) {
#   ret <- .Call( "testOpenMP", n.threads, PACKAGE="PhyloFun" )
#   return( ret )
# }
