mutProb <- function( comp.anno, brnch.ln, mutProbTbls ) {
  ret <- .Call( "mutationProbability", comp.anno, brnch.ln, mutProbTbls, 4,
               PACKAGE="PhyloFun" )
  return( ret )
}

condProbTbl <- function( branchLength, annos, stringifiedAnnotations,
                        annotsMutationProbTables, mutTblLengthColIndx ) {
# print( branchLength )
# print( annos )
print( class( stringifiedAnnotations ) )
# print( annotsMutationProbTables )
# print( mutTblLengthColIndx )
# print( unknownAnnot  )
  ret <- .Call( "conditionalProbabilityTable", branchLength, annos,
               stringifiedAnnotations, annotsMutationProbTables,
               mutTblLengthColIndx, PACKAGE="PhyloFun" )
  return( ret )
}

# mutProb <- function( RmutationProbTbls, RgoTerm, RbranchLength ) {
#   ret <- .Call( "mutationProbability", RmutationProbTbls, RgoTerm, RbranchLength,
#     PACKAGE="PhyloFun" )
#   return( ret )
# }
# 
# testOpenMP <- function( n.threads ) {
#   ret <- .Call( "testOpenMP", n.threads, PACKAGE="PhyloFun" )
#   return( ret )
# }
# 
# condProbsTbls <- function( unique.edge.lengths, annotations, stringified.annotations, annots.mutation.prob.table.list ) {
#   ret <- .Call( "jusuf", unique.edge.lengths,
#                annotations, stringified.annotations,
#                annots.mutation.prob.table.list, 5, 1, 'unknown', 20,
#                PACKAGE="PhyloFun" )
#   return( ret )
# }
