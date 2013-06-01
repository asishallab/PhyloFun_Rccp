require( PhyloFun )
load( "~/workspace/PhyloFun/data/p_mutation_tables_R_image.bin" )
load( "~/cptsTest.RData" )
st.p <- system.time( cpts.p <- condProbsTbls( test.data$uniq.branch.lengths,
                                             test.data$annos.gt$biological_process,
                                             test.data$annos.as.strs$biological_process,
                                             GO.TERM.MUTATION.PROBABILITIES.SEQUENCE.DISTANCE
                                             , 4, 10 ) )
