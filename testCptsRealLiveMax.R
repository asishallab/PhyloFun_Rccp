require( PhyloFun )
load( "~/workspace/PhyloFun/data/p_mutation_tables_R_image.bin" )
load( "./data/cptsTestMax.RData" )
st.p <- system.time( cpts.p <- condProbsTbls( test.data.max$uniq.branch.lengths,
                                             test.data.max$anno.space$biological_process,
                                             test.data.max$annos.as.strs$biological_process,
                                             GO.TERM.MUTATION.PROBABILITIES.SEQUENCE.DISTANCE
                                             , 4, 10 ) )
