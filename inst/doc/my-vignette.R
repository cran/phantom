## ----eval=F--------------------------------------------------------------
#  # load in data set
#  time.course.data = load.data(path_to_input_file/input_file_name.csv)

## ----eval=F--------------------------------------------------------------
#  # load in geneset
#  geneset = load.geneset(path_to_geneset_file/geneset.gmt)

## ----eval=F--------------------------------------------------------------
#  # print the names of genesets in reactome.geneset
#  g.names = geneset.names(geneset)

## ----eval=F--------------------------------------------------------------
#  # load in the demo data in phantom package
#  data("time.course.data")
#  
#  # load in the user downloaded REACTOME geneset (GMT file)
#  reactome.geneset = load.geneset(path_to_geneset_file/REACTOME.geneset.gmt)
#  
#  # run individual gene set mode phantom analysis and store the result in an object
#  obj = run.phantom(data = time.course.data, geneset_list = reactome.geneset,
#  query_geneset='REACTOME_ANTIVIRAL_MECHANISM_BY_IFN_STIMULATED_GENES', ncluster = 2, nsample = 1000)

## ----eval=F--------------------------------------------------------------
#  # load in the demo data in phantom package
#  data("time.course.data")
#  
#  # load in the user downloaded KEGGE geneset (GMT file)
#  kegg.geneset = load.geneset(path_to_geneset_file/kegg.geneset.gmt)
#  
#  # run batch mode phantom analysis and store the result in an object
#  
#  obj = run.phantom.batch(data = time.course.data, geneset_list = kegg.geneset,
#                    maxncluster = 5, nsample = 1000, report_pval = 0.05, report_nmin = 5,
#                    output_dir = file.path(getwd(),'/phantom_result'))

