Dytidriver is a method to identify the driver genes through involving the gene dysregulated expression, tissue-specific expression and variation frequency into the human functional interaction network (e.g. human FIN). It considers the fact that driver genes have impact on expression of their downstream genes, they likely interact with each other to form functional modules and those modules should tend to be expressed similarly in the same tissue. 


Requirements to run the Dytidriver:

1. set the work space of R in the ./source directory.

2. load the needed data with following codes.
GenNames=read.csv('/data/realGPL570.csv',blank.lines.skip = TRUE,header = TRUE)

load('/data/cancers/Lung.Rdata')

load('/data/cancers/new_InfluenceGraph.Rdata')

load('/data/New_tissue_matrix.Rdata')

3. source the tissue_ECC.R file and only_tissue.R file

4. Run the following codes:

ECC_tissue_matrix=main_function(patMutMatrix,patOutMatrix,influenceGraph,all_tissue_matrix)

lung_dytidriver=Add_tissue(patMutMatrix,ECC_tissue_matrix)
