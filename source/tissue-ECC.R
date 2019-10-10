GenNames=read.csv('/Users/jane/Desktop/folder/driver_mutation/network/DawnRank/final_examination/tissue_specified/realGPL570.csv',blank.lines.skip = TRUE,header = TRUE)
main_function=function(patMutMatrix,patOutMatrix,influenceGraph,all_tissue_matrix)
{pre_result=preprocess_func(patMutMatrix,patOutMatrix,influenceGraph)
##xl=.buildBipartiteGraph_1(pre_result$patMut,pre_result$patOut,pre_result$infGraph)
ECC_tissue_matrix=compute_ECC(influenceGraph,pre_result$inf,all_tissue_matrix,both = 1)
#ECC_tissue_matrix=compute_ECC(influenceGraph,xl,all_tissue_matrix,both = 1)
ECC_tissue_matrix
}
preprocess_func=function(patMutMatrix,patOutMatrix,influenceGraph)
{
  diag(influenceGraph)=0
  OutInter=intersect(row.names(influenceGraph),colnames(patOutMatrix))
  MutInter=intersect(colnames(influenceGraph),colnames(patMutMatrix))
  patMutMatrix=patMutMatrix[,MutInter]
  patOutMatrix=patOutMatrix[,OutInter]
  pats=intersect(row.names(patMutMatrix),row.names(patOutMatrix))
  patMutMatrix=patMutMatrix[pats,]
  patOutMatrix=patOutMatrix[pats,]
  influenceGraph=influenceGraph[OutInter,MutInter]
  list(patMut=patMutMatrix,patOut=patOutMatrix,infGraph=influenceGraph)
}
.buildBipartiteGraph_1=function(patMut,patOut,infG)
{
  BiG=matrix(0,nrow = ncol(patOut),ncol = ncol(patMut))
  row.names(BiG)=colnames(patOut)
  colnames(BiG)=colnames(patMut)
  
  Candidate_Driver=colnames(patMut)
  for(i in 1:length(Candidate_Driver))
  {
    RelateOut=row.names(infG)[which(infG[,Candidate_Driver[i]]==1)]
    RelatePat=which(patMut[,Candidate_Driver[i]]==1)
    if(length(RelatePat)==0)
    {next}
    for(j in 1:length(RelatePat))
    {
      pat_events=intersect(RelateOut,colnames(patOut)[which(patOut[RelatePat[j],]==1)])
      if(length(pat_events)>0)
      {BiG[pat_events,Candidate_Driver[i]]=1
      }
    }
  }
  count_mutated <- rowSums(BiG)
  count_outlier <- colSums(BiG)
  del_row <- which(count_mutated == 0)
  del_col <- which(count_outlier == 0)
  if (length(del_row) > 0) {
    BiG <- BiG[-del_row, ]
  }
  if (length(del_col) > 0) {
    BiG <- BiG[, -del_col]
  }    
  BiG
}

compute_ECC=function(influenceGraph,xl,all_tissue_matrix,both=1)
{
  MutMatrix=Mutated_Matrix_Degree(influenceGraph,xl)
  tissue_matrix=Mutated_Matrix_TissueandInf(MutMatrix,all_tissue_matrix,both =1)
  row_zero=which(rowSums(tissue_matrix)==0)
  if(length(row_zero)>0)
  {tissue_matrix=tissue_matrix[-row_zero,-row_zero]
  }
  ECC_tissue_matrix=matrix(0,nrow = nrow(tissue_matrix),ncol = ncol(tissue_matrix),dimnames = list(row.names(tissue_matrix),colnames(tissue_matrix)))
  for(i in 1:nrow(tissue_matrix))
  {
    Gen1=which(tissue_matrix[i,]!=0)
    k=i+1
    if(k<nrow(tissue_matrix))
    {for(j in k:nrow(tissue_matrix))
    {
      Gen2=which(tissue_matrix[j,]!=0)
      inter=intersect(Gen1,Gen2)
      if(length(inter)>0)
      {
        Gen1_intersum=sum(tissue_matrix[i,inter])
        Gen2_intersum=sum(tissue_matrix[j,inter])
        Gen1_sum=sum(tissue_matrix[i,])
        Gen2_sum=sum(tissue_matrix[j,])
        ECC_value=(Gen1_intersum+Gen2_intersum)/min(Gen1_sum,Gen2_sum)
        ECC_tissue_matrix[i,j]=ECC_value
        ECC_tissue_matrix[j,i]=ECC_value
      }
    }
    }
  }
  
  ECC_tissue_matrix

}
####################################################################
Mutated_Matrix_Degree=function(influenceGraph,xl)
{
  MutMatrix=matrix(0,nrow=ncol(xl),ncol=ncol(xl),dimnames=list(colnames(xl),colnames(xl)))
  MutMatrix[intersect(colnames(influenceGraph),colnames(xl)),intersect(colnames(influenceGraph),colnames(xl))]=influenceGraph[intersect(colnames(influenceGraph),colnames(xl)),intersect(colnames(influenceGraph),colnames(xl))]
  diag(MutMatrix)=0
  return(MutMatrix)
}
Mutated_Matrix_TissueandInf=function(MutMatrix,all_tissue_matrix,both=1)
{
  if(both==1)
  { 
    tissues_combine=combine_tissue(all_tissue_matrix,MutMatrix)##pcc
    aa=which(tissues_combine!=0)
    bb=which(MutMatrix==1)
    combine_tissueandinf=intersect(aa,bb)
    #browser()
    real_Mut=matrix(0,nrow = nrow(MutMatrix),ncol=ncol(MutMatrix),dimnames = list(row.names(MutMatrix),colnames(MutMatrix)))
    real_Mut[combine_tissueandinf]=tissues_combine[combine_tissueandinf]
  }
  else if(both==0)
  {
    tissues_combine=combine_tissue(all_tissue_matrix,MutMatrix)
    real_Mut=tissues_combine
  }
  return(real_Mut)
}
combine_tissue=function(all_tissue_matrix,MutMatrix)
{
  #GBM_cere_ppc_matrix=Mut_tissue_pcc(MutMatrix,GenNames,all_tissue_matrix$GBM_cere_matrix)
  #GBM_brain_ppc_matrix=Mut_tissue_pcc(MutMatrix,GenNames,all_tissue_matrix$GBM_brain_matrix)
  #aa=which(abs(GBM_brain_ppc_matrix)<0.3)
  #GBM_brain_ppc_matrix[aa]=0
  # bb=which(abs(GBM_cere_ppc_matrix)<0.3)
  # GBM_cere_ppc_matrix[bb]=0
  #tissues_combine=0.5*(abs(GBM_brain_ppc_matrix))+0.5*(abs(GBM_cere_ppc_matrix))####with 35 genes can be mapped to drivers
  
#   Lung_lung_ppc_matrix=Mut_tissue_pcc(MutMatrix,GenNames,all_tissue_matrix$Lung_lung_matrix)
#   aa=which(abs(Lung_lung_ppc_matrix)<0.3)
#   Lung_lung_ppc_matrix[aa]=0
#   tissues_combine=abs(Lung_lung_ppc_matrix)
  
  
  Breast_ov_ppc_matrix=Mut_tissue_pcc(MutMatrix,GenNames,all_tissue_matrix$Breast_ov_matrix)
  Breast_pro_ppc_matrix=Mut_tissue_pcc(MutMatrix,GenNames,all_tissue_matrix$Breast_pro_matrix)
  aa=which(abs(Breast_ov_ppc_matrix)<0.3)
  Breast_ov_ppc_matrix[aa]=0
  bb=which(abs(Breast_pro_ppc_matrix)<0.3)
  Breast_pro_ppc_matrix[bb]=0
  tissues_combine=0.5*(abs(Breast_ov_ppc_matrix))+0.5*(abs(Breast_pro_ppc_matrix))
  
  
  # Ovary_ov_ppc_matrix=Mut_tissue_pcc(MutMatrix,GenNames,all_tissue_matrix$Ovary_ovary_matrix)
  # aa=which(abs(Ovary_ov_ppc_matrix)<0.3)
  # Ovary_ov_ppc_matrix[aa]=0
  #tissues_combine=abs(Ovary_ov_ppc_matrix)
  
#   Prostate_pro_ppc_matrix=Mut_tissue_pcc(MutMatrix,GenNames,all_tissue_matrix$Breast_pro_matrix)
#   aa=which(abs(Prostate_pro_ppc_matrix)<0.3)
#   Prostate_pro_ppc_matrix[aa]=0
#   tissues_combine=abs(Prostate_pro_ppc_matrix)
  return(tissues_combine)
}

Mut_tissue_pcc=function(MutMatrix,GenNames,initial_tissue_matrix)
{
  Gen_name=intersect(row.names(MutMatrix),GenNames[,2])
  Mut_tissue_matrix=MutMatrix[Gen_name,Gen_name]
  expName=GenNames[match(Gen_name,GenNames[,2]),1]
  
  tissue_matrix=cor(t(initial_tissue_matrix[expName,]))
  dimnames(tissue_matrix)=list(Gen_name,Gen_name)
  NewMutmatrix=matrix(0,nrow = nrow(MutMatrix),ncol=ncol(MutMatrix),dimnames = list(row.names(MutMatrix),colnames(MutMatrix)))
  
  NewMutmatrix[Gen_name,Gen_name]=tissue_matrix[Gen_name,Gen_name]
  diag(NewMutmatrix)=0
  return(NewMutmatrix)
}