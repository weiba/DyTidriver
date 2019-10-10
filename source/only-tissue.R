Add_tissue=function(patMutMatrix,ECC_tissue_matrix)
{
patMutMatrix=patMutMatrix[,row.names(ECC_tissue_matrix)]
pat_tissue=c()
del_p=c()
  for(p in 1:nrow(patMutMatrix))
{
  mut_names=colnames(patMutMatrix)[which(patMutMatrix[p,]!=0)]
  if(length(mut_names)==0)
  {
    del_p=cbind(del_p,p)
    next
  }
  tissue_mut=matrix(0,nrow = ncol(patMutMatrix),ncol = 1,dimnames = list(colnames(patMutMatrix),'GenName'))
  
    for(j in 1:length(mut_names))
    {
      tissue_mut[mut_names[j],1]=sum(ECC_tissue_matrix[mut_names[j],])
    }
  pat_tissue=cbind(pat_tissue,tissue_mut)
  
  }
 if(length(del_p)>0)
{dimnames(pat_tissue)=list(colnames(patMutMatrix),row.names(patMutMatrix)[-del_p])
}
# else if(length(del_p)==0)
# {
#   dimnames(pat_tissue)=list(colnames(patMutMatrix),row.names(patMutMatrix))
# }

xx=rowSums(pat_tissue)
tissue_rank=matrix(xx,nrow=length(xx),dimnames=list(names(xx),'GenName'))
tissue_rank1=tissue_rank
tissue_rank1=tissue_rank[order(tissue_rank[,1],decreasing = T),]

return(tissue_rank1)
}