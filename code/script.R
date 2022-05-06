####Rhythmicity analysis####===========================================
rm(list=ls())
WD = "/path_to_unzipped_folder"
setwd(WD)
data = read.csv(paste0(WD,'/data/Example_data.csv'), row.names = 1)
clinical = read.csv(paste0(WD,'/data/Example_clinical.csv'),row.names = 1)

psychosis.index = which(clinical$Diagnosis == "PSYCHOSIS")
control.index = which(clinical$Diagnosis == "CONTROL")

n = nrow(data)
genes = row.names(data)

sub.clinical.psychosis = clinical[psychosis.index,]
sub.data.psychosis = data[,match(sub.clinical.psychosis$pair,colnames(data))]
tod.psychosis = sub.clinical.psychosis$CorrectedTOD
names(tod.psychosis) = sub.clinical.psychosis$pair

sub.clinical.control = clinical[control.index,]
sub.data.control = data[,match(sub.clinical.control$pair,colnames(data))]
tod.control = sub.clinical.control$CorrectedTOD
names(tod.control) = sub.clinical.control$pair

all(names(tod.psychosis) == colnames(sub.data.psychosis))
all(names(tod.control) == colnames(sub.data.control))


### 1. Rhythmicity analysis on psychosis subjects--------
library('minpack.lm')
source(paste0(WD,'/code/fitSinCurve.R'))
observed_p = data.frame(genes=genes,A=numeric(n), phase=numeric(n), offset=numeric(n), peak=numeric(n), R2=numeric(n))
set.seed(15213)
for (i in 1:n) {
  out = fitSinCurve(xx=tod.psychosis,observed=as.numeric(sub.data.psychosis[i,]))
  observed_p[i,-1] = c(out$A, out$phase, out$offset, out$peak, out$R2)
  if(i%%1000==0) print(i)
}

#Scatterplots
source(paste0(WD,'/code/circadianDrawing_axis.R'))
plotFolder = paste0(WD,'/output/CorePlots')
dir.create(plotFolder, recursive = T)
setwd(plotFolder)
core_genes = c("ENSG00000133794", "ENSG00000134852", "ENSG00000008405",
               "ENSG00000121671", "ENSG00000126368", "ENSG00000174738",
               "ENSG00000179094", "ENSG00000132326", "ENSG00000049246")
#corresponding symbols: c("ARNTL","CLOCK","CRY1", "CRY2","NR1D1","NR1D2","PER1", "PER2","PER3")
for(i in 1:length(core_genes)){
  agene = core_genes[i]
  fileName = paste0("psychosis_",agene,'.pdf')
  pdf(fileName)
  circadianDrawing(tod=tod.psychosis, expr=unlist(sub.data.psychosis[agene,]), apar=observed_p[observed_p$genes == agene,],
                   labels=rep(1, ncol(sub.data.psychosis)), specInfo='psychosis')
  dev.off()
}

#Permutation to derive p-values
nullFolder = paste0(WD,'/output/Rhythmicity_null_psychosis')
dir.create(nullFolder, recursive = T)
setwd(nullFolder)
groupName = 'psychosis'
library(doParallel)
registerDoParallel()
B = 10 #we use 1000 permutations in the paper, but to save time and computation power, the results in this folder reflect 10 permutations
result = foreach(b = 1:B) %dopar% {
  #for(b in 1:B){
  print(b)	
  library(minpack.lm)
  source(paste0(WD,'/code/fitSinCurve.R'))
  
  null_pare = data.frame(A=numeric(n), phase=numeric(n), offset=numeric(n), peak=numeric(n), R2=numeric(n))
  null_para_file = paste('null_',groupName,'_',b,'.rdata',sep='')
  
  set.seed(b)
  shuffleTOD = sample (tod.psychosis)
  
  for (i in 1:n) {
    out = fitSinCurve(xx=shuffleTOD,observed=unlist(sub.data.psychosis[i,]))
    null_pare[i,] = c(out$A, out$phase, out$offset, out$peak, out$R2)
  }		
  save(null_pare,file=null_para_file)	
}
null_pare_A = matrix(0,n,B)
null_pare_phase = matrix(0,n,B)
null_pare_offset = matrix(0,n,B)
null_pare_peak = matrix(0,n,B)
null_pare_R2 = matrix(0,n,B)

for(b in 1:B){
  print(b)
  file11 = paste('null_',groupName,'_',b,'.rdata',sep='')
  
  load(file11)
  null_pare_A[,b] = null_pare$A
  null_pare_phase[,b] = null_pare$phase
  null_pare_offset[,b] = null_pare$offset
  null_pare_peak[,b] = null_pare$peak
  null_pare_R2[,b] = null_pare$R2		
}
null_para = list(null_para_A=null_pare_A, null_para_phase=null_pare_phase, null_para_offset=null_pare_offset, 
                  null_para_peak=null_pare_peak, null_para_R2=null_pare_R2)
null_para_file = paste('null_',groupName,'.rdata',sep='')
save(null_para,file=null_para_file)
para_R2_pool = c(observed_p$R2,null_para$null_para_R2)
R2Rank_para = 1 - (rank(para_R2_pool)[1:length(observed_p$R2)] - 0.5)/length(para_R2_pool)
observed_p$pvalue = R2Rank_para
observed_p$qvalue = p.adjust(observed_p$pvalue, 'BH')
write.csv(observed_p,file =  paste0(WD,'/output/observed_para_psychosis.csv'))

### 2. Rhythmicity analysis on control subjects----------------------
library('minpack.lm')
source(paste0(WD,'/code/fitSinCurve.R'))
observed_c = data.frame(genes=genes,A=numeric(n), phase=numeric(n), offset=numeric(n), peak=numeric(n), R2=numeric(n))
set.seed(15213)
for (i in 1:n) {
  out = fitSinCurve(xx=tod.control,observed=as.numeric(sub.data.control[i,]))
  observed_c[i,-1] = c(out$A, out$phase, out$offset, out$peak, out$R2)
  if(i%%1000==0) print(i)
}

#Scatterplots
source(paste0(WD,'/code/circadianDrawing_axis.R'))
plotFolder = paste0(WD,'/output/CorePlots')
dir.create(plotFolder, recursive = T)
setwd(plotFolder)
core_genes = c("ENSG00000133794", "ENSG00000134852", "ENSG00000008405",
               "ENSG00000121671", "ENSG00000126368", "ENSG00000174738",
               "ENSG00000179094", "ENSG00000132326", "ENSG00000049246")
#corresponding symbols: c("ARNTL","CLOCK","CRY1", "CRY2","NR1D1","NR1D2","PER1", "PER2","PER3")
for(i in 1:length(core_genes)){
  agene = core_genes[i]
  fileName = paste0("control_",agene,'.pdf')
  pdf(fileName)
  circadianDrawing(tod=tod.control, expr=unlist(sub.data.control[agene,]), apar=observed_c[observed_c$genes == agene,],
                   labels=rep(1, ncol(sub.data.control)), specInfo='control')
  dev.off()
}

#Permutation to derive p-values
nullFolder = paste0(WD,'/output/Rhythmicity_null_control')
dir.create(nullFolder, recursive = T)
setwd(nullFolder)
groupName = 'control'
library(doParallel)
registerDoParallel()
B = 10 #we use 1000 permutations in the paper, but to save time and computation power, the results in this folder reflect 10 permutations
result = foreach(b = 1:B) %dopar% {
  #for(b in 1:B){
  print(b)	
  library(minpack.lm)
  source(paste0(WD,'/code/fitSinCurve.R'))
  
  null_pare = data.frame(A=numeric(n), phase=numeric(n), offset=numeric(n), peak=numeric(n), R2=numeric(n))
  null_para_file = paste('null_',groupName,'_',b,'.rdata',sep='')
  
  set.seed(b)
  shuffleTOD = sample (tod.control)
  
  for (i in 1:n) {
    out = fitSinCurve(xx=shuffleTOD,observed=unlist(sub.data.control[i,]))
    null_pare[i,] = c(out$A, out$phase, out$offset, out$peak, out$R2)
  }		
  save(null_pare,file=null_para_file)	
}
null_pare_A = matrix(0,n,B)
null_pare_phase = matrix(0,n,B)
null_pare_offset = matrix(0,n,B)
null_pare_peak = matrix(0,n,B)
null_pare_R2 = matrix(0,n,B)

for(b in 1:B){
  print(b)
  file11 = paste('null_',groupName,'_',b,'.rdata',sep='')
  
  load(file11)
  null_pare_A[,b] = null_pare$A
  null_pare_phase[,b] = null_pare$phase
  null_pare_offset[,b] = null_pare$offset
  null_pare_peak[,b] = null_pare$peak
  null_pare_R2[,b] = null_pare$R2		
}
null_para = list(null_para_A=null_pare_A, null_para_phase=null_pare_phase, null_para_offset=null_pare_offset, 
                 null_para_peak=null_pare_peak, null_para_R2=null_pare_R2)
null_para_file = paste('null_',groupName,'.rdata',sep='')
save(null_para,file=null_para_file)
para_R2_pool = c(observed_c$R2,null_para$null_para_R2)
R2Rank_para = 1 - (rank(para_R2_pool)[1:length(observed_c$R2)] - 0.5)/length(para_R2_pool)
observed_c$pvalue = R2Rank_para
observed_c$qvalue = p.adjust(observed_c$pvalue, 'BH')
write.csv(observed_c,file =  paste0(WD,'/output/observed_para_control.csv'))


## 3. Change of rhythmicity analysis-------------------
tod = c(tod.psychosis, tod.control)

nullFolder = paste0(WD,"/output/ChangeRhyth_null_psychosis_control")
dir.create(nullFolder, recursive = T)
setwd(nullFolder)
registerDoParallel()
B = 10
result = foreach(b = 1:B) %dopar% {
  print(b)	
  library(minpack.lm)
  source(paste0(WD,'/code/fitSinCurve.R'))
  
  set.seed(b)
  index1 = sample(1:length(tod), length(tod.psychosis))
  index2 = setdiff(1:length(tod), index1)
  
  shuffleTOD1 = tod[index1]
  shuffleTOD2 = tod[index2]
  
  n = nrow(data)
  null_pare1 = null_pare2 = data.frame(A=numeric(n), phase=numeric(n), offset=numeric(n), peak=numeric(n), R2=numeric(n))
  
  for (i in 1:n) {
    out1 = fitSinCurve(xx=shuffleTOD1,observed=unlist(sub.data.psychosis[i,]))
    null_pare1[i,] = c(out1$A, out1$phase, out1$offset, out1$peak, out1$R2)
    
    out2 = fitSinCurve(xx=shuffleTOD2,observed=unlist(sub.data.control[i,]))
    null_pare2[i,] = c(out2$A, out2$phase, out2$offset, out2$peak, out2$R2)
  }		
  save(null_pare1,null_pare2,file=paste('null_1psychosis_2control_',b,'.rdata',sep=''))	
}

observe_ctrl = read.csv(paste0(WD,"/output/observed_para_control.csv"),row.names = 1)
observe_case = read.csv(paste0(WD,"/output/observed_para_psychosis.csv"),row.names = 1)

R2change = observe_case$R2 - observe_ctrl$R2 #gain case
Ashift = abs(abs(observe_case$A) - abs(observe_ctrl$A))
phaseshift0 = observe_case$phase - observe_ctrl$phase
phaseshift = pmin(abs(phaseshift0), 24 - abs(phaseshift0))
intshift = abs(observe_case$offset - observe_ctrl$offset)

R2changeNULL = R2change
AshiftNULL = Ashift
phaseshiftNULL = phaseshift
intshiftNULL = intshift

for(b in 1:B){
  print(b)
  load(paste('null_1psychosis_2control_',b,'.rdata',sep=''))
  null_para_ctrl = null_pare2
  null_para_case = null_pare1
  
  null_para_case$R2[which(null_para_case$R2 < 0)] =0
  null_para_ctrl$R2[which(null_para_ctrl$R2 < 0)] =0
  
  aR2change = null_para_case$R2 - null_para_ctrl$R2
  aAshift = abs(abs(null_para_case$A) - abs(null_para_ctrl$A))
  aphaseshift = pmin(abs(null_para_case$phase - null_para_ctrl$phase),24 - abs(null_para_case$phase - null_para_ctrl$phase))
  aintshift = abs(null_para_case$offset - null_para_ctrl$offset)
  
  R2changeNULL = c(R2changeNULL,aR2change)
  AshiftNULL = c(AshiftNULL,aAshift)
  phaseshiftNULL = c(phaseshiftNULL,aphaseshift)
  intshiftNULL = c(intshiftNULL,aintshift)
}
R2gainPvalue = 1-rank(R2changeNULL)[1:n]/length(R2changeNULL) ## R2 psychosis > healthy
R2losePvalue = rank(R2changeNULL)[1:n]/length(R2changeNULL) ## R2 psychosis < healthy
AshiftPvalue = 1 - rank(AshiftNULL)[1:n]/length(AshiftNULL)
phasePvalue = 1 - rank(phaseshiftNULL)[1:n]/length(phaseshiftNULL)
intshiftPvalue = 1 - rank(intshiftNULL)[1:n]/length(intshiftNULL)
rthmicChange = data.frame(genes=observe_ctrl$gene, 
                           R2losePvalue=R2losePvalue,
                           R2gainPvalue=R2gainPvalue,
                           AshiftPvalue=AshiftPvalue,
                           phasePvalue=phasePvalue,
                           intshiftPvalue=intshiftPvalue)
write.csv(rthmicChange,paste0(WD,"/output/ChangeRhyth_psychosisVScontrol.csv"))


#### DE analysis ####============================================================
rm(list=ls())
WD = "/path_to_example_data_and_code_folder"
setwd(WD)
data = read.csv(paste0(WD,'/data/Example_data.csv'), row.names = 1)
clinical = read.csv(paste0(WD,'/data/Example_clinical.csv'),row.names = 1)
all(clinical$pair == colnames(data))

# #Day subjects and night subjects used in paper are selected by the index.d and index.n below. 
# #Here we demonstrate our analysis with all subjects i.e., data and clinical
# index.d = which(clinical$CorrectedTOD>=0 & clinical$CorrectedTOD<12)#day subjects
# index.n = which(clinical$CorrectedTOD<0|clinical$CorrectedTOD>=12)#night subjects
# 
# data.d = data[,index.d]
# clinical.d = clinical[index.d,]
# data.n = data[,index.n]
# clinical.n = clinical[index.n,]

dir=paste0(WD,"/output/DE_output")
dir.create(dir,recursive = T)
setwd(dir)

diagnosis = factor(clinical$Diagnosis,levels = c("CONTROL","PSYCHOSIS"))
coefficient = clinical[,c("Sex.Label","Age","Race.Label","RIN_C")] #candidate confounding variables

VariableListOne = combn(colnames(coefficient),1)
VariableListTwo = combn(colnames(coefficient),2)

source(paste0(WD,"/code/bestSelection.R"))
set.seed(15213)
GeneIndeces = 1:100 #demonstrate with the first 100 genes
library(snowfall) #if running on server, use parallel computing to speed up process by allowing multipe cpus e.g.,24
sfInit(parallel=TRUE,type="SOCK", cpus=1) 
sfExportAll()
indLM=sfLapply(GeneIndeces,function(x) try(bestModelSelection_LM(x,as.matrix(data),covariables3 = coefficient,
                                                                  Diagnosis = diagnosis,VariableListOne,VariableListTwo),silent=TRUE ) ) 

sfStop()
save(indLM,file='LM_true.RData')


#Run permutation test to correct for biased pval from LRT
null_dir = paste0(WD,"/output/DE_output/DE_null")
dir.create(null_dir,recursive = T)
setwd(null_dir)
B=10 #we use 100 permutations in the paper, but to save time and computation power, the results in this folder reflect 10 permutations
sfInit(parallel=TRUE,type="SOCK", cpus=1) 
for(b in 1:B){
  print(b)
  set.seed(b)
  data.b = data[,sample(ncol(data))]
  result.b = replicate(nrow(data.b),list())
  names(result.b) = rownames(data.b)
  sfExportAll()
  result.b=sfLapply(1:nrow(data.b),function(x) bestModelSelection_LM(x,data.b,covariables3 = coefficient,Diagnosis = diagnosis,
                                                                      VariableListOne,VariableListTwo)) 
  afile = paste('result_DE_null',b,'.rdata',sep='')
  save(result.b,file=afile)
}
sfStop()

setwd(null_dir) 
LM.NULL=c()
for(i in 1:B){
  print(i)
  aNull=get(load(paste("result_DE_null", i,".rdata", sep = "")))
  LM.NULL[[i]]<-aNull
}

unlistLM_null=unlist(LM.NULL)
names.pval=grep("lrt.pvalue", names(unlistLM_null))
pval.null=unname(as.numeric(unlistLM_null[names.pval]))

indLM=get(load(paste0(WD,"/output/DE_output/LM_true.RData")))
unlistLM_true=unlist(indLM)
names.pval=grep("lrt.pvalue", names(unlistLM_true))
pval=unname(as.numeric(unlistLM_true[names.pval]))


es.diagnosis = c()
st.diagnosis = c()
for( i in 1:length(indLM)){
  es.diagnosis[i] = indLM[[i]]$coefficients[2,1]
  st.diagnosis[i] = indLM[[i]]$coefficients[2,2]
}

p.pool = c(pval,pval.null)
p.corr = rank(p.pool)[1:length(pval)]/length(p.pool)
res=cbind(es.diagnosis,st.diagnosis, pval, p.corr)
genes=row.names(data)[1:100]
row.names(res) = make.names(genes, unique = TRUE)
colnames(res) = c("coefficient","sd","obs.p","corrected.p")

res=data.frame(res)
res$p.BH=p.adjust(res$corrected.p, "BH")
res_Sorted=res[order(res$corrected.p, decreasing = FALSE),]
head(res_Sorted)
write.csv(res_Sorted,paste0(WD,"/output/DE_output/DE_psychosisVScontrol.csv"))

