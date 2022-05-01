bestModelSelection_LM <- function(aGeneIndex,expressionMatrix,covariables3,Diagnosis,VariableListOne,VariableListTwo){
  curGeneExpression = as.numeric(expressionMatrix[aGeneIndex,])
  
  ## mixed model with adjusting for covariates
  basicFormula = update(curGeneExpression ~ Diagnosis,.~.)
  basicFormula1 = update(basicFormula,.~.+unlist(covariables3[,aVariable]))
  basicFormula2 = update(basicFormula,.~.+unlist(covariables3[,aVariable1]) + unlist(covariables3[,aVariable2]))
  
  bestNullFormula = update(basicFormula,.~.-Diagnosis)
  fm0 = lm(basicFormula)
  bestfm = fm0
  bestBIC = BIC(fm0)
  bestBICVariables = NULL
  
  ## adjusting for one covariate
  for(i in 1:length(VariableListOne)){
    aVariable = VariableListOne[i]
    fm1 = lm(basicFormula1)
    fm1BIC = tryCatch(BIC(fm1), error=function(e) e)
    
    if(!sum(class(fm1BIC)=="numeric")){
      cat(i)
      cat("th, one variable")
      cat("\n")
      cat("gene index: ")
      cat(aGeneIndex)
      cat("\n")
      next
    }
    
    if(fm1BIC<bestBIC){
      bestBIC = fm1BIC
      bestBICVariables = aVariable
      bestfm = fm1
    }
  }
  
  ## mixed model with adjusting for two covariates
  for(i in 1:dim(VariableListTwo)[2]){
    aVariable1 = VariableListTwo[1,i]
    aVariable2 = VariableListTwo[2,i]
    fm2 = lm(basicFormula2)
    
    fm2BIC = tryCatch(BIC(fm2), error=function(e) e)
    if(!sum(class(fm2BIC)=="numeric")){
      cat(i)
      cat("th, two variable")
      cat("\n")
      cat("gene index: ")
      cat(aGeneIndex)
      cat("\n")
      next
    }
    if(fm2BIC<bestBIC){
      bestBIC = fm2BIC
      bestBICVariables = VariableListTwo[,i]
      bestfm = fm2 	
    }	
  }
  
  if(length(bestBICVariables)==0){
    bestNullFormula = update(bestNullFormula,.~.)		
  } else if(length(bestBICVariables)==1){
    bestNullFormula = update(bestNullFormula,.~.+unlist(covariables3[,bestBICVariables]))		
  } else if(length(bestBICVariables)==2){
    bestNullFormula = update(bestNullFormula,.~.+unlist(covariables3[,bestBICVariables[1]])+unlist(covariables3[,bestBICVariables[2]]))		
  } else {
    stop('length of best BIC variables error')
  }
  
  
  
  bestNullfm = lm(bestNullFormula)
  lrt.pvalue = anova(bestfm,bestNullfm)[2,'Pr(>F)']
  
  coefficients = summary(bestfm)$coefficients
  
  
  # if(is.null(bestBICVariables)){
  #  residuals = curGeneExpression 
  #} else if(length(bestBICVariables)==1){
  # residuals = curGeneExpression - covariables3[,bestBICVariables]*coefficients['covariables3[, aVariable]','Estimate']
  #} else if(length(bestBICVariables)==2){
  #residuals = curGeneExpression - as.matrix(covariables3[,bestBICVariables])%*%as.matrix(coefficients[c('covariables3[, aVariable1]','covariables3[, aVariable2]'),'Estimate'])
  #}	
  result = list(coefficients=coefficients,bestBICVariables=bestBICVariables,lrt.pvalue=lrt.pvalue) #residuals=residuals
  return(result)
}

