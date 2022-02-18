
###########
# differential Expression
diffExpr <- function(data, targets){
  # getting pval
  design <- model.matrix(~targets$Type)
  colnames(design) <- levels(targets$Type)

  fit <- lmFit(data, design)
  fit <- eBayes(fit)

  # getting sd
  beta <- fit$coefficients[,2]
  pval <- fit$p.value[,2]
  #sd <- ((sqrt(fit$s2.post)) * (fit$stdev.unscaled))[,2]

  gene_list <- data.frame(log2FC = beta, pval = pval, df = fit$df.total)

  return(gene_list)
}

###########
# ECI

bootstrap.pval <- function(data, ES){
  B = 1000
  N = sum(data < ES)
  M = sum(data < 0 )
  alpha1 = pnorm(-2*qnorm(N/B) + qnorm(M/B))
  if(ES < 0) alpha1 = 1 - alpha1
  alpha = 2*alpha1
  return(alpha)
}


# geneExpr1 <-  sample1
# geneExpr2 <-  sample2
# targets1 <- targetsS1
# targets2<- targetsS2
# alphaU = 0.05
ECIbootstrapTest <- function(geneExpr1,geneExpr2,targets1,targets2){

  # differential gene expression
  gene_list1 <- diffExpr(geneExpr1,targets1)
  gene_list2 <- diffExpr(geneExpr2,targets2)

  #get beta values and pvalues from both studies
  beta1 = gene_list1$log2FC
  beta2 = gene_list2$log2FC
  pval1 = gene_list1$pval
  pval2 = gene_list2$pval
  names(beta1) = rownames(gene_list1)
  names(beta2) = rownames(gene_list2)
  names(pval1) = rownames(gene_list1)
  names(pval2) = rownames(gene_list2)

  # get eci for all genes
  eci <- getECI(beta1,beta2,pval1,pval2)

  n <- 10
  len <- dim(gene_list1)[1]
  num1 <- dim(targets1)[1]
  num2 <- dim(targets2)[1]

  group1 = which(targets1[,2] == "tumor")
  group2 = which(targets1[,2] == "normal")
  group3 = which(targets2[,2] == "tumor")
  group4 = which(targets2[,2] == "normal")
  #print(paste(length(group1),length(group2),length(group3),length(group4)))

  # bootstrap
  bootstrap <- matrix(NA, ncol = n, nrow = len)
  for(i in 1:n){
    # if (i %% 10 == 0){
    #  print(i)
    # }
    #Bootstrap sampling

    geneExpr1new <- matrix(NA, ncol = num1, nrow = len)
    geneExpr2new <- matrix(NA, ncol = num2, nrow = len)
    for(j in 1:len){
      new1 = c(sample(group1,length(group1),replace = TRUE),sample(group2,length(group2),replace = TRUE))
      new2 = c(sample(group3,length(group3),replace = TRUE),sample(group4,length(group4),replace = TRUE))

      geneExpr1new[j,] = geneExpr1[j,new1]
      geneExpr2new[j,] = geneExpr2[j,new2]
    }

    #diff gene expression
    gene_list_B1 = diffExpr(geneExpr1new,targets1)
    gene_list_B2 = diffExpr(geneExpr2new,targets2)

    #ECI
    beta_B1 = gene_list_B1$log2FC
    beta_B2 = gene_list_B2$log2FC
    pval_B1 = gene_list_B1$pval
    pval_B2 = gene_list_B2$pval
    names(beta_B1) = rownames(gene_list_B1)
    names(beta_B2) = rownames(gene_list_B2)
    names(pval_B1) = rownames(gene_list_B1)
    names(pval_B2) = rownames(gene_list_B2)

    #get eci for all genes
    bootstrap[,i] <- getECI(beta_B1,beta_B2,pval_B1,pval_B2)
  }

  #confidence intervals
  CI <- matrix(NA,ncol = 2, nrow = len)
  pval <- c()
  for(i in 1:len){
    pval[i] <- bootstrap.pval(bootstrap[i,], eci[i])
  }

  result <- data.frame(ECI = eci, p_value = pval)
  rownames(result) <- rownames(gene_list1)

  return(result)
}

