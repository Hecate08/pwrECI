#############################
# function to calculate power in a range




ECI_power <- function(alphaU = 0.05,
                      control_data_sd1 = 1,
                      control_data_sd2 = 1,
                      case_data_sd1 = 1,
                      case_data_sd2 = 1,
                      meanC1 = 3,
                      meanC2 = 3,
                      meanDiff1 = 1,
                      meanDiff2 = 1,
                      seed = 646,
                      ngenes = 1000,
                      sizeG = c(10,20),
                      updateProgress = NULL,
                      unbalanced = FALSE,
                      progressMonitor= NULL){

  print("ECI power calculation starts")

  # function
  sim_expr <- function(ncon,ncase,ngenes, meanC, meanDiff, control_data_sd, case_data_sd){
    case <- matrix(NA, ncol = ncase, nrow = ngenes)
    control <- matrix(NA, ncol = ncon, nrow = ngenes)
    for(i in 1:ngenes){
      # n random variable from truncated normal distribution for control
      control[i,] = truncnorm::rtruncnorm(ncon,a=0,b=Inf,mean = meanC, sd = control_data_sd)
      # n random variable from truncated normal distribution for case
      case[i,] = truncnorm::rtruncnorm(ncase,a=0,b=Inf, mean = meanC + meanDiff, sd = case_data_sd)
    }
    sample <- cbind(control,case)
    colnames(sample) <- seq(1,(ncon+ncase),1)
    rownames(sample) <- seq(1,ngenes,1)
    return(sample)
  }
  #print("loaded function")

  # variables
  if(unbalanced){
    print("unbalanced")
  } else {
    print("balanced")
  }
  print(sizeG)
  N <- 0
  if(unbalanced){
    N <- dim(sizeG)[1]
  } else {
    N <- length(sizeG)
  }
  #print("set values for balanced/unbalanced")

  if (is.function(updateProgress)) {
    text <- paste0(0,"/",N)
    updateProgress(detail = text)
  }

  start_time1 <- Sys.time()
  # main part
  time <- 0
  Power <- c()
  for(g in 1:N){

    start_time <- Sys.time()

    #############################
    # parameters
    set.seed(seed)

    # group information
    ncon1 = ncase1 = ncon2 = ncase2 = 0
    if(unbalanced){
      ncon1 = ncase1 = sizeG[g,1]
      ncon2 = ncase2 = sizeG[g,2]
    } else {
      ncon1 = ncase1 = ncon2 = ncase2 = sizeG[g]
    }
    ID1 <- seq(1,(ncon1+ncase1),1)
    ID2 <- seq(1,(ncon2+ncase2),1)
    Type1 <- c(rep("normal",ncon1),rep("tumor",ncase1))
    Type2 <- c(rep("normal",ncon2),rep("tumor",ncase2))
    targetsS1 <- data.frame(ID = ID1, Type = Type1)
    targetsS2 <- data.frame(ID = ID2, Type = Type2)

    #print(paste(ncon1,ncon2))

    #############################
    # main part
    #############################

    # simulate data
    sample1 <- sim_expr(ncon1,ncase1,ngenes, meanC1, meanDiff1, control_data_sd1, case_data_sd1)
    sample2 <- sim_expr(ncon2,ncase2,ngenes, meanC2, meanDiff2, control_data_sd2, case_data_sd2)


    # ECI permutation tests
    S1vsS2 <- ECIbootstrapTest(sample1,sample2,
                               targetsS1,targetsS2)

    #############################
    # power calculation
    #############################

    Power[g] <- sum(S1vsS2$p_value < alphaU)/dim(S1vsS2)[1]

    print(g)
    end_time <- Sys.time()
    if(g == 1) time <- end_time - start_time
    #if (rendered_by_shiny) shiny::incProgress(1/N)
    if(is.function(updateProgress)) {
      text <- paste0(g,"/",N," sim. done ","Estimated time: ",round(time*N,2)," min")
      updateProgress(value = 1/N, detail = text)
    }

    if(is.function(progressMonitor)) progressMonitor(g)
  }
  end_time1 <- Sys.time(); print(end_time1 - start_time1)
  return(Power)
}
