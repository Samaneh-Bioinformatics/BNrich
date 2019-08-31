#' @title Download necessary files to start BNrich
#' @description Download necessary files to start BNrich
#' @importFrom utils download.file
#' @return A list contain mapkG, PathName_final and pathway.id
#' mapkG is a list contains imported 187 preprocessed signaling pathways
#' PathName_final is a data.frame includes names and IDs of all 187 pathways
#' pathway.id is a character vector of pathways IDs
#'
#' @export
#'
#' @examples
#' start_files()

start_files <- function(){
  destfile <- "./R/BNrich-start.rda"
  fileURL <-
    "https://github.com/Samaneh-Bioinformatics/BNrich-RData/raw/master/BNrich-start.rda"
  if (!file.exists(destfile)) {
    print("please be patient, the files are downloading...")
    #download.file(fileURL ,destfile,method="auto")
    files <- c(fileURL)
    oldw <- getOption("warn")
    options(warn = -1)
    for (file in files) {
      tryCatch(download.file(file, destfile, method="auto"),
      error = function(e) print(paste(file, 'did not work out')))
    }
    options(warn = oldw)
   }
  return(destfile)
 }


#' @title Simplification networks -- applied to unifying nodes
#' @description Unifying nodes based imported signaling pathways and GE data
#' @param dataH A data frame contains (healthy) control objects data
#' @param dataD A data frame contains disease objects data
#' @param mapkG A list contains imported 187 signaling pathways
#' @param pathway.id A vector contains 187 KEEG pathway IDs
#' @importFrom graph nodes removeNode edgeMatrix
#' @return A list contain data_h,data_d,mapkG1 and pathway.id1
#'
#' @export
#'
#' @examples
#' unify_path()

unify_path <- function(dataH=dataH,dataD=dataD,mapkG = mapkG,pathway.id = pathway.id){
  NOD <- lapply(mapkG,nodes)
  NOD <- lapply(NOD,as.vector)
  data_h <- list()
  data_d <- list()
  Diff <- list()
  for (i in seq_along(mapkG)) {
    data_h[[i]] <- matrix()
    data_h[[i]] <- subset(dataH,rownames(dataH) %in% NOD[[i]])
    data_h[[i]] <- as.data.frame(t(data_h[[i]]))
    rownames(data_h[[i]]) <- NULL
    Diff[[i]] <- setdiff(NOD[[i]],colnames(data_h[[i]]))
    mapkG[[i]]=removeNode(Diff[[i]],mapkG[[i]])
    data_d[[i]] <- matrix()
    data_d[[i]] <- subset(dataD,rownames(dataD) %in% NOD[[i]])
    data_d[[i]] <- as.data.frame(t(data_d[[i]]))
    rownames(data_d[[i]]) <- NULL
  }

  mapkG1 <- mapkG
  pathway.id1 <- pathway.id
  for (i in length(mapkG1):1) {
    if(ncol(edgeMatrix(mapkG1[[i]]))<5){
      mapkG1 <- mapkG1[-i]
      data_h <- data_h[-i]
      data_d <- data_d[-i]
      pathway.id1 <- pathway.id1[-i]
    }}
  for (i in seq_along(mapkG1)) {
    data_h[[i]] <- data_h[[i]][,order(names(data_h[[i]]))]
    data_d[[i]] <- data_d[[i]][,order(names(data_d[[i]]))]
  }
  unify_results <- list("data_h"= data_h,"data_d"=data_d,
                        "mapkG1"= mapkG1,"pathway.id1" =pathway.id1)
  return(unify_results)
}

#' @title Construct Bayesian networks structures
#' @description Construct BNs structures using unified signaling pathways
#' @param mapkG1 A list contains unified signaling pathways
#' @importFrom bnlearn as.bn
#' @return A list contains Bayesian networks structures
#'
#' @export
#'
#' @examples
#' BN_struct(mapkG1)

BN_struct <- function(mapkG1){
  BN=list()
  BN <- lapply(mapkG1,as.bn)
  return(BN)
}

#' @title LASSO regression
#' @description LASSO regression – second step of simplification of BNs structures
#' @param BN A list of Bayesian networks achieved by BN_struct function
#' @param data_h A list contains data frames related to control objects
#' @param data_d A list contains data frames related to disease objects
#' @importFrom bnlearn nodes drop.arc
#' @importFrom glmnet cv.glmnet
#' @importFrom stats coef
#' @return A list contains two lists.BN_H and BN_D are simplified BNs
#'
#' @export
#'
#' @examples
#' LASSO_BN(BN,data_h,data_d)

LASSO_BN <- function(BN,data_h,data_d){
  oldw <- getOption("warn")
  options(warn = -1)

  BN_H <- BN
  BN_D <- BN
  for(k in seq_along(BN)){
    for(i in seq_along(BN[[k]]$nodes)) {
      if(length((BN_H[[k]]$nodes)[[i]][3]$parents) > 1) {
        Parents_h <- as.matrix(data_h[[k]][BN_H[[k]]$nodes[[i]][3]$parents])
        respons_h <- as.matrix(data_h[[k]][bnlearn::nodes(BN_H[[k]])[i]])
        Parents_d <- as.matrix(data_d[[k]][BN_D[[k]]$nodes[[i]][3]$parents])
        respons_d <- as.matrix(data_d[[k]][bnlearn::nodes(BN_D[[k]])[i]])
        cvfit_h = cv.glmnet(Parents_h,respons_h,grouped=FALSE,parallel=TRUE)
        cvfit_d = cv.glmnet(Parents_d,respons_d,grouped=FALSE,parallel=TRUE)
        for(j in nrow(coef(cvfit_h,s=cvfit_h$lambda.min)):1){
          if((coef(cvfit_h, s=cvfit_h$lambda.min)[j] == 0)&& (coef(cvfit_d, s=cvfit_d$lambda.min)[j] == 0)){
            BN_H[[k]] <- drop.arc(BN_H[[k]],from =row.names(coef(cvfit_h,s=cvfit_h$lambda.min))[j],to=bnlearn::nodes(BN_H[[k]])[i])
            BN_D[[k]] <- drop.arc(BN_D[[k]],from =row.names(coef(cvfit_d,s=cvfit_d$lambda.min))[j],to=bnlearn::nodes(BN_D[[k]])[i])
          }}}}}
  LASSO_results <- list("BN_H" = BN_H, "BN_D" = BN_D)

  options(warn = oldw)
  return(LASSO_results)}

#' @title Estimate parameters of BNs in control and disease states
#' @description Estimate parameters of BNs in control and disease states
#' @param BN_H A list contains simplified BNs structures for control objects
#' @param BN_D A list contains simplified BNs structures for disease objects
#' @param data_h A list contains data frames related to control objects for any BN
#' @param data_d A list contains data frames related to disease objects for any BN
#' @importFrom bnlearn bn.fit
#' @importFrom stats coef
#' @return A listcontains four lists BN_h, BN_d, coef_h and coef_d
#'
#' @export
#'
#' @examples
#' esti_par(BN_H,BN_D,data_h,data_d)


esti_par <- function(BN_H,BN_D,data_h,data_d){
  BN_h <- list()
  BN_d <- list()
  coef_h <- list()
  coef_d <- list()
  for (i in seq_along(BN_H)) {
    BN_h[[i]] <- bn.fit(BN_H[[i]],data_h[[i]])
    BN_d[[i]] <- bn.fit(BN_D[[i]],data_d[[i]])
    coef_h[[i]] <- coef(BN_h[[i]])
    coef_d[[i]] <- coef(BN_d[[i]])
  }
  esti_results <- list("BN_h"=BN_h,"BN_d"=BN_d,"coef_h"=coef_h,"coef_d"=coef_d)
  return(esti_results)
}

#' @title Estimate variance-covariance matrixes for any parameters of BNs
#' @description Estimate variance-covariance matrixes for any parameters of
#' @param data_h A list contains data frames related to control objects for any BN
#' @param coef_h A lists of parameters of BN_h achieved
#' @param BN_h A list of BNs learned by control objects data
#' @param data_d A list contains data frames related to disease objects for any BN
#' @param coef_d A lists of parameters of BN_d
#' @param BN_d A list of BNs learned by disease objects data
#' @importFrom corpcor pseudoinverse
#' @return A listcontains two lists var_mat_Bh and var_mat_Bd
#'
#' @export
#'
#' @examples
#' var_mat(data_h,coef_h,BN_h,data_d,coef_d,BN_d)


var_mat <- function(data_h,coef_h,BN_h,data_d,coef_d,BN_d){
  X_h <- list()
  for (k in seq_along(data_h)) {
    X_h[[k]] <- list()
    for(i in seq_along(colnames(data_h[[k]]))){
      X_h[[k]][[i]] <- matrix(nrow = nrow(data_h[[k]]),ncol = length(coef_h[[k]][[i]]))
      X_h[[k]][[i]][,1] <- as.vector(rep(1,nrow(data_h[[k]])))
      if(1<length(coef_h[[k]][[i]])){
        for (j in 2:length(coef_h[[k]][[i]])){
          X_h[[k]][[i]][,j] <- data_h[[k]][,names(coef_h[[k]][[i]][j])]}
      }}}
  var_mat_Bh <- list()
  for (k in seq_along(data_h)) {
    var_mat_Bh[[k]] <- list()
    for(i in seq_along(colnames(data_h[[k]]))){
      var_mat_Bh[[k]][[i]] <- BN_h[[k]][[i]]$sd*pseudoinverse((t(X_h[[k]][[i]]))%*%(X_h[[k]][[i]]))}}
  X_d <- list()
  for (k in seq_along(data_d)) {
    X_d[[k]] <- list()
    for(i in seq_along(colnames(data_d[[k]]))){
      X_d[[k]][[i]] <- matrix(nrow = nrow(data_d[[k]]),ncol = length(coef_d[[k]][[i]]))
      X_d[[k]][[i]][,1] <- as.vector(rep(1,nrow(data_d[[k]])))
      if(1<length(coef_d[[k]][[i]])){
        for (j in 2:length(coef_d[[k]][[i]])){
          X_d[[k]][[i]][,j] <- data_d[[k]][,names(coef_d[[k]][[i]][j])]}
      }}}
  var_mat_Bd <- list()
  for (k in seq_along(data_d)) {
    var_mat_Bd[[k]] <- list()
    for(i in seq_along(colnames(data_d[[k]]))){
      var_mat_Bd[[k]][[i]] <- BN_d[[k]][[i]]$sd*pseudoinverse((t(X_d[[k]][[i]]))%*%(X_d[[k]][[i]]))}}
  var_mat_results<- list("var_mat_Bh"=var_mat_Bh,"var_mat_Bd"=var_mat_Bd)
  return(var_mat_results)
}


#' @title Testing the equality regression coefficients
#' @description t-test for equality the corresponging parameters in any BN
#' @param data_h A list contains data frames related to control objects for any BN
#' @param coef_h A list contains parameters of BN_h
#' @param BN_h A list contains BNs learned by control objects data
#' @param data_d A list contains data frames related to disease objects for any BN
#' @param coef_d A list contains parameters of BN_d
#' @param BN_d A list contains BNs learned by disease objects data
#' @param var_mat_Bh A list contains covariance matrixes for any node of BN_h
#' @param var_mat_Bd A list contains covariance matrixes for any node of BN_d
#' @param pathway.id1 A vector contains modified KEEG pathway IDs
#' @importFrom stats pt p.adjust complete.cases
#' @return A data frame contains T-test results for all parameters in final BNs
#' @export
#'
#' @examples
#' parm_Ttest(data_h,coef_h,BN_h,data_d,coef_d,BN_d,var_mat_Bh,var_mat_Bd, pathway.id1)

parm_Ttest <- function(data_h,coef_h,BN_h,data_d,coef_d,BN_d,var_mat_Bh,var_mat_Bd, pathway.id1){
  t.test2 <- function(m1,m2,s1,s2,n1,n2){
    se <- sqrt((s1^2/(n1)) + (s2^2/(n2)))
    # welch-satterthwaite df
    df <- se^4/((((s1^2/n1)^2)/(n1-1))+(((s2^2/n2)^2)/(n2-1)))
    t <- (m1-m2)/se
    dat <- c(m1-m2, se, t, 2*pt(-abs(t),df))
    names(dat) <- c("Difference of means", "Std Error", "t", "p-value")
    return(dat)
  }
  Ttest_results <- data.frame()
  options(stringsAsFactors = FALSE)
  for (k in seq_along(data_h)) {
    n1= nrow(data_d[[k]])
    n2= nrow(data_h[[k]])
    for(i in seq_along(colnames(data_h[[k]]))){
      To <- names(BN_h[[k]][i])
      for (j in seq_along(coef_h[[k]][[i]])){
        m1=coef_d[[k]][[i]][j]
        m2=coef_h[[k]][[i]][j]
        s1= sqrt(var_mat_Bd[[k]][[i]][j,j])
        s2= sqrt(var_mat_Bh[[k]][[i]][j,j])
        if(j>1){
          From <- names(m2)
        }
        else {
          From <- "intercept"}
        T <- t.test2(m1,m2,s1,s2,n1,n2)
        Pval <- T["p-value"]
        Ttest_results <- rbind(Ttest_results,c(From, To,k,pathway.id1[k],Pval,m1,m2))
      }}}
  colnames(Ttest_results) <- c("From","To","pathway.number","pathwayID","Pval","coefficient in disease","coefficient in control")
  Ttest_results$fdr <- p.adjust(Ttest_results$Pval,method = "fdr")
  Ttest_results$pathway.number <-as.numeric(Ttest_results$pathway.number)
  Ttest_results <- Ttest_results[complete.cases(Ttest_results),]
  return(Ttest_results)
}

#' @title Analysis of significant final BNs
#' @description Fisher's exact test applied to PEA on final BNs
#'
#' @param Ttest_results A data frame contains T-test results for all parameters
#' @param fdr.value A numeric threshold to determine significant parameters
#' @param pathway.id1 A vector contains modified KEEG pathway IDs
#' @param PathName_final A data frame contains is IDs and names of KEEG pathways
#'
#' @importFrom stats fisher.test p.adjust
#' @return A data frame contains fisher test results for any final pathways
#'
#' @export
#'
#' @examples
#' BNrich(Ttest_results,pathway.id1,PathName_final)

BNrich <- function(Ttest_results,fdr.value = 0.05,pathway.id1,PathName_final){
  Ttest_results <- Ttest_results[order(Ttest_results$pathway.number),]
  BNrich_results=data.frame()
  b=1
  for (d in seq_len(max(Ttest_results$pathway.number))) {
    a=0
    f=1
    while(b<=nrow(Ttest_results) & Ttest_results$pathway.number[b]==d)
    {
      if(Ttest_results$fdr[b]<fdr.value) {a=a+1}
      f=f+1
      b=b+1
    }
    BNrich_results[d,1]=f-1
    BNrich_results[d,2]=a
  }
  N=sum(BNrich_results[,1])
  M=sum(BNrich_results[,2])

  for (d in seq_len(max(Ttest_results$pathway.number))){
    df <- matrix(c(BNrich_results[d,2],(M-BNrich_results[d,2]),(BNrich_results[d,1]-BNrich_results[d,2])
                   ,(N-M-BNrich_results[d,1]+BNrich_results[d,2])) ,nrow = 2)
    fisher <- fisher.test(df)
    BNrich_results[d,3]=fisher$p.value
    BNrich_results[d,4]=pathway.id1[d]
  }
  colnames(BNrich_results) <- c("nk","mk","p.value","ID")
  BNrich_results$fdr <- p.adjust(BNrich_results$p.value, method = "fdr")
  BNrich_results <- BNrich_results[,c("ID","p.value","fdr")]
  BNrich_results <- merge(BNrich_results, PathName_final, by = "ID")
  colnames(BNrich_results) <- c("pathwayID","p.value","fdr","pathway.number","Name")
  BNrich_results <- BNrich_results[order(BNrich_results$fdr),]
  return(BNrich_results)
}
