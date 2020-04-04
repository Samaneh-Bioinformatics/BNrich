#' @title Download data file
#' @description Download necessary data file to start BNrich
#' @param destfile A directory in user's own computer for save preprocessed data file
#' @param verbose A logical argument to show verbose results
#' @importFrom utils download.file
#' @return A list contain mapkG, PathName_final and pathway.id. The mapkG is a list contains imported 187 preprocessed signaling pathways, PathName_final is a data.frame includes names and IDs of all 187 pathways and pathway.id is a character vector of pathways IDs
#'
#' @export
#'
#' @examples
#' \donttest{
#' destfile = tempfile("files", fileext = ".rda")
#' files <- fetch_data_file()
#' load(destfile)
#' }
fetch_data_file <- function(destfile, verbose = FALSE){
  fileURL <-
    "https://github.com/Samaneh-Bioinformatics/RData/raw/master/BNrich-start.rda"
  if (!file.exists(destfile)) {
    if (verbose) {
      print("please be patient, the files are downloading...")
    }
    files <- c(fileURL)
    for (file in files) {
      tryCatch(utils::download.file(file, destfile, method="libcurl"),
               error = function(e) print(paste(file, 'did not work out')))
    }}}

#' @title Simplification networks -- applied to unifying nodes
#' @description Unifying nodes based imported signaling pathways and GE data
#' @param dataH A data frame contains (healthy) control objects data
#' @param dataD A data frame contains disease objects data
#' @param MapkG A list contains imported 187 signaling pathways
#' @param Pathway.id A vector contains 187 KEEG pathway IDs
#' @importFrom graph nodes removeNode edgeMatrix
#' @return A list contain data_h,data_d,mapkG1 and pathway.id1
#'
#' @export
#'
#' @examples
#' #All the 187 preprocessed signaling pathways can be entered in analysis by fetch_data_file().
#' #But here you enter a subset of those pathways to see how this package works.
#' files <- system.file("extdata", "test_files_to_start.RData", package = "BNrich", mustWork = TRUE)
#' load(files)
#' Data <- system.file("extdata", "Test_DATA.RData", package = "BNrich", mustWork = TRUE)
#' load(Data)
#' uni_Result <- unify_path(dataH, dataD, MapkG = sub_mapkG, Pathway.id = path.id)

unify_path <- function(dataH, dataD, MapkG, Pathway.id){
  NOD <- lapply(MapkG, graph::nodes)
  NOD <- lapply(NOD, as.vector)
  data_h <- list()
  data_d <- list()
  Diff <- list()
  for (i in seq_along(MapkG)) {
    data_h[[i]] <- matrix()
    data_h[[i]] <- subset(dataH, rownames(dataH) %in% NOD[[i]])
    data_h[[i]] <- as.data.frame(t(data_h[[i]]))
    rownames(data_h[[i]]) <- NULL
    Diff[[i]] <- setdiff(NOD[[i]], colnames(data_h[[i]]))
    MapkG[[i]]=removeNode(Diff[[i]], MapkG[[i]])
    data_d[[i]] <- matrix()
    data_d[[i]] <- subset(dataD, rownames(dataD) %in% NOD[[i]])
    data_d[[i]] <- as.data.frame(t(data_d[[i]]))
    rownames(data_d[[i]]) <- NULL
  }
  mapkG1 <- MapkG
  pathway.id1 <- Pathway.id
  for (i in length(mapkG1):1) {
    if(ncol(edgeMatrix(mapkG1[[i]]))<5){
      mapkG1 <- mapkG1[-i]
      data_h <- data_h[-i]
      data_d <- data_d[-i]
      pathway.id1 <- pathway.id1[-i]
    }}
  for (i in seq_along(mapkG1)) {
    data_h[[i]] <- data_h[[i]][, order(names(data_h[[i]]))]
    data_d[[i]] <- data_d[[i]][, order(names(data_d[[i]]))]
  }
  unify_results <- list("data_h"= data_h,"data_d"=data_d,"mapkG1"= mapkG1,"pathway.id1" =pathway.id1)
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
#' #All the 187 preprocessed signaling pathways can be entered in analysis by fetch_data_file().
#' #But here you enter a subset of those pathways to see how this package works.
#' files <- system.file("extdata", "test_files_to_start.RData", package = "BNrich", mustWork = TRUE)
#' load(files)
#' Data <- system.file("extdata", "Test_DATA.RData", package = "BNrich", mustWork = TRUE)
#' load(Data)
#' uni_Result <- unify_path(dataH, dataD, MapkG = sub_mapkG, Pathway.id = path.id)
#' M1 <- uni_Result$mapkG1
#' BN <- BN_struct(M1)
BN_struct <- function(mapkG1){
  BN <- list()
  BN <- lapply(mapkG1, bnlearn::as.bn)
  return(BN)
}

#' @title LASSO regression
#' @description LASSO regression – second step of simplification of BNs structures
#' @param BN A list of Bayesian networks achieved by BN_struct function
#' @param data_h A list contains data frames related to control objects
#' @param data_d A list contains data frames related to disease objects
#' @importFrom bnlearn drop.arc
#' @importFrom glmnet cv.glmnet
#' @importFrom stats coef
#' @return A list contains two lists.BN_h and BN_d are simplified BNs
#'
#' @export
#'
#' @examples
#' #All the 187 preprocessed signaling pathways can be entered in analysis by fetch_data_file().
#' #But here you enter a subset of those pathways to see how this package works.
#' files <- system.file("extdata", "test_files_to_start.RData", package = "BNrich", mustWork = TRUE)
#' load(files)
#' Data <- system.file("extdata", "Test_DATA.RData", package = "BNrich", mustWork = TRUE)
#' load(Data)
#' uni_Result <- unify_path(dataH, dataD, MapkG = sub_mapkG, Pathway.id = path.id)
#' M1 <- uni_Result$mapkG1
#' BN <- BN_struct(M1)
#' data_h1 <- uni_Result$data_h
#' data_d1 <- uni_Result$data_d
#' LASSO_Result <- LASSO_BN(BN = BN , data_h = data_h1 , data_d = data_d1)
LASSO_BN <- function(BN, data_h, data_d){
  BN_h <- BN
  BN_d <- BN
  for(k in seq_along(BN)){
    for(i in seq_along(BN[[k]]$nodes)) {
      if(length((BN_h[[k]]$nodes)[[i]][3]$parents) > 1) {
        Parents_h <- as.matrix(data_h[[k]][BN_h[[k]]$nodes[[i]][3]$parents])
        respons_h <- as.matrix(data_h[[k]][bnlearn::nodes(BN_h[[k]])[i]])
        Parents_d <- as.matrix(data_d[[k]][BN_d[[k]]$nodes[[i]][3]$parents])
        respons_d <- as.matrix(data_d[[k]][bnlearn::nodes(BN_d[[k]])[i]])
        cvfit_h = cv.glmnet(Parents_h, respons_h, grouped = FALSE)
        cvfit_d = cv.glmnet(Parents_d, respons_d, grouped = FALSE)
        for(j in nrow(coef(cvfit_h,s=cvfit_h$lambda.min)):1){
          if((coef(cvfit_h, s=cvfit_h$lambda.min)[j] == 0)&& (coef(cvfit_d, s=cvfit_d$lambda.min)[j] == 0)){
            BN_h[[k]] <- drop.arc(BN_h[[k]],from = row.names(coef(cvfit_h, s=cvfit_h$lambda.min))[j],to = bnlearn::nodes(BN_h[[k]])[i])
            BN_d[[k]] <- drop.arc(BN_d[[k]],from = row.names(coef(cvfit_d, s=cvfit_d$lambda.min))[j],to = bnlearn::nodes(BN_d[[k]])[i])
          }}}}}
  LASSO_results <- list("BN_h" = BN_h, "BN_d" = BN_d)
  return(LASSO_results)}

#' @title Estimate parameters of BNs in control and disease states
#' @description Estimate parameters of BNs in control and disease states
#' @param BN_H A list contains simplified BNs structures for control objects
#' @param BN_D A list contains simplified BNs structures for disease objects
#' @param data_h A list contains data frames related to control objects for any BN
#' @param data_d A list contains data frames related to disease objects for any BN
#' @importFrom bnlearn bn.fit
#' @importFrom stats coef
#' @return A listcontains four lists BNs_h, BNs_d, coef_h and coef_d
#'
#' @export
#'
#' @examples
#' #All the 187 preprocessed signaling pathways can be entered in analysis by fetch_data_file().
#' #But here you enter a subset of those pathways to see how this package works.
#' files <- system.file("extdata", "test_files_to_start.RData", package = "BNrich", mustWork = TRUE)
#' load(files)
#' Data <- system.file("extdata", "Test_DATA.RData", package = "BNrich", mustWork = TRUE)
#' load(Data)
#' uni_Result <- unify_path(dataH, dataD, MapkG = sub_mapkG, Pathway.id = path.id)
#' M1 <- uni_Result$mapkG1
#' BN <- BN_struct(M1)
#' data_h1 <- uni_Result$data_h
#' data_d1 <- uni_Result$data_d
#' LASSO_Result <- LASSO_BN(BN = BN , data_h = data_h1 , data_d = data_d1)
#' BN_h1 <- LASSO_Result$BN_h
#' BN_d1 <- LASSO_Result$BN_d
#' esti_result <- esti_par(BN_H = BN_h1, BN_D = BN_d1, data_h = data_h1, data_d = data_d1)
esti_par <- function(BN_H, BN_D, data_h, data_d){
  BNs_h <- list()
  BNs_d <- list()
  coef_h <- list()
  coef_d <- list()
  for (i in seq_along(BN_H)) {
    BNs_h[[i]] <- bn.fit(BN_H[[i]], data_h[[i]])
    BNs_d[[i]] <- bn.fit(BN_D[[i]], data_d[[i]])
    coef_h[[i]] <- coef(BNs_h[[i]])
    coef_d[[i]] <- coef(BNs_d[[i]])
  }
  esti_results <- list("BNs_h"=BNs_h,"BNs_d"=BNs_d,"coef_h"=coef_h,"coef_d"=coef_d)
  return(esti_results)
}

#' @title Estimate variance-covariance matrixes for any parameters of BNs
#' @description Estimate variance-covariance matrixes for any parameters of
#' @param Data_h A list contains data frames related to control objects for any BN
#' @param coef_H A lists of parameters of BN_h achieved
#' @param BNs_H A list of BNs learned by control objects data
#' @param Data_d A list contains data frames related to disease objects for any BN
#' @param coef_D A lists of parameters of BN_d
#' @param BNs_D A list of BNs learned by disease objects data
#' @importFrom corpcor pseudoinverse
#' @return A listcontains two lists var_mat_Bh and var_mat_Bd
#'
#' @export
#'
#' @examples
#' #All the 187 preprocessed signaling pathways can be entered in analysis by fetch_data_file().
#' #But here you enter a subset of those pathways to see how this package works.
#' files <- system.file("extdata", "test_files_to_start.RData", package = "BNrich", mustWork = TRUE)
#' load(files)
#' Data <- system.file("extdata", "Test_DATA.RData", package = "BNrich", mustWork = TRUE)
#' load(Data)
#' uni_Result <- unify_path(dataH, dataD, MapkG = sub_mapkG, Pathway.id = path.id)
#' M1 <- uni_Result$mapkG1
#' BN <- BN_struct(M1)
#' data_h1 <- uni_Result$data_h
#' data_d1 <- uni_Result$data_d
#' LASSO_Result <- LASSO_BN(BN = BN , data_h = data_h1 , data_d = data_d1)
#' BN_h1 <- LASSO_Result$BN_h
#' BN_d1 <- LASSO_Result$BN_d
#' esti_result <- esti_par(BN_H = BN_h1, BN_D = BN_d1, data_h = data_h1, data_d = data_d1)
#' BNs_H <- esti_result$BNs_h
#' BNs_D <- esti_result$BNs_d
#' coef_h <- esti_result$coef_h
#' coef_d <- esti_result$coef_d
#' var_result <- var_mat(data_h1, coef_h,  BNs_H, data_d1, coef_d, BNs_D)
var_mat <- function(Data_h, coef_H, BNs_H, Data_d, coef_D, BNs_D){
  X_h <- list()
  for (k in seq_along(Data_h)) {
    X_h[[k]] <- list()
    for(i in seq_along(colnames(Data_h[[k]]))){
      X_h[[k]][[i]] <- matrix(nrow = nrow(Data_h[[k]]),ncol = length(coef_H[[k]][[i]]))
      X_h[[k]][[i]][,1] <- as.vector(rep(1,nrow(Data_h[[k]])))
      if(1<length(coef_H[[k]][[i]])){
        for (j in 2:length(coef_H[[k]][[i]])){
          X_h[[k]][[i]][,j] <- Data_h[[k]][,names(coef_H[[k]][[i]][j])]}
      }}}
  var_mat_Bh <- list()
  for (k in seq_along(Data_h)) {
    var_mat_Bh[[k]] <- list()
    for(i in seq_along(colnames(Data_h[[k]]))){
      var_mat_Bh[[k]][[i]] <- (BNs_H[[k]][[i]]$sd)^2*pseudoinverse((t(X_h[[k]][[i]]))%*%(X_h[[k]][[i]]))}}
  X_d <- list()
  for (k in seq_along(Data_h)) {
    X_d[[k]] <- list()
    for(i in seq_along(colnames(Data_h[[k]]))){
      X_d[[k]][[i]] <- matrix(nrow = nrow(Data_h[[k]]),ncol = length(coef_D[[k]][[i]]))
      X_d[[k]][[i]][,1] <- as.vector(rep(1,nrow(Data_h[[k]])))
      if(1<length(coef_D[[k]][[i]])){
        for (j in 2:length(coef_D[[k]][[i]])){
          X_d[[k]][[i]][,j] <- Data_h[[k]][,names(coef_D[[k]][[i]][j])]}
      }}}
  var_mat_Bd <- list()
  for (k in seq_along(Data_h)) {
    var_mat_Bd[[k]] <- list()
    for(i in seq_along(colnames(Data_h[[k]]))){
      var_mat_Bd[[k]][[i]] <- (BNs_D[[k]][[i]]$sd)^2*pseudoinverse((t(X_d[[k]][[i]]))%*%(X_d[[k]][[i]]))}}
  var_mat_results<- list("var_mat_Bh"=var_mat_Bh,"var_mat_Bd"=var_mat_Bd)
  return(var_mat_results)
}

#' @title Testing the equality regression coefficients
#' @description t-test for equality the corresponging parameters in any BN
#' @param Data_h A list contains data frames related to control objects for any BN
#' @param coef_H A list contains parameters of BN_h
#' @param BNs_H A list contains BNs learned by control objects data
#' @param Data_d A list contains data frames related to disease objects for any BN
#' @param coef_D A list contains parameters of BN_d
#' @param BNs_D A list contains BNs learned by disease objects data
#' @param Var_mat_Bh A list contains covariance matrixes for any node of BN_h
#' @param Var_mat_Bd A list contains covariance matrixes for any node of BN_d
#' @param Pathway.id1 A vector contains modified KEEG pathway IDs
#' @importFrom stats pt p.adjust complete.cases
#' @return A data frame contains T-test results for all parameters in final BNs
#'
#' @export
#'
#' @examples
#' #All the 187 preprocessed signaling pathways can be entered in analysis by fetch_data_file().
#' #But here you enter a subset of those pathways to see how this package works.
#' files <- system.file("extdata", "test_files_to_start.RData", package = "BNrich", mustWork = TRUE)
#' load(files)
#' Data <- system.file("extdata", "Test_DATA.RData", package = "BNrich", mustWork = TRUE)
#' load(Data)
#' uni_Result <- unify_path(dataH, dataD, MapkG = sub_mapkG, Pathway.id = path.id)
#' M1 <- uni_Result$mapkG1
#' BN <- BN_struct(M1)
#' data_h1 <- uni_Result$data_h
#' data_d1 <- uni_Result$data_d
#' LASSO_Result <- LASSO_BN(BN = BN , data_h = data_h1 , data_d = data_d1)
#' BN_h1 <- LASSO_Result$BN_h
#' BN_d1 <- LASSO_Result$BN_d
#' esti_result <- esti_par(BN_H = BN_h1, BN_D = BN_d1, data_h = data_h1, data_d = data_d1)
#' BNs_H <- esti_result$BNs_h
#' BNs_D <- esti_result$BNs_d
#' coef_h <- esti_result$coef_h
#' coef_d <- esti_result$coef_d
#' var_result <- var_mat(data_h1, coef_h,  BNs_H, data_d1, coef_d, BNs_D)
#' Var_H = var_result$var_mat_Bh
#' Var_D = var_result$var_mat_Bd
#' path.id1 <- uni_Result$pathway.id1
#' Ttest_result <- parm_Ttest(data_h1, coef_h, BNs_H, data_d1, coef_d, BNs_D, Var_H, Var_D, path.id1)
parm_Ttest <- function(Data_h, coef_H, BNs_H, Data_d, coef_D, BNs_D, Var_mat_Bh , Var_mat_Bd, Pathway.id1){
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
  for (k in seq_along(Data_h)) {
    n1= nrow(Data_d[[k]])
    n2= nrow(Data_h[[k]])
    for(i in seq_along(colnames(Data_h[[k]]))){
      To <- names(BNs_H[[k]][i])
      for (j in seq_along(coef_H[[k]][[i]])){
        m1=coef_D[[k]][[i]][j]
        m2=coef_H[[k]][[i]][j]
        s1= sqrt(Var_mat_Bd[[k]][[i]][j,j])
        s2= sqrt(Var_mat_Bh[[k]][[i]][j,j])
        if(j>1){
          From <- names(m2)
        }
        else {
          From <- "intercept"}
        T <- t.test2(m1,m2,s1,s2,n1,n2)
        Pval <- T["p-value"]
        Ttest_results <- rbind(Ttest_results,c(From, To,k,Pathway.id1[k],Pval,m1,m2))
      }}}
  colnames(Ttest_results) <- c("From","To","pathway.number","pathwayID","Pval","coefficient in disease","coefficient in control")
  Ttest_results$fdr <- p.adjust(Ttest_results$Pval,method = "fdr")
  Ttest_results$pathway.number <-as.numeric(Ttest_results$pathway.number)
  Ttest_results <- Ttest_results[complete.cases(Ttest_results),]
  return(Ttest_results)
}

#' @title Analysis of significant final BNs
#' @description Fisher's exact test applied to PEA on final BNs
#' @param Ttest_Results A data frame contains T-test results for all parameters
#' @param Pathway.id1 A vector contains modified KEEG pathway IDs
#' @param PathName_Final A data frame contains is IDs and names of KEEG pathways
#' @param fdr.value A numeric threshold to determine significant parameters
#' @importFrom stats fisher.test p.adjust
#' @return A data frame contains fisher test results for any final pathways
#'
#' @export
#'
#' @examples
#' #All the 187 preprocessed signaling pathways can be entered in analysis by fetch_data_file().
#' #But here you enter a subset of those pathways to see how this package works.
#' files <- system.file("extdata", "test_files_to_start.RData", package = "BNrich", mustWork = TRUE)
#' load(files)
#' Data <- system.file("extdata", "Test_DATA.RData", package = "BNrich", mustWork = TRUE)
#' load(Data)
#' uni_Result <- unify_path(dataH, dataD, MapkG = sub_mapkG, Pathway.id = path.id)
#' M1 <- uni_Result$mapkG1
#' BN <- BN_struct(M1)
#' data_h1 <- uni_Result$data_h
#' data_d1 <- uni_Result$data_d
#' LASSO_Result <- LASSO_BN(BN = BN , data_h = data_h1 , data_d = data_d1)
#' BN_h1 <- LASSO_Result$BN_h
#' BN_d1 <- LASSO_Result$BN_d
#' esti_result <- esti_par(BN_H = BN_h1, BN_D = BN_d1, data_h = data_h1, data_d = data_d1)
#' BNs_H <- esti_result$BNs_h
#' BNs_D <- esti_result$BNs_d
#' coef_h <- esti_result$coef_h
#' coef_d <- esti_result$coef_d
#' var_result <- var_mat(data_h1, coef_h,  BNs_H, data_d1, coef_d, BNs_D)
#' Var_H = var_result$var_mat_Bh
#' Var_D = var_result$var_mat_Bd
#' path.id1 <- uni_Result$pathway.id1
#' Ttest_result <- parm_Ttest(data_h1, coef_h, BNs_H, data_d1, coef_d, BNs_D, Var_H, Var_D, path.id1)
#' BNrich_result <- BNrich(Ttest_result, path.id1, Path.Name)
BNrich <- function(Ttest_Results, Pathway.id1, PathName_Final, fdr.value = 0.05){
  Ttest_results <- Ttest_Results[order(Ttest_Results$pathway.number),]
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
    BNrich_results[d,4]=Pathway.id1[d]
  }
  colnames(BNrich_results) <- c("nk","mk","p.value","ID")
  BNrich_results$fdr <- p.adjust(BNrich_results$p.value, method = "fdr")
  BNrich_results <- BNrich_results[,c("ID","p.value","fdr")]
  BNrich_results <- merge(BNrich_results, PathName_Final, by = "ID")
  colnames(BNrich_results) <- c("pathwayID","p.value","fdr","pathway.number","Name")
  BNrich_results <- BNrich_results[order(BNrich_results$fdr),]
  return(BNrich_results)
}
