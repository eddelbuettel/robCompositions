#' EM-based replacement of rounded zeros in compositional data
#' 
#' Parametric replacement of rounded zeros for compositional data using
#' classical and robust methods based on ilr coordinates with a special
#' choice of balances.
#' 
#' Statistical analysis of compositional data including zeros runs into
#' problems, because log-ratios cannot be applied.  Usually, rounded zeros are
#' considerer as missing not at random missing values.
#' 
#' The algorithm iteratively imputes parts with rounded zeros whereas in each
#' step (1) compositional data are expressed in pivot coordinates (2) tobit regression is
#' applied (3) the rounded zeros are replaced by the expected values (4) the
#' corresponding inverse ilr mapping is applied. After all parts are
#' imputed, the algorithm starts again until the imputations do not change.
#' 
#' @aliases imputeBDL imputeBDLs impRZilr print.replaced checkData adjustImputed
#' @param x data.frame or matrix
#' @param maxit maximum number of iterations
#' @param eps convergency criteria 
#' @param method either "lm", "lmrob" or "pls"
#' @param dl Detection limit for each variable. zero for variables with
#' variables that have no detection limit problems.
#' @param variation, if TRUE those predictors are chosen in each step, who's variation is lowest to the predictor.
#' @param nPred, if determined and variation equals TRUE, it fixes the number of predictors 
#' @param nComp if determined, it fixes the number of pls components. If
#' \dQuote{boot}, the number of pls components are estimated using a
#' bootstraped cross validation approach.
#' @param bruteforce sets imputed values above the detection limit to the
#' detection limit. Replacement above the detection limit are only exeptionally
#' occur due to numerical instabilities. The default is FALSE!
#' @param noisemethod adding noise to imputed values. Experimental
#' @param noise TRUE to activate noise (experimental)
#' @param R number of bootstrap samples for the determination of pls
#' components. Only important for method \dQuote{pls}.
#' @param correction normal or density
#' @param verbose additional print output during calculations.
#' @importFrom cvTools cvFit 
#' @importFrom zCompositions multRepl
#' @import pls
#' @return \item{x }{imputed data} \item{criteria }{change between last and
#' second last iteration} \item{iter }{number of iterations} \item{maxit
#' }{maximum number of iterations} \item{wind}{index of zeros}
#' \item{nComp}{number of components for method pls} \item{method}{chosen
#' method}
#' @author Matthias Templ, method subPLS from Jiajia Chen
#' @references
#' Martin-Fernandez, J.A., Hron, K., Templ, M., Filzmoser, P., Palarea-Albaladejo, J. (2012)
#' Model-based replacement of rounded zeros in compositional data: Classical and robust approaches.
#' \emph{Computational Statistics and Data Analysis}, 56 (9), 2688-2704.
#'
#' Templ, M., Hron, K., Filzmoser, P., Gardlo, A. (2016).
#' Imputation of rounded zeros for high-dimensional compositional data.
#' \emph{Chemometrics and Intelligent Laboratory Systems}, 155, 183-190.
#' 
#' Chen, J., Zhang, X., Hron, K., Templ, M., Li, S. (2018). 
#' Regression imputation with Q-mode clustering for rounded zero replacement in high-dimensional compositional data. 
#' \emph{Journal of Applied Statistics}, 45 (11), 2067-2080.
#' @seealso \code{\link{imputeUDLs}}
#' @keywords manip multivariate
#' @export
#' @importFrom MASS rlm
#' @examples
#' 
#' p <- 10
#' n <- 50
#' k <- 2
#' T <- matrix(rnorm(n*k), ncol=k)
#' B <- matrix(runif(p*k,-1,1),ncol=k)
#' X <- T %*% t(B)
#' E <-  matrix(rnorm(n*p, 0,0.1), ncol=p)
#' XE <- X + E
#' data <- data.frame(pivotCoordInv(XE))
#' col <- ncol(data)
#' row <- nrow(data)
#' DL <- matrix(rep(0),ncol=col,nrow=1)
#' for(j in seq(1,col,2))
#' {DL[j] <- quantile(data[,j],probs=0.06,na.rm=FALSE)}
#' 
#' for(j in 1:col){        
#'   data[data[,j]<DL[j],j] <- 0
#' }
#' \dontrun{
#' # under dontrun because of long exectution time
#' imp <- imputeBDL(data,dl=DL,maxit=10,eps=0.1,R=10,method="subPLS")
#' imp
#' imp <- imputeBDL(data,dl=DL,maxit=10,eps=0.1,R=10,method="pls", variation = FALSE)
#' imp
#' imp <- imputeBDL(data,dl=DL,maxit=10,eps=0.1,R=10,method="lm")
#' imp
#' imp <- imputeBDL(data,dl=DL,maxit=10,eps=0.1,R=10,method="lmrob")
#' imp
#' 
#' data(mcad)
#' ## generate rounded zeros artificially:
#' x <- mcad
#' x <- x[1:25, 2:ncol(x)]
#' dl <- apply(x, 2, quantile, 0.1)
#' for(i in seq(1, ncol(x), 2)){
#'   x[x[,i] < dl[i], i] <- 0
#' } 
#' ni <- sum(x==0, na.rm=TRUE) 
#' ni/(ncol(x)*nrow(x)) * 100
#' dl[seq(2, ncol(x), 2)] <- 0
#' replaced_lm <- imputeBDL(x, dl=dl, eps=1, method="lm",  
#'   verbose=FALSE, R=50, variation=TRUE)$x
#' replaced_lmrob <- imputeBDL(x, dl=dl, eps=1, method="lmrob",  
#'   verbose=FALSE, R=50, variation=TRUE)$x
#' replaced_plsfull <- imputeBDL(x, dl=dl, eps=1, 
#'   method="pls", verbose=FALSE, R=50, 
#'   variation=FALSE)$x 
#' }
#' 
#' 
#' 
`imputeBDL` <-
  function(x, maxit=10, eps=0.1, method="subPLS", 
           dl=rep(0.05, ncol(x)), variation=TRUE,	nPred=NULL, 
           nComp = "boot", bruteforce=FALSE,  
           noisemethod="residuals", noise=FALSE, R=10, 
           correction="normal", verbose=FALSE){
    
    ## check if all is fine:
    # check if values are in (0, dl[i]):
    checkDL <- function(x, dl, indexNA){
      check <- logical(ncol(x))
      for(i in 1:ncol(x)){
        critvals <- x[indexNA[,i],i] > dl[i]
        check[i] <- any(critvals)
        if(check[i]){ 
          x[which(x[indexNA[,i],i] > dl[i]),i] <- dl[i]
        }
      }
      if(any(check) & verbose){
        message(paste(sum(critvals), "/", sum(w), "replaced values has been corrected"))      
      }
      x
    }
    
    ## check if data are fine
    checkData(x, dl)
    
    ## some specific checks
    stopifnot((method %in% c("lm", "MM", "lmrob", "pls", "subPLS")))
    if(method=="pls" & ncol(x)<5) stop("too less variables/parts for method pls")
    if(!(correction %in% c("normal","density"))){
      stop("correction method must be normal or density")
    }
    if(method == "pls" & variation){
      stop("if variation is TRUE then pls is not supported.")
    }
    
    ## how to deal with user input on nComp
    if(is.null(nComp)){
      pre <- FALSE
      nC <- NULL
    } else if(nComp=="boot"){
      nC <- integer(ncol(x))
      pre <- TRUE
    } else if(length(nComp) == ncol(x)){
      nC <- nComp
      pre <- FALSE
    } else  {
      pre <- FALSE	
    }
    
    
    ## store some important values
    n <- nrow(x) 
    d <- ncol(x)
    rs <- rowSums(x, na.rm = TRUE)
    
    if(method == "subPLS"){
      w <- x == 0
      indexFinalCheck <- x == 0
      mulzero <- zCompositions::multRepl(x,label=0,dl=dl)
      dd <- as.dist(robCompositions::variation(mulzero))
      myClusterFunction <- function(x) {
        if (!requireNamespace("fpc", quietly = TRUE)) {
          stop("Package 'fpc' is required for this function. Please install it.")
        }
        if (requireNamespace("fpc", quietly = TRUE)) {
          res <- fpc::pamk(x)
        }
        return(res)
      }
      g <- myClusterFunction(dd)$nc
      pos <- as.matrix(cutree(hclust(dd, method="ward.D"),g))
      indNA <- apply(x,2,function(x){any(x==0)})
      gz <- intersect(order(-table(pos[indNA])),which(table(pos)>1))
      
      fun.PLS <- function(x,data,dl,pos,b,R)
      {
        data_b <- data[,pos==b]
        data_nb <- data[,pos!=b]
        col1 <- ncol(data_b)
        row <- nrow(data_b)
        u <- cbind(cenLR(data_b)$x,cenLR(data_nb)$x)
        u1 <- u[,1:col1]
        n <- bootnComp(u[,-(1:col1),drop=FALSE],as.matrix(u1),R,plotting=FALSE)$res
        # form <- paste(paste(colnames(u1), collapse = "+"), "~", paste(colnames(u), collapse = "+"))
        # form <- as.formula(form)
        # pls <- mvr(form, data=cbind(u1, u), method="simpls", 
        #            ncomp=n) 
        pls <- mvr(as.matrix(u1)~ as.matrix(u[,-(1:col1)]),ncomp=n,method="simpls")
        # pls <- mvr(u1~ u[,-(1:col1)],ncomp=n,method="simpls")
        mean <- matrix(predict(pls,ncomp=n),ncol=col1)
        sigma <- cov(u1-mean)
        sig_b <- x[,pos==b]==0
        cz <- which(colSums(sig_b)!=0)
        rz <- which(rowSums(sig_b)!=0)
        lj <- dl[pos==b]
        for(j in cz)
        {
          linjie <- as.matrix(cenLR(cbind(rep(lj[j],row),data_b[,-j,drop=FALSE]))$x[,1,drop=FALSE])
          u1[sig_b[,j],j] <- mean[sig_b[,j],j]-sqrt(sigma[j,j])*(dnorm((linjie[sig_b[,j]]-mean[sig_b[,j],j])/sqrt(sigma[j,j]))/pnorm((linjie[sig_b[,j]]-mean[sig_b[,j],j])/sqrt(sigma[j,j])))
        }
        for(i in rz)
        {if(rowSums(sig_b)[i]<col1)
        {u1[i,!sig_b[i,]] <- (sum(!sig_b[i,])*u1[i,!sig_b[i,]]-rep(sum(u1[i,]),sum(!sig_b[i,])))/(sum(!sig_b[i,]))}
        }
        chabu <- adjustImputed(cenLRinv(u1)[rowSums(sig_b)<col1,],x[rowSums(sig_b)<col1,pos==b],sig_b[rowSums(sig_b)<col1,])
        data[rowSums(sig_b)<col1,pos==b] <- chabu
        return(data)
      }
      
      it <- 1
      criteria <- 99999999
      impute <- mulzero
      while(it <= maxit & criteria >= eps){
        xold <- impute
        for (m in 1:length(gz))
        {impute <- fun.PLS(x,impute,dl,pos,gz[m],R)}
        it <- it+1
        criteria <- sum(((xold-impute)/impute)^2, na.rm = TRUE)
      }
      impute <- checkDL(impute, dl, indexFinalCheck)
      res <- list(x=impute, criteria=criteria, iter=it,
                  maxit=maxit, wind=w, nComp=nC, nPred=nPred,
                  variation=variation,
                  method=method, dl=dl)
      class(res) <- "replaced"
      return(res)
    }
    
    # if(method == "pls_clust"){
    #   indexFinalCheck <- x == 0
    #   mulzero <- zCompositions::multRepl(x,label=0,dl=dl)
    #   dd <- as.dist(robCompositions::variation(mulzero))
    # myClusterFunction <- function(x) {
    #   if (!requireNamespace("fpc", quietly = TRUE)) {
    #     stop("Package 'fpc' is required for this function. Please install it.")
    #   }
    #   res <- fpc::pamk(x)
    #   return(res)
    # }
    # g <- myClusterFunction(dd)$nc
    #   pos <- as.matrix(cutree(hclust(dd, method="ward.D"),g))
    #   indNA <- apply(x,2,function(x){any(x==0)})
    #   gz <- intersect(order(-table(pos[indNA])),which(table(pos)>1))
    #   
    #   fun.PLS <- function(x,data,dl,pos,b,R)
    #   {
    #     data_b <- data[,pos==b]
    #     data_nb <- data[,pos!=b]
    #     col1 <- ncol(data_b)
    #     row <- nrow(data_b)
    #     u <- cbind(cenLR(data_b)$x.clr, cenLR(data_nb)$x.clr)
    #     u1 <- u[,1:col1]
    #     n <- bootnComp(u[,-(1:col1),drop=FALSE],u1,R,plotting=FALSE)$res
    #     pls <- mvr(u1~ u[,-(1:col1)],ncomp=n,method="simpls")
    #     mean <- matrix(predict(pls,ncomp=n),ncol=col1)
    #     sigma <- cov(u1-mean)
    #     sig_b <- x[,pos==b]==0
    #     cz <- which(colSums(sig_b)!=0) 
    #     rz <- which(rowSums(sig_b)!=0)
    #     lj <- dl[pos==b]
    #     for(j in cz)
    #     {
    #       linjie <- clr(cbind(rep(lj[j],row),data_b[,-j,drop=FALSE]))[,1,drop=FALSE]
    #       u1[sig_b[,j],j] <- mean[sig_b[,j],j]-sqrt(sigma[j,j])*(dnorm((linjie[sig_b[,j]]-mean[sig_b[,j],j])/sqrt(sigma[j,j]))/pnorm((linjie[sig_b[,j]]-mean[sig_b[,j],j])/sqrt(sigma[j,j])))
    #     }
    #     for(i in rz)
    #     {if(rowSums(sig_b)[i]<col1)
    #     {u1[i,!sig_b[i,]] <- (sum(!sig_b[i,])*u1[i,!sig_b[i,]]-rep(sum(u1[i,]),sum(!sig_b[i,])))/(sum(!sig_b[i,]))}
    #     }
    #     chabu <- adjustImputed(clrInv(u1)[rowSums(sig_b)<col1,],x[rowSums(sig_b)<col1,pos==b],sig_b[rowSums(sig_b)<col1,])
    #     data[rowSums(sig_b)<col1,pos==b] <- chabu
    #     return(data)
    #   }
    #   
    #   it <- 1
    #   criteria <- 99999999
    #   impute <- mulzero
    #   while(it <= maxit & criteria >= eps){
    #     xold <- impute
    #     for (m in 1:length(gz))
    #     {impute <- fun.PLS(x,impute,dl,pos,gz[m],R)}
    #     it <- it+1
    #     criteria <- sum(((xold-impute)/impute)^2, na.rm = TRUE)
    #   }
    #   impute <- checkDL(impute, dl, indexFinalCheck)
    #   
    #   res <- list(x=impute, criteria=criteria, iter=it, 
    #               maxit=maxit, wind=w, nComp=nC, nPred=nPred,
    #               variation=variation,
    #               method=method, dl=dl)
    #   class(res) <- "replaced"
    #   return(res)
    # }
    
    ## values below the detection limit are rounded zeros to be imputed.
    ## checkData() validated this above; previously its thresholding result was
    ## discarded (return value dropped), so below-DL values that were not exactly
    ## 0 were silently skipped. Apply the thresholding here so they are imputed.
    for(i in seq_len(ncol(x))) x[x[, i] < dl[i], i] <- 0
    ## set zeros to NA for easier handling
    x[x == 0] <- NA
    x[x < 0] <- NA
    indexFinalCheck <- is.na(x)
    
    ## sort variables of x based on 
    ## decreasing number of zeros in the variables
    cn <- colnames(x)
    wcol <- - abs(apply(x, 2, function(x) sum(is.na(x))))
    o <- order(wcol)
    x <- x[,o]
    if(verbose) cat("number of variables with zeros:\n", sum(wcol != 0))
    ## --> now work in revised order of variables
    ## dl must also be in correct order
    dlordered <- dl[o]
    
    #################
    ## index of missings / non-missings
    w <- is.na(x)
    wn <- !is.na(x)
    xcheck <- x
    w2 <- is.na(x)
    
    ## initialisation
    indNA <- apply(x, 2, function(x){any(is.na(x))})
    #    print(indNA)
    for(i in 1:length(dl)){
      ind <- is.na(x[,i])
      #		if(length(ind) > 0) x[ind,i] <- dl[i]*runif(sum(ind),1/3,2/3)
      if(length(ind) > 0) x[ind,i] <- dlordered[i] *2/3
    }
    xOrig <- x
    
    ## nPred if variation == TRUE and cv for number of predictors
    ## evaluate the vector containing the amount of predictors
    if(!is.null(nPred) & length(nPred) == 1){
      nPred <- rep(nPred, ncol(x))
    }
    if(!is.null(nPred) & length(nPred) > 1){
      stop("nPred must be NULL or a vector of length 1.")
    }
    if(is.null(nPred) & variation){
      ptmcv <- proc.time()
      if(verbose) cat("\n cross validation to estimate number of predictors\n")
      ii <- 1
      if(verbose) pb <- txtProgressBar(min = 0, max = sum(indNA), style = 3)
      nPred <- numeric(nrow(x)) 
      rtmspe <- NULL
      for(i in which(indNA)){
        xneworder <- cbind(x[, i, drop=FALSE], x[, -i, drop=FALSE]) 
        rv <- variation(x, method = "Pairwise")[1,]
        cve <- numeric()
        for(np in seq(3, min(c(27,ncol(x),floor(nrow(x)/2))), 3)){
          s <- sort(rv)[np]
          cols <- which(rv <= s)[1:np]
          xn <- xneworder[, cols]
          xilr <- data.frame(pivotCoord(xn))
          colnames(xilr)[1] <- "Y"
          call <- call(method, formula = Y ~ .)
          # perform cross-validation
          cve[np] <- suppressWarnings(cvFit(call, data = xilr, 
                                            y = xilr$Y, 
                                            cost = cvTools::rtmspe,
                                            K = 5, R = 1, costArgs = list(trim = 0.1), seed = 1234)$cv)
        }
        nPred[i] <- which.min(cve)
        ## update progress bar
        if(verbose) setTxtProgressBar(pb, ii); ii <- ii + 1
      }  
      if(verbose) close(pb)
      ptmcv <- proc.time() - ptmcv
    }
    
    ###############################
    ###   start the iteration   ###
    ###############################
    
    if(verbose) cat("\n start the iteration:")
    it <- 1; criteria <- 99999999
    while(it <= maxit & criteria >= eps){
      if(verbose) cat("\n iteration", it, "; criteria =", criteria)	
      xold <- x  
      for(i in which(indNA)){
        if(verbose) cat("\n replacement on (sorted) part", i)
        ## sort data columns first
        ## i'th positions must be first
        xneworder <- cbind(x[, i, drop=FALSE], x[, -i, drop=FALSE]) 
        ## if based on variation matrix:
        if(variation){
          orig <- xneworder  
          rv <- variation(x, method = "Pairwise")[1,]
          s <- sort(rv)[nPred[i]]
          cols <- which(rv <= s)[1:nPred[i]]
          xneworder <- xneworder[, cols]
        }
        ## factors for preserving abs values:
        # fac <- multis(xneworder)
        ## detection limit in ilr-space
        forphi <- cbind(rep(dlordered[i], n), xneworder[,-1,drop=FALSE])
        if(any(is.na(forphi))) break()
        phi <- pivotCoord(forphi)[,1]
        #		part <- cbind(x[,i,drop=FALSE], x[,-i,drop=FALSE])
        xneworder[xneworder < 2*.Machine$double.eps] <- 2*.Machine$double.eps
        xilr <- data.frame(pivotCoord(xneworder))
        #        c1 <- colnames(xilr)[1]					
        #        colnames(xilr)[1] <- "V1"	
        response <- as.matrix(xilr[, 1, drop=FALSE])
        predictors <- as.matrix(xilr[, -1, drop=FALSE])
        if(method=="lm"){ 
          reg1 <- lm(response ~ predictors)
          yhat <- predict(reg1, newdata=data.frame(predictors))
        } else if(method=="MM" | method=="lmrob"){
          reg1 <- MASS::rlm(response ~ predictors, method="MM",maxit = 100)#rlm(V1 ~ ., data=xilr2, method="MM",maxit = 100)
          yhat <- predict(reg1, newdata=data.frame(predictors))
        } else if(method=="pls"){
          if(it == 1 & pre){ ## evaluate ncomp.
            nC[i] <- bootnComp(xilr[, 2:ncol(xilr), drop=FALSE], y=xilr[, 1], R, 
                               plotting=FALSE)$res #$res2
          }
          if(verbose) cat("   ;   ncomp:", nC[i])
          reg1 <- mvr(as.matrix(response) ~ as.matrix(predictors), ncomp=nC[i], method="simpls")
          yhat <- predict(reg1, newdata=data.frame(predictors), ncomp=nC[i])
        }
        
        #		s <- sqrt(sum(reg1$res^2)/abs(nrow(xilr)-ncol(xilr))) ## quick and dirty: abs()
        s <- sqrt(sum(reg1$res^2)/nrow(xilr)) 
        yhat <- as.numeric(yhat)
        ex <- (phi - yhat)/s 
        if(correction=="normal"){
          yhat2sel <- ifelse(dnorm(ex[w[, i]]) > .Machine$double.eps,
                             yhat[w[, i]] - s*dnorm(ex[w[, i]])/pnorm(ex[w[, i]]),
                             yhat[w[, i]])
        } else if(correction=="density"){
          stop("correction equals density is no longer be supported")
          # den <- density(ex[w[,i]])
          # distr <- sROC::kCDF(ex[w,i])
        }
        if(any(is.na(yhat)) || any(is.infinite(yhat))) stop("Problems in ilr because of infinite or NA estimates")
        # check if we are under the DL:
        if(any(yhat2sel >= phi[w[, i]])){
          yhat2sel <- ifelse(yhat2sel > phi[w[, i]], phi[w[, i]], yhat2sel)
        }
        xilr[w[, i], 1] <- yhat2sel
        xinv <- pivotCoordInv(xilr)
        ## if variation:
        if(variation == TRUE){
          xneworder <- adjustImputed(xinv, xneworder, w2[, cols])
          orig[, cols] <- xneworder #* fac
          xinv <- orig
        }
        ## reordering of xOrig
        if(i %in% 2:(d-1)){
          xinv <- cbind(xinv[,2:i], xinv[,c(1,(i+1):d)])
        }
        if(i == d){
          xinv <- cbind(xinv[,2:d], xinv[,1])
        }
        if(!variation){ 
          x <- adjustImputed(xinv, xOrig, w2)
        } else {
          x <- xinv
        }
        #		x <- adjust3(xinv, xOrig, w2) 
        #		## quick and dirty:
        #		x[!w] <- xOrig[!w]
      }
      
      it <- it + 1
      criteria <- sum( ((xold - x)/x)^2, na.rm=TRUE) ## DIRTY: (na.rm=TRUE)
      if(verbose & criteria != 0) cat("\n iteration", it, "; criteria =", criteria)
    }
    
    #### add random error ###
    if(noise){
      for(i in which(indNA)){
        if(verbose) cat("\n add noise on variable", i)
        
        # add error terms
        inderr <- w[,i]
        if(noisemethod == "residuals") {
          error <- sample(residuals( reg1 )[inderr], 
                          size=wcol[i], replace=TRUE)
          reg1$res[inderr] <- error
        } else {
          mu <- median(residuals( reg1 )[inderr])
          sigma <- mad(residuals( reg1 )[inderr])
          error <- rnorm(wcol[i], mean=mu, sd=sigma)
          reg1$res[inderr] <- error		   
        }
        # return realizations
        yhat[inderr] <- yhat[inderr] + error
        
        
        s <- sqrt(sum(reg1$res^2)/nrow(xilr)) ## quick and dirty: abs()
        ex <- (phi - yhat)/s 
        yhat2sel <- ifelse(dnorm(ex[w[, i]]) > .Machine$double.eps,
                           yhat[w[, i]] - s*dnorm(ex[w[, i]])/pnorm(ex[w[, i]]),
                           yhat[w[, i]])
        if(any(is.na(yhat)) || any(is.infinite(yhat))) stop("Problems in ilr because of infinite or NA estimates")
        # check if we are under the DL:
        if(any(yhat2sel >= phi[w[, i]])){
          yhat2sel <- ifelse(yhat2sel > phi[w[, i]], phi[w[, i]], yhat2sel)
        }
        xilr[w[, i], 1] <- yhat2sel
        xinv <- pivotCoordInv(xilr) 
        ## reordering of xOrig
        if(i %in% 2:(d-1)){
          xinv <- cbind(xinv[,2:i], xinv[,c(1,(i+1):d)])
        }
        if(i == d){
          xinv <- cbind(xinv[,2:d], xinv[,1])
        }
        
        #	   x <- adjust2(xinv, xOrig, w) 
        #	   ## quick and dirty:
        #	   x[!w] <- xOrig[!w]
      }
    }
    ### end add random error ###
    #   x <- adjust3(x, xOrig, w)
    #   x[!w] <- xOrig[!w] 
    x <- x[,order(o)] ## checked: reordering is OK!
    w <- w[,order(o)] ## keep wind (zero index) aligned with returned original-order x
    colnames(x) <- cn
    
    if(verbose){
      message(paste(sum(w)), "values below detection limit has been imputed \n below the corresponding detection limits")
    }
    
    
    x <- checkDL(x, dl, indexFinalCheck)
    
    res <- list(x=x, criteria=criteria, iter=it, 
                maxit=maxit, wind=w, nComp=nC, nPred=nPred,
                variation=variation,
                method=method, dl=dl)
    class(res) <- "replaced"
    return(res)
  }

#' @rdname imputeBDL
#' @export
#' @param xImp imputed data set
#' @param xOrig original data set
#' @param wind index matrix of rounded zeros
adjustImputed <- function(xImp, xOrig, wind){
  ## aim: 
  ## (1) ratios must be preserved
  ## (2) do not change original values
  ## (3) adapt imputations
  xneu  <- xImp
  s1 <- rowSums(xOrig, na.rm = TRUE)
  ## per row: consider rowsums of imputed data
  sumPrevious <- sumAfter <- numeric(nrow(xImp))
  for (i in 1:nrow(xImp)) {
    if(any(wind[i, ]) & !all(wind[i,])){
      sumPrevious[i] <- sum(xOrig[i, !wind[i, ]]) 
      sumAfter[i] <- sum(xImp[i, !wind[i,]])
    } else{ 
      sumPrevious[i] <- sumAfter[i] <- 1
    }
  }
  # how much is rowsum increased by imputation:
  fac <- sumPrevious/sumAfter
  
  #    # decrese rowsums of orig.
  #    s1[i] <- s1[i]/fac
  #  }
  ## for non-zeros overwrite them:
  xneu[!wind] <- xOrig[!wind]
  ## adjust zeros:
  for(i in 1:nrow(xImp)){
    if(any(wind[i,])){
      xneu[i,wind[i,]] <- fac[i]*xneu[i,wind[i,]]
    }
  }
  return(xneu)
}

#' @rdname imputeBDL
#' @export
`checkData` <- function(x, dl){
  if(any(is.na(x))) stop("your data includes missing values.\n Use impKNNa() or impCoda() to impute them first.")
  if( is.vector(x) ) stop("x must be a matrix or data frame")
  ## check if only numeric variables are in x:
  cl <- lapply(x, class)
  if(!all(cl %in% "numeric")) stop("some of your variables are not of class numeric.")
  if( length(dl) < ncol(x)) stop(paste("dl has to be a vector of ", ncol(x)))
  if(any(is.na(x))) stop("missing values are not allowed. \n Use impKNNa or impCoda to impute them first.")
  
  #     pre <- TRUE
  #      if(length(nComp) != ncol(x) & nComp!="boot") stop("nComp must be NULL, boot or of length ncol(x)")
  #    } else if(nComp == "boot"){#
  #		pre <- TRUE
  #	} else {
  #		pre <- FALSE
  #	}
  #################
  ## zeros to NA:
  # check if values are in (0, dl[i]):
  check <- logical(ncol(x))
  for(i in 1:ncol(x)){
    #      check[i] <- any(x[,i] < dl[i] & x[,i] != 0)
    x[x[,i] < dl[i],i] <- 0
  }
  #    if(any(check)){warning("values below detection limit have been set to zero and will be imputed")}
  check2 <- any(x < 0)
  if(check2){warning("values below 0 set have been set to zero and will be imputed")}
  x[x == 0] <- NA
  x[x < 0] <- NA
  indexFinalCheck <- is.na(x)
  
  ## check if rows consists of only zeros:
  checkRows <- unlist(apply(x, 1, function(x) all(is.na(x))))
  if(any(checkRows)){ 
    w <- which(checkRows)
    cat("\n--------\n")
    message("Rows with only zeros are not allowed")
    message("Remove this rows before running the algorithm")
    cat("\n--------\n")      
    stop(paste("Following rows with only zeros:", w))
  }  
  
  ## check if cols consists of only zeros:
  checkCols <- unlist(apply(x, 2, function(x) all(is.na(x))))
  if(any(checkCols)){ 
    w <- which(checkCols)
    cat("\n--------\n")
    message("Cols with only zeros are not allowed")
    message("Remove this columns before running the algorithm")
    cat("\n--------\n")      
    stop(paste("\n Following cols with only zeros:", colnames(x)[w]))
  }  
  
  indNA <- apply(x, 2, function(x){any(is.na(x))})
  
  ## check if for any variable with zeros,
  ## the detection limit should be larger than 0:
  if(any(dl[indNA]==0)){
    w <- which(dl[indNA]==0)
    invalidCol <- colnames(x)[w]
    for(i in 1:length(invalidCol)){
      cat("-------\n")
      cat(paste("Error: variable/part", invalidCol[i], 
                "has detection limit 0 but includes zeros"))
      cat("\n-------\n")
    }
    stop(paste("Set detection limits larger than 0 for variables/parts \n including zeros"))
  }
  
}

#' @rdname imputeBDL
#' @method print replaced  
#' @export
#' @param ... further arguments passed through the print function
print.replaced <- function(x, ...){
  message(paste("\n", sum(x$w), "values below detection limit were imputed \n below their corresponding detection limits.\n"))
}

#' Bootstrap to find optimal number of components
#' 
#' Combined bootstrap and cross validation procedure to find optimal number of
#' PLS components
#' 
#' Heavily used internally in function impRZilr.
#' 
#' @param X predictors as a matrix
#' @param y response
#' @param R number of bootstrap replicates
#' @param plotting if TRUE, a diagnostic plot is drawn for each bootstrap
#' replicate
#' @return Including other information in a list, the optimal number of
#' components
#' @author Matthias Templ
#' @seealso \code{\link{impRZilr}}
#' @keywords manip
#' @export
#' @examples
#' 
#' ## we refer to impRZilr()
#' 
bootnComp <- function(X, y, R=99, plotting=FALSE){
  ind <- 1:nrow(X)
  d <- matrix(, ncol=R, nrow=nrow(X))#nrow(X))
  nc <- integer(R)
  for(i in 1:R){
    bootind <- sample(ind)
    #    XX <- X
    #    yy <- y
    if(is.null(ncol(y))){
      ds <- cbind(X[bootind,], as.numeric(y[bootind]))
      colnames(ds)[ncol(ds)] <- "V1"
      res1 <- mvr(V1~., data=data.frame(ds), method="simpls", 
                  validation="CV")
    } else {
      ds <- cbind(X[bootind,], as.matrix(y[bootind,])) 
      form <- paste(paste(colnames(ds[, (ncol(X)+1):ncol(ds)]), collapse = "+"), "~", paste(colnames(ds[, 1:ncol(X)]), collapse = "+"))
      form <- as.formula(form)
      res1 <- mvr(form, data=data.frame(ds), method="simpls", 
                  validation="CV") 
    }
    
    
    d[1:res1$ncomp, i] <- res1$validation$PRESS
    nc[i] <- which.min(res1$validation$PRESS)
    #    d[1:reg1$ncomp,i] <- as.numeric(apply(reg1$validation$pred, 3, 
    #                                  function(x) sum(((y - x)^2)) ) )
  }
  d <- na.omit(d)
  sdev <- apply(d, 1, sd, na.rm=TRUE)
  means <- apply(d, 1, mean, na.rm=TRUE)
  mi <- which.min(means)
  r <- round(ncol(X)/20)
  mi2 <- which.min(means[r:length(means)])+r-1
  w <- means < min(means) + sdev[mi]
  threshold <- min(means) + sdev[mi]
  sdev <- sdev
  means <- means
  mi <- mi
  means2 <- means
  means2[!w] <- 999999999999999
  res <- which.min(means2)
  mi3 <- which.max(w)
  #  minsd <- means - sdev > means[mi]
  #  check <- means
  #  check[!minsd] <- 99999999
  if(plotting) plot(means, type="l")
  #  res <- which.min(check)
  list(res3=mi3, res2=mi2, res=res, means=means)
}
