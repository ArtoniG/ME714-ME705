## *** Diagnostiic analysis for linear mixed models fitted via library nlme 

## *** Authors: Francisco Marcelo M. Rocha, Juvencio S. Nobre and Julio M. Singer

#####################################################################################################################
## This function generates the following diagnostic plots for gaussian linear mixed models fitted via function lme ##
## in library nlme.                                                                                                ##
#####################################################################################################################

# plot 1: Standardized marginal residuals versus fitted values and corresponding histogram
# plot 2: Standardized Mahalanobis distance versus unit indices
# plot 3: Chi-squared QQ plot for Mahalanobis distance
# plot 4: tandardized Lesaffre-Verbeke measure versus unit indices
# plot 5: Standardized conditional residuals versus fitted values and corresponding histogram
# plot 6: Normal QQ plot and histogram for standardized least confounded conditional residuals
# plot 7: Cook's Conditional distance versus observation indices
# plot 8: Cook's Conditional distance 1 (D1i) versus observation indices
# plot 9: Cook's Conditional distance 2 (D2i) versus observation indices
# plot 10: Cook's Conditional distance 3 (D3i) versus observation indices
# plot 11: Generalized marginal leverage (L1) versus unit indices
# plot 12: Generalized leverage for random components (L2) versus unit indices
# plot 13: Generalized joint leverage (L) versus unit indices
# plot 14: Generalized joint leverage [L1i(jj)] versus observation indices
# plot 15: Generalized joint leverage [L2i(jj)] versus observation indices
# plot 16: Generalized joint leverage [L3i(jj)] versus observation indices

###################################################################################################################
## 1. The data must be arranged in the LONG format, with units indexed by a numerical (not necessarily equally   ##
##    spaced) variable. If the function "groupedData()" is used, the unit index must be transformed to numerical.## 
## 2. The horizontal axis contains either the unit indices or the observation indices                            ##
## 3. The outliers are labelled in the format unit.obs (e.g. 13.5 denotes the fifth observation of unit          ##
##    labelled 13)                                                                                               ##
###################################################################################################################

## An example is

## dataset<-groupedData(response ~  time|id, data=dataset_label)
        # Construct a groupedData object
        # Note that the grouping variable must be a factor
## model1<-lme(response ~ time, random = ~time|id, na.action=na.omit, data=dataset)  
        # use of the function lme in nlme to fit a mixed model 
## resid1<-residdiag.nlme(model1,limit=2,plotid=1:16)
        # the computed quantities will be saved in the object labelled resid1
        # the limit option indicates cutpoints for the residual plots
        # plotid indicates the required residual plots
## names(resid1)  
        # produces the labels of the quantities saved in resid1
## resid1$least.confounded.residuals   
        # produces the least confounded standardized conditional residuals


###############################################################################################################
## *** References:                                                                                           ##
##                                                                                                           ##
##  - Nobre, J.S. and Singer, J.M. (2007). Residual Analysis for Linear Mixed Models.                        ##
##        Biometrical Journal, 49, 1-13.                                                                     ##
##  - Nobre, J.S. and Singer, J.M. (2011). Leverage analysis for linear mixed models.                        ##
##        Journal of Applied Statistics, 38, 1063-1072.                                                      ##
##  - Singer, J.M., Nobre, J.S. and Rocha, F.M.M. (2014). Diagnostic and treatment in linear mixed models.   ##
##        Submitted.                                                                                         ##
##  - Pinheiro, J.C. and Bates, D.M. (2000). Mixed-effects models in S and S-plus. 1st edition.              ##
##        New York: Springer                                                                                 ##
##  - Scheipl, F., Greven, S. and Kuechenhoff, H. (2008) Size and power of tests for a zero random effect    ##
##       variance or polynomial regression in additive and linear mixed models.                              ##
##       Computational Statistics & Data Analysis, 52, 3283-3299.                                            ##
###############################################################################################################

residdiag.nlme = function (fit, limit,plotid=NULL,option=0) {
require(MASS) 
require(Matrix)
require(car)


####################################################################################

# Function for extracting square root of a matrix

sqrt.matrix <- function(mat) {              
  mat <- as.matrix(mat)  # new line of code
  singular_dec<-svd(mat,LINPACK=F)
  U<-singular_dec$u
  V<-singular_dec$v
  D<-diag(singular_dec$d)
  sqrtmatrix<-U%*%sqrt(D)%*%t(V)
  #  return(list(sqrt=sqrtmatrix))
}

#####################################################################################

## This function extracts various objects of the function lme

extract.lmeDesign2 <- function(m){
    start.level = 1
    data <- if(any(!complete.cases(m$data))){
      warning("Removing incomplete cases from supplied data.") 
      m$data[complete.cases(m$data),]
    } else m$data
    grps <- nlme::getGroups(m)
    n <- length(grps)
    X <- list()
    grp.dims <- m$dims$ncol
    Zt <- model.matrix(m$modelStruct$reStruct, data)
    cov <- as.matrix(m$modelStruct$reStruct)
    i.col <- 1
    n.levels <- length(m$groups)
    Z <- matrix(0, n, 0)
    if (start.level <= n.levels) {
      for (i in 1:(n.levels - start.level + 1)) {
        if(length(levels(m$groups[[n.levels-i+1]]))!=1)
        {
          X[[1]] <- model.matrix(~m$groups[[n.levels - i +
                                              1]] - 1, 
                                 contrasts.arg = c("contr.treatment",
                                                   "contr.treatment"))
        }
        else X[[1]]<-matrix(1, n, 1)
        X[[2]] <- as.matrix(Zt[, i.col:(i.col + grp.dims[i] -
                                          1)])
        i.col <- i.col + grp.dims[i]
        Z <- cbind(mgcv::tensor.prod.model.matrix(X),Z)
      }
      Vr <- matrix(0, ncol(Z), ncol(Z))
      start <- 1
      for (i in 1:(n.levels - start.level + 1)) {
        k <- n.levels - i + 1
        for (j in 1:m$dims$ngrps[i]) {
          stop <- start + ncol(cov[[k]]) - 1
          Vr[ncol(Z)+1-(stop:start),ncol(Z)+1-(stop:start)] <- cov[[k]]
          start <- stop + 1
        }
      }
    }
    X <- if(class(m$call$fixed) == "name" &&  !is.null(m$data$X)){
      m$data$X
    } else   {
      model.matrix(formula(eval(m$call$fixed)),data)
    }
    y<-as.vector(matrix(m$residuals, ncol=NCOL(m$residuals))[,NCOL(m$residuals)] + 
                   matrix(m$fitted, ncol=NCOL(m$fitted))[,NCOL(m$fitted)])
    return(list(
      Vr=Vr, #Cov(RanEf)/Var(Error)
      X=X,
      Z=Z,
      sigmasq=m$sigma^2,
      lambda=unique(diag(Vr)),
      y=y,
      k=n.levels
    )
    )
}

#####################################################################################
#####################################################################################

## Extracting information from lme fitted model and dataset 
  
  data.fit<-extract.lmeDesign2(fit)
  data<-fit$data
  y<- data.fit$y
  X<- data.fit$X
  N<-length(y)                                    # Number of observations
 id<- as.numeric(names(getResponse(fit)))
 subject<-as.numeric(unique(id))
  n<-length(as.numeric(names(table(id))))         # Number of units
  vecni<-(table(id))                              # Vector with number of observations per unit
  p<-ncol(X)                                      # Number of fixed parameters

  obs=numeric()

for(i in 1:n)
 {
  obs=append(obs,1:vecni[i])                       # Observation labels for each unit
 }

#####################################################################################
  
## Construction of the Z matrix

mataux<-model.matrix(fit$modelStruct$reStruct, data)
mataux<- as.data.frame(cbind(mataux,id))
#Z<- as.matrix((subset(split(mataux,id==min(id),drop=T)$`TRUE`,select=-id)))



for(i in (as.numeric(unique(id)))){ 
  ifelse( i==min(as.numeric(unique(id))), 
          Z<- as.matrix((subset(split(mataux,id==min(id),drop=T)$`TRUE`,select=-id))),
          Z<-as.matrix(bdiag(Z,as.matrix(subset(split(mataux,id==i,drop=T)$`TRUE`,select=-id)))))
}

#####################################################################################

## Estimate of the Gamma (Gam) matrix (random effects covariance matrix)
  
 g<-  getVarCov(fit,type="random.effects")

#####################################################################################

## Estimate of random effects covariance matrix for each unit

 q<-dim(g)[1]                           # Total number of random effects
 Gam<-as.matrix(kronecker(diag(length(as.numeric(unique(id)))),g))

#####################################################################################
  
# Estimate of the covariance matrix of conditional errors (homoskedastic conditional independence model)

R <- getVarCov(fit,type="conditional",individual=1)[[1]]
for (i in 2:length(as.numeric(unique(id)))){
  R <- as.matrix(bdiag(R,getVarCov(fit,type="conditional",individual=i)[[1]] ) )
  
}

#####################################################################################

# Construction of covariance matrix of Y

  V<- (Z%*%Gam%*%t(Z)) + R
  iV<-ginv(V)                        # inverse of V

#####################################################################################

# Construction of the Q matrix 

  varbeta<-ginv((t(X)%*%iV%*%X))
  Q<-(iV-iV%*%X%*%(varbeta)%*%t(X)%*%iV )  

#####################################################################################

# EBLUE and EBLUP
  
  eblue<-as.vector(fixef(fit))
  eblup<-Gam%*%t(Z)%*%iV%*%(y-X%*%eblue)

#####################################################################################

## Residual analysis

  
  predm<-X%*%eblue                   # Predicted values for expected response
  predi<-X%*%eblue+Z%*%eblup         # Predicted values for units
  resm<-(y-predm)                    # Marginal residuals
  resc<-(y-predi)                    # Conditional residuals

#####################################################################################

## Variance of marginal residuals
  
  var.resm<-V-X%*%ginv(t(X)%*%iV%*%X)%*%t(X) 

#####################################################################################

## Standardized marginal residuals 
  
resmp<-matrix(0,N,1)
auxni=as.vector(vecni)
for (t in 1:n){ 
  li<- sum(vecni[1:t-1])+1
  ls<- sum(vecni[1:t])
  if(vecni[t]==1){
    
    auxr2 <- solve(sqrt(var.resm[li:ls,li:ls]))
    resmpi<-(auxr2)%*%resm[li:ls]
    resmp[li:ls,]<- resmpi
  }
  else
  {  
    auxr1 <- sqrt.matrix(var.resm[li:ls,li:ls])
    auxr2 <- solve(sqrt.matrix(var.resm[li:ls,li:ls]))
    resmpi<- auxr2%*%resm[li:ls]
    resmp[li:ls,]<- resmpi
  } 
}

#####################################################################################

## Variance of conditional residuals
 
  var.resc<- R%*%Q%*%R

#####################################################################################
    
## Standardized conditional residuals
  
rescp<-matrix(0,N,1)

for (t in 1:n){ 
  li<- sum(vecni[1:t-1])+1
  ls<- sum(vecni[1:t])
  if(vecni[t]==1){
    
    auxr2 <- solve(sqrt(var.resc[li:ls,li:ls]))
    rescpi<-(auxr2)%*%resc[li:ls]
    rescp[li:ls,]<- rescpi
  }
  else
  {  
    
    auxr2 <- solve(sqrt.matrix(var.resc[li:ls,li:ls]))
    rescpi<- auxr2%*%resc[li:ls]
    rescp[li:ls,]<- rescpi
  } 
}

#####################################################################################
    
## Mahalanobis's distance 

  aux=Gam%*%t(Z)%*%Q%*%Z%*%Gam

  qm<-q-1
  dm<-matrix(0,n,1)

  for(j in 1:n) 
  {
    if(q==1)
    {
      gbi<-aux[j,j]
      eblupi<-eblup[(q*j-qm):(q*j)]
      dmi<- t(eblupi)%*%ginv(gbi)%*%eblupi
      dm[j]<-dmi
    }
    else
    {
      gbi<-aux[(q*j-qm):(q*j),(q*j-qm):(q*j)]
      eblupi<-eblup[(q*j-qm):(q*j)]
      dmi<- t(eblupi)%*%ginv(gbi)%*%eblupi
      dm[j]<-dmi
     }
  }
  dmp<-dm/sum(dm)

################################################################################################

## Standardized Lesaffre and Verbeke's measure

lesverb<- rep(0,n) 
auxni=as.vector(vecni)
for (t in 1:n){ 
  li<- sum(vecni[1:t-1])+1
  ls<- sum(vecni[1:t])
  if(vecni[t]==1){
    
    auxr2 <- solve(sqrt(var.resm[li:ls,li:ls]))
    Ri<-(auxr2)%*%resm[li:ls]
    auxt<- diag(vecni[t])-Ri%*%t(Ri)
    lesverb[t]<- sum(diag(auxt%*%t(auxt)))
  }
  else
  {  
    
    auxr2 <- solve(sqrt.matrix(var.resm[li:ls,li:ls]))
    Ri<- auxr2%*%resm[li:ls]
    auxt<- diag(vecni[t])-Ri%*%t(Ri)
    lesverb[t]<- sum(diag(auxt%*%t(auxt)))
  } 
}
lesverbp<- lesverb/sum(lesverb)

########################################################################################
  
## Least confounded residuals

R.half<- sqrt.matrix(R)
auxqn<-eigen(( R.half %*% Q %*% R.half), symmetric = T, only.values = FALSE) 
lt<-t(sqrt(solve(diag((auxqn$values[1:(N-p)])))) %*% t(auxqn$vectors[1:(N-p),1:(N-p)]) %*% solve(sqrt.matrix(R[1:(N-p),1:(N-p)])) )
var.resmcp<- lt %*% var.resc[1:(N-p),1:(N-p)] %*% t(lt)
#resmc<- (lt %*% resc[1:(N-p)] )
resmcp<- (lt %*% resc[1:(N-p)] )/sqrt(diag(var.resmcp))
#resmc<- (lt %*% resc[1:(N-p)] )/sqrt(s2)

################################################################################

## Cook's  Distance

In<-diag(rep(1,N))               
Dc=matrix(c(rep(0,3*N)),N,3)
Dccond<-rep(0,N) 
Dcor<-rep(0,N) 
k<-(n-1)*q +p
for (i in 1:N){
  Ui<-In[,i] 
  auxi<- solve(crossprod(Ui,crossprod(t(Q),Ui)), crossprod(Ui,crossprod(t(Q),y)) )    
  aux2i<- solve( crossprod(X,crossprod(t(iV),X)), crossprod(X, (crossprod(t(iV),crossprod(t(Ui),auxi))) ) )  
  aux3i<- crossprod(t(Gam),crossprod(Z,crossprod(t(Q),crossprod(t(Ui),auxi)))) 
  Dc[i,1]<- crossprod(aux2i,crossprod(X,crossprod(t(X),aux2i))  )/k
  Dc[i,2]<- crossprod(aux3i,crossprod(Z,crossprod(t(Z),aux3i)))/k 
  Dc[i,3]<- crossprod(2*aux2i,crossprod(X,crossprod(t(Z),aux3i)))/k  
  Dcor[i]<- crossprod(aux2i,crossprod(X,crossprod(t(iV),crossprod(t(X),aux2i))))/p 
}
DC1<-Dc[,1];DC2<-Dc[,2];DC3<-Dc[,3]
Dccond<- DC1+DC2+DC3  

################################################################################

## Leverage

Covb<-fit$varFix

L1<-X%*%Covb%*%t(X)%*%iV                                                    # Generalized leverage marginal matrix
L1d<-diag(L1)                                                               # Leverage for observations
L1i<- as.numeric(as.matrix(aggregate(L1d, by=list(id), FUN=mean)[2]))
trace.L1d<- sum(L1d)
L<-L1+Z%*%Gam%*%t(Z)%*%Q                                                    # Generalized joint leverage matrix
Ld<-diag(L)                                                                 # Joint leverage latrix for observations
Li<- as.numeric(as.matrix(aggregate(Ld, by=list(id), FUN=mean)[2]))  
L2<-Z%*%Gam%*%t(Z)                                                          # Generalized leverage matrix for the random effect component
L2d<-diag(L2)
L2i<- as.numeric(as.matrix(aggregate(L2d, by=list(id), FUN=mean)[2]))  
trace.L2d<- sum(L2d)  

########################################################################################

## QQplot2

qqPlot2<-function(x, distribution="norm", ..., ylab=deparse(substitute(x)),
                  xlab=paste(distribution, "quantiles"), main=NULL, las=par("las"),
                  envelope=.95,  
                  col=palette()[1], col.lines=palette()[2], lwd=2, pch=1, cex=par("cex"),
                  cex.lab=par("cex.lab"),cex.axis=par("cex.axis"), 
                  line=c("quartiles", "robust", "none"), 
                  labels = if(!is.null(names(x))) names(x) else seq(along=x),
                  id.method = "y", 
                  id.n = if(id.method[1]=="identify") Inf else 0,
                  id.cex=1, id.col=palette()[1], grid=TRUE)
{
  line <- match.arg(line)
  good <- !is.na(x)
  ord <- order(x[good])
  ord.x <- x[good][ord]
  ord.lab <- labels[good][ord]
  q.function <- eval(parse(text=paste("q", distribution, sep="")))
  d.function <- eval(parse(text=paste("d", distribution, sep="")))
  n <- length(ord.x)
  P <- ppoints(n)
  z <- q.function(P, ...)
  plot(z, ord.x, type="n", xlab=xlab, ylab=ylab, main=main, las=las,cex.lab=cex.lab, cex.axis=cex.axis)
  if(grid){
    grid(lty=1, equilogs=FALSE)
    box()}
  points(z, ord.x, col=col, pch=pch, cex=cex)
  if (line == "quartiles" || line == "none"){
    Q.x <- quantile(ord.x, c(.25,.75))
    Q.z <- q.function(c(.25,.75), ...)
    b <- (Q.x[2] - Q.x[1])/(Q.z[2] - Q.z[1])
    a <- Q.x[1] - b*Q.z[1]
    abline(a, b, col=col.lines, lwd=lwd)
  }
  if (line=="robust") {
    coef <- coef(rlm(ord.x ~ z))
    a <- coef[1]
    b <- coef[2]
    abline(a, b)
  }
  conf <- if (envelope == FALSE) .95 else envelope
  zz <- qnorm(1 - (1 - conf)/2)
  SE <- (b/d.function(z, ...))*sqrt(P*(1 - P)/n)
  fit.value <- a + b*z
  upper <- fit.value + zz*SE
  lower <- fit.value - zz*SE
  if (envelope != FALSE) {
    lines(z, upper, lty=2, lwd=lwd, col=col.lines)
    lines(z, lower, lty=2, lwd=lwd, col=col.lines)
  }
  showLabels(z, ord.x, labels=ord.lab,
             id.method = id.method, id.n = id.n, id.cex=id.cex, id.col=id.col)
}


## Diagnostic plots
  
#######################################################################################  

  plotg = function(plotid){
    #cat("\n To select the graphic use plotid \n
#1- Standardized marginal residuals versus fitted values and corresponding histogram
#2- Standardized Mahalanobis distance versus unit indices
#3- Chi-squared QQ plot for Mahalanobis distance
#4- tandardized Lesaffre-Verbeke measure versus unit indices
#5- Standardized conditional residuals versus fitted values and corresponding histogram
#6- Normal QQ plot and histogram for standardized least confounded conditional residuals
#7- Cook's Conditional distance versus observation indices
#8- Cook's Conditional distance 1 (D1i) versus observation indices
#9- Cook's Conditional distance 2 (D2i) versus observation indices
#10- Cook's Conditional distance 3 (D3i) versus observation indices
#11- Generalized marginal leverage (L1) versus unit indices
#12- Generalized leverage for random components (L2) versus unit indices
#13- Generalized joint leverage (L) versus unit indices
#14- Generalized joint leverage [L1i(jj)] versus observation indices
#15- Generalized joint leverage [L2i(jj)] versus observation indices
#16- Generalized joint leverage [L3i(jj)] versus observation indices
#\n")
    cat("\n Graph plotting", plotid)
    
    if(plotid==1)
    {
      par(mfrow=c(1,2), mar=c(14, 5, 1, 2))
      plot(predm, resmp, xlab=expression(paste("Marginal fitted values")),
      ylab=expression(paste("Stand. marginal residuals")), pch=20, cex=1.0, cex.lab=1.5, cex.axis=1.3, ylim=c(-1.3*max(abs(range(resmp))),1.3*max(abs(range(resmp)))))
      abline(h=0, lty=2)
      abline(h=-limit, lty=2)
      abline(h=limit, lty=2)
      index=which(abs(resmp)>limit)
      if(length(index)>0)
      {
      text(predm[index], resmp[index], paste(id[index], obs[index], sep="."), adj=c(1,-.5), cex=.8, font=2)
      }      
      hist(resmp, freq=F, main="", xlab=expression(paste("Standardized marginal residuals")), cex=0.9, cex.lab=1.5, cex.axis=1.3)
      
    }
    if(plotid==2)
    {
      par(mfrow=c(1,1),mar=c(2.5, 2.5, 2.5, 2.5))
      plot(dmp, ylab=expression(paste("Stand. Mahalanobis dist.")), xlab="Unit index",
           pch=20, ylim=c(0,2*max(dmp)), cex=1.0, cex.lab=1.5, cex.axis=1.3)
      abline(h=2*mean(dmp), lty=2)
      index=which(dmp>2*mean(dmp))
      index1<-subject[index]
      if(length(index)>0)
      {
        text(index,dmp[index], index1, adj=c(1,-.5), cex=.8, font=2)
      }
     }
    if(plotid==3)
    {
      par(mfrow=c(1,1),mar=c(2.5, 2.5, 2.5, 2.5))
      quant.chisq<-qqPlot2(dm, distribution='chisq', df=q, pch=20, cex=1.0,cex.lab=1.5,cex.axis=1.3,
      ylab=expression(paste("Mahalanobis distance")),xlab="Chi-squared quantiles")
    }
    if(plotid==4)
    {
      par(mfrow=c(1,1),mar=c(2.5, 2.5, 2.5, 2.5))
      plot(lesverbp,ylab=expression(paste("Standardized Lesaffre-Verbeke measure")),
           xlab="Unit index", cex=1.0, cex.lab=1.5, cex.axis=1.3, pch=20, ylim=c(0,2*max(abs(range(lesverbp)))))
      abline(h=2*mean(lesverbp),lty=2)
      index=which(lesverbp>2*mean(lesverbp))
      index1<-subject[index]
      if(length(index)>0)
      {
        text(index, lesverbp[index], index1, adj=c(1,-.5), cex=.8, font=2)
      }
    }
    if(plotid==5)
    {
#      par(mfrow=c(1,2),mar=c(14, 5, 1, 2))
      par(mfrow=c(1,2))
       plot(predi, rescp, xlab=expression(paste("Valores preditos")),
           cex=1.0, cex.lab=1.5, cex.axis=1.3, ylab=expression(paste("Res. cond. padronizados")), pch=20, ylim=c(-1.3*max(abs(range(rescp))),1.3*max(abs(range(rescp)))))
      abline(h=0,lty=2)
      abline(h=limit,lty=2)
      abline(h=-limit,lty=2)
      index=which(abs(rescp)>limit)
      index1<-subject[index]
      if(length(index)>0)
      {
        text(predi[index], rescp[index], paste(id[index], obs[index],sep="."), adj=c(1,-.5), cex=.8, font=2)
      }
      hist(rescp, freq=F, main="", ylab="densidade",xlab=expression(paste("Resíduos condicionais padronizados")), cex=1.0, cex.lab=1.5, cex.axis=1.3)
    }
#    if(plotid==6)
#    {
#      par(mfrow=c(1,2),mar=c(14, 5, 1, 2))
#      qqPlot2(resmcp, ylab="Standardized least confounded residuals", xlab="N(0,1) quantiles", pch=20, cex=0.75, cex.lab=1.5,cex.axis=1.3)      
#      hist(resmcp, freq=F, xlab="Standardized least confounded residuals", main="", cex=1.0, cex.lab=1.5, cex.axis=1.3, pch=20)
#    }
    if(plotid==6)
    {
    if (option == 0)
    {
#      par(mfrow=c(2,2),mar=c(14, 5, 1, 2))
       par(fin = c(2, 1), mfrow = c(2, 2), mar = c(2.5, 2.5, 2.5, 2.5))
       plot(resmcp,ylim=c(min(-3,min(resmcp)),max(3,max(resmcp))),xlab="índice",ylab="resíduo de confundimento mínimo padronizado",cex.axis=1.4,cex.lab=1.1)
       abline(-2,0,lty=2)
       abline(2,0,lty=2)
       abline(0,0,lty=2)
       #plot(predm,resmcp,ylim=c(min(-3,min(resmcp)),max(3,max(resmcp))),xlab="valor ajustado",ylab="res?duo de confundimento m?nimo padronizado")
       boxplot(resmcp, freq=F, ylab="resíduo de confundimento mínimo padronizado", main="", cex=1.0, cex.lab=1.1, cex.axis=1.4, pch=20, horizontal = TRUE)
       hist(resmcp, freq=F, xlab="resíduo de confundimento mínimo padronizado", ylab="densidade", main="", cex=1.0, cex.lab=1.1, cex.axis=1.4, pch=20)
       qqPlot2(resmcp, ylab="resíduo de confundimento mínimo padronizado", xlab="quantis da N(0,1)", pch=20, cex=0.75, cex.lab=1.1,cex.axis=1.4)      
       }
       else
       {
       par(mfrow=c(1,1))
       qqPlot2(resmcp, ylab="resíduo de confundimento mínimo padronizado", xlab="quantis da N(0,1)", pch=20, cex=0.75, cex.lab=1.1,cex.axis=1.4)      
       }
    }

    if(plotid==7)
    {
      par(mfrow=c(1,1))
      plot(Dccond,pch=16,xlab="Observation index",ylab="Cook's conditional distance")
      abline(h=(2*mean(Dccond)),lty=2)                   ### Os pontos de corte abaixo sao definidos de forma descritiva!!!
      index=which(Dccond>(2*mean(Dccond)))
      index1<-id[index]
      if(length(index1)>0)
      {
        text(index, Dccond[index], paste(id[index], obs[index],sep="."), adj=c(1,-.5), cex=.8, font=2)
      }
    }
    if(plotid==8)
    {
      par(mfrow=c(1,1))
      plot(DC1,pch=16,xlab="Observation index",ylab="Cook's conditional distance Di1",ylim=c(0,max(DC1)))
      abline(h=(2*mean(DC1)),lty=2)
      index=which(DC1>(2*mean(DC1)))
      if(length(index)>0)
      {
        text(index, DC1[index], paste(id[index], obs[index],sep="."), adj=c(1,-.5), cex=.8, font=2)
      }
    }
    if(plotid==9)
    {
      par(mfrow=c(1,1))
      plot(DC2,pch=16,xlab="Observation index",ylab="Cook's conditional distance D2i",ylim=c(0,max(DC2)))
      abline(h=(2*mean(DC2)),lty=2)
      index=which(DC2>(2*mean(DC2)))
      if(length(index)>0)
      {
        text(index, DC2[index], paste(id[index], obs[index],sep="."), adj=c(1,-.5), cex=.8, font=2)
      }
    }
    if(plotid==10)
    {
      
      plot(abs(DC3),pch=16,xlab="Observation index",ylab="Cook's conditional distance D3i",ylim=c(0,(1.5*max(DC3)) ))
      abline(h=(2*mean(abs(DC3))),lty=2)
      index=which(abs(DC3)>(2*mean(abs(DC3))))
      if(length(index)>0)
      {
        text(index, abs(DC3)[index], paste(id[index], obs[index],sep="."), adj=c(1,-.5), cex=.8, font=2)
      }
    }
    if(plotid==11)
    {
      par(mfrow=c(1,1))
      plot(L1i,pch=16,xlab="Unit index",ylab="Generalized marginal leverage L1i: ",ylim=c(0,2*max(L1i)))
      abline(h=((2*p)/n),lty=2)
      index=which(L1i>((2*p)/n))    
      if(length(index)>0)
      {
        text(index, L1i[index], paste(subject[index], obs[index],sep="."), adj=c(1,-.5), cex=.8, font=2)
      }
    }
    if(plotid==12)
    {
      par(mfrow=c(1,1))
      plot(L2i,pch=16,xlab="Unit index",ylab="L2i",ylim=c(0,(2*max(L2i))))
      abline(h=(2*mean(L2i)),lty=2)
      index=which(L2i>((2*mean(L2i))))    
      if(length(index)>0)
      {
        text(index, L2i[index], paste(subject[index], obs[index],sep="."), adj=c(1,-.5), cex=.8, font=2)
      }
    }
    if(plotid==13)
    {
      par(mfrow=c(1,1))
      plot(Li,pch=16,xlab="Unit index",ylab="Generalized joint leverage Li",ylim=c(0,(2*max(Li))))
      abline(h=((2*mean(Li))),lty=2)
      index=which(Li>((2*mean(Li))))    
      if(length(index)>0)
      {
        text(index, Li[index], paste(subject[index], obs[index],sep="."), adj=c(1,-.5), cex=.8, font=2)
      }
    }
    if(plotid==14)
    {
      par(mfrow=c(1,1))
      plot(L1d,pch=16,xlab="Observation index",ylab="Generalized joint leverage L1i(jj)",ylim=c(0,(1.5*max(L1d))))
      abline(h=((2*p)/n),lty=2)
      index=which(L1d>((2*p)/n))
      if(length(index)>0)
      {
        text(index, L1d[index], paste(id[index], obs[index],sep="."), adj=c(1,-.5), cex=.8, font=2)
      }
    }
    if(plotid==15)
    {
      par(mfrow=c(1,1))
      plot(L2d,pch=16,xlab="Observation index",ylab="Generalized joint leverage L2i(jj)",ylim=c(0,(1.5*max(L2d))))
      abline(h=((2*mean(L2d))),lty=2)
      index=which(L2d>((2*mean(L2d))))
      if(length(index)>0)
      {
        text(index, L2d[index], paste(id[index], obs[index],sep="."), adj=c(1,-.5), cex=.8, font=2)
      }
    }
    if(plotid==16)
    {
      par(mfrow=c(1,1))
      plot(Ld,pch=16,xlab="Observation index",ylab="Generalized joint leverage L3i(jj)",ylim=c(0,(1.5*max(Ld))))
      abline(h=(2*mean(Ld)),lty=2)
      index=which(Ld>(2*mean(Ld)))
      if(length(index)>0)
      {
        text(index, Ld[index], paste(id[index], obs[index],sep="."), adj=c(1,-.5), cex=.8, font=2)
      }
    }
        
    
  }
 
############################################################################################
 
  if (is.null(plotid)) {
    cat("\n To choose plot, select plotid \n
1- Standardized marginal residuals versus fitted values and corresponding histogram
2- Standardized Mahalanobis distance versus unit indices
3- Chi-squared QQ plot for Mahalanobis distance
4- tandardized Lesaffre-Verbeke measure versus unit indices
5- Standardized conditional residuals versus fitted values and corresponding histogram
6- Normal QQ plot and histogram for standardized least confounded conditional residuals
7- Cook's Conditional distance versus observation indices
8- Cook's Conditional distance 1 (D1i) versus observation indices
9- Cook's Conditional distance 2 (D2i) versus observation indices
10- Cook's Conditional distance 3 (D3i) versus observation indices
11- Generalized marginal leverage (L1) versus unit indices
12- Generalized leverage for random components (L2) versus unit indices
13- Generalized joint leverage (L) versus unit indices
14- Generalized joint leverage [L1i(jj)] versus observation indices
15- Generalized joint leverage [L2i(jj)] versus observation indices
16- Generalized joint leverage [L3i(jj)] versus observation indices
        \n")
    return (1);  
  }
  

# Obtain plots

  for (g in plotid) {
   plotg(g)
    cat("\n Press ENTER to continue...")
    readline()
  }
  
useful.results <- list(
     std.marginal.residuals=cbind(Subject=id,Predicted=as.numeric(resmp)),
     std.conditional.residuals=cbind(Subject=id,Predicted=as.numeric(rescp)),
     std.mahalanobis.distance=cbind(Subject=as.numeric(unique(id)),std.md=as.numeric(dmp)),
     std.lesaffreverbeke.measure=cbind(Subject=as.numeric(unique(id)),std.LV.m=lesverbp),
     least.confounded.residuals=cbind(l.c.r=as.numeric(resmcp)),
     cook.conditional.distance=cbind(Subject=id,DCCond=Dccond),
     DC1=cbind(Subject=id,DC1=Dc[,1]),
     DC2=cbind(Subject=id,DC2=Dc[,2]),
     DC3=cbind(Subject=id,DC3=Dc[,3]),
     L1=cbind(Subject=subject,L1=L1i),
     L2=cbind(Subject=subject,L1=L2i),
     L=cbind(Subject=subject,L=Li),
     L1ijj=cbind(Subject=id,L1ijj=L1d),
     L2ijj=cbind(Subject=id,L2ijj=L2d),
     L3ijj=cbind(Subject=id,L3ijj=Ld)
     )

}

