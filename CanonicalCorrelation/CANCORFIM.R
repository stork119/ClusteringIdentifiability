# Implementation from internet.
# License GPL
#### CANCOR FIM #### 

CANCOR.FIM <- function(){
  # Description :
  # Class CANCOR.FIM inherits from class CANCOR.
  # Class implements function that calculate canonical correaltion for input
  # being Fisher Information Matrix
  # Parameters  :
  #
  
  ### MAIN CODE
  cancorFIM  <-function(covin, xin, yin){

    # functions calculates cononical correlations from the FIM
    # Parameters:
    # covin -- FIM,
    # xin   -- indices from matrix covin that denoting set A,
    # yin   -- indices from matrix covin that denoting set B

    eps=10^{-6} # defining  numerical zero
    
    #SVD of submatrices AA and BB
    E1<-svd(covin[xin,xin,drop=F],LINPACK=TRUE)
    E2<-svd(covin[yin,yin,drop=F],LINPACK=TRUE)
    
    # finding actual dimensions of the matrices
    pin=ncol(covin)
    p1=sum(E1$d>eps)
    p2=sum(E2$d>eps)
    p=p1+p2
    
    
    # converting the FIM accordingly 
    covxx=t(E1$v[,1:p1])%*%covin[xin,xin]%*%(E1$v[,1:p1])
    covyy=t(E2$v[,1:p2])%*%covin[yin,yin]%*%(E2$v[,1:p2])
    covxy=t(E1$v[,1:p1])%*%covin[xin,yin]%*%(E2$v[,1:p2])
    covyx=t(E2$v[,1:p2])%*%covin[yin,xin]%*%(E1$v[,1:p1])
    cov=matrix(0,p,p)
    cov[1:p1,1:p1]=covxx
    cov[p1+1:p2,p1+1:p2]=covyy
    cov[1:p1,p1+1:p2]=covxy
    cov[p1+1:p2,1:p1]=covyx
    x=1:p1
    y=p1+1:p2
    
    
    # conventional calculation of the canonical correations
    if(is.null(colnames(cov))) 
      rownames(cov) <- colnames(cov) <- as.character(1:p)
    y <- (1:p)[y]
    eigx <- svd(cov[x,x,drop=F],LINPACK=TRUE)
    eigy <- svd(cov[y,y,drop=F],LINPACK=TRUE)
    rownames(eigx$v) <- rownames(cov)[x]
    rownames(eigy$v) <- rownames(cov)[y]
    sx <- eigx$v %*% diag(1/sqrt(eigx$d),p1,p1) %*% t(eigx$v)
    sy <- eigy$v %*% diag(1/sqrt(eigy$d),p2,p2) %*% t(eigy$v)
    A <- crossprod(sy %*% cov[y,x] %*% sx)
    B <- crossprod(sx %*% cov[x,y] %*% sy)
    eigA <- svd(A,LINPACK=TRUE)
    eigB <- svd(B,LINPACK=TRUE)
    m <- min(length(x),length(y))
    xcoef <- sx%*%eigA$v
    xcoef <- sign(xcoef[1,1]) * xcoef
    ycoef <- sy%*%eigB$v
    ycoef <- sign(ycoef[1,1]) * ycoef
    colnames(xcoef) <- rep("",ncol(xcoef))
    colnames(ycoef) <- rep("",ncol(ycoef))
    colnames(ycoef)[1:m] <- colnames(xcoef)[1:m] <- paste("Can",1:m,sep=".")
    sd <- diag(1/sqrt(diag(cov)[c(x,y)]))
    rownames(sd) <- colnames(sd) <- colnames(cov)[c(x,y)]
    xcontr <- t(xcoef[,1:m,drop=F]) %*% cov[x,,drop=F] %*% sd
    ycontr <- t(ycoef[,1:m,drop=F]) %*% cov[y,,drop=F] %*% sd
    list(cor=sqrt(eigA$d[1:m]), xcoef=xcoef, ycoef=ycoef,
         xcontr=t(xcontr), ycontr=t(ycontr))
    
  }
  
  get_information_cancor <- function(m,
                                               c1,
                                               c2) {
    ###TODO - set as class parameter
    infcontrol=0.99999
    ccor = cancorFIM(m, c1, c2)
    mean(-log(1-min(infcontrol^2,ccor$cor*ccor$cor)))
  }
  
  get_correlations_cancor <- function(m,
                                                c1,
                                                c2) {
    # Function returns vector of cannonical correlations between components c1 and c2 of the matrix m
    # this is interface function for easier use of CanCor_FIM
    ###TODO - set as class parameter
    one=0.9999999999 #numerical 1
    ccor = cancorFIM(m,c1,c2)
    ccor=pmin(as.vector(ccor$cor),one)
    ccor
  }
  
  val <- CANCOR("FIM",
                cancor                  = cancorFIM,
                get_information_cancor  = get_information_cancor,
                get_correlations_cancor = get_correlations_cancor
  )
  
  class(val) <- append(class(val),
                       "CANCORFIM")
  return(val)
}
