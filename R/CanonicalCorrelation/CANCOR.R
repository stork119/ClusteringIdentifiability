#### CANCOR #### 

source(paste(dirname(parent.frame(2)$ofile),"/CANCORSM.R", sep = ""))
source(paste(dirname(parent.frame(2)$ofile),"/CANCORFIM.R", sep = ""))


CANCOR <- function(type = NULL,
                   cancor                  = NULL,
                   get_information_cancor  = NULL,
                   get_correlations_cancor = NULL){
  # Description :
  # Abstract class consisitng declaration of functions for calculating canonical correaltions 
  # Parameters  :
  #
  
  ### MAIN CODE
  
  val <- list(type                    = type,
              cancor                  = cancor,
              get_information_cancor  = get_information_cancor,
              get_correlations_cancor = get_correlations_cancor)
  class(val) <- "CANCOR"
  val
}
print.CANCOR <- function(x,...){
  cat("CANCOR :")
  print(x$type)
}

info.CANCOR <- function(x,...){
  cat("CANCOR :")
  print(x$type)
}

