#### CANCOR SM #### 
CANCOR.SM <- function(){
  # Description :
  # Class CANCOR.SM inherits from class CANCOR.
  # Class implements function that calculate canonical correaltion for input
  # being Sensitvity Matrix
  # Parameters  :
  #
 
  ### MAIN CODE
  
  get_information_cancor <- function(m,
                                     c1,
                                     c2) {
    ###TODO - set as class parameter
    infcontrol=0.99999
    ccor = cancor(m[,c1], m[,c2], FALSE, FALSE)
    mean(-log(1-min(infcontrol^2,ccor$cor*ccor$cor)))
  }
  
  get_correlations_cancor <- function(m,
                                               c1,
                                               c2) {
    ccor <- cancor(m[,c1], m[,c2], FALSE, FALSE)
    ccor$cor
  }
  
  val <- CANCOR("SM",
                cancor                  = cancor,
                get_information_cancor  = get_information_cancor,
                get_correlations_cancor = get_correlations_cancor
  )
  
  class(val) <- append(class(val),
                       "CANCORSM")
  return(val)
}
