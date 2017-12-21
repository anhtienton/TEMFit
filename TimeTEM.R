TimeTEM <- function(Vmax, KM, C){
  source("fileReader.R",chdir = T)
  library("deSolve");

  
# Code description
# We solve the following system of equations
# Blac_p  :: Concentration of betalactam in periplasm
# Blac_m  :: Concentration of betalactam in the media
# Pep     :: Concentration of peptydoglycan
# TEM     :: Concentration of TEM-1

parameters <- c(kcatPBP=kcatPBP, EPBP=EPBP,CPep=Cpep,KampPBP=KampPBP, KMPBP=KMPBP, D=D, Blac_m=C, Vmax_e=Vmax_e, n=n, K_e=K_e, Vmax=Vmax, KM=KM)
state <- c(Blac_p = Blac_p_init, Pep_init = Cpep)

  TEMrun<-function(t, state, parameters) {
    with(as.list(c(state, parameters)),{
      # rate of change
      dBlac_p=D*(Blac_m - Blac_p) - (Vmax*Blac_p/(KM + Blac_p)) - (kcatPBP*EPBP*CPep)/((1+(Blac_p/KampPBP))*KMPBP+CPep) - ((Vmax_e)*(Blac_p^n))/((K_e^n)+(Blac_p^n));
      dPep=(kcatPBP*EPBP*CPep)/((1+(Blac_p/KampPBP))*KMPBP+CPep);
  
      # return the rate of change
      list(c(dBlac_p, dPep))
    }) # end with(as.list ... +}
  }

times <- seq(0, t_double, by = 1)

out <- deSolve::ode(y = state, times = times, func = TEMrun, parms = parameters)

Cp_us=tail(out[,2],1)
output=c(Cp_us);

}

printRates <- function(){
  
  parameters <- c(kcatPBP=kcatPBP, EPBP=EPBP,CPep=Cpep,KampPBP=KampPBP, KMPBP=KMPBP, D=D, Blac_m=C, Vmax_e=Vmax_e, n=n, K_e=K_e, Vmax=Vmax, KM=KM)
  state <- c(Blac_p = Blac_p_init, Pep_init = Cpep)
  
  dRate=unname(parameters[6]*(parameters[7] - state[1])) #dRate=D*(Blac_m - Blac_p)   
  
  ppbiRate=unname((parameters[1]*parameters[2]*parameters[3])/((1+(state[1]/parameters[4]))*parameters[5]+parameters[3])) #pbpiRate=(kcatPBP*EPBP*CPep)/((1+(Blac_p/KampPBP))*KMPBP+CPep)
  
  write(paste("Diffusion rate:",dRate,"\n"),fileToWrite,append=TRUE)
  write(paste("PBP Inhibition rate:",ppbiRate,"\n"),fileToWrite,append=TRUE)
}
