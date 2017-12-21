args = commandArgs(trailingOnly=TRUE)

if (length(args)!=4) {
  stop("At least one argument must be supplied. Provide option for efflux, and file location for cell specification, enzyme specification, and dosage specification.\n", call.=FALSE)
}

efflux <-  as.numeric(args[1])
cellLocation <- args[2]
enzymeLocation <- args[3]
dosageLocation <- args[4]

cat("Checking packages...\n")
list.of.packages <- c("deSolve")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages))
  install.packages(new.packages,repos = "http://cran.us.r-project.org")
  
cat("Checking files...\n")
source("fileReader.R",chdir = T)
  
cat("Solving differential...\n")
source("TimeTEM.R")  ## I call the differential solver here to be used for the efflux part

fileName <- paste("output",Sys.Date(),".txt")
fileToWrite <- file(fileName,"w")  

if (efflux==0) {

  kcatPBP=0.45          # Enzyme turnover rate of PBPs (1/s)
  EPBP=0.14;            # Concentration of folded and active PBPs (uM)
  Cpep=930;             # Concentration of peptidoglycan in the periplasm  (uM)
  KMPBP=0.042*1e6;      # Michaelis constant of PBPs  (uM)
  KampPBP=0.02;         # Dissociation constant of beta-lactam from PBPs (1/uM)
  D=0.102;              # Diffusion coefficient (1/s)
  t_exp=2*60*60        # Experimental time (seconds)
  
  Cp_mut=0.5*(C - KM -(Vmax/D) + sqrt( (-C+ KM + (Vmax/D))*(-C+ KM + (Vmax/D)) +4*C*KM));
  printRates()
  Cp_WT=0.5*(C - KM_WT - (Vmax_WT/D) + sqrt((-C+KM_WT+ (Vmax_WT/D))^2 +4*C*KM_WT));
  
  F_WT=(kcatPBP*EPBP*Cpep)/((1+(Cp_WT/KampPBP))*KMPBP+Cpep)
  F_mut=(kcatPBP*EPBP*Cpep)/((1+(Cp_mut/KampPBP))*KMPBP+Cpep)
  rel_fit=(F_mut - F_WT)*t_exp
  
} else {
  
  t_exp=2*60*60        # Experimental time (seconds)
  cat("Calculating rate of change...","\n")
  Cp_mut=TimeTEM(Vmax = Vmax, KM = KM, C = C);
  printRates()
  Cp_WT=TimeTEM(Vmax = Vmax_WT, KM =KM_WT, C = C)
  
  cat("Calculating relative fitness...","\n")
  F_WT=(kcatPBP*EPBP*Cpep)/((1+(Cp_WT/KampPBP))*KMPBP+Cpep)
  F_mut=(kcatPBP*EPBP*Cpep)/((1+(Cp_mut/KampPBP))*KMPBP+Cpep)
  rel_fit=(F_mut - F_WT)*t_exp
  
}

hRate=unname((Vmax*Cp_mut/(KM + Cp_mut))) #hRate=(Vmax*Blac_p/(KM + Blac_p))
eRate=unname(((Vmax_e)*(Cp_mut^n))/((K_e^n)+(Cp_mut^n))) #eRate=((Vmax_e)*(Blac_p^n))/((K_e^n)+(Blac_p^n))

write(paste("Hydrolysis rate:",hRate,"\n"),fileToWrite,append=TRUE)
write(paste("Efflux rate:",eRate,"\n"),fileToWrite,append=TRUE)

write(paste("Concentration of beta-lactam in the periplasm:",Cp_mut,"\n"),file = fileToWrite,append=TRUE)
write(paste("relative fitness of the mutant to WT:",rel_fit,"\n"),fileToWrite,append=TRUE)  

cat("Done.\n")