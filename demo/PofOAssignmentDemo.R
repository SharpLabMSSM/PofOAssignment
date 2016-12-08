#' Demo of how to run the PofOAssignment pipeline.

###load the package
library("PofOAssignment")

base.dir= find.package("PofOAssignment")
### phased genotype vcf file with full path
phasedVcf=paste0(base.dir, "/inst/extdata/PhasedData/TestPhasedData.vcf")
### pedigree file with full path
ped=paste0(base.dir, "/inst/extdata/PhasedData/TestPed.txt")

### Trio family data ###
cat("\n ########## Trio PofOrigin Assignment process ######### \n")

cat("\n Trio:Processing Step1 .....\n")
### Step1: phased data to detect Mendel Error with 2 cores 
TrioMendelError(phasedVcf,ped,2)

cat("\n Trio:Processing Step2 .....\n")
### Step2: Sliding window to fix phase switch error with default 4 core
phasedVcf1=paste0(base.dir,"/inst/extdata/MendelErrorFixed/TestPhasedData.MendelErrorFixed.vcf")
TrioPhaseCorrect(phasedVcf1,ped)

cat("\n Trio:Processing Step3 .....\n")
### Step3 : Recover missing data in slidind window at het sites in offsprings before sliding window
phasedVcf2=paste0(base.dir,"/inst/extdata/PhaseSwitchSlidingWindow/TestPhasedData.MendelErrorFixed.PhaseSwitch.vcf")
TrioRecoverMissing(phasedVcf1,phasedVcf2,ped)
cat("\n Trio: Finished assigning parent of origin to sites in offsprings !\n")


### Mother-Offsprings Duo family data ###

cat("\n ########## Duo:Mother-Offspring PofOrigin Assignment process ######### \n")
cat("\n Duo-Mother-Offspring:Processing Step1 .....\n")
### Step1: phased data to detect Mendel Error with 2 cores
DuoMotherOffspringMendelError(phasedVcf,ped,2)

cat("\n Duo-Mother-Offspring:Processing Step2 .....\n")
### Step2: Sliding window to fix phase switch error with default 4 core
phasedVcf1=paste0(base.dir,"/inst/extdata/MendelErrorFixed/MotherOffsprings/TestPhasedData.MendelErrorFixed.vcf")
DuoMotherOffspringPhaseCorrect(phasedVcf1,ped)

cat("\n Duo-Mother-Offspring:Processing Step3 .....\n")
### Step3 : Recover missing data in slidind window at het sites in offsprings before sliding window
phasedVcf2=paste0(base.dir,"/inst/extdata/PhaseSwitchSlidingWindow/MotherOffsprings/TestPhasedData.MendelErrorFixed.PhaseSwitch.vcf")
DuoMotherOffspringRecoverMissing(phasedVcf1,phasedVcf2,ped)
cat("\n Duo-Mother-Offspring: Finished assigning parent of origin to sites in offsprings .....\n")


### Father-Offsprings Duo family data ###
cat("\n ########## Duo:Father-Offspring PofOrigin Assignment process ######### \n")
cat("\n Duo-Father-Offspring:Processing Step1 .....\n")
### Step1: phased data to detect Mendel Error with 2 cores
DuoFatherOffspringMendelError(phasedVcf,ped,2)

cat("\n Duo-Father-Offspring:Processing Step2 .....\n")
### Step2: Sliding window to fix phase switch error with default 4 core
phasedVcf1=paste0(base.dir,"/inst/extdata/MendelErrorFixed/FatherOffsprings/TestPhasedData.MendelErrorFixed.vcf")
DuoFatherOffspringPhaseCorrect(phasedVcf1,ped)

cat("\n Duo-Father-Offspring:Processing Step3 .....\n")
### Step3 : Recover missing data in slidind window at het sites in offsprings before sliding window
phasedVcf2=paste0(base.dir,"/inst/extdata/PhaseSwitchSlidingWindow/FatherOffsprings/TestPhasedData.MendelErrorFixed.PhaseSwitch.vcf")
DuoFatherOffspringRecoverMissing(phasedVcf1,phasedVcf2,ped)
cat("\n Duo-Father-Offspring: Finished assigning parent of origin to sites in offsprings .....\n")

