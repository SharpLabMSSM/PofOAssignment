#' @title
#' Step 3 of 3 : Recover het sites in offsprings set to missing during sliding window (step2) in Trios
#' @description
#' Compare offspring genotypes in imputed QC data and phase corrected data at het sites in offspring and if is missing in phase corrected data, use Medelian inheritance rule.
#' Takes in phased vcf files of trio familyd data before sliding window and after sliding windo, and pedigree file.
#' @param vcf1 Phased genotype vcf file for family before sliding window step (Step2) but after Mendel error (Step1)
#' @param vcf2 Phased genotype vcf file (after phasing and QC and before sliding window step)
#' @param ped Tab delimited Pedigree file with four columns(FamilyID, IndividualID, FatherID, MotherID) with NO header.
#' @examples
#' TrioRecoverMissing("/my/dir/MendelErrorFixed/myphased.MendelErrorFixed.vcf","/my/dir/seSwitchSlidingWindown/myphased.MendelErrorFixed.PhaseSwitch.vcf","/my/dir/data/myped.txt")
#' @return NULL
#' @export
#' @note
#' Automatically creates dir named 'PhaseSwitchMissingFix' one level above the data dir for the result file of this process
#' @author Bharati Jadhav
TrioRecoverMissing = function(vcf1,vcf2, ped){
  mainTic = Sys.time()
  ###checking arguments
  if (missing(vcf1))
      stop("Need to specify the name of phased vcf with absolute path.")

  if (missing(vcf2))
    stop("Need to specify the name of phased vcf with absolute path.")

  if (missing(ped))
      stop("Need to specify the name of pedigree file with absolute path.")

  inDir2=dirname(vcf2)
  ### creating output dir
  outDir=paste(dirname(inDir2),"PhaseSwitchMissingFix",sep="/")

  inVcfFile1=vcf1
  inVcfFile2=vcf2
  pedFile=ped

  print(paste("Input vcf file1 is : ",inVcfFile1))
  print(paste("Input vcf file2 is : ",inVcfFile2))
  print(paste("Input pedigree file is : ",pedFile))
  print(paste("Output directory is : ",outDir))

  ###create output dir if not exists
  if(dir.exists(outDir)){
    print(paste0("Output dir ",outDir," already exists!"))
  }else{
    ###create dir
    print(paste0("Output dir ",outDir," does not exists, creating ..... "))
    dir.create(outDir,showWarnings=FALSE,mode = "0755")
  }

  ###creating output files
  #\\. matches a dot (instead of . which matches any symbol) and $ to remove extension at the end of string
  coreName=sub("\\.[[:alnum:]]+$", "", basename(as.character(inVcfFile2)))

  statFile=paste(outDir,paste0(coreName,".TrioPhaseSwitchFillMissStat.txt"),sep="/")
  outVcfFile=paste(outDir,paste0(coreName,".FillMiss.vcf"),sep="/")

  ###adding header to stat file
  ln=paste("chr","TrioId(Child|Dad|Mom)","hetCt_before","hetCt_after","%ofIncreaseInHetSite",sep=" ")
  cat(paste(ln,"\n"),file=statFile, append=TRUE)

  print("Reading vcf before sliding window .....")

  ###traditional file reading approach an alternative to readVcf but takes 4X time than readVcf()
  ###check if vcf file exists
  if(file.access(inVcfFile1, mode = 0) != 0){
      stop(paste("File:", inVcfFile1, "does not appear to exist!"))
  }
  ###check if file is readable
  if(file.access(inVcfFile1, mode = 4) != 0){
      stop(paste("File:", inVcfFile1, "appears to exist but is not readable!"))
  }

  ### if everything is fine with vcf then vcf file read file
  com=paste("grep ^##",inVcfFile1, " | wc -l")
  skipLines <- as.numeric(system(command=com, intern=TRUE))
  #skipLines=12
  names <- scan(inVcfFile1,what = character(),skip=skipLines, nlines=1)
  sourceVcf = read.table(inVcfFile1)
  colnames(sourceVcf) = names

  print("Reading Pedigree file .....")
  ###check if ped file exists
  if(file.access(pedFile, mode = 0) != 0){
      stop(paste("File:", pedFile, "does not appear to exist!"))
  }

  ###check if file is readable
  if(file.access(pedFile, mode = 4) != 0){
      stop(paste("File:", pedFile, "appears to exist but is not readable!"))
  }
  ### if everything is fine withh ped file read file
  ###tab delimited file without header line and having 4 columns, Family ID, Individual ID, Paternal ID, Maternal ID without header
  trioped = read.table(pedFile)
  print("Finished Reading Pedigree file! ")
  colnames(trioped) = c("FID","IID","PATID","MATID")
  trioped1= subset(trioped,PATID !=0 & MATID !=0 )

  #find the mother offspring column index for filtered imputed data
  iidxSource=match(trioped1$IID, colnames(sourceVcf))
  midxSource=match(trioped1$MATID,colnames(sourceVcf))
  fidxSource=match(trioped1$PATID,colnames(sourceVcf))
  matchtrioSource=na.omit(cbind(iidxSource,fidxSource,midxSource))
  print(paste("No. of full trios in source vcf: ",nrow(matchtrioSource)))
  sourceOrdered=matrix(0, nrow=nrow(sourceVcf), ncol = 9, byrow=T)
  sourceOrdered=sourceVcf[,1:9]
  print("Ordering source vcf in pedigree order [child,father,mother]  .....")
  for(i in 1:nrow(matchtrioSource))
  {
      trio= sourceVcf[,c(1:9, matchtrioSource[i,1],matchtrioSource[i,2],matchtrioSource[i,3])]
      trio[,10:12]=apply(trio[,10:12],2, substr,1,3)

      sourceOrdered = cbind(sourceOrdered,trio[,10])
      colnames(sourceOrdered)[ncol(sourceOrdered)]=colnames(trio)[10]
      #add mother's genotypes
      sourceOrdered = cbind(sourceOrdered,trio[,11])
      colnames(sourceOrdered)[ncol(sourceOrdered)]=colnames(trio)[11]
      #add mother's genotypes
      sourceOrdered = cbind(sourceOrdered,trio[,12])
      colnames(sourceOrdered)[ncol(sourceOrdered)]=colnames(trio)[12]

  }
  rm(sourceVcf)
  rm(trio)
  print("Reading Target vcf file :file after sliding window .....")

  ###check if vcf file exists
  if(file.access(inVcfFile2, mode = 0) != 0){
      stop(paste("File:", inVcfFile2, "does not appear to exist!"))
  }
  ###check if file is readable
  if(file.access(inVcfFile2, mode = 4) != 0){
      stop(paste("File:", inVcfFile2, "appears to exist but is not readable!"))
  }

  ### if everything is fine with vcf then vcf file read file
  com=paste("grep ^##",inVcfFile2, " | wc -l")
  skipLines <- as.numeric(system(command=com, intern=TRUE))
  #skipLines <- 0
  names2 <- scan(inVcfFile2,what = character(),skip=skipLines, nlines=1)
  targetVcf = read.table(inVcfFile2)
  colnames(targetVcf) = names2
  #find the mother offspring column index for phase corrected data
  iidxTarget=match(trioped1$IID, colnames(targetVcf))
  fidxTarget=match(trioped1$PATID, colnames(targetVcf))
  midxTarget=match(trioped1$MATID,colnames(targetVcf))
  matchtrioTarget=na.omit(cbind(iidxTarget,fidxTarget,midxTarget))
  print(paste("No. of full trios in target vcf : ",nrow(matchtrioTarget)))
  #now compare each mothee-offspring genotypes in imputed, filtered data and phase corrected data(using sliding window algorithm).
  #if any snp is set to missing in phase corrected data check that site in imputed filter data and if child is heterozygous and
  #mother is homozygous fill the missing genotype with these genotypes
  fixedTarget=matrix(0, nrow=nrow(targetVcf), ncol = 9, byrow=T)
  fixedTarget=targetVcf[,1:9]
  print("Fixing the het sites in child which are set to missing during sliding window .....")
  for(i in 1:nrow(matchtrioTarget))
  {
      #get the duo genotype information from the ordered source vcf
      hetCt_before=0
      hetCt_after =0
      fixrows=0
      trioSO= sourceOrdered[,c(1:9,matchtrioTarget[i,1],matchtrioTarget[i,2],matchtrioTarget[i,3])]
      trioTarget=targetVcf[,c(1:9, matchtrioTarget[i,1],matchtrioTarget[i,2],matchtrioTarget[i,3])]
      #beforeCt=sum(duoTarget[,10] == ".|.")
      hetCt_before=sum(trioTarget[,10] == "0|1" |trioTarget[,10] == "1|0" )

      fixrows=which((trioTarget[,10] == ".|.")&(trioSO[,10] == "1|0" | trioSO[,10] == "0|1" ) & (trioSO[,11] == "1|1" | trioSO[,11] == "0|0" | trioSO[,12] == "1|1" | trioSO[,12] == "0|0"))

      ### find the het rows where the child genotype could be 0|1
      fixrows01=which(((trioSO[fixrows,12] == "0|0" )&(trioSO[fixrows,11] == "1|1" | trioSO[fixrows,11] == "0|1" | trioSO[fixrows,11] == "1|0")) | ((trioSO[fixrows,12] == "0|1" | trioSO[fixrows,12] == "1|0")& (trioSO[fixrows,11] == "1|1")))
      ### find the het rows where the child genotype could be 1|0
      fixrows10=which(((trioSO[fixrows,12] == "1|1") & (trioSO[fixrows,11] == "0|0" | trioSO[fixrows,11] == "0|1" | trioSO[fixrows,11] == "1|0"))|((trioSO[fixrows,12] == "0|1" | trioSO[fixrows,12] == "1|0") & (trioSO[fixrows,11] == "0|0")))
      ###get the subset of het site row where child's genotype could be 0|1 and set it to 0|1
      trioTarget[fixrows[fixrows01],10] = "0|1"
      ###get the subset of het site row where child's genotype could be 1|0 and set it to 1|0
      trioTarget[fixrows[fixrows10],10] = "1|0"
      hetCt_after=sum(trioTarget[,10] == "0|1" |trioTarget[,10] == "1|0" )

      ln=paste(trioTarget[1,1],as.character(paste(colnames(trioTarget)[10],colnames(trioTarget)[11],colnames(trioTarget)[12],sep=":")),hetCt_before,hetCt_after,round((((hetCt_after-hetCt_before)*100)/hetCt_before),2),sep=" ")

      cat(paste(ln,"\n"),file=statFile, append=TRUE)

      fixedTarget=cbind(fixedTarget,trioTarget[,10])
      colnames(fixedTarget)[ncol(fixedTarget)]=colnames(trioTarget)[10]

      fixedTarget=cbind(fixedTarget,trioTarget[,11])
      colnames(fixedTarget)[ncol(fixedTarget)]=colnames(trioTarget)[11]

      fixedTarget=cbind(fixedTarget,trioTarget[,12])
      colnames(fixedTarget)[ncol(fixedTarget)]=colnames(trioTarget)[12]

  }

  write.table(fixedTarget,file=outVcfFile,row.names=F,col.names=T, quote=F, sep="\t")

  mainToc=Sys.time()
  print("Total processing time")
  print(mainToc-mainTic)
  rm(list = ls())
}###function end
