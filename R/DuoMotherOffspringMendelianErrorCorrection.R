#' @title
#' Step 1 of 3 : Mendel Error detection in mother-offspring Duos and set it to missing
#' @description
#' Scan each mother-offspring duo for Mendel errors and set such sites to missing in all family members using parallel computing
#' Takes in phased vcf files of duo family data, pedigree file and number of cores (optional, default is 4 cores)
#' @param vcf Phased genotype vcf file for mother-offspring duos with full path [mandatory]
#' @param ped Tab delimited Pedigree file with full path having four columns (FamilyID,IndividualID,FatherID,MotherID) with NO header and 0 missing for missing family member [mandatory]
#' @param n Number of cores to use with parallel processing (default 4) [optional]
#' @examples
#' DuoMotherOffspringMendelError("/my/dir/data/myphased.vcf","/my/dir/data/myped.txt",2)
#' DuoMotherOffspringMendelError("/my/dir/data/myphased.vcf","/my/dir/data/myped.txt")
#' @note
#' Automatically creates dir named 'MendelErrorFixed/MotherOffsprings' one level above the data dir for the result file of this process
#' @export
#' @author Bharati Jadhav
DuoMotherOffspringMendelError = function(vcf, ped, ncore){
  mainTic = Sys.time()
  ###checking arguments
  if (missing(vcf))
      stop("Need to specify the name of phased vcf with absolute path.")

  if (missing(ped))
      stop("Need to specify the name of pedigree file with absolute path.")
  if(missing(ncore)) {
      ncore=4
      print("Number of cores not specified using 4 cores by default!")
  }


  inDir=dirname(vcf)
  ###going one level back and creating 'MendelErrorFixed' dir
  outDir=paste(dirname(inDir),"MendelErrorFixed","MotherOffsprings",sep="/")

  inVcfFile=vcf
  pedFile=ped

  print(paste("Input vcf file is : ",inVcfFile))
  print(paste("Input pedigree file is : ",pedFile))
  print(paste("Output directory is : ",outDir))

  ###create output dir if not exists
  if(dir.exists(outDir)){
      print(paste0("Output dir ",outDir," already exists, no need to create one!"))
  }else{
    ###create dir
      print(paste0("Output dir ",outDir," does not exists, creating ..... "))
      dir.create(outDir,showWarnings=FALSE,mode = "0755")
  }

  ###creating output files
  #\\. matches a dot (instead of . which matches any symbol) and $ to remove extension at the end of string
  coreName=sub("\\.[[:alnum:]]+$", "", basename(as.character(inVcfFile)))

  statFile=paste(outDir,paste0(coreName,".MendelErrorStat.txt"),sep="/")
  outVcfFile=paste(outDir,paste0(coreName,".MendelErrorFixed.vcf"),sep="/")
  library(doParallel)
  cl <- makeCluster(ncore,outfile=statFile)
  registerDoParallel(cl)
  tic = Sys.time()

  ###adding header to stat file
  cat(paste("TrioId(Child|Dad|Mom)","chr","pos","snpId","TrioGT",sep="\t"),"\n")
  print("Reading vcf file .....")

  ###check if vcf file exists
  if(file.access(inVcfFile, mode = 0) != 0){
      stop(paste("File:", inVcfFile, "does not appear to exist!"))
  }

  ###check if file is readable
  if(file.access(inVcfFile, mode = 4) != 0){
      stop(paste("File:", inVcfFile, "appears to exist but is not readable!"))
  }

  ### if everything is fine withh ped file read file
  com=paste("grep ^##",inVcfFile, " | wc -l")
  skipLines <- as.numeric(system(command=com, intern=TRUE))
  #skipLines=12
  names <- scan(inVcfFile,what = character(),skip=skipLines, nlines=1)
  duovcf = read.table(inVcfFile)
  colnames(duovcf) = names
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

  duoped = read.table(pedFile)
  print(paste("Finished Reading Pedigree file! ","Read total",nrow(duoped),"lines",sep=" "))
  colnames(duoped) = c("FID","IID","PATID","MATID")
  ###keeping only full mother-offspring duo
  duoped1= subset(duoped,PATID ==0 & MATID !=0 )
  ###find the column index where IID, PATID and MATID matched in the vcf file
  iidx=match(duoped1$IID, colnames(duovcf))
  midx=match(duoped1$MATID, colnames(duovcf))
  matchduo=na.omit(cbind(iidx,midx))
  print(paste("No. of full mother-offsprings duos : ",nrow(matchduo)))
  if(nrow(matchduo) < 1){
      stop("No duo samples to process!")
  }

  pofo = matrix(0, nrow=nrow(duovcf), ncol = 9, byrow=T)
  pofo=duovcf[,1:9]

  print("Finding Mendel error in mother-offsprings duos....")
  pofo1=pofo
  tic = Sys.time()

  #########################################################################################
  #Check Medelian Error in each duo and set the sites missing in both mothe rand offspring
  #########################################################################################

  #for each duo
  pofo=foreach(i = 1:nrow(matchduo), .combine='cbind',.inorder=TRUE) %dopar% {

      duo= duovcf[,c(1:9, matchduo[i,1],matchduo[i,2])]
      duo[,10:11]=apply(duo[,10:11],2, substr,1,3)

      errCt=0
      ###for duo
      errCt=sum((duo[,10] == "0|0" &  duo[,11] == "1|1") | (duo[,10] == "1|1" &  duo[,11] == "0|0"))
      ln= which((duo[,10] == "0|0" &  duo[,11] == "1|1") | (duo[,10] == "1|1" &  duo[,11] == "0|0") )
      if(length(ln)>0){
          duo[ln,10]=".|."
          duo[ln,11]=".|."
      }#if
      cat(paste(duo[1,1],as.character(paste(colnames(duo)[10],colnames(duo)[11],sep=":")),errCt,sep=" "),"\n")
      result = cbind(duo[,10],duo[,11])
      colnames(result)=c(colnames(duo)[10],colnames(duo)[11])

      return(result)
  }#dopar
  print ("Parallel duo processing time")
  toc=Sys.time()
  print(toc-tic)
  finalVcf=cbind(duovcf[,1:9],pofo)

  write.table(finalVcf,file=outVcfFile,row.names=F,col.names=T, quote=F, sep="\t")

  mainToc=Sys.time()
  print("Total processing time")
  print(mainToc-mainTic)
  stopCluster(cl)
  rm(list=ls())
}###end of the function

