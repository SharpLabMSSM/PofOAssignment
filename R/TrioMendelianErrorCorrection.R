#' @title
#' Step 1 of 3 : Mendel Error detection in Trios and set it to missing
#' @description
#' Scan each trio for Mendel errors and set such sites to missing in all family members using parallel computing.
#' Takes in phased vcf files of trio familyd data, pedigree file and number of cores (optional, default is 4 cores)
#' @author Bharati Jadhav
#' @param vcf Phased genotype vcf file for Trio family with full path [mandatory]
#' @param ped Tab delimited Pedigree file with full path having four columns (FamilyID,IndividualID,FatherID,MotherID) with NO header and 0 missing for missing family member [mandatory]
#' @param n Number of cores to use with parallel processing (default 4) [optional]
#' @return NULL
#' @examples
#' TrioMendelError("/my/dir/data/myphased.vcf","/my/dir/data/myped.txt",2)
#' TrioMendelError("/my/dir/data/myphased.vcf","/my/dir/data/myped.txt")
#' @export
#' @note
#' Automatically creates dir named 'MendelErrorFixed' one level above the data dir for the result file of this process

TrioMendelError = function(vcf, ped, ncore){
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
  outDir=paste(dirname(inDir),"MendelErrorFixed",sep="/")
  
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
  names <- scan(inVcfFile,what = character(),skip=skipLines, nlines=1)
  triovcf = read.table(inVcfFile)
  colnames(triovcf) = names
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
  print(paste("Finished Reading Pedigree file! ","Read total",nrow(trioped),"lines",sep=" "))
  colnames(trioped) = c("FID","IID","PATID","MATID")
  ###keeping only full trios
  trioped1= subset(trioped,PATID !=0 & MATID !=0 )
  print(paste("Trioped1 : ",nrow(trioped1)))
  ###find the column index where IID, PATID and MATID matched in the vcf file
  iidx=match(trioped1$IID, colnames(triovcf))
  midx=match(trioped1$MATID, colnames(triovcf))
  fidx=match(trioped1$PATID, colnames(triovcf))
  matchtrio=na.omit(cbind(iidx,fidx,midx))
  #matchtrio=cbind(iidx,fidx,midx)
  print(paste("No. of full trios : ",nrow(matchtrio)))
  #########################################################################################
  #Check Medelian Error in each trio and set the sites missing in both mothe rand offspring
  ########################################################################################
  pofo = matrix(0, nrow=nrow(triovcf), ncol = 9, byrow=T)
  pofo=triovcf[,1:9]
  print("Mendelian error for parent-offsprings trios....")
  pofo1=pofo
  tic = Sys.time()
  ###start each trio paralleley
  pofo=foreach(i = 1:nrow(matchtrio), .combine='cbind',.inorder=TRUE) %dopar% {
    trio= triovcf[,c(1:9, matchtrio[i,1],matchtrio[i,2],matchtrio[i,3])]
    trio[,10:12]=apply(trio[,10:12],2, substr,1,3)
    #error count in this ith duo
    errCt=0
    lines=0
    for(j in 1:nrow(trio))
    {
      ### for trios
      if((trio[j,11] == "0|0" &  trio[j,12] == "0|0" & trio[j,10] != "0|0") |(trio[j,11] == "1|1" &  trio[j,12] == "1|1" & trio[j,10] != "1|1")){
        
        cat(paste(as.character(paste(colnames(trio)[10],colnames(trio)[11],colnames(trio)[12],sep=":")),trio[j,1],trio[j,2],trio[j,3],paste(trio[j,10],trio[j,11],trio[j,12],sep=":"),sep="\t"),"\n")
        trio[j,10] = ".|."
        trio[j,11] = ".|."
        trio[j,12] = ".|."
        errCt=errCt+1
      }else if((trio[j,11] == "0|0" &  trio[j,12] == "1|1" & (trio[j,10] == "0|0" |trio[j,10] == "1|1"))|
                 (trio[j,11] == "1|1" &  trio[j,12] == "0|0" & (trio[j,10] == "0|0" | trio[j,10] == "1|1" ))){
        
        cat(paste(as.character(paste(colnames(trio)[10],colnames(trio)[11],colnames(trio)[12],sep=":")),trio[j,1],trio[j,2],trio[j,3],paste(trio[j,10],trio[j,11],trio[j,12],sep=":"),sep="\t"),"\n")
        trio[j,10] = ".|."
        trio[j,11] = ".|."
        trio[j,12] = ".|."
        errCt=errCt+1
      }else if((trio[j,11] == "0|0" & (trio[j,12] == "0|1" |trio[j,12] == "1|0" ) & trio[j,10] == "1|1" ) |
                 (trio[j,11] == "1|1" & (trio[j,12] == "0|1" |trio[j,12] == "1|0" ) & trio[j,10] == "0|0" )){
        
        cat(paste(as.character(paste(colnames(trio)[10],colnames(trio)[11],colnames(trio)[12],sep=":")),trio[j,1],trio[j,2],trio[j,3],paste(trio[j,10],trio[j,11],trio[j,12],sep=":"),sep="\t"),"\n")
        trio[j,10] = ".|."
        trio[j,11] = ".|."
        trio[j,12] = ".|."
        errCt=errCt+1
      }else if((trio[j,12] == "0|0" & (trio[j,11] == "0|1" |trio[j,11] == "1|0" ) & trio[j,10] == "1|1" ) |
                 (trio[j,12] == "1|1" & (trio[j,11] == "0|1" |trio[j,11] == "1|0" ) & trio[j,10] == "0|0" )){
        
        cat(paste(as.character(paste(colnames(trio)[10],colnames(trio)[11],colnames(trio)[12],sep=":")),trio[j,1],trio[j,2],trio[j,3],paste(trio[j,10],trio[j,11],trio[j,12],sep=":"),sep="\t"),"\n")
        trio[j,10] = ".|."
        trio[j,11] = ".|."
        trio[j,12] = ".|."
        errCt=errCt+1
      }
      
    }#for row
    result = cbind(trio[,10],trio[,11],trio[,12])
    colnames(result)=c(colnames(trio)[10],colnames(trio)[11],colnames(trio)[12])
    return(result)
  }#dopar
  toc=Sys.time()
  print(toc-tic)
  finalVcf=cbind(triovcf[,1:9],pofo)
  write.table(finalVcf,file=outVcfFile,row.names=F,col.names=T, quote=F, sep="\t")
  toc=Sys.time()
  print(toc-tic)
  stopCluster(cl)
  print(paste("Finished processing : ", inVcfFile ))
  #rm(list=setdiff(ls(),"mainTic"))
  mainToc=Sys.time()
  print ("Total processing time")
  print(mainToc-mainTic)
  #rm(mainTic,mainToc)
   rm(list=ls())
  return(NULL)
}##function end
