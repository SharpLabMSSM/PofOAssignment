#' @title
#' Step 2 of 3: Sliding window algorithm to detect and fix phase switch error in mother-offspring duos
#' @description
#' Scan each mother-offpsring duo for to detect and fix phase switch error in offspring using sliding window (size 100 and step 50 sites) using parallel computing
#' Takes in phased vcf files of trio familyd data, pedigree file and number of cores (optional, default is 4 cores)
#' @param vcf Mendel error fixed Phased genotype vcf file for Duo family with full path [mandatory]
#' @param ped Tab delimited Pedigree file with full path having four columns (FamilyID,IndividualID,FatherID,MotherID) with NO header and 0 missing for missing family member [mandatory]
#' @param n Number of cores to use with parallel processing (default 4) [optional]
#' @examples
#' DuoMotherOffspringPhaseCorrect("/my/dir/MendelErrorFixed/MotherOffsprings/myphased.MendelErrorFixed.vcf","/my/dir/data/myped.txt",2)
#' DuoMotherOffspringPhaseCorrect("/my/dir/MendelErrorFixed/MotherOffsprings/myphased.MendelErrorFixed.vcf","/my/dir/data/myped.txt")
#' @author Bharati Jadhav
#' @export
#' @note
#' Automatically creates dir named 'PhaseSwitchSlidingWindow/MotherOffsprings' one level above the data dir for the result file of this process
DuoMotherOffspringPhaseCorrect = function(vcf, ped, ncore){
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
    outDir=paste(dirname(dirname(inDir)),"PhaseSwitchSlidingWindow","MotherOffsprings",sep="/") 
    inVcfFile=vcf
    pedFile=ped

    print(paste("Input vcf file is : ",inVcfFile))
    print(paste("Input pedigree file is : ",pedFile))
    print(paste("Output directory is : ",outDir))

    ###create output dir if not exists
    if(dir.exists(outDir)){
        print(paste0("Output dir ",outDir," already exists!"))
    }else{
        ###create dir
        print(paste0("Output dir ",outDir," does not exists, creating ..... "))
        dir.create(outDir,showWarnings=FALSE,mode = "0755",recursive=TRUE)
    }
    ###creating output files
    #\\. matches a dot (instead of . which matches any symbol) and $ to remove extension at the end of string
    coreName=sub("\\.[[:alnum:]]+$", "", basename(as.character(inVcfFile)))

    statFile=paste(outDir,paste0(coreName,".PhaseSwitchErrorStat.txt"),sep="/")
    outVcfFile=paste(outDir,paste0(coreName,".PhaseSwitch.vcf"),sep="/")

    library(doParallel)
    cl <- makeCluster(16,outfile=statFile)
    registerDoParallel(cl)
    tic = Sys.time()

    ###adding header to stat file
    cat(paste("chr","TrioId","N_MissSite","%MissSite","PossibleSwitches",sep=" "),"\n")

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
    #skipLines <-0
    names <- scan(inVcfFile,what = character(),skip=skipLines, nlines=1)
    duovcf = read.table(inVcfFile)
    colnames(duovcf) = names
    print(paste("Finished Reading vcf file! ","Read total",nrow(duovcf),"lines",sep=" "))

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
    colnames(duoped) = c("FID","IID","PATID","MATID")

    #extract only those rows that have child, father and mother info
    duoped1= subset(duoped,MATID !=0 & PATID ==0 )
    iidx=match(duoped1$IID, colnames(duovcf))
    midx=match(duoped1$MATID, colnames(duovcf))
    matchduo=na.omit(cbind(iidx,midx))
    print(paste("No. of mother offspring duos : ",nrow(matchduo)))

    pofo = matrix(0, nrow=nrow(duovcf), ncol = 9, byrow=T)
    pofo=duovcf[,1:9]

    #################################################
    #Check switches and fix it using sliding window
    #################################################
    print("Finding Phase switch error for mother-offsprings pairs....")
    pofo1=pofo
    tic = Sys.time()
    #open and read the duos which has switch error
    print("Finding Phase switch error for mother-offsprings duos....")
    pofo=foreach(i = 1:nrow(matchduo), .combine='cbind',.inorder=TRUE) %dopar% {

        op=matrix(0,nrow=nrow(duovcf),ncol=2)
        duo= duovcf[,c(1:9, matchduo[i,1],matchduo[i,2])]
        duo[,10:11]=apply(duo[,10:11],2, substr,1,3)
        wSize=100
        step=50
        ct=0
        errCt=0
        haplo1 = matrix(0, nrow=1, ncol = 3, byrow=T)
        colnames(haplo1)=c("WindowNo","C1MA1","C1MA2")
        haplo2 = matrix(0, nrow=1, ncol = 3, byrow=T)
        colnames(haplo2)=c("WindowNo","C2MA1","C2MA2")
        from=wSize
        startIdx=1
        endIdx=wSize
        to=nrow(duo)+(wSize-1)
        idx=1

        while(startIdx <= nrow(duovcf))
        {

              if(endIdx < nrow(duo)){
                  rows=c(startIdx:endIdx)
                  r=rows[1:step]
              }else if(endIdx > nrow(duo)){
                  rows=c(startIdx:nrow(duo))
                  wSize=length(rows)
                  r=rows
                  step=length(r)
              }#else if

              CHA1=substr(duo[rows,10],1,1)
              CHA2=substr(duo[rows,10],3,3)
              MatA1=substr(duo[rows,11],1,1)
              MatA2=substr(duo[rows,11],3,3)

              A1Winscore = c(startIdx,sum(CHA1==MatA1),sum(CHA1==MatA2))

              haplo1=rbind(haplo1,A1Winscore)

              if(length(which(A1Winscore[2:3] == wSize)) > 0){
                    colstat=as.vector(which(A1Winscore[2:3] == wSize))
                    colstat=colstat+1
              }else {
                    colstat=0

              }
              A1WFilter=as.matrix(haplo1[idx+1,c(1,colstat)])
              if(nrow(A1WFilter)>2){
                  #As this is mother-offspring duo we do not have Father's alleles, following was commented, and considred condition if CA1 matches with both alleles of mothers MA1 and MA2 ,
                  #for example CA1=0 amd MA1=0 and MA2=0 for whole wondow
                  srd1 = ifelse(length(grep("MA", rownames(A1WFilter))) > 1,"C1MA1","NA")
              }else if(nrow(A1WFilter)==2){
                  srd1=rownames(A1WFilter)[2]
              }else {srd1="NA"}

              A2Winscore = c(startIdx,sum(CHA2==MatA1),sum(CHA2==MatA2))

              haplo2=rbind(haplo2,A2Winscore)
              if(length(which(A2Winscore[2:3] == wSize)) > 0){
                  colstat=as.vector(which(A2Winscore[2:3] == wSize))
                  colstat=colstat+1
              }else { colstat=0  }

              A2WFilter=as.matrix(haplo2[idx+1,c(1,colstat)])

              if(nrow(A2WFilter)>2){
                  srd2 = ifelse(length(grep("MA", rownames(A2WFilter))) > 1,"C2MA2","NA")
              }else if(nrow(A2WFilter)==2){
                  srd2=rownames(A2WFilter)[2]
              }else {srd2="NA"}

              if((srd1 == "NA") & (srd2 == "C2MA1" | srd2 == "C2MA2" )){
                  op[r,2]=paste("DA",CHA1[1:step],sep=":")
                  op[r,1]= paste("MA",CHA2[1:step],sep=":")
              }else if((srd1 == "C1MA1" | srd1 == "C1MA2") & ( srd2 == "NA")){
                  op[r,1]= paste("MA",CHA1[1:step],sep=":")
                  op[r,2]=paste("DA",CHA2[1:step],sep=":")
              }else if((srd1 == "C1MA1" | srd1 == "C1MA2") & ( srd2 == "C2MA1" | srd2 == "C2MA2")){
                  op[r,1]= paste("MA",CHA1[1:step],sep=":")
                  op[r,2]=paste("DA",CHA2[1:step],sep=":")
              }else{ op[r,]= "." }

              if(srd1 == "NA" & srd2 =="NA"){
                  ct = ct+length(r)
                  errCt=errCt+1
              }

              startIdx= startIdx + step
              endIdx = startIdx+(wSize-1)
              idx=idx+1
      } #while  window
      op[,1]=ifelse(op[,1] != ".",substr(op[,1],4,4),".")
      op[,2]=ifelse(op[,2] != ".",substr(op[,2],4,4),".")

      cat(paste(duovcf[1,1],as.character(paste(colnames(duo)[10],colnames(duo)[11],sep=":")),ct,(ct/nrow(duovcf))*100,errCt,sep=" "),"\n")
      result = paste(op[,1],op[,2],sep="|")
      result = cbind(result,duo[,11])
      colnames(result)=c(colnames(duo)[10],colnames(duo)[11])
      return(result)
  } #for  all duos
  print("Finished Sliding window scan!")
  print("Parallel processing duo time")
  toc=Sys.time()
  print(toc-tic)
  finalVcf=cbind(pofo1,pofo)
  write.table(finalVcf,file=outVcfFile,row.names=F,col.names=T, quote=F, sep="\t")
  print("Total processing time")
  mainToc=Sys.time()
  print(mainToc-mainTic)
  stopCluster(cl)
  rm(list=ls())
}### end of the function
