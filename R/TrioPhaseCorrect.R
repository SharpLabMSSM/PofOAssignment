#' @title
#' Step 2 of 3: Sliding window algorithm to detect and fix phase switch error
#' @description
#' Scan each trio for to detect and fix phase switch error in offspring using sliding window (size 100 and step 50 sites) using parallel computing.
#' Takes in phased vcf files of trio familyd data, pedigree file and number of cores (optional, default is 4 cores)
#' @param vcf Mendel error fixed Phased genotype vcf file for Trio family with full path [mandatory]
#' @param ped Tab delimited Pedigree file with full path having four columns (FamilyID,IndividualID,FatherID,MotherID) with NO header and 0 missing for missing family member [mandatory]
#' @param n Number of cores to use with parallel processing (default 4) [optional]
#' @examples
#' TrioPhaseCorrect("/my/dir/MendelErrorFixed/myphased.MendelErrorFixed.vcf","/my/dir/data/myped.txt",2)
#' TrioPhaseCorrect("/my/dir/MendelErrorFixed/myphased.MendelErrorFixed.vcf","/my/dir/data/myped.txt")
#' @author Bharati Jadhav
#' @export
#' @note
#' Automatically creates dir named 'PhaseSwitchSlidingWindow' one level above the data dir for the result file of this process
TrioPhaseCorrect = function(vcf, ped, ncore){
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
    outDir=paste(dirname(inDir),"PhaseSwitchSlidingWindow",sep="/")

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
        dir.create(outDir,showWarnings=FALSE,mode = "0755")
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
    #triovcf = read.vcfR(inVcfFile)
    #print(paste("Finished Reading vcf file! ","Read total",nrow(triovcf),"lines",sep=" "))

    ###traditional file reading approach an alternative to readVcf but takes 4X time than readVcf()
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

    #################################################
    #Check switches and fix it using sliding window
    #################################################

    pofo = matrix(0, nrow=nrow(triovcf), ncol = 9, byrow=T)
    pofo=triovcf[,1:9]

    print("Finding Phase switch error for parent-offsprings trios....")
    pofo1=pofo
    tic = Sys.time()

    pofo=foreach(i = 1:nrow(matchtrio), .combine='cbind',.inorder=TRUE) %dopar% {

        op=matrix(0,nrow=nrow(triovcf),ncol=2)
        trio= triovcf[,c(1:9, matchtrio[i,1],matchtrio[i,2],matchtrio[i,3] )]
        trio[,10:12]=apply(trio[,10:12],2, substr,1,3)
        wSize=100
        step=50
        ct=0
        errCt=0
        haplo1 = matrix(0, nrow=1, ncol = 5, byrow=T)
        colnames(haplo1)=c("WindowNo","C1MA1","C1MA2", "C1DA1","C1DA2")
        haplo2 = matrix(0, nrow=1, ncol = 5, byrow=T)
        colnames(haplo2)=c("WindowNo","C2MA1","C2MA2","C2DA1","C2DA2")
        from=wSize
        startIdx=1
        endIdx=wSize
        to=nrow(trio)+(wSize-1)
        idx=1

        while(startIdx <= nrow(triovcf))
        {

            if(endIdx < nrow(trio)){
                rows=c(startIdx:endIdx)
                r=rows[1:step]
            }else if(endIdx > nrow(trio)){
                rows=c(startIdx:nrow(trio))
                wSize=length(rows)
                r=rows
                step=length(r)
            } #else if

            CHA1=substr(trio[rows,10],1,1)
            CHA2=substr(trio[rows,10],3,3)
            PatA1=substr(trio[rows,11],1,1)
            PatA2=substr(trio[rows,11],3,3)
            MatA1=substr(trio[rows,12],1,1)
            MatA2=substr(trio[rows,12],3,3)

            A1Winscore = c(startIdx,sum(CHA1==MatA1),sum(CHA1==MatA2),sum(CHA1==PatA1),sum(CHA1==PatA2))
            haplo1=rbind(haplo1,A1Winscore)

            if(length(which(A1Winscore[2:5] == wSize)) > 0){
                colstat=as.vector(which(A1Winscore[2:5] == wSize))
                colstat=colstat+1
            }else { colstat=0 }
                A1WFilter=as.matrix(haplo1[idx+1,c(1,colstat)])

            A2Winscore = c(startIdx,sum(CHA2==MatA1),sum(CHA2==MatA2), sum(CHA2==PatA1),sum(CHA2==PatA2))
            haplo2=rbind(haplo2,A2Winscore)
            if(length(which(A2Winscore[2:5] == wSize)) > 0){
                colstat=as.vector(which(A2Winscore[2:5] == wSize))
                colstat=colstat+1
            }else { colstat=0 }
            A2WFilter=as.matrix(haplo2[idx+1,c(1,colstat)])

            ###if child's first 100% matches with one or both allele of mom or one or both alleles of dad or one allele of mom and one allele of dad

            #if((nrow(A1WFilter)>2) & (nrow(A2WFilter)>=2)){
            if((nrow(A1WFilter)>2)){

                if(nrow(A2WFilter)>=2){

                    ### first check if child's allele1 maches with one or both alleles of mmom and child's 2 allele macthes with one or both alleles of data

                    if((length(grep("MA", rownames(A1WFilter))) >= 1) & (length(grep("DA", rownames(A2WFilter))) >= 1)){
                          srd1 = "C1MA"
                          srd2 = "C2DA"
                    }
                    else if((length(grep("DA", rownames(A1WFilter))) >= 1) & (length(grep("MA", rownames(A2WFilter))) >= 1)){
                          ### if child's allele1 maches with one or both alleles from Dad and child's allele2 matched with one or more alleles of Mom
                          srd1 = "C1DA"
                          srd2 = "C2MA"
                    }else{
                          srd1="NA"
                          srd2="NA"
                    }
                ### does not match any of the allele of dad
                }else{
                    srd1="NA"
                    srd2="NA"
                }

          }else if(nrow(A1WFilter) == 2){
              ### if child's allele1 matches with only one allele from mom or dad
              srd1=ifelse(grepl("MA",rownames(A1WFilter)[2]),"C1MA",ifelse(grepl("DA",rownames(A1WFilter)[2]),"C1DA","NA"))
              #### check if child's allele two matches with one or two alles of parents
              if(nrow(A2WFilter)>2){
                  ### if child's allele1 is coming from mom then second must be match one one or both alleles of dad
                  if(grepl("MA", srd1)){
                        srd2 = ifelse(length(grep("DA", rownames(A2WFilter))) >= 1,"C2DA","NA")
                  }else if(grepl("DA", srd1)){
                        ### if child's allele1 coming from Dad then second must be match one one or both alleles of mom
                        srd2 = ifelse(length(grep("MA", rownames(A2WFilter))) >= 1,"C2MA","NA")
                  }
              }else if(nrow(A2WFilter)==2){
                  ### child's allele2 matches with only
                  srd2=ifelse(grepl("MA",rownames(A2WFilter)[2]),"C2MA",ifelse(grepl("DA",rownames(A2WFilter)[2]),"C2DA","NA"))
              }else {srd2="NA"}
          }else{
                #child's allele1 does not match any of the alleles of aprents
                srd1="NA"
                srd2="NA"
          }
          ### set up the child's allele order where 1st will be coming from mom and 2nd will be coming from dad
          if((srd1== "C1DA") & (srd2 == "C2MA")){
              ### if C1 is coming from dad  then put that allele second and allele from mom first
              op[r,2]=paste("DA",CHA1[1:step],sep=":")
              op[r,1]= paste("MA",CHA2[1:step],sep=":")
          }else if((srd1 == "C1MA") & (srd2 == "C2DA")){
              ### if C1 is coming from mom then put that allele first and allele from  dad second
              op[r,1]= paste("MA",CHA1[1:step],sep=":")
              op[r,2]=paste("DA",CHA2[1:step],sep=":")
          }else{
            ### otherwise set both alleles to missing
              op[r,]= "."
          }

          if(srd1 == "NA" | srd2 =="NA"){
              ###count number of missing sites
              ct = ct+length(r)
              ###possible switch
              errCt=errCt+1
          }

          startIdx= startIdx + step
          endIdx = startIdx+(wSize-1)
          idx=idx+1
        } #while
        op[,1]=ifelse(op[,1] != ".",substr(op[,1],4,4),".")
        op[,2]=ifelse(op[,2] != ".",substr(op[,2],4,4),".")

        cat(paste(triovcf[1,1],as.character(paste(colnames(trio)[10],colnames(trio)[11],colnames(trio)[12],sep=":")),ct,(ct/nrow(triovcf))*100,"%",errCt,sep=" "),"\n")

        result = paste(op[,1],op[,2],sep="|")
        result = cbind(result,trio[,11])
        result = cbind(result,trio[,12])
        colnames(result)=c(colnames(trio)[10],colnames(trio)[11],colnames(trio)[12])
        return(result)
  } #par for loop
  print("Finished Sliding window scan!")
  toc=Sys.time()
  print(toc-tic)


  finalVcf=cbind(pofo1,pofo)
  write.table(finalVcf,file=outVcfFile,row.names=F,col.names=T, quote=F, sep="\t")

  stopCluster(cl)
  mainToc=Sys.time()
  print ("Total processing time")
  print(mainToc-mainTic)
  rm(list=ls())
}##function end

