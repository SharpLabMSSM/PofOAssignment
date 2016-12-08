#' @title
#' Step 3 of 3 : Recover het sites in offsprings set to missing during sliding window (step2) in father-offspring Duos
#' @description
#' Compare offspring's genotypes in imputed QC data and phase corrected data at het sites in offspring where mothe ris homozygous and if is missing in phase corrected data, use Medelian inheritance rule.
#' Takes in phased vcf files of family data before sliding window and after sliding windo, and pedigree file.
#' @param vcf1 Phased genotype vcf file for family before sliding window step (Step2) but after Mendel error (Step1)
#' @param vcf2 Phased genotype vcf file (after phasing and QC and before sliding window step)
#' @param ped Tab delimited Pedigree file with four columns(FamilyID, IndividualID, FatherID, MotherID) with NO header.
#' @return NULL
#' @examples
#' DuoMotherOffspringRecoverMissing("/my/dir/MendelErrorFixed/FatherOffsprings/myphased.MendelErrorFixed.vcf","/my/dir/PhaseSwitchSlidingWindown/FatherOffsprings/myphased.MendelErrorFixed.PhaseSwitch.vcf","/my/dir/data/myped.txt")
#' @export
#' @note
#' Automatically creates dir named 'PhaseSwitchMissingFix/FatherOffsprings' one level above the data dir for the result file of this process
#' @author Bharati Jadhav
DuoFatherOffspringRecoverMissing = function(vcf1,vcf2, ped){
    mainTic = Sys.time()
    ###checking arguments
    if (missing(vcf1))
        stop("Need to specify the name of phased vcf with absolute path.")

    if (missing(vcf2))
        stop("Need to specify the name of phased vcf with absolute path.")

    if (missing(ped))
        stop("Need to specify the name of pedigree file with absolute path.")

    inDir2=dirname(vcf2)
    outDir=paste(dirname(dirname(inDir2)),"PhaseSwitchMissingFix","FatherOffsprings",sep="/")
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

    statFile=paste(outDir,paste0(coreName,".DuoPhaseSwitchFillMissStat.txt"),sep="/")
    outVcfFile=paste(outDir,paste0(coreName,".FillMiss.vcf"),sep="/")
    ###adding header to stat file
    ln=paste("chr","DuoId(Child|Mom)","hetCt_before","hetCt_after","%ofIncreaseInHetSite",sep=" ")
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
    #skipLines <- 12
    names <- scan(inVcfFile1,what = character(),skip=skipLines, nlines=1)
    sourceVcf = read.table(inVcfFile1)
    colnames(sourceVcf) = names


    ###check if ped file exists
    if(file.access(pedFile, mode = 0) != 0){
        stop(paste("File:", pedFile, "does not appear to exist!"))
    }

    ###check if file is readable
    if(file.access(pedFile, mode = 4) != 0){
        stop(paste("File:", pedFile, "appears to exist but is not readable!"))
    }
    ### if everything is fine with ped file read file
    print("Reading Pedigree file .....")

    duoped = read.table(pedFile)

    colnames(duoped) = c("FID","IID","PATID","MATID")
    duoped1= subset(duoped,MATID ==0 & PATID !=0 )

    #find the mother offspring column index for filtered imputed data
    iidxSource=match(duoped1$IID, colnames(sourceVcf))

    fidxSource=match(duoped1$PATID, colnames(sourceVcf))
    matchduoSource=na.omit(cbind(iidxSource,fidxSource))
    print(paste("No. of full trios in source vcf: ",nrow(matchduoSource)))

    sourceOrdered=matrix(0, nrow=nrow(sourceVcf), ncol = 9, byrow=T)
    sourceOrdered=sourceVcf[,1:9]
    print("Ordering source vcf in pedigree order [child,Father]  .....")
    for(i in 1:nrow(matchduoSource))
    {
        duo= sourceVcf[,c(1:9, matchduoSource[i,1],matchduoSource[i,2])]
        duo[,10:11]=apply(duo[,10:11],2, substr,1,3)

        sourceOrdered = cbind(sourceOrdered,duo[,10])
        colnames(sourceOrdered)[ncol(sourceOrdered)]=colnames(duo)[10]
        #add father's genotypes
        sourceOrdered = cbind(sourceOrdered,duo[,11])
        colnames(sourceOrdered)[ncol(sourceOrdered)]=colnames(duo)[11]
    }
    rm(sourceVcf)
    rm(duo)

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
    #skipLines <-0
    names2 <- scan(inVcfFile2,what = character(),skip=skipLines, nlines=1)
    targetVcf = read.table(inVcfFile2)
    colnames(targetVcf) = names2

    #find the mother offspring column index for phase corrected data
    iidxTarget=match(duoped1$IID, colnames(targetVcf))
    fidxTarget=match(duoped1$PATID, colnames(targetVcf))
    matchduoTarget=na.omit(cbind(iidxTarget,fidxTarget))

    #now compare each father-offspring genotypes in imputed, filtered data and phase corrected data(using sliding window algorithm).
    #if any snp is set to missing in phase corrected data check that site in imputed filter data and if child is heterozygous and
    #mother is homozygous fill the missing genotype with these genotypes
    fixedTarget=matrix(0, nrow=nrow(targetVcf), ncol = 9, byrow=T)
    fixedTarget=targetVcf[,1:9]
    for(i in 1:nrow(matchduoTarget))
    {
        #get the duo genotype information from the ordered source vcf
        #hetCt=0
	hetCtBefore=0
	hetCtAfter=0
        fixrows=0
        duoSO= sourceOrdered[,c(1:9, matchduoTarget[i,1],matchduoTarget[i,2])]
        duoTarget=targetVcf[,c(1:9, matchduoTarget[i,1],matchduoTarget[i,2])]
        #beforeCt=sum(duoTarget[,10] == ".|.")
        hetCtBefore=sum(duoTarget[,10] == "0|1" |duoTarget[,10] == "1|0" )
        fixrows=which((duoTarget[,10] == ".|.") & ((duoSO[,10] == "1|0" ) | (duoSO[,10] == "0|1" )) & ((duoSO[,11] == "1|1" ) | (duoSO[,11] == "0|0" )))
        #duoTarget[fixrows,10]= ifelse((duoSo[fixrows,11] == "0|0") & ((duoSo[fixrows,10] == "0|1")|(duoSo[fixrows,10] == "1|0")),"0|1",(ifelse((duoSo[fixrows,11] == "1|1") & ((duoSo[fixrows,10] == "0|1")|(duoSo[fixrows,10] ==
        duoTarget[fixrows,10]= ifelse(duoSO[fixrows,11] == "0|0","1|0",ifelse(duoSO[fixrows,11] == "1|1","0|1",".|."))
        hetCtAfter=sum(duoTarget[,10] == "0|1" |duoTarget[,10] == "1|0" )

        ln=paste(duoTarget[1,1],as.character(paste(colnames(duoTarget)[10],colnames(duoTarget)[11],sep=":")),hetCtBefore,hetCtAfter,round((((hetCtAfter-hetCtBefore)*100)/hetCtBefore),2),sep=" ")

        cat(paste(ln,"\n"),file=statFile, append=TRUE)

        fixedTarget=cbind(fixedTarget,duoTarget[,10])
        colnames(fixedTarget)[ncol(fixedTarget)]=colnames(duoTarget)[10]

        fixedTarget=cbind(fixedTarget,duoTarget[,11])
        colnames(fixedTarget)[ncol(fixedTarget)]=colnames(duoTarget)[11]

    }###for

    write.table(fixedTarget,file=outVcfFile,row.names=F,col.names=T, quote=F, sep="\t")

    mainToc=Sys.time()
    print(mainToc-mainTic)

    rm(list = ls())
}### end of function

