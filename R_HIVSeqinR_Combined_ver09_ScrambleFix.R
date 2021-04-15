###########################
#As long as this following section is user-verified to be ready, press CTRL-A, then CTRL-R to run all lines in R
###########################
rm(list=ls())
closeAllConnections() 
graphics.off() 
cat("\014") #clean console
MyWD <- getwd()
setwd(MyWD)
### Dependencies ###
#1. install library(Biostrings)
#2. install library(muscle)
#3. install blast+ suite, create HXB2 library installed in directory:
MyBlastnDir <- "/Users/guin/" #change as appropriate
#Follow blast+ suite installation instructions at
#https://www.ncbi.nlm.nih.gov/books/NBK279671/
#To create HXB2 library in blast+ suite, use commend:
# makeblastdb -in HXB2.fasta -parse_seqids -dbtype nucl 
# Make sure to rename R_HXB2.fasta to HXB2.fasta after copying file into dir where blast+ suite is installed
#4. Primer Settings:  Currently Optimized for these two 2nd round PCR primers. NGS library prep method must be shearing and adaptor-ligation based so that post de novo assembled contig is flanked by both primers
Primer2ndF_HXB2 <- "GCGCCCGAACAGGGACCTGAAAGCGAAAG"
Primer2ndR_HXB2 <- "TAAGCCTCAATAAAGCTTGCCTTGAGTGC"  #3UTRi "CTAGGGAACCCACTGCTTAAGCCT" or U5-547R "TAAGCCTCAATAAAGCTTGCCTTGAGTGC"
#5. Copy and paste all raw *.seq files in folder RAW_FASTA.  No dashes, *, or IUPAC mixture symbols are allowed (eg. N, R, S, Y etc)
##  Optional: Copy and paste all raw *.csv files in folder RAW_QC
#Final check:  Dir should have four files and at least one folder:
#(1) HIV Genome RefSeq ver04 AA.fasta
#(2) HIV Genome RefSeq ver04 NUC.fasta
#(3) R_HXB2.fasta
#(4) R_HIVSeqinR_Combined_ver04.R #This script
#(5) ./RAW_FASTA/ #folder containing fasta file(s) of post de novo assembled contigs
#(6) ./RAW_QC/ #optional, only if using the Massachusetts General Hospital CCDB Core service


#The following sections are automated and do not require user edits unless users wish to adjust specific settings

###########################
#R01_FASTApaste_fasta2csv_ver05.R
###########################

#Objective:  (1) Paste all FASTA files in one folder into a single file (2) Output as FASTA and csv
#Folder should have nothing beside this script and individual fasta files
#Folder must be named RAW_FASTA

#Step 00. Setup R environment
StartSysTime <- Sys.time()
rm(list= ls()[!(ls() %in% c("MyBlastnDir","StartSysTime","Primer2ndF_HXB2","Primer2ndR_HXB2"))])
closeAllConnections() 
graphics.off() 
cat("\014") #send CTRL+L to R, clean console
MyWD <- getwd()
setwd(MyWD)
library(Biostrings)

#Step 01.  User Settings: List all files
AllFiles <- list.files("./RAW_FASTA")
Remove <- c("Output","FASTApaste_fasta2csv_ver04.R","Output_Concatenate.fasta")
AllFiles <- AllFiles[!AllFiles %in% Remove]
Output_FASTA <- "Output_Concatenate.fasta"
MyAllFASTA <- DNAStringSet()


#Step 02.  Loop through each file in the dir
for (n in 1:length(AllFiles)){  
  #Step01. Read fasta
  InputFASTAFileName <- AllFiles[n]
  MyFASTA <- readDNAStringSet(paste("./RAW_FASTA/",InputFASTAFileName,sep=""))
  #Step02.  Append to DNAStringSet object
  MyAllFASTA <- append(MyAllFASTA,MyFASTA)
}

writeXStringSet(MyAllFASTA,Output_FASTA)


#Step 03. Convert FASTA to csv

input <- readLines("Output_Concatenate.fasta")
output <- file("Output_Concatenate.csv","w")

currentSeq <- 0
newLine <- 0

for(i in 1:length(input)) {
  if(strtrim(input[i], 1) == ">") {
    if(currentSeq == 0) {
      writeLines(gsub(" ","",gsub(">","",paste(input[i],"\t"))), output, sep="")
      currentSeq <- currentSeq + 1
    } else {
      writeLines(gsub(" ","",gsub(">","",paste("\n",input[i],"\t", sep=""))), output, sep="")
    }
  } else {
    writeLines(paste(input[i]), output, sep="")
  }
}

close(output)


###############################
#R02_RC_ver03.R
###############################

#Objective:  Reverse Complement of a FASTA file

#Step01.  User settings
MyInputFASTA <- "Output_Concatenate.fasta"

#Step 02.  Reverse Complement
MyFASTA_DNAStringSet <- readDNAStringSet(MyInputFASTA)
Output_fasta <- paste(gsub(".fasta","",ignore.case=TRUE,MyInputFASTA),"_RC.fasta",sep="")
MyFASTA_DNAStringSet_RC <- reverseComplement(MyFASTA_DNAStringSet)
writeXStringSet(MyFASTA_DNAStringSet_RC,Output_fasta)

#Step 03.  Convert FASTA to CSV
MyOutputDF <- data.frame()
MySEQID <- c()
MySEQRAW <- c()
for (k in 1:length(MyFASTA_DNAStringSet_RC)){
  MySEQID <- append(MySEQID,names(MyFASTA_DNAStringSet_RC[k]))
  MySEQRAW <- append(MySEQRAW,as.character(MyFASTA_DNAStringSet_RC[[k]]))
}
MyOutputDF <- data.frame(cbind(MySEQID,MySEQRAW))
colnames(MyOutputDF) <- c("SEQID","SEQRAWRC")
write.csv(MyOutputDF,paste(gsub(".fasta","",Output_fasta),".csv",sep=""),row.names=FALSE)


###############################

#Run terminal from R studio #Miracle Fix

#MyWD <- getwd() 
#setwd(MyWD)
#system(paste("cd ",MyWD,sep=""))
MyBlastCommand <- paste(
  "blastn -query Output_Concatenate.fasta -db "
  ,MyBlastnDir
  ,"HXB2.fasta -num_alignments 1 -reward 1 -penalty -1 -gapopen 2 -gapextend 1 -out Output_Blastn_HXB2MEGA28_tabdelim.txt -outfmt \"6 qseqid qlen sseqid sgi slen qstart qend sstart send evalue bitscore length pident nident btop stitle sstrand\""
  ,sep=""
)
system(MyBlastCommand)


###############################
#R03b_MyPipelineQC_ver05_QC123update.R
###############################


MyQCPipeline <- function () {
  #Step 00.  Prep R
  rm(list= ls()[!(ls() %in% c("MyBlastnDir","StartSysTime","Primer2ndF_HXB2","Primer2ndR_HXB2"))])
  closeAllConnections() 
  graphics.off() 
  cat("\014") #send CTRL+L to R, clean console
  MyWD <- getwd()
  setwd(MyWD)
  
  #Step 01.  User settings  ##CHANGE_ME###
  Current_QCpipeline_version <- "ver05"
  #QC_Rscript_filename <- "MyPipelineQC_ver03_subfolder.R"
  My_Min_Acceptable_Coverage <- 10 #<1000 (subsetting >=1000) #20171207update to min
  My_QC01_PercMinAcceptCovergeInLen <- 0.8
  My_QC02_MinCoverageThreshold <- 3 #20171207update to min
  My_QCPercMix_Thresholds <- c(0.2,0.1,0.05,0.02,0.01) #e.g. if there is >=20%, or 10%..etc any mixture per position, it is poor quality
  My_QCPercBias_Cutoff <- 0.8 #ie if >=0.8 mismatch base is biased, then it's poor quality
  My_QC03_MaxBiasAbsoluteAllowed <- 0.1 #20171207 Changed to %
  Output_df_filename <- paste("Output_QC_",Current_QCpipeline_version,sep="")
  
  
  #Step 02.  Get list of file names in folder
  FileList_fin <- list.files("./RAW_QC")
  #FileList_fin <- subset(FileList_raw,FileList_raw!=QC_Rscript_filename)
  #Step 02i. Optional.  Manual: Ensure File names are formatted correctly
  #write(FileList_fin,"Output_FileListManualCheck.txt")
  
  
  
  #Step 03.  Final summary output setup, one element per sample
  Output01a_CurrentFile <- c()  #"88fF4_HIV_275761w_BE8_coverage_of_total_216298_reads.csv"
  Output01b_CurrentSample <- c() #88fF4
  Output02a_CurrentContig <- c() #216298
  Output02a1_CurrentContigSimple <- c()
  Output02b_TotalSmallReads <- c() #p1
  Output04a_PositionsLessThenMinCoverageAccepted <- c() #115
  Output04b_TotalContigLength <- c() #669
  Output04c_PercPositionGoodCoverage <- c()
  Output04d_QC01_Depth <- c()
  Output04e_MinCoverage <- c()
  Output04f_MaxCoverage <- c()
  Output04g_MedianCoverage <- c()
  Output04h_QC02_MinCoverage <- c()
  Output05a_MaxFreqMismatch <- c() #0.5247496
  Output05b_PositionsMoreThenMinCoverageAccepted <- c()
  Output06a_Threshold20_NonBias <- c() #2
  Output07a_Threshold10_NonBias <- c()
  Output08a_Threshold05_NonBias <- c()
  Output09a_Threshold02_NonBias <- c()
  Output10a_Threshold01_NonBias <- c()
  Output06b_Threshold20_Bias <- c() #2
  Output07b_Threshold10_Bias <- c()
  Output08b_Threshold05_Bias <- c()
  Output09b_Threshold02_Bias <- c()
  Output10b_Threshold01_Bias <- c()
  Output11_QC03_MixturesSGA <- c()
  
  
  
  #Step 04.  Loop through each file
  for (i in 1:length(FileList_fin)) {
    #Step 05a.  Read raw QC file & Replace all NAs as zeros & Replace all white spaces
    CurrentSample <- FileList_fin[i]
    CurrentDF <- read.delim(paste("./RAW_QC/",FileList_fin[i],sep="")) #colnames(CurrentDF)
    CurrentDF[is.na(CurrentDF)] <- 0
    CurrentDF$Contig <- gsub(" ","",CurrentDF$Contig)
    #Step 05b.  Add Read1+Read2 columns for each nucleotide A (col6+10) >> 14, T (col7+11) >> 15, C (col8+12), G (col9+13)
    CurrentDF[,14] <- CurrentDF[,6]+CurrentDF[,10]  #Read1 + Read2:  A+A
    CurrentDF[,15] <- CurrentDF[,7]+CurrentDF[,11]  #Read1 + Read2:  T+T
    CurrentDF[,16] <- CurrentDF[,8]+CurrentDF[,12]  #Read1 + Read2:  C+C
    CurrentDF[,17] <- CurrentDF[,9]+CurrentDF[,13]  #Read1 + Read2:  G+G
    #Step 06.  Loop through each contig:  Output 1-10
    ListAllContigs <- unique(CurrentDF$Contig)
    ########  If there were zero reads ######
    if (length(ListAllContigs)==0){
      Output01a_CurrentFile <- append(Output01a_CurrentFile,CurrentSample)
      Output01b_CurrentSample <- append(Output01b_CurrentSample,unlist(strsplit(CurrentSample,"_"))[1])
      Output02a_CurrentContig <- append(Output02a_CurrentContig,"ZeroReads")
      Output02a1_CurrentContigSimple <- append(Output02a1_CurrentContigSimple,"ZeroReads")
      Output02b_TotalSmallReads <- append(Output02b_TotalSmallReads,0)
      Output04a_PositionsLessThenMinCoverageAccepted <- append(Output04a_PositionsLessThenMinCoverageAccepted,NA)
      Output04b_TotalContigLength <- append(Output04b_TotalContigLength,NA)
      Output04c_PercPositionGoodCoverage <- append(Output04c_PercPositionGoodCoverage,NA)
      Output04d_QC01_Depth <- append(Output04d_QC01_Depth,"QC01_Depth_Failed")
      Output04e_MinCoverage <- append(Output04e_MinCoverage,NA)
      Output04f_MaxCoverage <- append(Output04f_MaxCoverage,NA)
      Output04g_MedianCoverage <- append(Output04g_MedianCoverage,NA)
      Output04h_QC02_MinCoverage <- append(Output04h_QC02_MinCoverage,"QC02_MinCoverage_Failed")
      Output05a_MaxFreqMismatch <- append(Output05a_MaxFreqMismatch,NA)
      Output05b_PositionsMoreThenMinCoverageAccepted <- append(Output05b_PositionsMoreThenMinCoverageAccepted,NA)
      Output06a_Threshold20_NonBias <- append(Output06a_Threshold20_NonBias,NA)
      Output07a_Threshold10_NonBias <- append(Output07a_Threshold10_NonBias,NA)
      Output08a_Threshold05_NonBias <- append(Output08a_Threshold05_NonBias,NA)
      Output09a_Threshold02_NonBias <- append(Output09a_Threshold02_NonBias,NA)
      Output10a_Threshold01_NonBias <- append(Output10a_Threshold01_NonBias,NA)
      Output06b_Threshold20_Bias <- append(Output06b_Threshold20_Bias,NA)
      Output07b_Threshold10_Bias <- append(Output07b_Threshold10_Bias,NA)
      Output08b_Threshold05_Bias <- append(Output08b_Threshold05_Bias,NA)
      Output09b_Threshold02_Bias <- append(Output09b_Threshold02_Bias,NA)
      Output10b_Threshold01_Bias <- append(Output10b_Threshold01_Bias,NA)
    } else {
      for (k in 1:length(unique(CurrentDF$Contig))){
        #Step 07.  Export File Name (1a), Sample Name (1b), and TotalSmallReads on Filename (2)
        Output01a_CurrentFile <- append(Output01a_CurrentFile,CurrentSample)
        Output01b_CurrentSample <- append(Output01b_CurrentSample,unlist(strsplit(CurrentSample,"_"))[1])
        Output02b_TotalSmallReads <- append(Output02b_TotalSmallReads,as.numeric(unlist(strsplit(CurrentSample,"_"))[8]))
        #Step 08.  Calculate (number of base positions with less than accaptable coverage)/(total bp length of sequence)
        ############# PASS if QC01_Depth>=1000 more than 0.8 of its base positions ###############
        CurrentDF_Contig <- subset(CurrentDF,CurrentDF$Contig==ListAllContigs[k])
        Output02a_CurrentContig <- append(Output02a_CurrentContig,CurrentDF_Contig$Contig[1])
        Output02a1_CurrentContigSimple <- append(Output02a1_CurrentContigSimple,unlist(strsplit(CurrentDF_Contig$Contig[1],"_"))[4])
        CurrentDF_CoverageMinExcl <- subset(CurrentDF_Contig,CurrentDF_Contig[,4]>=My_Min_Acceptable_Coverage)  #Col4: Number of perfect matches
        Output04a_PositionsLessThenMinCoverageAccepted <- append(Output04a_PositionsLessThenMinCoverageAccepted,nrow(CurrentDF_Contig)-nrow(CurrentDF_CoverageMinExcl))
        Output04b_TotalContigLength <- append(Output04b_TotalContigLength,unlist(strsplit(CurrentDF_Contig$Contig[1],"_"))[3])
        PercPositionGoodCoverage <- nrow(CurrentDF_CoverageMinExcl)/nrow(CurrentDF_Contig)
        Output04c_PercPositionGoodCoverage <- append(Output04c_PercPositionGoodCoverage,PercPositionGoodCoverage)
        if(PercPositionGoodCoverage>=My_QC01_PercMinAcceptCovergeInLen){
          Output04d_QC01_Depth <- append(Output04d_QC01_Depth,"QC01_Depth_Passed")
        } else {
          Output04d_QC01_Depth <- append(Output04d_QC01_Depth,"QC01_Depth_Failed")
        }
        ############# Update 20170203 #################################################
        ############# PASS if QC02_3rd-to-MinCoverage>=10 ##############################################
        Output04e_MinCoverage <- append(Output04e_MinCoverage,min(CurrentDF_Contig[,4]+CurrentDF_Contig[,5]))
        Output04f_MaxCoverage <- append(Output04f_MaxCoverage,max(CurrentDF_Contig[,4]+CurrentDF_Contig[,5]))
        Output04g_MedianCoverage <- append(Output04g_MedianCoverage,median(CurrentDF_Contig[,4]+CurrentDF_Contig[,5]))
        #20170203. If any position has <10 reads coverage, fail.  First and last position occasional fail. 
        #Solution:  sort the coverage vector "CurrentDF_Contig4+5", and the third number has to be >=10
        if((sort(CurrentDF_Contig[,4]+CurrentDF_Contig[,5])[3])>=My_QC02_MinCoverageThreshold){   
          Output04h_QC02_MinCoverage <- append(Output04h_QC02_MinCoverage,"QC02_MinCoverage_Passed")    
        } else {
          Output04h_QC02_MinCoverage <- append(Output04h_QC02_MinCoverage,"QC02_MinCoverage_Failed")    
        }
        #####################################CurrentDF_CoverageMinExcl##############################################################
        #Step 09. ONLY in positions passing QC01_depth >=1000:  Calculate ###Perc### Mismatch & Replace all NAs as zeros for EACH base positions
        ######## The reason I chose ONLY positions passing QC01_depth was because... it makes no sense to look at %mismatch when coverage was 10#######
        #######  This part of the pipeline is used to control +1 template (ie. not-real-SGA)... therefore I want to get rid of sequencing quality bias ######
        if (nrow(CurrentDF_CoverageMinExcl)==0){
          Output05a_MaxFreqMismatch <- append(Output05a_MaxFreqMismatch,NA)
          Output05b_PositionsMoreThenMinCoverageAccepted <- append(Output05b_PositionsMoreThenMinCoverageAccepted,NA)
          Output06a_Threshold20_NonBias <- append(Output06a_Threshold20_NonBias,NA)
          Output07a_Threshold10_NonBias <- append(Output07a_Threshold10_NonBias,NA)
          Output08a_Threshold05_NonBias <- append(Output08a_Threshold05_NonBias,NA)
          Output09a_Threshold02_NonBias <- append(Output09a_Threshold02_NonBias,NA)
          Output10a_Threshold01_NonBias <- append(Output10a_Threshold01_NonBias,NA)
          Output06b_Threshold20_Bias <- append(Output06b_Threshold20_Bias,NA)
          Output07b_Threshold10_Bias <- append(Output07b_Threshold10_Bias,NA)
          Output08b_Threshold05_Bias <- append(Output08b_Threshold05_Bias,NA)
          Output09b_Threshold02_Bias <- append(Output09b_Threshold02_Bias,NA)
          Output10b_Threshold01_Bias <- append(Output10b_Threshold01_Bias,NA)
        } else {
          Mismatch <- CurrentDF_CoverageMinExcl[,5]/(CurrentDF_CoverageMinExcl[,4]+CurrentDF_CoverageMinExcl[,5]) ## %Mismatch as a vector##
          Mismatch[is.na(Mismatch)] <- 0
          Output05a_MaxFreqMismatch <- append(Output05a_MaxFreqMismatch,max(Mismatch))
          Output05b_PositionsMoreThenMinCoverageAccepted <- append(Output05b_PositionsMoreThenMinCoverageAccepted,nrow(CurrentDF_CoverageMinExcl))
          #Step 10a.  Count all positions >=X% mismatch; regardless of bias/non-bias nucleotide ##User_Defined##
          Mismatch_Threshold_TF <- Mismatch >= My_QCPercMix_Thresholds[1] 
          Output06a_Threshold20_NonBias <- append(Output06a_Threshold20_NonBias,sum(Mismatch_Threshold_TF == TRUE))
          Mismatch_Threshold_TF <- Mismatch >= My_QCPercMix_Thresholds[2] 
          Output07a_Threshold10_NonBias <- append(Output07a_Threshold10_NonBias,sum(Mismatch_Threshold_TF == TRUE))
          Mismatch_Threshold_TF <- Mismatch >= My_QCPercMix_Thresholds[3] 
          Output08a_Threshold05_NonBias <- append(Output08a_Threshold05_NonBias,sum(Mismatch_Threshold_TF == TRUE))
          Mismatch_Threshold_TF <- Mismatch >= My_QCPercMix_Thresholds[4] 
          Output09a_Threshold02_NonBias <- append(Output09a_Threshold02_NonBias,sum(Mismatch_Threshold_TF == TRUE))
          Mismatch_Threshold_TF <- Mismatch >= My_QCPercMix_Thresholds[5] 
          Output10a_Threshold01_NonBias <- append(Output10a_Threshold01_NonBias,sum(Mismatch_Threshold_TF == TRUE))
          #Step10b.  Count all positions >=X% (eg. 5%) mismatch AND bias >=X% (eg. 80%)
          Mismatch_Threshold_Bias_TF <- Mismatch >= My_QCPercMix_Thresholds[1] & 
            (CurrentDF_CoverageMinExcl[,14]/CurrentDF_CoverageMinExcl[,5] >= My_QCPercBias_Cutoff | 
               CurrentDF_CoverageMinExcl[,15]/CurrentDF_CoverageMinExcl[,5] >= My_QCPercBias_Cutoff |
               CurrentDF_CoverageMinExcl[,16]/CurrentDF_CoverageMinExcl[,5] >= My_QCPercBias_Cutoff |
               CurrentDF_CoverageMinExcl[,17]/CurrentDF_CoverageMinExcl[,5] >= My_QCPercBias_Cutoff)
          Output06b_Threshold20_Bias <- append(Output06b_Threshold20_Bias,sum(Mismatch_Threshold_Bias_TF == TRUE))
          Mismatch_Threshold_Bias_TF <- Mismatch >= My_QCPercMix_Thresholds[2] & 
            (CurrentDF_CoverageMinExcl[,14]/CurrentDF_CoverageMinExcl[,5] >= My_QCPercBias_Cutoff | 
               CurrentDF_CoverageMinExcl[,15]/CurrentDF_CoverageMinExcl[,5] >= My_QCPercBias_Cutoff |
               CurrentDF_CoverageMinExcl[,16]/CurrentDF_CoverageMinExcl[,5] >= My_QCPercBias_Cutoff |
               CurrentDF_CoverageMinExcl[,17]/CurrentDF_CoverageMinExcl[,5] >= My_QCPercBias_Cutoff)
          Output07b_Threshold10_Bias <- append(Output07b_Threshold10_Bias,sum(Mismatch_Threshold_Bias_TF == TRUE))
          Mismatch_Threshold_Bias_TF <- Mismatch >= My_QCPercMix_Thresholds[3] & 
            (CurrentDF_CoverageMinExcl[,14]/CurrentDF_CoverageMinExcl[,5] >= My_QCPercBias_Cutoff | 
               CurrentDF_CoverageMinExcl[,15]/CurrentDF_CoverageMinExcl[,5] >= My_QCPercBias_Cutoff |
               CurrentDF_CoverageMinExcl[,16]/CurrentDF_CoverageMinExcl[,5] >= My_QCPercBias_Cutoff |
               CurrentDF_CoverageMinExcl[,17]/CurrentDF_CoverageMinExcl[,5] >= My_QCPercBias_Cutoff)
          Output08b_Threshold05_Bias <- append(Output08b_Threshold05_Bias,sum(Mismatch_Threshold_Bias_TF == TRUE))
          Mismatch_Threshold_Bias_TF <- Mismatch >= My_QCPercMix_Thresholds[4] & 
            (CurrentDF_CoverageMinExcl[,14]/CurrentDF_CoverageMinExcl[,5] >= My_QCPercBias_Cutoff | 
               CurrentDF_CoverageMinExcl[,15]/CurrentDF_CoverageMinExcl[,5] >= My_QCPercBias_Cutoff |
               CurrentDF_CoverageMinExcl[,16]/CurrentDF_CoverageMinExcl[,5] >= My_QCPercBias_Cutoff |
               CurrentDF_CoverageMinExcl[,17]/CurrentDF_CoverageMinExcl[,5] >= My_QCPercBias_Cutoff)
          Output09b_Threshold02_Bias <- append(Output09b_Threshold02_Bias,sum(Mismatch_Threshold_Bias_TF == TRUE))
          Mismatch_Threshold_Bias_TF <- Mismatch >= My_QCPercMix_Thresholds[5] & 
            (CurrentDF_CoverageMinExcl[,14]/CurrentDF_CoverageMinExcl[,5] >= My_QCPercBias_Cutoff | 
               CurrentDF_CoverageMinExcl[,15]/CurrentDF_CoverageMinExcl[,5] >= My_QCPercBias_Cutoff |
               CurrentDF_CoverageMinExcl[,16]/CurrentDF_CoverageMinExcl[,5] >= My_QCPercBias_Cutoff |
               CurrentDF_CoverageMinExcl[,17]/CurrentDF_CoverageMinExcl[,5] >= My_QCPercBias_Cutoff)
          Output10b_Threshold01_Bias <- append(Output10b_Threshold01_Bias,sum(Mismatch_Threshold_Bias_TF == TRUE))
        }
      }
    }
  } 
  
  ########### QC03_MixturesSGA ########################
  
  #Criteria:  If there were more than 10 positions (absolute count) that has >=5% Mixtures AND which were biased, fail it.
  #20171207: Change to %: If there were over or equal to 10% of contig length that has >=5% Biased base, fail it.
  Output11_QC03_MixturesSGA <- (Output08b_Threshold05_Bias/as.numeric(Output04b_TotalContigLength))>=My_QC03_MaxBiasAbsoluteAllowed
  Output11_QC03_MixturesSGA[is.na(Output11_QC03_MixturesSGA)] <- "QC03_MixturesSGA_Failed"
  Output11_QC03_MixturesSGA <- gsub("TRUE","QC03_MixturesSGA_Failed",gsub("FALSE","QC03_MixturesSGA_Passed",Output11_QC03_MixturesSGA))
  
  ########### QC 1 + 2 + 3 ############################
  
  Output12_QC123Comb <- 
    substr(Output04d_QC01_Depth,nchar(Output04d_QC01_Depth)-5,nchar(Output04d_QC01_Depth)-2) == "Pass" &
    substr(Output04h_QC02_MinCoverage,nchar(Output04h_QC02_MinCoverage)-5,nchar(Output04h_QC02_MinCoverage)-2) == "Pass" &
    substr(Output11_QC03_MixturesSGA,nchar(Output11_QC03_MixturesSGA)-5,nchar(Output11_QC03_MixturesSGA)-2) == "Pass"
  Output12_QC123Comb <- gsub("TRUE","QC04_123Comb_Passed",gsub("FALSE","QC04_123Comb_Failed",Output12_QC123Comb))
  #summary(as.factor(Output12_QC123Comb)) #tally - #passed vs #failed
  
  
  
  #####################################################
  
  mydf <- data.frame(
    Output01a_CurrentFile,
    Output01b_CurrentSample,
    Output02a_CurrentContig,
    Output02a1_CurrentContigSimple,
    Output02b_TotalSmallReads,
    ###### QC01 Depth ######
    Output04a_PositionsLessThenMinCoverageAccepted,
    Output04b_TotalContigLength,
    Output04c_PercPositionGoodCoverage,
    Output04d_QC01_Depth,
    ##### QC02 MinCoverage ######
    Output04e_MinCoverage,
    Output04f_MaxCoverage,
    Output04g_MedianCoverage,
    Output04h_QC02_MinCoverage,
    ##### QC03 MismatchSGA ######
    Output05a_MaxFreqMismatch,
    Output05b_PositionsMoreThenMinCoverageAccepted,
    Output06a_Threshold20_NonBias,
    Output06b_Threshold20_Bias,
    Output07a_Threshold10_NonBias,
    Output07b_Threshold10_Bias, 
    Output08a_Threshold05_NonBias,
    Output08b_Threshold05_Bias, 
    Output09a_Threshold02_NonBias, 
    Output09b_Threshold02_Bias, 
    Output10a_Threshold01_NonBias, 
    Output10b_Threshold01_Bias,
    Output11_QC03_MixturesSGA,
    Output12_QC123Comb
  )
  
  MyNrows <- length(Output11_QC03_MixturesSGA)
  write.csv(mydf,file=paste(Output_df_filename,".csv",sep=""),row.names=FALSE)
  

}

#Check to see the RAW_QC exists
if (any(list.files()=="RAW_QC") == FALSE) {
  #Skip this# 
} else {
  #run the QCPipeline
  MyQCPipeline()
}


#############################################
#R04_BigSummary_ver45_PSC5DEFECT.R Part I
#############################################

#Objective: General summary table for each MGH DNA core output contig


#Step 00. Setup R #PC@work R version 3.2.2
rm(list= ls()[!(ls() %in% c("MyBlastnDir","StartSysTime","Primer2ndF_HXB2","Primer2ndR_HXB2"))])
closeAllConnections() 
graphics.off() 
cat("\014") #send CTRL+L to R, clean console
Sys.setenv(TZ="EST")
MyWD <- getwd()
setwd(MyWD)
library(Biostrings) #PC@work Biostrings_2.38.1
library(muscle) #PC@work muscle_3.12.0
#library(ape)
#sessionInfo()

#Step 01. User Settings
Input00_RAWFASTA <- "Output_Concatenate.csv"
Input02_BlastHXB2 <- "Output_Blastn_HXB2MEGA28_tabdelim.txt"
Input03_RC <- "Output_Concatenate_RC.csv"
MyBigSumVer <- "ver45"

#Step 02. Read ALL Input Files as DF
DF00_RAWFASTA <- read.delim(Input00_RAWFASTA,header=FALSE,col.names=c("SEQID","SEQRAW"))
#DF01_BlastNT <- read.delim(Input01_BlastNT,header=FALSE)
#colnames(DF01_BlastNT) <- c("SEQID", "qlen", "sseqid", "sgi", "slen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "length", "pident", "nident", "btop", "stitle", "sstrand")
DF02_BlastHXB2 <- read.delim(Input02_BlastHXB2,header=FALSE)
colnames(DF02_BlastHXB2) <- c("SEQID", "qlen", "sseqid", "sgi", "slen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "length", "pident", "nident", "btop", "stitle", "sstrand")
DF03_RC <- read.csv(Input03_RC,header=TRUE)

#Step 03. Joint RAW and RC
DF_FINAL <- merge(DF00_RAWFASTA,DF03_RC,by="SEQID",all=TRUE)

####################################################

###Make unique identifier per contig
MyUniqSampleContigID <- c()
SEQID <- as.character(DF_FINAL$SEQID)
for (n in 1:length(SEQID)){
  MyUniqSampleID <- unlist(strsplit(SEQID[n],"_"))[1]
  ### If SEQID can be splitted into 5 parts, look for the part that begins with p; if none, returns "check"
  if (length(unlist(strsplit(SEQID[n],"_")))==5) {
    if (substring(unlist(strsplit(SEQID[n],"_"))[5],1,1)=="p"){
      MyUniqSampleContig <- unlist(strsplit(SEQID[n],"_"))[5]
    } else {
      MyUniqSampleContig <- "check"
    }
    ### If SEQID can be splitted into 6 parts, look for the part that begins with p; if none, returns the 6th element
  } else if (length(unlist(strsplit(SEQID[n],"_")))==6){
    if (substring(unlist(strsplit(SEQID[n],"_"))[5],1,1)=="p"){
      MyUniqSampleContig <- unlist(strsplit(SEQID[n],"_"))[5]
    } else if (substring(unlist(strsplit(SEQID[n],"_"))[6],1,1)=="p"){
      MyUniqSampleContig <- unlist(strsplit(SEQID[n],"_"))[6]
    } else {
      MyUniqSampleContig <- unlist(strsplit(SEQID[n],"_"))[6]
    }
    ### If SEQID can be splitted into 7 parts, look for the part that begins with p; if none, returns "check"
  } else if (length(unlist(strsplit(SEQID[n],"_")))==7){
    if (substring(unlist(strsplit(SEQID[n],"_"))[5],1,1)=="p"){
      MyUniqSampleContig <- unlist(strsplit(SEQID[n],"_"))[5]
    } else if (substring(unlist(strsplit(SEQID[n],"_"))[6],1,1)=="p"){
      MyUniqSampleContig <- unlist(strsplit(SEQID[n],"_"))[6]
    } else if (substring(unlist(strsplit(SEQID[n],"_"))[7],1,1)=="p") {
      MyUniqSampleContig <- unlist(strsplit(SEQID[n],"_"))[7]
    } else {
      MyUniqSampleContig <- unlist(strsplit(SEQID[n],"_"))[6]
    }
  } else {
    MyUniqSampleContig <- "check"
  }
  SampleContigID <- paste(MyUniqSampleID,"_",MyUniqSampleContig,sep="")
  MyUniqSampleContigID <- append(MyUniqSampleContigID,SampleContigID)
}
DF_FINAL <- cbind(DF_FINAL,MyUniqSampleContigID)

#View(DF_FINAL)

#####################################################################################

#Analysis01.  Human or HIV?  

#01b. Produce T/F vector 
MyHIVTF <- c()
#01c. Retrieve unique SEQID from DF00_RAWFASTA as character vector
MyUniqSEQID <- c()
for(i in 1:nrow(DF_FINAL)){
  MyUniqSEQID <- append(MyUniqSEQID,as.character(DF_FINAL$SEQID[i]))
}
#01d. For each uniq SEQID in DF_FINAL, return True if there is ###at least one HIV in this particular contig####
MyHIVTF_FINAL <- logical()
MyAlignedLength <- c()
for (j in 1:length(MyUniqSEQID)){
  MyDF_Pre <- subset(DF02_BlastHXB2,DF02_BlastHXB2$SEQID == MyUniqSEQID[j])
  #Remove any mappings that falls within 1-638 or 9632-9719
  MyDF <- subset(MyDF_Pre,(MyDF_Pre$sstart>=638 & MyDF_Pre$sstart<=9632)|(MyDF_Pre$send>=638 & MyDF_Pre$send<=9632)|(MyDF_Pre$sstart<=638 & MyDF_Pre$send>=9632)|(MyDF_Pre$send<=638 & MyDF_Pre$sstart>=9632))
  #If BlastHXB2.MEGA has an entry
  if(nrow(MyDF)>=1){
    #Remove any overlapped mappings, take the longer one
    MyOverlapTF <- logical()
    for (b in 1:nrow(MyDF)){
      MyStartCoordinate <- MyDF$qstart[b]
      MyEndCoordinate <- MyDF$qend[b]
      #Is there another row that maps to a bigger range than me (this row)?
      MyTestDF <- subset(MyDF,MyDF$qstart<=MyStartCoordinate&MyDF$qend>=MyEndCoordinate)
      if(nrow(MyTestDF)==1){
        MyOverlapTF <- append(MyOverlapTF,FALSE) #I have no overlap
      } else if (nrow(MyTestDF)>1){
        MyOverlapTF <- append(MyOverlapTF,TRUE) #I have overlap
      } else {
        MyOverlapTF <- append(MyOverlapTF,NA)
      }
    }
    MyDF <- cbind(MyDF,MyOverlapTF)
    MyDF <- subset(MyDF,MyDF$MyOverlapTF==FALSE)
    #View(MyDF_Pre)
    #View(MyDF_Pre2)
    #View(MyDF)
    # Check to see the length of contig aligned 
    LengthAligned <- 0
    for (a in 1:nrow(MyDF)){
      LengthAligned <- LengthAligned + abs(MyDF$qstart[a]-MyDF$qend[a])+1
    }
    if (LengthAligned/MyDF$qlen[1]>=0.8){
      MyHIVTF_FINAL <- append(MyHIVTF_FINAL,TRUE) #return TRUE if there is any True in a logical vector
      MyAlignedLength <- append(MyAlignedLength,LengthAligned)
    } else {
      MyHIVTF_FINAL <- append(MyHIVTF_FINAL,FALSE)
      MyAlignedLength <- append(MyAlignedLength,LengthAligned)
    }
    #If BlastHXB2.MEGA has no entry
  } else {
    MyHIVTF_FINAL <- append(MyHIVTF_FINAL,FALSE)
    MyAlignedLength <- append(MyAlignedLength,0)
  }
}
DF_FINAL <- cbind(DF_FINAL,MyHIVTF_FINAL)
#View(DF_FINAL)

####################################################################

#Analysis02. Plus or minus? NCBI HXB2


### Important:  Must subset DF02_BlastHXB2 to sstart-send NOT-FALL-BETWEEN 1-638 ###

#02a. Produce MyStrand+- vector of whether the Hxb2 match is plus ("Plus") or minus ("True") or plus&minus ("Mix")
MyStrand <- DF02_BlastHXB2[,"sstrand"]
#02b. For each uniq SEQID in DF_FINAL, evaluate whether it is pure-plus, pure-minus, or plus&minus-mix based on NCBI-HXB2 blast+
MyStrandBlastHXB2_FINAL <- c()
MyHXB2MapStart <- c()
MyHXB2MapEnd <- c()
for (j in 1:length(MyUniqSEQID)){
  MyDF <- subset(DF02_BlastHXB2,DF02_BlastHXB2$SEQID == MyUniqSEQID[j])
  MyDF_638 <- subset(MyDF, MyDF$sstart > 622 & MyDF$send > 622) #637 left out 1stPCR 623
  MyDF_638 <- MyDF_638[order(MyDF_638$qstart,MyDF_638$qend),]
  MyStrandVector <- as.character(MyDF_638$sstrand)
  if (length(MyStrandVector)==1){
    MyStrandBlastHXB2_FINAL <- append(MyStrandBlastHXB2_FINAL,MyStrandVector[1]) 
    MyHXB2MapStart <- append(MyHXB2MapStart,as.character(MyDF_638$sstart)[1])
    MyHXB2MapEnd <- append(MyHXB2MapEnd,as.character(MyDF_638$send)[1])
  } else if (length(MyStrandVector)==0){
    MyStrandBlastHXB2_FINAL <- append(MyStrandBlastHXB2_FINAL,"NoMapping")
    MyHXB2MapStart <- append(MyHXB2MapStart,"NoMapping")
    MyHXB2MapEnd <- append(MyHXB2MapEnd,"NoMapping")
  } else if (length(MyStrandVector)>1){
    if (length(unique(MyStrandVector))==1){  #test for pure-plus and pure-minus
      #20180409 Add "Scramble-but-not-mix"
      #20191216 fix due to S5215
      #"Plus" case:
      if (MyStrandVector[1]=="plus"){
        MyTestScramble <- is.unsorted(MyDF_638$sstart)
        if (MyTestScramble==TRUE){
          MyStrandBlastHXB2_FINAL <- append(MyStrandBlastHXB2_FINAL,"plusScramble")
        } else {
          MyStrandBlastHXB2_FINAL <- append(MyStrandBlastHXB2_FINAL,MyStrandVector[1])
        }
        #"minus" case:  
      } else if (MyStrandVector[1]=="minus") {
        MyTestScramble <- is.unsorted(rev(MyDF_638$send))
        if (MyTestScramble==TRUE){
          MyStrandBlastHXB2_FINAL <- append(MyStrandBlastHXB2_FINAL,"minusScramble")
        } else {
          MyStrandBlastHXB2_FINAL <- append(MyStrandBlastHXB2_FINAL,MyStrandVector[1])
        }
      } else {
        MyStrandBlastHXB2_FINAL <- append(MyStrandBlastHXB2_FINAL,"Check")
      }
      MyHXB2MapStart <- append(MyHXB2MapStart,paste(length(MyStrandVector),"UniqMatches",sep=""))
      MyHXB2MapEnd <- append(MyHXB2MapEnd,paste(length(MyStrandVector),"UniqMatches",sep=""))
    } else {
      MyStrandBlastHXB2_FINAL <- append(MyStrandBlastHXB2_FINAL,"mix")
      MyHXB2MapStart <- append(MyHXB2MapStart,paste(length(MyStrandVector),"UniqMatches",sep=""))
      MyHXB2MapEnd <- append(MyHXB2MapEnd,paste(length(MyStrandVector),"UniqMatches",sep=""))
    }
  } else {
    MyStrandBlastHXB2_FINAL <- append(MyStrandBlastHXB2_FINAL,"Check")
    MyHXB2MapStart <- append(MyHXB2MapStart,NA)
    MyHXB2MapEnd <- append(MyHXB2MapEnd,NA)
  }
}
DF_FINAL <- cbind(DF_FINAL,MyStrandBlastHXB2_FINAL)
DF_FINAL <- cbind(DF_FINAL,MyHXB2MapStart)
DF_FINAL <- cbind(DF_FINAL,MyHXB2MapEnd)
# View(DF_FINAL)

#02c. Output FinalSEQ  ### IF MIX+- OR NoMapping with NBCI_HXB2, always output the RAW_MGH SEQ ###
MyFinalSEQ <- c()
for (k in 1:nrow(DF_FINAL)){
  if(DF_FINAL$MyStrandBlastHXB2_FINAL[k]=="minus"){
    MyFinalSEQ <- append(MyFinalSEQ,as.character(DF_FINAL$SEQRAWRC[k]))
  } else if (DF_FINAL$MyStrandBlastHXB2_FINAL[k]=="plus"){
    MyFinalSEQ <- append(MyFinalSEQ,as.character(DF_FINAL$SEQRAW[k]))
  } else if (DF_FINAL$MyStrandBlastHXB2_FINAL[k]=="mix"){
    MyFinalSEQ <- append(MyFinalSEQ,as.character(DF_FINAL$SEQRAW[k]))
  } else if (DF_FINAL$MyStrandBlastHXB2_FINAL[k]=="NoMapping"){
    MyFinalSEQ <- append(MyFinalSEQ,as.character(DF_FINAL$SEQRAW[k]))
    #Added 20180409 for Scramble-not-mix
    #20191216 ScrambleFix
  } else if (DF_FINAL$MyStrandBlastHXB2_FINAL[k]=="plusScramble") {
    MyFinalSEQ <- append(MyFinalSEQ,as.character(DF_FINAL$SEQRAW[k]))
  } else if (DF_FINAL$MyStrandBlastHXB2_FINAL[k]=="minusScramble") {
    MyFinalSEQ <- append(MyFinalSEQ,as.character(DF_FINAL$SEQRAWRC[k]))
  } else if (DF_FINAL$MyStrandBlastHXB2_FINAL[k]=="Check") {
    MyFinalSEQ <- append(MyFinalSEQ,as.character(DF_FINAL$SEQRAW[k]))
  } else {
    MyFinalSEQ <- append(MyFinalSEQ,as.character(DF_FINAL$SEQRAW[k]))
  }
}
DF_FINAL <- cbind(DF_FINAL,MyFinalSEQ)
#02d. Length
MyLenFinal <- c()
for (k in 1:nrow(DF_FINAL)){
  MyLenFinal <- append(MyLenFinal,nchar(MyFinalSEQ[k]))
}
DF_FINAL <- cbind(DF_FINAL,MyLenFinal)
DF_FINAL <- cbind(DF_FINAL,MyAlignedLength)


#Add len info to Unique SampleID_Contig_Len
MyUniqSampleContigID_Len <- mapply(paste,DF_FINAL$MyUniqSampleContigID,"_",DF_FINAL$MyLenFinal,sep="")
DF_FINAL <- cbind(DF_FINAL,MyUniqSampleContigID_Len)

#View(DF_FINAL)

############## HYPERMUT #########################

MyDF_Hypermut <- subset(DF_FINAL,DF_FINAL$MyHIVTF_FINAL=="TRUE" & DF_FINAL$MyLenFinal>=8000 & (DF_FINAL$MyStrandBlastHXB2_FINAL=="plus" | DF_FINAL$MyStrandBlastHXB2_FINAL=="minus"))

#View(MyDF_Hypermut)

#If there is no row after subsetting, fill in all fields with NA and move on
if (nrow(MyDF_Hypermut)==0) {
  print("No sequence made it to MyDF_Hypermut, 5Psi and gag check; create empty table")
  DF04_Hypermut <- data.frame(a=as.character(DF_FINAL$SEQID),b="NA",c="NA",d="NA",e="NA",f="NA",g="NA",h="NA",i="NA",j="NA",k="NA",l="NA")
  colnames(DF04_Hypermut) <- c("SEQID","Muts@MutSite","Outof@MutSite","Muts@ControlSites","Outof@ControlSites","RateRatio","FisherExactPvalue","MyHypermutTF","MyPrimerTF","MyPsiDELrelative","MyPsiINrelatiive","MyGagATG")
  #View(DF04_Hypermut)
  DF_FINAL <- merge(DF_FINAL,DF04_Hypermut,by="SEQID",all=TRUE)  #no omission; all rows retained
} else {
  #do nothing, move onto in-house Los Alamos Hypermut below
}
#View(DF_FINAL)

#03ab. Export FinalSEQ that are (1) HIV AND (2) >=8000bp as FASTA + Add HXB2-9719 AND (3) Genome orientation must be pure-plus or pure-minus, no "mix" no "no mapping" (4) Must only have 1UniqMapping in HXB2-blast & must map to one of the primer positions FLv1 and FLv2:  623,638,769,9632,9686
MyHIVDNAStringSetFinal <- DNAStringSet()
for (k in 1:nrow(MyDF_Hypermut)){
  MySEQ <- DNAStringSet(MyDF_Hypermut$MyFinalSEQ[k])
  names(MySEQ) <- MyDF_Hypermut$SEQID[k]
  MyHIVDNAStringSetFinal <- append(MyHIVDNAStringSetFinal,MySEQ)
  #print(MyDF_Hypermut$SEQID[k])
}
### Add HXB2 to top ###
HXB29719 <- readDNAStringSet("R_HXB2.fasta")
MyHIVDNAStringSetFinal <- append(MyHIVDNAStringSetFinal,HXB29719,0)
### Write to file ###
writeXStringSet(MyHIVDNAStringSetFinal,"Output_HIV_Over8000_toHypermut_ver05fixscramble.fasta")
length(MyHIVDNAStringSetFinal)

#03c. Align with muscle
#Step 01. User Settings
InputFASTAFileName <- "Output_HIV_Over8000_toHypermut_ver05fixscramble.fasta"
#Step 02. Muscle Multiple Sequence Alignment
MyFASTA <- readDNAStringSet(InputFASTAFileName)
MyNUCAlign <- muscle(MyFASTA
                     ,gapopen=-1000) #20180301 Riddhima
writeXStringSet(DNAStringSet(MyNUCAlign),paste(InputFASTAFileName,"_Raligned.fasta",sep=""))

## MANUALLY REMOVE WEIRD ALIGNMENTS HERE ##

## AUTOTRIM ##
InputFASTAFileNameRaligned <- paste(InputFASTAFileName,"_Raligned.fasta",sep="")
MyFASTA_InputFASTAFileNameRaligned <- readDNAStringSet(InputFASTAFileNameRaligned)
#Step1. Index Beginning
MyCoordinatesF_List <- vmatchPattern(Primer2ndF_HXB2,MyFASTA_InputFASTAFileNameRaligned[1],max.mismatch=1,fixed=TRUE)
MyCoordinatesF_Unlist <- unlist(startIndex(MyCoordinatesF_List))[1]
#Step2. Index End
MyCoordinatesR_List <- vmatchPattern(Primer2ndR_HXB2,MyFASTA_InputFASTAFileNameRaligned[1],max.mismatch=1,fixed=TRUE)
MyCoordinatesR_Unlist <- unlist(startIndex(MyCoordinatesR_List))[length(unlist(startIndex(MyCoordinatesR_List)))] + nchar(Primer2ndR_HXB2) - 1
#Step3. Trim
MyFASTA_InputFASTAFileNameRaligned_Trimmed <- subseq(MyFASTA_InputFASTAFileNameRaligned,start=MyCoordinatesF_Unlist,end=MyCoordinatesR_Unlist)
writeXStringSet(MyFASTA_InputFASTAFileNameRaligned_Trimmed,paste(InputFASTAFileNameRaligned,"_GLtrimmed.fasta",sep=""))
#Step4. Trim_psigag
MyCoordinatesGag_List <- vmatchPattern("GGTGCGAGAGCGTCAGTATT",MyFASTA_InputFASTAFileNameRaligned_Trimmed[1],max.mismatch=0,fixed=TRUE)
MyCoordinatesGag_Unlist <- unlist(startIndex(MyCoordinatesGag_List))[1] - 1
MyFASTA_InputFASTAFileNameRaligned_Trimmed_PsiGag <- subseq(MyFASTA_InputFASTAFileNameRaligned_Trimmed,start=1,end=MyCoordinatesGag_Unlist)
writeXStringSet(MyFASTA_InputFASTAFileNameRaligned_Trimmed_PsiGag,paste(InputFASTAFileNameRaligned,"_GLtrimmed_psigag.fasta",sep=""))

#####################
#R_Hypermut_ver01.R
#####################

MyInhouseHypermut <- function() {
  #STEP1.  LOAD THE FASTA FILE and CONVERT into DF
  MyFASTA <- readDNAStringSet("Output_HIV_Over8000_toHypermut_ver05fixscramble.fasta_Raligned.fasta_GLtrimmed.fasta")
  MyOutput01_SampleID <- c()
  MyOutput02_RawSeq <- c()
  for (i in 1:length(MyFASTA)){
    MyOutput01_SampleID <- append(MyOutput01_SampleID,names(MyFASTA[i]))
    MyOutput02_RawSeq <- append(MyOutput02_RawSeq,as.character(MyFASTA[i]))
  }
  
  #STEP2.  Index REF for 3G and 3F sites
  MyREF01_Length <- nchar(MyOutput02_RawSeq[1])
  MyResult <- ""
  for (j in 1:MyREF01_Length){
    #2a. Count and Output category per-base-position: C=Control, N=not-G, G=3G (GGD), F=3F (GAD)
    #Site1.  Is it G?
    if (substring(MyOutput02_RawSeq[1],j,j) == "G"){
      #Site2.  Is it -?  If so, check next position until it is not -.  Then see if it is G or A.
      k=1
      while(substring(MyOutput02_RawSeq[1],j+k,j+k) == "-"){k <- k+1}
      if (substring(MyOutput02_RawSeq[1],j+k,j+k) == "G" | substring(MyOutput02_RawSeq[1],j+k,j+k) == "A"){
        #Site3.  Is it -?  If so, check next position until it is not -.  Then see if it is !C.
        n=1
        while(substring(MyOutput02_RawSeq[1],j+k+n,j+k+n) == "-"){n <- n+1}
        if (substring(MyOutput02_RawSeq[1],j+k,j+k) == "G" & substring(MyOutput02_RawSeq[1],j+k+n,j+k+n) != "C"){
          MyResult <- paste(MyResult,"G",sep="")
        } else if (substring(MyOutput02_RawSeq[1],j+k,j+k) == "A" & substring(MyOutput02_RawSeq[1],j+k+n,j+k+n) != "C"){
          MyResult <- paste(MyResult,"F",sep="")
        } else {MyResult <- paste(MyResult,"C",sep="")}
      } else if (substring(MyOutput02_RawSeq[1],j+k,j+k) == "-"){
        MyResult <- paste(MyResult,"N",sep="")
      } else {MyResult <- paste(MyResult,"C",sep="")}
    } else {MyResult <- paste(MyResult,"N",sep="")}
  }
  write(MyResult,file="Output_Hypermut01_REFCategory.txt")
  MyREFIndex_G <- unlist(gregexpr(pattern="G",MyResult))
  MyREFIndex_F <- unlist(gregexpr(pattern="F",MyResult))
  MyREFIndex_C <- unlist(gregexpr(pattern="C",MyResult))
  MyREFIndex_N <- unlist(gregexpr(pattern="N",MyResult))
  
  #STEP3.  COUNT NUMBERS of 3G and 3F ## PLus ## revert A>G in original
  MyOutput03a1_APOBEC3G2A_count <- c()
  MyOutput03b1_APOBEC3F2A_count <- c()
  MyOutput03c1_APOBEC3FGsum2A_count <- c()
  MyOutput03d1_Control_count <- c()
  MyRevertedSEQ_A2G <- c()
  MyQuery_GsumIndex_AllSamples <- c()
  MyQuery_FsumIndex_AllSamples <- c()
  MyQuery_GFsumIndex_AllSamples <- c()
  
  for (m in 1:length(MyOutput02_RawSeq)){
    MyQuery_Gindex <- c()
    MyQuery_Findex <- c()
    #3a. APOBEC 3G Count
    g=0
    for (o in 1:length(MyREFIndex_G)){
      if(substring(MyOutput02_RawSeq[m],MyREFIndex_G[o],MyREFIndex_G[o]) == "A"){
        g <- g+1
        MyQuery_Gindex <- append(MyQuery_Gindex,MyREFIndex_G[o])
      }
    } 
    MyOutput03a1_APOBEC3G2A_count <- append(MyOutput03a1_APOBEC3G2A_count,g)
    if (m!=1){
      MyQuery_GsumIndex_AllSamples <- append(MyQuery_GsumIndex_AllSamples,MyQuery_Gindex)
    }
    #3b. APOBEC 3F Count
    f=0
    for (p in 1:length(MyREFIndex_F)){
      if(substring(MyOutput02_RawSeq[m],MyREFIndex_F[p],MyREFIndex_F[p]) == "A"){
        f <- f+1
        MyQuery_Findex <- append(MyQuery_Findex,MyREFIndex_F[p])
      }
    } 
    MyOutput03b1_APOBEC3F2A_count <- append(MyOutput03b1_APOBEC3F2A_count,f)
    MyOutput03c1_APOBEC3FGsum2A_count <- append(MyOutput03c1_APOBEC3FGsum2A_count,f+g)
    if (m!=1){
      MyQuery_FsumIndex_AllSamples <- append(MyQuery_FsumIndex_AllSamples,MyQuery_Findex)
    }
    #3c. Control Count
    c=0
    for (q in 1:length(MyREFIndex_C)){
      if(substring(MyOutput02_RawSeq[m],MyREFIndex_C[q],MyREFIndex_C[q]) == "A"){c <- c+1}
    } 
    MyOutput03d1_Control_count <- append(MyOutput03d1_Control_count,c)
    #3d. Revert A2G in hypermuted positions
    if (length(append(MyQuery_Gindex,MyQuery_Findex))==0){
      MyRevertedSEQ_A2G <- append(MyRevertedSEQ_A2G,MyOutput02_RawSeq[m])
    } else {
      #REVERT
      MyCurrentSeq <- MyOutput02_RawSeq[m]
      MyQuery_GFsumIndex <- sort(append(MyQuery_Gindex,MyQuery_Findex))
      for (x in 1:length(MyQuery_GFsumIndex)){
        Start <- MyQuery_GFsumIndex[x]
        End <- MyQuery_GFsumIndex[x]
        substr(MyCurrentSeq,Start,End) <- "G"
      }
      MyRevertedSEQ_A2G <- append(MyRevertedSEQ_A2G,MyCurrentSeq)
      #3e. Genome Distribution append
      MyQuery_GFsumIndex_AllSamples <- append(MyQuery_GFsumIndex_AllSamples,MyQuery_GFsumIndex)
    }
  }
  MyFASTA_Fixed <- DNAStringSet(MyRevertedSEQ_A2G)
  names(MyFASTA_Fixed) <- names(MyFASTA) 
  writeXStringSet(MyFASTA_Fixed,file="Output_Hypermut02_MyFASTA_Reverted.fasta")
  
  
  #STEP3.1 Summary data: NOT
  MyOutput03a2_APOBEC3G2A_NOTcount <- length(MyREFIndex_G)-MyOutput03a1_APOBEC3G2A_count
  MyOutput03b2_APOBEC3F2A_NOTcount <- length(MyREFIndex_F)-MyOutput03b1_APOBEC3F2A_count
  MyOutput03c2_APOBEC3FGsum2A_NOTcount <- length(MyREFIndex_G)+length(MyREFIndex_F)-MyOutput03c1_APOBEC3FGsum2A_count
  MyOutput03d2_Control_NOTcount <- length(MyREFIndex_C)-MyOutput03d1_Control_count
  #STEP3.2 Summary data: TOTAL
  MyOutput03a3_APOBEC3G2A_TOTALcount <- rep(length(MyREFIndex_G),length(MyOutput01_SampleID))
  MyOutput03b3_APOBEC3F2A_TOTALcount <- rep(length(MyREFIndex_F),length(MyOutput01_SampleID))
  MyOutput03c3_APOBEC3FGsum2A_TOTALcount <- rep(length(MyREFIndex_G)+length(MyREFIndex_F),length(MyOutput01_SampleID))
  MyOutput03d3_Control_TOTALcount <- rep(length(MyREFIndex_C),length(MyOutput01_SampleID))
  #STEP3.3 Summary data:  PERCENT HypermutedG>A/TotalG
  MyOutput03a4_APOBEC3G2A_PERC <- MyOutput03a1_APOBEC3G2A_count/length(MyREFIndex_G)*100
  MyOutput03b4_APOBEC3F2A_PERC <- MyOutput03b1_APOBEC3F2A_count/length(MyREFIndex_F)*100
  MyOutput03c4_APOBEC3FGsum2A_PERC <- MyOutput03c1_APOBEC3FGsum2A_count/(length(MyREFIndex_G)+length(MyREFIndex_F))*100
  MyOutput03d4_Control_PERC <- MyOutput03d1_Control_count/length(MyREFIndex_C)*100
  #STEP3.4 Summary data:  PERCENT Contribution G sites / total mutated sites
  MyOutput03a5_APOBEC3G2A_PERC <- MyOutput03a1_APOBEC3G2A_count/(MyOutput03a1_APOBEC3G2A_count+MyOutput03b1_APOBEC3F2A_count)*100
  MyOutput03b5_APOBEC3F2A_PERC <- MyOutput03b1_APOBEC3F2A_count/(MyOutput03a1_APOBEC3G2A_count+MyOutput03b1_APOBEC3F2A_count)*100
  
  
  #STEP 4.  Reproduce WebHypermut 2.0 tool:  FISHER'S
  MyFisherVector <- c()
  for (m in 1:length(MyOutput02_RawSeq)){
    Q1 <- MyOutput03c1_APOBEC3FGsum2A_count[m]
    Q2 <- MyOutput03c2_APOBEC3FGsum2A_NOTcount[m]
    Q3 <- MyOutput03d1_Control_count[m]
    Q4 <- MyOutput03d2_Control_NOTcount[m]
    table <- rbind(c(Q1,Q2),c(Q3,Q4))
    MyFisher <- fisher.test(table,alternative="greater")
    MyFisherVector <- append(MyFisherVector,MyFisher$p.value)
  }
  
  #STEP 6.  Export DF
  MyDF <- as.data.frame(cbind(
    MyOutput01_SampleID
    , MyOutput02_RawSeq
    , MyOutput03a1_APOBEC3G2A_count
    , MyOutput03a2_APOBEC3G2A_NOTcount
    , MyOutput03a3_APOBEC3G2A_TOTALcount
    , MyOutput03a4_APOBEC3G2A_PERC
    , MyOutput03a5_APOBEC3G2A_PERC
    , MyOutput03b1_APOBEC3F2A_count
    , MyOutput03b2_APOBEC3F2A_NOTcount
    , MyOutput03b3_APOBEC3F2A_TOTALcount
    , MyOutput03b4_APOBEC3F2A_PERC
    , MyOutput03b5_APOBEC3F2A_PERC
    , MyOutput03c1_APOBEC3FGsum2A_count
    , MyOutput03c2_APOBEC3FGsum2A_NOTcount
    , MyOutput03c3_APOBEC3FGsum2A_TOTALcount
    , MyOutput03c4_APOBEC3FGsum2A_PERC
    , MyOutput03d1_Control_count
    , MyOutput03d2_Control_NOTcount
    , MyOutput03d3_Control_TOTALcount
    , MyOutput03d4_Control_PERC
    , MyRevertedSEQ_A2G
    , as.numeric(MyFisherVector)
  ))
  write.csv(MyDF,file="Output_Hypermut03_MyDF_Summary.csv")
}

MyInhouseHypermut()
DF04_Hypermut <- read.csv("Output_Hypermut03_MyDF_Summary.csv",header=TRUE)[-1,-1]
colnames(DF04_Hypermut) <- c("SEQID"
                             , "MyOutput02_RawSeq"
                             , "MyOutput03a1_APOBEC3G2A_count"
                             , "MyOutput03a2_APOBEC3G2A_NOTcount"
                             , "MyOutput03a3_APOBEC3G2A_TOTALcount"
                             , "MyOutput03a4_APOBEC3G2A_PERC"
                             , "MyOutput03a5_APOBEC3G2A_PERC"
                             , "MyOutput03b1_APOBEC3F2A_count"
                             , "MyOutput03b2_APOBEC3F2A_NOTcount"
                             , "MyOutput03b3_APOBEC3F2A_TOTALcount"
                             , "MyOutput03b4_APOBEC3F2A_PERC"
                             , "MyOutput03b5_APOBEC3F2A_PERC"
                             , "MyOutput03c1_APOBEC3FGsum2A_count"
                             , "MyOutput03c2_APOBEC3FGsum2A_NOTcount"
                             , "MyOutput03c3_APOBEC3FGsum2A_TOTALcount"
                             , "MyOutput03c4_APOBEC3FGsum2A_PERC"
                             , "MyOutput03d1_Control_count"
                             , "MyOutput03d2_Control_NOTcount"
                             , "MyOutput03d3_Control_TOTALcount"
                             , "MyOutput03d4_Control_PERC"
                             , "MyRevertedSEQ_A2G"
                             , "FisherExactPvalue")
MyHypermutTF <- DF04_Hypermut$FisherExactPvalue<=0.05
DF04_Hypermut <- cbind(DF04_Hypermut,MyHypermutTF)
DF_FINAL <- merge(DF_FINAL,DF04_Hypermut,by="SEQID",all=TRUE)  #no omission; all rows retained



################################################################
### 5' CHECK ### 
##################################################################

MyPsiInput <- readDNAStringSet("Output_HIV_Over8000_toHypermut_ver05fixscramble.fasta_Raligned.fasta_GLtrimmed_psigag.fasta")
HXB2gapPostAlignment <- as.numeric(nchar(as.character(MyPsiInput[1]))-nchar(gsub("-","",as.character(MyPsiInput[1]))))
HXBPsi_Aligned <- as.character(MyPsiInput[1])
HXBPsiLen <- nchar(gsub("-","",MyPsiInput[1]))
MyPsiInput <- MyPsiInput[2:length(MyPsiInput)]
SEQID <- names(MyPsiInput)
SEQ <- as.character(MyPsiInput)

### Primer: Yes/No within 100bp mismatchPattern returns TRUE with 5bp max.mismatch allowance
#Primer2ndF: U5-638F: 5'-GCGCCCGAACAGGGACYTGAAARCGAAAG-3'  len=29
#Primer2ndR: U5-547R: 5'-GCACTCAAGGCAAGCTTTATTGAGGCTTA-3' RC: TAAGCCTCAATAAAGCTTGCCTTGAGTGC
Primer2ndF <- DNAString("GCGCCCGAACAGGGACYTGAAARCGAAAG")  #fixed=FALSE allows Y to match with C and T
MyPrimerTF <- c()
for(n in 1:length(SEQ)) {
  #Step1. Yes/No Primer
  MyMatchCount <- length(matchPattern(Primer2ndF,DNAString(substring(SEQ[n],1,100)),max.mismatch=5,fixed=FALSE))
  if (MyMatchCount==1){ 
    MyPrimerTF <- append(MyPrimerTF,"YesPrimer")
  } else if (MyMatchCount==0) {
    MyPrimerTF <- append(MyPrimerTF,"NoPrimer")
  } else {
    MyPrimerTF <- append(MyPrimerTF,"CheckPrimer")
  }
}

### Count Indel ###
MyPsiDELrelative <- c()
MyPsiINrelatiive <- c()
for (n in 1:length(SEQ)){
  CountDEL <- 0
  CountIN <- 0
  for(c in 1:nchar(HXBPsi_Aligned)){
    MyCurrentNUC_HXB2 <- substring(HXBPsi_Aligned,c,c)
    MyCurrentNUC_QUERY <- substring(SEQ[n],c,c)
    if(MyCurrentNUC_HXB2!="-"&MyCurrentNUC_QUERY=="-"){
      CountDEL <- CountDEL + 1
    } else if (MyCurrentNUC_HXB2=="-"&MyCurrentNUC_QUERY!="-"){
      CountIN <- CountIN + 1
    } else {
      CountDEL <- CountDEL + 0
      CountIN <- CountIN + 0
    }
  }
  MyPsiDELrelative <- append(MyPsiDELrelative,as.numeric(CountDEL))
  MyPsiINrelatiive <- append(MyPsiINrelatiive,CountIN) 
}

### Check Gag ATG ###
MyGagATG <- logical()
for (n in 1:length(SEQ)){
  MyGagTF <- substring(SEQ[n],nchar(SEQ[n])-2,nchar(SEQ[n]))=="ATG"
  MyGagATG <- append(MyGagATG,MyGagTF)
}

My5Prime_Final <- cbind(SEQID
                        ,MyPrimerTF
                        ,MyPsiDELrelative
                        ,MyPsiINrelatiive
                        ,MyGagATG
)
#View(My5Prime_Final)

DF_FINAL <- merge(DF_FINAL,My5Prime_Final,by="SEQID",all=TRUE)

#View(DF_FINAL)
write.csv(DF_FINAL,paste("Output_MyBigSummary_DF_FINAL_half.csv",sep=""))

#########################################################

### Subset for HIVSeqinR and convert to FASTA: Anything that is not hypermut and not scramble  ###

dir.create("./ToHIVSeqinR/")
MyHIVSeqinR_InputFileName <- "./ToHIVSeqinR/Output_toHIVSeqinR.fasta"
DF_HIVSeqinR <- subset(DF_FINAL,DF_FINAL$MyHypermutTF==FALSE & DF_FINAL$MyHIVTF_FINAL==TRUE & DF_FINAL$MyStrand!="mix")
#View(DF_HIVSeqinR)
if (nrow(DF_HIVSeqinR)>=1){
  for (n in 1:nrow(DF_HIVSeqinR)){
    #print(DF_HIVSeqinR$MyVerdict[n])
    #print(DF_HIVSeqinR$SEQID[n])
    #print(DF_HIVSeqinR$MyFinalSEQ[n])
    MyHIVSeqinR_Output_SEQID <- paste(">",DF_HIVSeqinR$SEQID[n],sep="")
    MyHIVSeqinR_Output_FinalSEQ <- as.character(DF_HIVSeqinR$MyFinalSEQ[n])
    write(MyHIVSeqinR_Output_SEQID,file=MyHIVSeqinR_InputFileName,append=TRUE)
    write(MyHIVSeqinR_Output_FinalSEQ,file=MyHIVSeqinR_InputFileName,append=TRUE)
  }
} else {
  print("No sequence made it to ToHIVSeqinR; create empty table")
  MyHIVSeqinRPivot <- data.frame(a=as.character(DF_FINAL$SEQID),b="NA",c="NA",d="NA",e="NA",f="NA",g="NA",h="NA",i="NA",j="NA",k="NA",l="NA",m="NA",n="NA",o="NA",p="NA",q="NA",r="NA",s="NA",t="NA",u="NA",v="NA",w="NA",x="NA",y="NA",z="NA")
  colnames(MyHIVSeqinRPivot) <- c("MyHIVSeqinR_OutputCSV_DF_Uniq","MyGAGPassFail","MyPOLPassFail","MyPRPassFail","MyRTPassFail","MyRNaseHPassFail","MyINTPassFail","MyENVPassFail","MyGP41PassFail","MyAnyATGFailed","MyGAGLenUntrimmed","MyPOLLenUntrimmed","MyPRLenUntrimmed","MyRTLenUntrimmed","MyRNaseHLenUntrimmed","MyINTLenUntrimmed","MyENVLenUntrimmed","MyGP41LenUntrimmed","MyGAGLenHXB2","MyPOLLenHXB2","MyPRLenHXB2","MyRTLenHXB2","MyRNaseHLenHXB2","MyINTLenHXB2","MyENVLenHXB2","MyGP41LenHXB2")
  #View(DF04_Hypermut)
  DF_FINAL <- merge(DF_FINAL,MyHIVSeqinRPivot,by.x="SEQID",by.y="MyHIVSeqinR_OutputCSV_DF_Uniq",all=TRUE)  #no omission; all rows retained
}
#View(DF_FINAL)


###################################
#R05_MyPipeline ver02.6_Misalign.R
###################################

#Step 00a.  R Start up
rm(list= ls()[!(ls() %in% c("MyBlastnDir","StartSysTime","Primer2ndF_HXB2","Primer2ndR_HXB2"))])
closeAllConnections()
graphics.off()
cat("\014") #send CTRL+L to R, clean console
Sys.setenv(TZ="EST")
MyWD <- getwd()
setwd(MyWD)
library(Biostrings) #PC@work Biostrings_2.38.1
library(muscle) #PC@work muscle_3.12.0


#Step00b. User input
current_pipeline_version <- "ver02.6"

#GLOBAL#
dir.create("./ToHIVSeqinR_DB/")

#RAM Efficiency#
AllFiles <- list.files("./ToHIVSeqinR")
for (z in 1:length(AllFiles)){
  rm(list= ls()[!(ls() %in% c("AllFiles","z","current_pipeline_version","StartSysTime","MyBlastnDir","Primer2ndF_HXB2","Primer2ndR_HXB2"))])
  MyInputFASTA <- paste("./ToHIVSeqinR/",AllFiles[z],sep="")
  MyTargetNuc <- readDNAStringSet(MyInputFASTA) 
  MyRefNuc_HXB2 <- readDNAStringSet("HIV Genome RefSeq ver04 NUC.fasta") ###CHANGE_ME###
  MyRefAA_HXB2 <- readAAStringSet("HIV Genome RefSeq ver04 AA.fasta") ###CHANGE_ME###
  #Step 01c. Muscle settings
  MyGapPenaltyNUC <- -5000 ###CHANGE_ME### enter interger from -1 to negative infinity; larger negative == more stringent 
  MuscleNUCLog <- paste("Output_HIVSeqinR_",current_pipeline_version,"_muscle_nuc_log.txt",sep="") ###CHANGE_ME###
  MyGapPenaltyAA <- -1 ###CHANGE_ME### 
  MuscleAALog <- paste("Output_HIVSeqinR_",current_pipeline_version,"_muscle_aa_log.txt",sep="") ###CHANGE_ME### 
  #Step 01d. IUPAC AA list
  IUPAC_AA <- "ABCDEFGHIKLMNPQRSTVWXYZ" #This list is meant for HXB2, and therefore has no stop (otherwise need to use \\*)
  #QC2b: Nuc %iden cutoff setting
  PerIden_Cutoff_Nuc <- 0.7 ###CHANGE_ME### this value was arbituary set on 20151130 after looking at HXB2 vs SEQ0002 A1 (ACH2.Well7B) and A2 (Nina's N3.Well-1G); changed to 0.76 on 20171127 for K0214_HIV_AB10_CONTIG_8713_p1
  #QC3:  AA Too Short %
  AA_TooShort <- 0.95 ###CHANGE_ME### this value was arbituary set on 20151130. It means "if observed AA is less than 50% of HXB2's AA length, flag as 'failed QC'
  AA_TooShort_Gag <- 0.97 #This allows for 15/53aa missing in p6
  ## AA_TooShort is 0.95 expect gag 0.994 to adjust for p6
  AA_TooLong <- 1.2 #added 20171205
  ##20161016 (1) 50% was way too low... 8E5 would pass because pol's stop 
  ##20161016 (2) after testing 90% and 95%, I decide to finalize it at 95% to flag EVERYTHING that is even a little suspicious... for manual checking
  ##20161016 (3) The end result should be compared to the user's manual check
  MyWorkingAln_Nuc_fasta <- paste("Output_HIVSeqinR_",current_pipeline_version,"_MyWorkingAln_Nuc_muscleAlign_full_untrimmed.fasta",sep="")  ###CHANGE_ME###
  MyWorkingAln_Nuc_muscleAlign <- paste("Output_HIVSeqinR_",current_pipeline_version,"_MyWorkingAln_Nuc_muscleAlign_trimmed.fasta",sep="") ###CHANGE_ME###
  MyWorkingAln_AA_muscleAlign <- paste("Output_HIVSeqinR_",current_pipeline_version,"_MyWorkingAln_AA_muscleAlign_oneb4stop.fasta",sep="") ###CHANGE_ME###
  MyWorkingAln_AA_AAStringSet_trimmed <- paste("Output_HIVSeqinR_",current_pipeline_version,"_MyWorkingAln_NUC2AA_muscleAlign_trimmed.fasta",sep="")
  ################################################################
  #Step 02.  Set up df for final data output
  
  Output_df <- paste("./ToHIVSeqinR_DB/Output_HIVSeqinR",current_pipeline_version,"_mydf",z,".csv",sep="") ###CHANGE_ME###
  Output_01_MyWorkingRefNuc_Name <- c() 
  length(Output_01_MyWorkingRefNuc_Name) #36 OK :)
  Output_01a_MyWorkingRefNuc_UnAlign_Seq <- c()
  length(Output_01a_MyWorkingRefNuc_UnAlign_Seq)
  Output_02_MyWorkingPreAlnNuc_Name <-c()
  length(Output_02_MyWorkingPreAlnNuc_Name) #36 OK :)
  Output_02a_QC0_TargetSeqATCG <- c()
  length(Output_02a_QC0_TargetSeqATCG)
  Output_03a_MyWorkingAln_NucRef <- c()
  length(Output_03a_MyWorkingAln_NucRef) #36 OK :)
  Output_03b_MyWorkingAln_NucTarget <- c()
  length(Output_03b_MyWorkingAln_NucTarget) #36 OK :)
  #MyWorkingAln_Nuc.fasta
  Output_03a1_MyWorkingAln_NucRef_Len <- c()
  length(Output_03a1_MyWorkingAln_NucRef_Len)  #36 OK :)
  Output_03a2_MyWorkingAln_NucTarget_Len <- c()
  length(Output_03a2_MyWorkingAln_NucTarget_Len) #36 OK :)
  Output_04a_MyIndex_start <- c()
  length(Output_04a_MyIndex_start) #36 OK :)
  Output_04b_MyIndex_end <- c()  #20151203 development stopped for vmatchPattern; ignore results
  length(Output_04b_MyIndex_end) #36 OK :)
  Output_04c_MyIndex_ATCG_start <- c()  #20151203 development stopped for vmatchPattern; ignore results
  length(Output_04c_MyIndex_ATCG_start) #36 OK :)
  Output_04d_MyIndex_ATCG_end <- c()
  length(Output_04d_MyIndex_ATCG_end) #36 OK :)
  Output_05a_QC1_lentrim1vs2 <- c() #QC1. trim1 and trim2 must be identical in seq len
  length(Output_05a_QC1_lentrim1vs2) 
  Output_05b_QC2a_PercIdentity <- c() #QC2. %iden #MyNote: It's up to your own judgement to decide which genomic region has %identity match that is way too low
  length(Output_05b_QC2a_PercIdentity)
  Output_05c_QC2b_PerIden_Cutoff_Nuc <- c()
  length(Output_05c_QC2b_PerIden_Cutoff_Nuc)
  Output_06a_MyFinal_AA <- c()
  length(Output_06a_MyFinal_AA)
  Output_06b_MyFinal_AA_len <- c()
  length(Output_06b_MyFinal_AA_len)
  Output_06b2_MyAA_Expected_len <- c()
  length(Output_06b2_MyAA_Expected_len)
  Output_06c_MyFinal_Nuc <- c() #This Nuc is the one-before-stop, no dashes
  length(Output_06c_MyFinal_Nuc)
  Output_06d_MyFinal_Nuc_len <- c()
  length(Output_06d_MyFinal_Nuc_len)
  Output_06e_MyFinal_Stop_4bp <- c()
  length(Output_06e_MyFinal_Stop_4bp)
  Output_06f_QC3_AA_TooShort <- c()
  length(Output_06f_QC3_AA_TooShort)
  Output_07a_AA_trimmed_HXB2 <- c()
  length(Output_07a_AA_trimmed_HXB2)
  Output_07b_AA_trimmed_sample <- c()
  length(Output_07b_AA_trimmed_sample)
  Output_08a_ATG <- c()
  length(Output_08a_ATG)
  Output_08b_QC4_ATG <- c()
  length (Output_08b_QC4_ATG)

  
  Store_TATexon1_NUC_nodash <- c()
  Store_TATexon2GP41_NUC_nodash <- c()
  Store_REVexon1_NUC_nodash <- c()
  Store_REVexon2GP41_NUC_nodash <- c()
  
  ################################################################
  #Step 03.  Build Loop Framework:  Loop through all items in RefNuc and TargetNuc files (load gene first, then load TargetNuc to align with it)
  for (i in 1:length(MyRefNuc_HXB2)){
    #print (MyRefNuc_HXB2[i]) #Loop retrieve each items in MyRefNuc_HXB2, one gene at a time
    #print (names(MyRefNuc_HXB2[i])) ###Print Progress###
    for (k in 1:length(MyTargetNuc)){
      Output_01_MyWorkingRefNuc_Name <- append(Output_01_MyWorkingRefNuc_Name,names(MyRefNuc_HXB2[i])) #write current working ref into df
      #print(MyTargetNuc[k])
      #print(names(MyTargetNuc[k])) ###Print Progress###
      Output_01a_MyWorkingRefNuc_UnAlign_Seq <- append(Output_01a_MyWorkingRefNuc_UnAlign_Seq,toString(MyRefNuc_HXB2[i]))
      Output_02_MyWorkingPreAlnNuc_Name <- append(Output_02_MyWorkingPreAlnNuc_Name,names(MyTargetNuc[k]))
      ##################QC0 for any non-ATCG characters in the Target Seq file#################################33
      QC0_NonA <- gsub("A","",toString(MyTargetNuc[k]))
      #nchar(QC0_NonA)
      QC0_NonAT <- gsub("T","",QC0_NonA)
      #nchar(QC0_NonAT)
      QC0_NonATC <- gsub("C","",QC0_NonAT)
      #nchar(QC_NonATC)
      QC0_NonATCG <- gsub("G","",QC0_NonATC)
      #nchar(QC0_NonATCG)
      if (nchar(QC0_NonATCG)==0){
        Output_02a_QC0_TargetSeqATCG <- append(Output_02a_QC0_TargetSeqATCG,"Passed QC0")
      } else {
        Output_02a_QC0_TargetSeqATCG <- append(Output_02a_QC0_TargetSeqATCG,"Failed QC0; contains non-ATCG characters")
      }
      #############################################################
      #Step 03a. PreAln: Combine Ref and Target into a DNAStringSet Object
      MyWorkingPreAln_Nuc <- DNAStringSet() 
      MyWorkingPreAln_Nuc <- append(MyWorkingPreAln_Nuc,MyRefNuc_HXB2[i])
      MyWorkingPreAln_Nuc <- append(MyWorkingPreAln_Nuc,MyTargetNuc[k])
      #print(class(MyWorkingPreAln_Nuc))
      #print(MyWorkingPreAln_Nuc) #OK
      #############################################################
      #Step 03b. Align them with muscle, save as appended
      MyWorkingAln_Nuc <- muscle(MyWorkingPreAln_Nuc,gapopen=MyGapPenaltyNUC,log=MuscleNUCLog,verbose=TRUE)
      #class(MyWorkingAln_Nuc) #"DNAMultipleAlignment"
      MyWorkingAln_Nuc_DNAStringSet <- DNAStringSet(MyWorkingAln_Nuc) #Convert "DNAMultipleAlignment" to "DNAStringSet" for fasta writing
      #print(MyWorkingAln_Nuc)#OK
      writeXStringSet(MyWorkingAln_Nuc_DNAStringSet,MyWorkingAln_Nuc_fasta,append=TRUE) #UNTRIMMED_NUC_ALIGNMENT(FULL_UNTRIMMED)
      Output_03a_MyWorkingAln_NucRef <- append(Output_03a_MyWorkingAln_NucRef,toString(MyWorkingAln_Nuc_DNAStringSet[1]))
      Output_03a1_MyWorkingAln_NucRef_Len <- append(Output_03a1_MyWorkingAln_NucRef_Len,width(MyWorkingAln_Nuc_DNAStringSet[1]))
      Output_03b_MyWorkingAln_NucTarget <- append(Output_03b_MyWorkingAln_NucTarget,toString(MyWorkingAln_Nuc_DNAStringSet[2]))
      Output_03a2_MyWorkingAln_NucTarget_Len <- append(Output_03a2_MyWorkingAln_NucTarget_Len,width(MyWorkingAln_Nuc_DNAStringSet[2]))
      ############################################################
      #Step 04a.  Index the alignment:  startIndex & endIndex of class ByPos_MIndex (Biostrings)
      ##20151203 vmatchPattern Indexing development has stopped because it doesn't work for some
      ##Therefore, ignore 04a and 04b output
      MyIndex <- vmatchPattern(toString(MyRefNuc_HXB2[i]),MyWorkingAln_Nuc_DNAStringSet[1],max.mismatch=0,min.mismatch=0,with.indels=TRUE)
      MyIndex_start <- as.numeric(unlist(startIndex(MyIndex)))
      #print(MyIndex_start) #if no match, function returns numeric(0)
      MyIndex_end <- as.numeric(unlist(endIndex(MyIndex)))
      #print(MyIndex_end) #if no match, function returns numeric(0)
      if (length(MyIndex_start)==0){
        MyIndex_start <- 0
      }
      if (length(MyIndex_end)==0){
        MyIndex_end <- 0
      }
      #20151129 is.numeric worked but T/F
      #20151129 the major problem is vmatchPattern is not finding a match eg. in A2 vs gag
      Output_04a_MyIndex_start <- append(Output_04a_MyIndex_start,MyIndex_start)
      Output_04b_MyIndex_end <- append(Output_04b_MyIndex_end,MyIndex_end)
      ############################################################
      #Step 04b.  Index the alignment: regular expression
      MyIndex_A <- unlist(gregexpr("A",MyWorkingAln_Nuc_DNAStringSet[1],ignore.case=TRUE))
      MyIndex_T <- unlist(gregexpr("T",MyWorkingAln_Nuc_DNAStringSet[1],ignore.case=TRUE))
      MyIndex_C <- unlist(gregexpr("C",MyWorkingAln_Nuc_DNAStringSet[1],ignore.case=TRUE))
      MyIndex_G <- unlist(gregexpr("G",MyWorkingAln_Nuc_DNAStringSet[1],ignore.case=TRUE))
      MyIndex_ATCG <- c(MyIndex_A,MyIndex_T,MyIndex_C,MyIndex_G) #combine into one vector
      #MyIndex_ATCG_sorted <- sort(MyIndex_ATCG) #sort alphabetically
      #if (MyIndex_A==-1|MyIndex_T==-1|MyIndex_C==-1|MyIndex_G==-1){
      #Output_04c_MyIndex_ATCG_sorted_start <- append(Output_04c_MyIndex_ATCG_sorted_start,0)
      #Output_04d_MyIndex_ATCG_sorted_end <- append(Output_04d_MyIndex_ATCG_sorted_end,0)	
      #} else {
      Output_04c_MyIndex_ATCG_start <- append(Output_04c_MyIndex_ATCG_start,min(MyIndex_ATCG))
      Output_04d_MyIndex_ATCG_end <- append(Output_04d_MyIndex_ATCG_end,max(MyIndex_ATCG))
      ############################################################
      #Step 05. QC %match in alignment including dashes
      MyWorkingAln_Nuc_DNAStringSet_trim <- DNAStringSet(MyWorkingAln_Nuc_DNAStringSet,start=min(MyIndex_ATCG),end=max(MyIndex_ATCG))
      writeXStringSet(MyWorkingAln_Nuc_DNAStringSet_trim,MyWorkingAln_Nuc_muscleAlign,append=TRUE) #TRIMMED_NUC_ALIGNED
      MyWorkingAln_Nuc_DNAStringSet_trim_1 <- unlist(strsplit(toString(MyWorkingAln_Nuc_DNAStringSet_trim[1]),"")) #strsplit is the same as seqinr 
      MyWorkingAln_Nuc_DNAStringSet_trim_2 <- unlist(strsplit(toString(MyWorkingAln_Nuc_DNAStringSet_trim[2]),"")) #Trim 1 is always HXB2; Trim 2 is my data
      #05a:  QC1. trim1 and trim2 must be identical in seq len
      if (length(MyWorkingAln_Nuc_DNAStringSet_trim_1)==length(MyWorkingAln_Nuc_DNAStringSet_trim_2)){
        Output_05a_QC1_lentrim1vs2 <- append(Output_05a_QC1_lentrim1vs2,"Passed QC1")
      } else {
        Output_05a_QC1_lentrim1vs2 <- append(Output_05a_QC1_lentrim1vs2,"Manual Check Required")
      }
      #05b:  QC2. %iden
      mismatch_count <- 0
      denominator <- length(MyWorkingAln_Nuc_DNAStringSet_trim_1)
      for (m in 1:length(MyWorkingAln_Nuc_DNAStringSet_trim_1)){
        if (MyWorkingAln_Nuc_DNAStringSet_trim_1[m]==MyWorkingAln_Nuc_DNAStringSet_trim_2[m]){
          mismatch_count <- mismatch_count + 0
        } else {
          mismatch_count <- mismatch_count + 1
        }
      }
      PercIdentity <- (denominator-mismatch_count)/denominator
      Output_05b_QC2a_PercIdentity <- append(Output_05b_QC2a_PercIdentity,PercIdentity)
      #05c. QC2b. %iden pass/fail
      if (PercIdentity>=PerIden_Cutoff_Nuc){
        Output_05c_QC2b_PerIden_Cutoff_Nuc <- append(Output_05c_QC2b_PerIden_Cutoff_Nuc,"Passed QC2")
      } else {
        Output_05c_QC2b_PerIden_Cutoff_Nuc <- append(Output_05c_QC2b_PerIden_Cutoff_Nuc,paste("Failed QC2; less than",PerIden_Cutoff_Nuc,"iden vs HXB2"))
      }
      ############################################################
      #Step 06. Trim start only; translate 'til stop codon
      ############################################################
      ############Exception 1:  exclude HXB2_FL_9719bp############
      if (names(MyRefNuc_HXB2[i])=="HXB2_K03455_FULL"){  
        Output_06a_MyFinal_AA <- append(Output_06a_MyFinal_AA,NA)
        Output_06b_MyFinal_AA_len <- append(Output_06b_MyFinal_AA_len,NA)
        Output_06b2_MyAA_Expected_len <- append(Output_06b2_MyAA_Expected_len,NA)
        Output_06c_MyFinal_Nuc <- append(Output_06c_MyFinal_Nuc,NA)
        Output_06d_MyFinal_Nuc_len <- append(Output_06d_MyFinal_Nuc_len,NA)
        Output_06e_MyFinal_Stop_4bp <- append(Output_06e_MyFinal_Stop_4bp,NA)
        Output_06f_QC3_AA_TooShort <- append(Output_06f_QC3_AA_TooShort,NA)
        Output_07a_AA_trimmed_HXB2 <- append(Output_07a_AA_trimmed_HXB2,NA)
        Output_07b_AA_trimmed_sample <- append(Output_07b_AA_trimmed_sample,NA)
        ###################################################################################################
        ###############Exception 2:  V3 has no stop codon, trim according to alignment#####################
      } else if (names(MyRefNuc_HXB2[i])=="HXB2_V3105_7110-7217"){  
        MyExpectedAALength <- 105/3
        MyWorking_V3NUC_char <- toString(MyWorkingAln_Nuc_DNAStringSet_trim[2])
        MyFinal_Nuc <- gsub("-","",MyWorking_V3NUC_char)
        MyFinal_AA_AAString <- translate(DNAString(MyFinal_Nuc))
        MyFinal_AA_character <- toString(MyFinal_AA_AAString)
        Output_06a_MyFinal_AA <- append(Output_06a_MyFinal_AA,MyFinal_AA_character)
        Output_06b_MyFinal_AA_len <- append(Output_06b_MyFinal_AA_len,nchar(MyFinal_AA_character))
        Output_06b2_MyAA_Expected_len <- append(Output_06b2_MyAA_Expected_len,MyExpectedAALength)
        Output_06c_MyFinal_Nuc <- append(Output_06c_MyFinal_Nuc,MyFinal_Nuc) #no dash
        Output_06d_MyFinal_Nuc_len <- append(Output_06d_MyFinal_Nuc_len,nchar(MyFinal_Nuc))
        Output_06e_MyFinal_Stop_4bp <- append(Output_06e_MyFinal_Stop_4bp,NA)
        MyActualAALength <- as.numeric(nchar(MyFinal_AA_character))
        if (MyActualAALength>=MyExpectedAALength*AA_TooLong){
          Output_06f_QC3_AA_TooShort <- append(Output_06f_QC3_AA_TooShort,"Failed QC3; AA_TooLong")
        } else if (MyActualAALength<=MyExpectedAALength*AA_TooShort){
          Output_06f_QC3_AA_TooShort <- append(Output_06f_QC3_AA_TooShort,"Failed QC3; AA_TooShort")
        } else if (MyActualAALength>MyExpectedAALength*AA_TooShort & MyActualAALength<MyExpectedAALength*AA_TooLong) {
          Output_06f_QC3_AA_TooShort <- append(Output_06f_QC3_AA_TooShort,"Passed QC3")
        } else {
          Output_06f_QC3_AA_TooShort <- append(Output_06f_QC3_AA_TooShort,"Failed QC3; Check")
        }
        #Step 06f. Export AA aligned with HXB2
        #(1) Convert my curent seq to AAStringSet
        #MyFinal_AA_character #My current AA as string, unaligned
        MyCurrentSeq_name <- names(MyTargetNuc[k]) #My current AA name
        names(MyFinal_AA_character) <- MyCurrentSeq_name
        MyFinal_AA_AAStringSet_named <- AAStringSet(MyFinal_AA_character)
        #(2) Load the HXB2 equivilent gene
        names(MyRefNuc_HXB2[i]) #My curernt gene: HXB2 region
        MyCurrentRefSeq_HXB2_AAStringSet <- MyRefAA_HXB2[names(MyRefNuc_HXB2[i])] #search for HXB2 RefSeq with the same name in the AARefSeq file
        #(3) Combine the two into the same AAStringSet file
        COMB_MyCurrentSeq_MyRefSeq_AAStringSet <- append(MyCurrentRefSeq_HXB2_AAStringSet,MyFinal_AA_AAStringSet_named)
        #(4) muscle alignment
        MyFinal_AA_Aligned <- muscle(COMB_MyCurrentSeq_MyRefSeq_AAStringSet,gapopen=MyGapPenaltyAA,log=MuscleAALog,verbose=TRUE)
        #(5) convert to AAStringSet
        MyFinal_AA_Aligned_AAStringSet <- AAStringSet(MyFinal_AA_Aligned)
        writeXStringSet(MyFinal_AA_Aligned_AAStringSet,MyWorkingAln_AA_muscleAlign,append=TRUE) #aligned_AA_oneb4stop
        #(6) since V3 has no stop, repeat AA alignment with Biostrings pairwiseAlignment() for output_07a and b
        COMB_MyCurrentSeq_MyRefSeq_AAStringSet[1]#HXB2 V3
        COMB_MyCurrentSeq_MyRefSeq_AAStringSet[2]#MySample V3
        MyWorkingV3_unaligned_AAStringSet_pairAlign <- pairwiseAlignment(COMB_MyCurrentSeq_MyRefSeq_AAStringSet[1],COMB_MyCurrentSeq_MyRefSeq_AAStringSet[2],type="global-local")
        MyFinal_AA_paired <- AAStringSet((pattern(MyWorkingV3_unaligned_AAStringSet_pairAlign)))
        MyFinal_AA_paired2 <- append(MyFinal_AA_paired,AAStringSet((subject(MyWorkingV3_unaligned_AAStringSet_pairAlign))))
        names(MyFinal_AA_paired2) <- c(names(COMB_MyCurrentSeq_MyRefSeq_AAStringSet[1]),names(COMB_MyCurrentSeq_MyRefSeq_AAStringSet[2]))
        writeXStringSet(MyFinal_AA_paired2,MyWorkingAln_AA_AAStringSet_trimmed,append=TRUE)
        #rm(MyFinal_AA_paired2)
        #(5) Save to df
        Output_07a_AA_trimmed_HXB2 <- append(Output_07a_AA_trimmed_HXB2,toString(pattern(MyWorkingV3_unaligned_AAStringSet_pairAlign)))
        Output_07b_AA_trimmed_sample <- append(Output_07b_AA_trimmed_sample,toString(subject(MyWorkingV3_unaligned_AAStringSet_pairAlign)))
        ##############################################################################################################
        ###############Exception 3: Tat exon 1##########################################################################
      } else if (names(MyRefNuc_HXB2[i])=="HXB2_tat-exon1_5831-6045") {   
        Temp_TATexon1_NUC_nodash <- gsub("-","",toString(MyWorkingAln_Nuc_DNAStringSet_trim[2]))
        names(Temp_TATexon1_NUC_nodash) <- names(MyWorkingAln_Nuc_DNAStringSet_trim[2])
        Store_TATexon1_NUC_nodash <- append(Store_TATexon1_NUC_nodash,Temp_TATexon1_NUC_nodash)
        #Store_TATexon1_NUC_nodash <- append(Store_TATexon1_NUC_nodash,gsub("-","",toString(MyWorkingAln_Nuc_DNAStringSet_trim[2])))
        #Store_TATexon1_NUC_nodash_names <- append(Store_TATexon1_NUC_nodash_names,names(MyWorkingAln_Nuc_DNAStringSet_trim[2]))
        Output_06a_MyFinal_AA <- append(Output_06a_MyFinal_AA,NA)
        Output_06b_MyFinal_AA_len <- append(Output_06b_MyFinal_AA_len,NA)
        Output_06b2_MyAA_Expected_len <- append(Output_06b2_MyAA_Expected_len,NA)
        Output_06c_MyFinal_Nuc <- append(Output_06c_MyFinal_Nuc,NA)
        Output_06d_MyFinal_Nuc_len <- append(Output_06d_MyFinal_Nuc_len,NA)
        Output_06e_MyFinal_Stop_4bp <- append(Output_06e_MyFinal_Stop_4bp,NA)
        Output_06f_QC3_AA_TooShort <- append(Output_06f_QC3_AA_TooShort,NA)
        Output_07a_AA_trimmed_HXB2 <- append(Output_07a_AA_trimmed_HXB2,NA)
        Output_07b_AA_trimmed_sample <- append(Output_07b_AA_trimmed_sample,NA)
        ##############################################################################################################
        ###############Exception 4: Tat exon 2##########################################################################
      } else if (names(MyRefNuc_HXB2[i])=="HXB2_tat-exon2_nuc_8379-8466") { 
        Temp_TATexon2GP41_NUC <- DNAStringSet(MyWorkingAln_Nuc_DNAStringSet,start=min(MyIndex_ATCG),end=max(MyIndex_ATCG)+60)
        Temp_TATexon2GP41_NUC_nodash <- gsub("-","",toString(Temp_TATexon2GP41_NUC[2]))
        names(Temp_TATexon2GP41_NUC_nodash) <- names(Temp_TATexon2GP41_NUC[2])
        Store_TATexon2GP41_NUC_nodash <- append(Store_TATexon2GP41_NUC_nodash,Temp_TATexon2GP41_NUC_nodash)
        #Store_TATexon2GP41_NUC_nodash <- append(Store_TATexon2GP41_NUC_nodash,gsub("-","",toString(Store_TATexon2GP41_NUC[1])))
        #Store_TATexon2GP41_NUC_nodash_names <- append(Store_TATexon2GP41_NUC_nodash_names,names(Store_TATexon2GP41_NUC[1]))
        Output_06a_MyFinal_AA <- append(Output_06a_MyFinal_AA,NA)
        Output_06b_MyFinal_AA_len <- append(Output_06b_MyFinal_AA_len,NA)
        Output_06b2_MyAA_Expected_len <- append(Output_06b2_MyAA_Expected_len,NA)
        Output_06c_MyFinal_Nuc <- append(Output_06c_MyFinal_Nuc,NA)
        Output_06d_MyFinal_Nuc_len <- append(Output_06d_MyFinal_Nuc_len,NA)
        Output_06e_MyFinal_Stop_4bp <- append(Output_06e_MyFinal_Stop_4bp,NA)
        Output_06f_QC3_AA_TooShort <- append(Output_06f_QC3_AA_TooShort,NA)
        Output_07a_AA_trimmed_HXB2 <- append(Output_07a_AA_trimmed_HXB2,NA)
        Output_07b_AA_trimmed_sample <- append(Output_07b_AA_trimmed_sample,NA)
        ##############################################################################################################
        ###############Exception 5: Rev exon 1##########################################################################
      } else if (names(MyRefNuc_HXB2[i])=="HXB2_rev-exon1_5970-6045") {  
        Temp_REVexon1_NUC_nodash <- gsub("-","",toString(MyWorkingAln_Nuc_DNAStringSet_trim[2]))
        names(Temp_REVexon1_NUC_nodash) <- names(MyWorkingAln_Nuc_DNAStringSet_trim[2])
        Store_REVexon1_NUC_nodash <- append(Store_REVexon1_NUC_nodash,Temp_REVexon1_NUC_nodash)
        #Store_REVexon1_NUC_nodash <- append(Store_REVexon1_NUC_nodash,gsub("-","",toString(MyWorkingAln_Nuc_DNAStringSet_trim[2])))
        #Store_REVexon1_NUC_nodash_names <- append(Store_REVexon1_NUC_nodash_names,names(MyWorkingAln_Nuc_DNAStringSet_trim[2]))
        Output_06a_MyFinal_AA <- append(Output_06a_MyFinal_AA,NA)
        Output_06b_MyFinal_AA_len <- append(Output_06b_MyFinal_AA_len,NA)
        Output_06b2_MyAA_Expected_len <- append(Output_06b2_MyAA_Expected_len,NA)
        Output_06c_MyFinal_Nuc <- append(Output_06c_MyFinal_Nuc,NA)
        Output_06d_MyFinal_Nuc_len <- append(Output_06d_MyFinal_Nuc_len,NA)
        Output_06e_MyFinal_Stop_4bp <- append(Output_06e_MyFinal_Stop_4bp,NA)
        Output_06f_QC3_AA_TooShort <- append(Output_06f_QC3_AA_TooShort,NA)
        Output_07a_AA_trimmed_HXB2 <- append(Output_07a_AA_trimmed_HXB2,NA)
        Output_07b_AA_trimmed_sample <- append(Output_07b_AA_trimmed_sample,NA)
        ##############################################################################################################
        ###############Exception 6: Rev exon 2##########################################################################
      } else if (names(MyRefNuc_HXB2[i])=="HXB2_rev-exon2_nuc_8379-8650") {  #Exception: Rev Exon 2
        Temp_REVexon2GP41_NUC <- DNAStringSet(MyWorkingAln_Nuc_DNAStringSet,start=min(MyIndex_ATCG),end=max(MyIndex_ATCG)+60)
        Temp_REVexon2GP41_NUC_nodash <- gsub("-","",toString(Temp_REVexon2GP41_NUC[2]))
        names(Temp_REVexon2GP41_NUC_nodash) <- names(Temp_REVexon2GP41_NUC[2])
        Store_REVexon2GP41_NUC_nodash <- append(Store_REVexon2GP41_NUC_nodash,Temp_REVexon2GP41_NUC_nodash)
        #Store_REVexon2GP41_NUC <- DNAStringSet(MyWorkingAln_Nuc_DNAStringSet,start=min(MyIndex_ATCG)+1,end=max(MyIndex_ATCG)+145)[2]
        #Store_REVexon2GP41_NUC_nodash <- gsub("-","",toString(Store_REVexon2GP41_NUC[1]))
        #Store_REVexon2GP41_NUC_nodash_names <- append(Store_REVexon2GP41_NUC_nodash_names,names(Store_REVexon2GP41_NUC[1]))
        Output_06a_MyFinal_AA <- append(Output_06a_MyFinal_AA,NA)
        Output_06b_MyFinal_AA_len <- append(Output_06b_MyFinal_AA_len,NA)
        Output_06b2_MyAA_Expected_len <- append(Output_06b2_MyAA_Expected_len,NA)
        Output_06c_MyFinal_Nuc <- append(Output_06c_MyFinal_Nuc,NA)
        Output_06d_MyFinal_Nuc_len <- append(Output_06d_MyFinal_Nuc_len,NA)
        Output_06e_MyFinal_Stop_4bp <- append(Output_06e_MyFinal_Stop_4bp,NA)
        Output_06f_QC3_AA_TooShort <- append(Output_06f_QC3_AA_TooShort,NA)
        Output_07a_AA_trimmed_HXB2 <- append(Output_07a_AA_trimmed_HXB2,NA)
        Output_07b_AA_trimmed_sample <- append(Output_07b_AA_trimmed_sample,NA)
        ##############################################################################################################
        ############### For the rest of the genes #######################################################################
      } else { 
        MyWorkingAln_Nuc_DNAStringSet_trim_begin <- DNAStringSet(MyWorkingAln_Nuc_DNAStringSet,start=min(MyIndex_ATCG))
        MyWorkingAln_Nuc_DNAStringSet_trim_begin_nodash <- gsub("-","",toString(MyWorkingAln_Nuc_DNAStringSet_trim_begin[2]))
        MyWorking_AA_AAString <- translate(DNAString(MyWorkingAln_Nuc_DNAStringSet_trim_begin_nodash))
        #Step 06a. Index first stop
        MyIndex_FirstStop <- unlist(gregexpr("\\*",MyWorking_AA_AAString,ignore.case=TRUE)) #class integer
        #Step 06b. Store the translation & len
        MyFinal_AA_AAString <- MyWorking_AA_AAString[1:MyIndex_FirstStop[1]-1] 
        MyFinal_AA_character <- toString(MyFinal_AA_AAString)
        #Step 06c. Store NUC corresponding to AA, length
        MyFinal_Nuc <- substr(MyWorkingAln_Nuc_DNAStringSet_trim_begin_nodash,1,nchar(MyFinal_AA_character)*3)
        if (nchar(MyFinal_AA_AAString)==0){
          Output_06c_MyFinal_Nuc <- append(Output_06c_MyFinal_Nuc,NA)
          Output_06d_MyFinal_Nuc_len <- append(Output_06d_MyFinal_Nuc_len,0)		    
        } else {
          Output_06c_MyFinal_Nuc <- append(Output_06c_MyFinal_Nuc,MyFinal_Nuc)
          Output_06d_MyFinal_Nuc_len <- append(Output_06d_MyFinal_Nuc_len,nchar(MyFinal_Nuc))
        }
        #Step 06d. Store the 4bp STOP codon
        MyFinal_Stop_4bp <- substr(MyWorkingAln_Nuc_DNAStringSet_trim_begin_nodash,nchar(MyFinal_AA_AAString)*3+1,nchar(MyFinal_AA_AAString)*3+4)
        Output_06e_MyFinal_Stop_4bp <- append(Output_06e_MyFinal_Stop_4bp,MyFinal_Stop_4bp)
        #Step 06e. QC3 Any seq 50% (?) shorter than expected would be flagged
        MyExpectedAALength <- width(MyRefAA_HXB2[names(MyRefNuc_HXB2[i]),])
        Output_06b2_MyAA_Expected_len <- append(Output_06b2_MyAA_Expected_len,MyExpectedAALength)
        MyActualAALength <- as.numeric(nchar(MyFinal_AA_character))
        #Step 06f. Export AA aligned oneb4stop
        if (nchar(MyFinal_AA_AAString)==0){
          #do not append to fasta
          Output_07a_AA_trimmed_HXB2 <- append(Output_07a_AA_trimmed_HXB2,"NA")
          Output_07b_AA_trimmed_sample <- append(Output_07b_AA_trimmed_sample,"NA")
        } else {
          #(1) Convert my curent seq to AAStringSet
          #MyFinal_AA_character #My current AA as string, unaligned
          MyCurrentSeq_name <- names(MyTargetNuc[k]) #My current AA name
          names(MyFinal_AA_character) <- MyCurrentSeq_name
          MyFinal_AA_AAStringSet_named <- AAStringSet(MyFinal_AA_character)
          #(2) Load the HXB2 equivilent gene
          names(MyRefNuc_HXB2[i]) #My curernt gene: HXB2 region
          MyCurrentRefSeq_HXB2_AAStringSet <- MyRefAA_HXB2[names(MyRefNuc_HXB2[i])] #search for HXB2 RefSeq with the same name in the AARefSeq file
          #(3) Combine the two into the same AAStringSet file
          COMB_MyCurrentSeq_MyRefSeq_AAStringSet <- append(MyCurrentRefSeq_HXB2_AAStringSet,MyFinal_AA_AAStringSet_named)
          #(4) muscle alignment
          MyFinal_AA_Aligned <- muscle(COMB_MyCurrentSeq_MyRefSeq_AAStringSet,gapopen=MyGapPenaltyAA,log=MuscleAALog,verbose=TRUE)
          #(5) convert to AAStringSet
          MyFinal_AA_Aligned_AAStringSet <- AAStringSet(MyFinal_AA_Aligned)
          writeXStringSet(MyFinal_AA_Aligned_AAStringSet,MyWorkingAln_AA_muscleAlign,append=TRUE) #aligned_AA_oneb4stop
          #Step 06g. Export AA sligned to HXB2 using pairwiseAlignment()
          #(1) Retrieve unaligned NUC from both HXB2 and sample, combine into one DNAStringSet object, name them
          class(MyWorkingAln_Nuc_DNAStringSet_trim_begin_nodash) #My Sample NUC; class=char
          class(toString(MyRefNuc_HXB2[i])) #HXB2 NUC; class=char
          MyWorking_unaligned_DNAStringSet <- append(DNAStringSet(toString(MyRefNuc_HXB2[i])),DNAStringSet(MyWorkingAln_Nuc_DNAStringSet_trim_begin_nodash))
          names(MyWorking_unaligned_DNAStringSet) <- c(names(MyRefNuc_HXB2[i]),names(MyTargetNuc[k]))
          #(2) Translate both nodash FASTA
          MyWorking_unaligned_AAStringSet <- translate(MyWorking_unaligned_DNAStringSet)
          #(3) Align the AA
          MyWorking_aligned_AAStringSet_pairAlign <- pairwiseAlignment(MyWorking_unaligned_AAStringSet[1],MyWorking_unaligned_AAStringSet[2],type="global-local")
          #(4) Convert back to AAStringSet
          MyFinal_AA_paired <- AAStringSet(pattern(MyWorking_aligned_AAStringSet_pairAlign))
          MyFinal_AA_paired2 <- append(MyFinal_AA_paired,AAStringSet(subject(MyWorking_aligned_AAStringSet_pairAlign)))
          names(MyFinal_AA_paired2) <- c(names(MyRefNuc_HXB2[i]),names(MyTargetNuc[k]))
          writeXStringSet(MyFinal_AA_paired2,MyWorkingAln_AA_AAStringSet_trimmed,append=TRUE)
          #(5) Save to df
          Output_07a_AA_trimmed_HXB2 <- append(Output_07a_AA_trimmed_HXB2,toString(pattern(MyWorking_aligned_AAStringSet_pairAlign)))
          Output_07b_AA_trimmed_sample <- append(Output_07b_AA_trimmed_sample,toString(subject(MyWorking_aligned_AAStringSet_pairAlign)))
        }
        #20171204 Moved here because polyprotein AA lengths should be taken from trimmed
        if (nchar(MyFinal_AA_AAString)==0){
          Output_06a_MyFinal_AA <- append(Output_06a_MyFinal_AA,NA)
          Output_06b_MyFinal_AA_len <- append(Output_06b_MyFinal_AA_len,0) 	
        } else {
          Output_06a_MyFinal_AA <- append(Output_06a_MyFinal_AA,MyFinal_AA_character)
          if (names(MyRefNuc_HXB2[i])=="HXB2_PR_2253-2549" | names(MyRefNuc_HXB2[i])=="HXB2_RT_2550-3869" | names(MyRefNuc_HXB2[i])=="HXB2_RNaseH_3870-4229" | names(MyRefNuc_HXB2[i])=="HXB2_PR-RT_2253-3869" | names(MyRefNuc_HXB2[i])=="HXB2_PR-RNaseH_2253-4229" | names(MyRefNuc_HXB2[i])=="HXB2_PR-INT_2253-5093" | names(MyRefNuc_HXB2[i])=="HXB2_RT-INT_3720-5093" | names(MyRefNuc_HXB2[i])=="HXB2_gp120_6315-7757"){
            MyAALength_TRIMMED <- nchar(gsub("-","",toString(subject(MyWorking_aligned_AAStringSet_pairAlign))))
            Output_06b_MyFinal_AA_len <- append(Output_06b_MyFinal_AA_len,MyAALength_TRIMMED)
          } else {
            Output_06b_MyFinal_AA_len <- append(Output_06b_MyFinal_AA_len,nchar(MyFinal_AA_character)) 
          }
        }  
        ########## QC 3 ###############
        ### for pol polyproteins, fail QC3 if QC2 fail because large deletions within the polyprotein doesn't mean there is a stop codon (eg S0527) ###
        if (names(MyRefNuc_HXB2[i])=="HXB2_PR_2253-2549" | names(MyRefNuc_HXB2[i])=="HXB2_RT_2550-3869" | names(MyRefNuc_HXB2[i])=="HXB2_RNaseH_3870-4229" | names(MyRefNuc_HXB2[i])=="HXB2_PR-RT_2253-3869" | names(MyRefNuc_HXB2[i])=="HXB2_PR-RNaseH_2253-4229" | names(MyRefNuc_HXB2[i])=="HXB2_PR-INT_2253-5093" | names(MyRefNuc_HXB2[i])=="HXB2_RT-INT_3720-5093" | names(MyRefNuc_HXB2[i])=="HXB2_gp120_6315-7757"){
          if (PercIdentity>=PerIden_Cutoff_Nuc){ #if it passed QC2, normal route
            if (MyAALength_TRIMMED>MyExpectedAALength*AA_TooShort & MyAALength_TRIMMED<MyExpectedAALength*AA_TooLong){
              Output_06f_QC3_AA_TooShort <- append(Output_06f_QC3_AA_TooShort,"Passed QC3")
            } else if (MyAALength_TRIMMED<=MyExpectedAALength*AA_TooShort) {
              Output_06f_QC3_AA_TooShort <- append(Output_06f_QC3_AA_TooShort,"Failed QC3; AA_TooShort")
            } else if (MyAALength_TRIMMED>=MyExpectedAALength*AA_TooLong) {
              Output_06f_QC3_AA_TooShort <- append(Output_06f_QC3_AA_TooShort,"Failed QC3; AA_TooLong")
            } else {
              Output_06f_QC3_AA_TooShort <- append(Output_06f_QC3_AA_TooShort,"Failed QC3; Check")
            }
          } else { #if it didn't pass QC2, fail QC3
            Output_06f_QC3_AA_TooShort <- append(Output_06f_QC3_AA_TooShort,"Failed QC3; polyprotein failing QC2")
          }
          ### adjust for 95% of gag p6:  p6=53aa; 95%=3aa; gag500aa-3aa=497aa; 497/500=0.994 ###
        } else if (names(MyRefNuc_HXB2[i])=="HXB2_gag_790-2289"){
          if (MyActualAALength>MyExpectedAALength*AA_TooShort_Gag & MyActualAALength<MyExpectedAALength*AA_TooLong){
            Output_06f_QC3_AA_TooShort <- append(Output_06f_QC3_AA_TooShort,"Passed QC3")
          } else if (MyActualAALength<=MyExpectedAALength*AA_TooShort_Gag) {
            Output_06f_QC3_AA_TooShort <- append(Output_06f_QC3_AA_TooShort,"Failed QC3; AA_TooShort")
          } else if (MyActualAALength>=MyExpectedAALength*AA_TooLong) {
            Output_06f_QC3_AA_TooShort <- append(Output_06f_QC3_AA_TooShort,"Failed QC3; AA_TooLong")
          } else {
            Output_06f_QC3_AA_TooShort <- append(Output_06f_QC3_AA_TooShort,"Failed QC3; Check")
          }
          ## Not pol's polyproteins nor gag nor gp120; including pol##
        } else {  
          if (MyActualAALength>MyExpectedAALength*AA_TooShort & MyActualAALength<MyExpectedAALength*AA_TooLong){
            Output_06f_QC3_AA_TooShort <- append(Output_06f_QC3_AA_TooShort,"Passed QC3")
          } else if (MyActualAALength<=MyExpectedAALength*AA_TooShort) {
            Output_06f_QC3_AA_TooShort <- append(Output_06f_QC3_AA_TooShort,"Failed QC3; AA_TooShort")
          } else if (MyActualAALength>=MyExpectedAALength*AA_TooLong) {
            Output_06f_QC3_AA_TooShort <- append(Output_06f_QC3_AA_TooShort,"Failed QC3; AA_TooLong")
          } else {
            Output_06f_QC3_AA_TooShort <- append(Output_06f_QC3_AA_TooShort,"Failed QC3; Check")
          }
        }
      }
    }
  }
  
  ############## TAT and REV exon1-2 slab################This is done after all the regions*samples were processed
  #Step 07.  Preparation & Review
  #Two files were loaded for this section:
  #1. MyRefNUC_TATREVslab_HXB2 #class(DNAStringSet)
  #2. MyRefAA_TATREVslab_HXB2 #class(AAStringSet)
  #Four vectors were exported from Step 6:
  #1. Store_TATexon1_NUC_nodash #class char
  #2. Store_TATexon2GP41_NUC_nodash #class char
  #3. Store_REVexon1_NUC_nodash #class char
  #4. Store_REVexon2GP41_NUC_nodash #class char
  HXB2_TAT_SLAB_NUC <- "ATGGAGCCAGTAGATCCTAGACTAGAGCCCTGGAAGCATCCAGGAAGTCAGCCTAAAACTGCTTGTACCAATTGCTATTGTAAAAAGTGTTGCTTTCATTGCCAAGTTTGTTTCATAACAAAAGCCTTAGGCATCTCCTATGGCAGGAAGAAGCGGAGACAGCGACGAAGAGCTCATCAGAACAGTCAGACTCATCAAGCTTCTCTATCAAAGCAACCCACCTCCCAACCCCGAGGGGACCCGACAGGCCCGAAGGAATAGAAGAAGAAGGTGGAGAGAGAGACAGAGACAGATCCATTCGAT"
  HXB2_REV_SLAB_NUC <- "ATGGCAGGAAGAAGCGGAGACAGCGACGAAGAGCTCATCAGAACAGTCAGACTCATCAAGCTTCTCTATCAAAGCAACCCACCTCCCAACCCCGAGGGGACCCGACAGGCCCGAAGGAATAGAAGAAGAAGGTGGAGAGAGAGACAGAGACAGATCCATTCGATTAGTGAACGGATCCTTGGCACTTATCTGGGACGATCTGCGGAGCCTGTGCCTCTTCAGCTACCACCGCTTGAGAGACTTACTCTTGATTGTAACGAGGATTGTGGAACTTCTGGGACGCAGGGGGTGGGAAGCCCTCAAATATTGGTGGAATCTCCTACAGTATTGGAGTCAGGAACTAAAGAA"
  HXB2_TAT_SLAB_AA <- "MEPVDPRLEPWKHPGSQPKTACTNCYCKKCCFHCQVCFITKALGISYGRKKRRQRRRAHQNSQTHQASLSKQPTSQPRGDPTGPKE*KKKVERETETDPFD"
  HXB2_REV_SLAB_AA <- "MAGRSGDSDEELIRTVRLIKLLYQSNPPPNPEGTRQARRNRRRRWRERQRQIHSISERILGTYLGRSAEPVPLQLPPLERLTLDCNEDCGTSGTQGVGSPQILVESPTVLESGTKE"
  #########################################################
  #Step 07. Initial QC
  if (
    length(Store_TATexon1_NUC_nodash)==
    length(Store_TATexon2GP41_NUC_nodash)&
    length(Store_TATexon2GP41_NUC_nodash)==
    length(Store_REVexon1_NUC_nodash)&
    length(Store_REVexon1_NUC_nodash)==
    length(Store_REVexon2GP41_NUC_nodash)
  ){
    #Step 08.  Slab TAT
    for (n in 1:length(Store_TATexon1_NUC_nodash)){
      #Step 08a.  Log information
      Output_01_MyWorkingRefNuc_Name <- append(Output_01_MyWorkingRefNuc_Name,"HXB2_TAT_slab")
      Output_01a_MyWorkingRefNuc_UnAlign_Seq <- append(Output_01a_MyWorkingRefNuc_UnAlign_Seq,HXB2_TAT_SLAB_NUC)
      Output_02_MyWorkingPreAlnNuc_Name <- append(Output_02_MyWorkingPreAlnNuc_Name,names(Store_TATexon1_NUC_nodash[n]))
      Output_02a_QC0_TargetSeqATCG <- append(Output_02a_QC0_TargetSeqATCG,NA)
      Output_03a_MyWorkingAln_NucRef <- append(Output_03a_MyWorkingAln_NucRef,NA)
      Output_03b_MyWorkingAln_NucTarget <- append(Output_03b_MyWorkingAln_NucTarget,NA)
      Output_03a1_MyWorkingAln_NucRef_Len <- append(Output_03a1_MyWorkingAln_NucRef_Len,NA)
      Output_03a2_MyWorkingAln_NucTarget_Len <- append(Output_03a2_MyWorkingAln_NucTarget_Len,NA)
      Output_04a_MyIndex_start <- append(Output_04a_MyIndex_start,NA)
      Output_04b_MyIndex_end <- append(Output_04b_MyIndex_end,NA)
      Output_04c_MyIndex_ATCG_start <- append(Output_04c_MyIndex_ATCG_start,NA)
      Output_04d_MyIndex_ATCG_end <- append(Output_04d_MyIndex_ATCG_end,NA)
      Output_05a_QC1_lentrim1vs2 <- append(Output_05a_QC1_lentrim1vs2,NA)
      Output_05b_QC2a_PercIdentity <- append(Output_05b_QC2a_PercIdentity,NA)
      Output_05c_QC2b_PerIden_Cutoff_Nuc <- append(Output_05c_QC2b_PerIden_Cutoff_Nuc,NA)
      #Step 08b. Slab exon1 & 2 as DNAStringSet
      TAT_slab_NUC <- paste(Store_TATexon1_NUC_nodash[n],Store_TATexon2GP41_NUC_nodash[n],sep="") #ACH-1 (A1)
      names(TAT_slab_NUC) <- names(Store_TATexon1_NUC_nodash[n])
      TAT_slab_NUC_DNAStringSet <- DNAStringSet(TAT_slab_NUC) #This is 30 nucleotides longer than the usual tat end
      #     #Step 08c. Grab HXB2 as DNAStringSet 
      #     DNAStringSet(HXB2_TAT_SLAB_NUC)
      #     #Step 08d. Combine HXB2-TAT and Slabbed-TAT into one DNAStringSet
      #     MyWorkingTATslab_DNAStringSet <- append(DNAStringSet(HXB2_TAT_SLAB_NUC),TAT_slab_NUC_DNAStringSet)
      #     #Step 08e. Muscle align them
      #     MyWorkingTATslab_muscle <- muscle(MyWorkingTATslab_DNAStringSet,gapopen=MyGapPenaltyNUC,log=MuscleNUCLog,verbose=TRUE)
      #     #Step 08f. Convert this back to DNAStringSet
      #     MyWorkingTATslab_muscle_DNAStringSet <- DNAStringSet(MyWorkingTATslab_muscle)
      #     #Step 08g. Aim to remove gaps, but in DNAStringSet format
      #     length(MyWorkingTATslab_DNAStringSet) #2
      #     nchar(gsub("-","",toString(MyWorkingTATslab_DNAStringSet[1]))) #HXB2
      #     nchar(gsub("-","",toString(MyWorkingTATslab_DNAStringSet[2])))
      #Step 08c. Translate
      TAT_slab_AA_AAStringSet <- translate(TAT_slab_NUC_DNAStringSet)
      #TAT_slab_AA_AAStringSet #OK, MEPV
      #nchar(TAT_slab_AA_AAStringSet) #101
      if (nchar(TAT_slab_AA_AAStringSet)>0){
        #Step 08d. Index first stop
        MyIndex_FirstStop <- unlist(gregexpr("\\*",TAT_slab_AA_AAStringSet,ignore.case=TRUE)) #class integer
        #what if there are no stops
        if (MyIndex_FirstStop[1]>=1) {
          MyIndex_FirstStop <- MyIndex_FirstStop[1]
        } else if (MyIndex_FirstStop[1]<1) {
          MyIndex_FirstStop <- nchar(TAT_slab_AA_AAStringSet)
        } else {
          MyIndex_FirstStop <- nchar(TAT_slab_AA_AAStringSet)
        }
        #MyIndex_FirstStop #aa87 for ACH-2, correct
        #Step 08e. Store the translation & len
        TAT_slab_AA_AAString <- toString(TAT_slab_AA_AAStringSet[1]) #convert to AAString
        MyFinal_AA_character <- substr(TAT_slab_AA_AAString,1,MyIndex_FirstStop[1]-1)
        #MyFinal_AA_character #correct
        #nchar(MyFinal_AA_character) #correct 86 for ACH-2, premature stop
        Output_06a_MyFinal_AA <- append(Output_06a_MyFinal_AA,MyFinal_AA_character)
        Output_06b_MyFinal_AA_len <- append(Output_06b_MyFinal_AA_len,nchar(MyFinal_AA_character))
        #Step 08f. Store NUC corresponding to AA, length
        MyFinal_Nuc <- substr(TAT_slab_NUC,1,nchar(MyFinal_AA_character)*3)
        #MyFinal_Nuc #Correct
        #nchar(MyFinal_Nuc) #correct 258
        Output_06c_MyFinal_Nuc <- append(Output_06c_MyFinal_Nuc,MyFinal_Nuc)
        Output_06d_MyFinal_Nuc_len <- append(Output_06d_MyFinal_Nuc_len,nchar(MyFinal_Nuc))
        #Step 08g. Store the 4bp STOP codon
        MyFinal_Stop_4bp <- substr(TAT_slab_NUC,nchar(MyFinal_AA_character)*3+1,nchar(MyFinal_AA_character)*3+4)
        #MyFinal_Stop_4bp #wrong; TAGA, expecting TAGT AHHHHHHHHHHh NO... this is ACH2's premature stop
        Output_06e_MyFinal_Stop_4bp <- append(Output_06e_MyFinal_Stop_4bp,MyFinal_Stop_4bp)
        #Step 08h. QC3 Any seq 50% (?) shorter than expected would be flagged
        MyExpectedAALength <- width(MyRefAA_HXB2[18])
        Output_06b2_MyAA_Expected_len <- append(Output_06b2_MyAA_Expected_len,MyExpectedAALength)
        MyActualAALength <- as.numeric(nchar(MyFinal_AA_character))
        if (MyActualAALength>MyExpectedAALength*AA_TooShort & MyActualAALength<MyExpectedAALength*AA_TooLong){
          Output_06f_QC3_AA_TooShort <- append(Output_06f_QC3_AA_TooShort,"Passed QC3")
        } else if (MyActualAALength<=MyExpectedAALength*AA_TooShort) {
          Output_06f_QC3_AA_TooShort <- append(Output_06f_QC3_AA_TooShort,"Failed QC3; AA_TooShort")
        } else if (MyActualAALength>=MyExpectedAALength*AA_TooLong) {
          Output_06f_QC3_AA_TooShort <- append(Output_06f_QC3_AA_TooShort,"Failed QC3; AA_TooLong")
        } else {
          Output_06f_QC3_AA_TooShort <- append(Output_06f_QC3_AA_TooShort,"Failed QC3; Check")
        }
        #Step 08i. Export AA aligned
        TAT_AA_w_HXB2_char <- c(HXB2_TAT_SLAB_AA,MyFinal_AA_character)
        names(TAT_AA_w_HXB2_char) <- c("HXB2_TAT_SLAB_AA",names(Store_TATexon1_NUC_nodash[n]))
        TAT_AA_w_HXB2_AAStringSet <- AAStringSet(TAT_AA_w_HXB2_char)
        TAT_AA_w_HXB2_muscle_Aln <- muscle(TAT_AA_w_HXB2_AAStringSet,gapopen=MyGapPenaltyAA,log=MuscleAALog,verbose=TRUE) #class AAMultipleAlignment
        TAT_AA_w_HXB2_muscle_Aln_AAStringSet <- AAStringSet(TAT_AA_w_HXB2_muscle_Aln)
        writeXStringSet(TAT_AA_w_HXB2_muscle_Aln_AAStringSet,MyWorkingAln_AA_muscleAlign,append=TRUE)
        #Step 08j. Export NUC aligned trimmed
        #(1) Grab HXB2 and mycurrent, combine into one DNAStringSet
        HXB2_TAT_SLAB_NUC #character
        names(HXB2_TAT_SLAB_NUC) <- "HXB2_TAT_SLAB_NUC"
        nchar(HXB2_TAT_SLAB_NUC) #303,correct
        HXB2_TAT_SLAB_NUC_DNAStringSet <- DNAStringSet(HXB2_TAT_SLAB_NUC)
        HXB2_TAT_SLAB_NUC_DNAStringSet_2 <- append(HXB2_TAT_SLAB_NUC_DNAStringSet,TAT_slab_NUC_DNAStringSet)
        #(2) Write both NUC and AA as fasta
        writeXStringSet(HXB2_TAT_SLAB_NUC_DNAStringSet_2,MyWorkingAln_Nuc_muscleAlign,append=TRUE) #aligned_NUC_trimmed
        #####################################################################################
        #Step 08k. Export AA aligned
        #(1) Grab unaligned NUC
        DNAStringSet(HXB2_TAT_SLAB_NUC)
        MyWorkingTATslab_unaligned_DNAStringSet <- append(DNAStringSet(HXB2_TAT_SLAB_NUC),TAT_slab_NUC_DNAStringSet)
        #(2) Translate both nodash fasta
        MyWorkingTATslab_unaligned_AAStringSet <- translate(MyWorkingTATslab_unaligned_DNAStringSet)
        #(3) Align the AAs
        MyWorkingTATslab_unaligned_AAStringSet_muscle <- muscle(MyWorkingTATslab_unaligned_AAStringSet,gapopen=MyGapPenaltyNUC,log=MuscleNUCLog,verbose=TRUE)
        ##########use Biostrings function pairwiseAlignment(), default settings##########################
        MyWorkingTATslab_unaligned_AAStringSet_pairAlign <- pairwiseAlignment(MyWorkingTATslab_unaligned_AAStringSet[1],MyWorkingTATslab_unaligned_AAStringSet[2],type="global-local")
        MyWorkingTATslab_unaligned_AAStringSet_pairAlign #good!  stop codons preserved
        class(MyWorkingTATslab_unaligned_AAStringSet_pairAlign)  #[1] "PairwiseAlignmentsSingleSubject"
        summary(MyWorkingTATslab_unaligned_AAStringSet_pairAlign)
        #writePairwiseAlignments(MyWorkingTATslab_unaligned_AAStringSet_pairAlign,file="test_TATpair.txt")
        #pairwiseAlignment() trim to HXB2 for me automatically...
        #(4) Convert back to AAStringSet
        #AAStringSet((pattern(MyWorkingTATslab_unaligned_AAStringSet_pairAlign)))
        #AAStringSet((subject(MyWorkingTATslab_unaligned_AAStringSet_pairAlign)))
        MyFinal_AA_paired <- AAStringSet((pattern(MyWorkingTATslab_unaligned_AAStringSet_pairAlign)))
        MyFinal_AA_paired2 <- append(MyFinal_AA_paired,AAStringSet((subject(MyWorkingTATslab_unaligned_AAStringSet_pairAlign))))
        names(MyFinal_AA_paired2) <- c("HXB2_TAT_SLAB_AA",names(Store_TATexon1_NUC_nodash[n]))
        writeXStringSet(MyFinal_AA_paired2,MyWorkingAln_AA_AAStringSet_trimmed,append=TRUE)
        #rm(MyFinal_AA_paired2)
        #(5) Save to df
        Output_07a_AA_trimmed_HXB2 <- append(Output_07a_AA_trimmed_HXB2,toString(pattern(MyWorkingTATslab_unaligned_AAStringSet_pairAlign)))
        Output_07b_AA_trimmed_sample <- append(Output_07b_AA_trimmed_sample,toString(subject(MyWorkingTATslab_unaligned_AAStringSet_pairAlign)))
        ###########################################
      } else {
        Output_06a_MyFinal_AA <- append(Output_06a_MyFinal_AA,NA)
        Output_06b_MyFinal_AA_len <- append(Output_06b_MyFinal_AA_len,0)
        Output_06c_MyFinal_Nuc <- append(Output_06c_MyFinal_Nuc,NA)
        Output_06d_MyFinal_Nuc_len <- append(Output_06d_MyFinal_Nuc_len,0)
        Output_06e_MyFinal_Stop_4bp <- append(Output_06e_MyFinal_Stop_4bp,NA)
        MyExpectedAALength <- width(MyRefAA_HXB2[18])
        Output_06b2_MyAA_Expected_len <- append(Output_06b2_MyAA_Expected_len,MyExpectedAALength)
        Output_06f_QC3_AA_TooShort <- append(Output_06f_QC3_AA_TooShort,"Failed QC3; AA_TooShort")
        Output_07a_AA_trimmed_HXB2 <- append(Output_07a_AA_trimmed_HXB2,NA)
        Output_07b_AA_trimmed_sample <- append(Output_07b_AA_trimmed_sample,NA)
      }
    }
    #Step 09.  Slab REV
    for (o in 1:length(Store_REVexon1_NUC_nodash)){
      #Step 09a.  Log information
      Output_01_MyWorkingRefNuc_Name <- append(Output_01_MyWorkingRefNuc_Name,"HXB2_REV_slab")
      Output_01a_MyWorkingRefNuc_UnAlign_Seq <- append(Output_01a_MyWorkingRefNuc_UnAlign_Seq,HXB2_REV_SLAB_NUC)
      Output_02_MyWorkingPreAlnNuc_Name <- append(Output_02_MyWorkingPreAlnNuc_Name,names(Store_REVexon1_NUC_nodash[o]))
      Output_02a_QC0_TargetSeqATCG <- append(Output_02a_QC0_TargetSeqATCG,NA)
      Output_03a_MyWorkingAln_NucRef <- append(Output_03a_MyWorkingAln_NucRef,NA)
      Output_03b_MyWorkingAln_NucTarget <- append(Output_03b_MyWorkingAln_NucTarget,NA)
      Output_03a1_MyWorkingAln_NucRef_Len <- append(Output_03a1_MyWorkingAln_NucRef_Len,NA)
      Output_03a2_MyWorkingAln_NucTarget_Len <- append(Output_03a2_MyWorkingAln_NucTarget_Len,NA)
      Output_04a_MyIndex_start <- append(Output_04a_MyIndex_start,NA)
      Output_04b_MyIndex_end <- append(Output_04b_MyIndex_end,NA)
      Output_04c_MyIndex_ATCG_start <- append(Output_04c_MyIndex_ATCG_start,NA)
      Output_04d_MyIndex_ATCG_end <- append(Output_04d_MyIndex_ATCG_end,NA)
      Output_05a_QC1_lentrim1vs2 <- append(Output_05a_QC1_lentrim1vs2,NA)
      Output_05b_QC2a_PercIdentity <- append(Output_05b_QC2a_PercIdentity,NA)
      Output_05c_QC2b_PerIden_Cutoff_Nuc <- append(Output_05c_QC2b_PerIden_Cutoff_Nuc,NA)
      #Step 09b.  Slab
      REV_slab_NUC <- paste(Store_REVexon1_NUC_nodash[o],Store_REVexon2GP41_NUC_nodash[o],sep="") #ACH-1 (A1)
      names(REV_slab_NUC) <- names(Store_REVexon1_NUC_nodash[o])
      REV_slab_NUC_DNAStringSet <- DNAStringSet(REV_slab_NUC)
      #Step 09c. Translate
      REV_slab_AA_AAStringSet <- translate(REV_slab_NUC_DNAStringSet)
      #Step 09d. Index first stop
      MyIndex_FirstStop <- unlist(gregexpr("\\*",REV_slab_AA_AAStringSet,ignore.case=TRUE)) #class integer
      #what if there is no stop codon
      if (MyIndex_FirstStop[1]>=1) {
        MyIndex_FirstStop <- MyIndex_FirstStop[1]
      } else if (MyIndex_FirstStop[1]<1) {
        MyIndex_FirstStop <- nchar(REV_slab_AA_AAStringSet)
      } else {
        MyIndex_FirstStop <- nchar(REV_slab_AA_AAStringSet)
      }
      #Step 09e. Store the translation & len
      REV_slab_AA_AAString <- toString(REV_slab_AA_AAStringSet[1]) #convert to AAString
      MyFinal_AA_character <- substr(REV_slab_AA_AAString,1,MyIndex_FirstStop[1]-1)
      if (nchar(MyFinal_AA_character)>0){
        Output_06a_MyFinal_AA <- append(Output_06a_MyFinal_AA,MyFinal_AA_character)
        Output_06b_MyFinal_AA_len <- append(Output_06b_MyFinal_AA_len,nchar(MyFinal_AA_character))
        #Step 09f. Store NUC corresponding to AA, length
        MyFinal_Nuc <- substr(REV_slab_NUC,1,nchar(MyFinal_AA_character)*3)
        Output_06c_MyFinal_Nuc <- append(Output_06c_MyFinal_Nuc,MyFinal_Nuc)
        Output_06d_MyFinal_Nuc_len <- append(Output_06d_MyFinal_Nuc_len,nchar(MyFinal_Nuc))
        #Step 09g. Store the 4bp STOP codon
        MyFinal_Stop_4bp <- substr(REV_slab_NUC,nchar(MyFinal_AA_character)*3+1,nchar(MyFinal_AA_character)*3+4)
        Output_06e_MyFinal_Stop_4bp <- append(Output_06e_MyFinal_Stop_4bp,MyFinal_Stop_4bp)
        #Step 09h. QC3 Any seq 50% (user settings) shorter than expected would be flagged
        MyExpectedAALength <- width(MyRefAA_HXB2[19])
        Output_06b2_MyAA_Expected_len <- append(Output_06b2_MyAA_Expected_len,MyExpectedAALength)
        MyActualAALength <- as.numeric(nchar(MyFinal_AA_character))
        if (MyActualAALength>MyExpectedAALength*AA_TooShort & MyActualAALength<MyExpectedAALength*AA_TooLong){
          Output_06f_QC3_AA_TooShort <- append(Output_06f_QC3_AA_TooShort,"Passed QC3")
        } else if (MyActualAALength<=MyExpectedAALength*AA_TooShort) {
          Output_06f_QC3_AA_TooShort <- append(Output_06f_QC3_AA_TooShort,"Failed QC3; AA_TooShort")
        } else if (MyActualAALength>=MyExpectedAALength*AA_TooLong) {
          Output_06f_QC3_AA_TooShort <- append(Output_06f_QC3_AA_TooShort,"Failed QC3; AA_TooLong")
        } else {
          Output_06f_QC3_AA_TooShort <- append(Output_06f_QC3_AA_TooShort,"Failed QC3; Check") 
        }
        #Step 09i. Export AA aligned; one-before-stop
        REV_AA_w_HXB2_char <- c(HXB2_REV_SLAB_AA,MyFinal_AA_character)
        names(REV_AA_w_HXB2_char) <- c("HXB2_REV_SLAB_AA",names(Store_REVexon1_NUC_nodash[o]))
        REV_AA_w_HXB2_AAStringSet <- AAStringSet(REV_AA_w_HXB2_char)
        REV_AA_w_HXB2_muscle_Aln <- muscle(REV_AA_w_HXB2_AAStringSet,gapopen=MyGapPenaltyAA,log=MuscleAALog,verbose=TRUE) #class AAMultipleAlignment
        REV_AA_w_HXB2_muscle_Aln_AAStringSet <- AAStringSet(REV_AA_w_HXB2_muscle_Aln)
        writeXStringSet(REV_AA_w_HXB2_muscle_Aln_AAStringSet,MyWorkingAln_AA_muscleAlign,append=TRUE)
        #Step 09j. Export NUC aligned
        #(1) Grab HXB2 and mycurrent, combine into one DNAStringSet
        HXB2_REV_SLAB_NUC #character
        #names(HXB2_REV_SLAB_NUC) <- "HXB2_REV_SLAB_NUC"
        HXB2_REV_SLAB_NUC_DNAStringSet <- DNAStringSet(HXB2_REV_SLAB_NUC)
        HXB2_REV_SLAB_NUC_DNAStringSet_2 <- append(HXB2_REV_SLAB_NUC_DNAStringSet,REV_slab_NUC_DNAStringSet)
        #(2) Write NUC as fasta
        writeXStringSet(HXB2_REV_SLAB_NUC_DNAStringSet_2,MyWorkingAln_Nuc_muscleAlign,append=TRUE) #aligned_NUC_trimmed
        #####################################################################################
        #Step 09k. Export AA aligned; trimmed to HXB2
        #(1) Grab unaligned NUC
        #DNAStringSet(HXB2_REV_SLAB_NUC)
        MyWorkingREVslab_unaligned_DNAStringSet <- append(DNAStringSet(HXB2_REV_SLAB_NUC),REV_slab_NUC_DNAStringSet)
        #(2) Translate both nodash fasta
        MyWorkingREVslab_unaligned_AAStringSet <- translate(MyWorkingREVslab_unaligned_DNAStringSet)
        #(3) Align the AAs
        #MyWorkingREVslab_unaligned_AAStringSet_muscle <- muscle(MyWorkingREVslab_unaligned_AAStringSet,gapopen=MyGapPenaltyNUC,log=MuscleNUCLog,verbose=TRUE)
        ######AHHHHHH muscle is producing wrong alignment! ignoring/skipping/eliminating stops ##############
        ##########use Biostrings function pairwiseAlignment(), default settings##########################
        MyWorkingREVslab_unaligned_AAStringSet_pairAlign <- pairwiseAlignment(MyWorkingREVslab_unaligned_AAStringSet[1],MyWorkingREVslab_unaligned_AAStringSet[2],type="global-local")
        #MyWorkingREVslab_unaligned_AAStringSet_pairAlign #good!  stop codons preserved
        #class(MyWorkingREVslab_unaligned_AAStringSet_pairAlign)  #[1] "PairwiseAlignmentsSingleSubject"
        #summary(MyWorkingREVslab_unaligned_AAStringSet_pairAlign)
        #writePairwiseAlignments(MyWorkingREVslab_unaligned_AAStringSet_pairAlign,file="test_REVpair.txt")
        #(4) Convert back to AAStringSet
        #AAStringSet((pattern(MyWorkingREVslab_unaligned_AAStringSet_pairAlign)))
        #AAStringSet((subject(MyWorkingREVslab_unaligned_AAStringSet_pairAlign)))
        MyFinal_AA_paired <- AAStringSet((pattern(MyWorkingREVslab_unaligned_AAStringSet_pairAlign)))
        MyFinal_AA_paired2 <- append(MyFinal_AA_paired,AAStringSet((subject(MyWorkingREVslab_unaligned_AAStringSet_pairAlign))))
        names(MyFinal_AA_paired2) <- c("HXB2_REV_SLAB_AA",names(Store_REVexon1_NUC_nodash[o]))
        writeXStringSet(MyFinal_AA_paired2,MyWorkingAln_AA_AAStringSet_trimmed,append=TRUE)
        #(5) Save to df
        Output_07a_AA_trimmed_HXB2 <- append(Output_07a_AA_trimmed_HXB2,toString(pattern(MyWorkingREVslab_unaligned_AAStringSet_pairAlign)))
        Output_07b_AA_trimmed_sample <- append(Output_07b_AA_trimmed_sample,toString(subject(MyWorkingREVslab_unaligned_AAStringSet_pairAlign)))
        ###########################################
      } else {
        Output_06a_MyFinal_AA <- append(Output_06a_MyFinal_AA,NA)
        Output_06b_MyFinal_AA_len <- append(Output_06b_MyFinal_AA_len,0)
        Output_06c_MyFinal_Nuc <- append(Output_06c_MyFinal_Nuc,NA)
        Output_06d_MyFinal_Nuc_len <- append(Output_06d_MyFinal_Nuc_len,0)
        Output_06e_MyFinal_Stop_4bp <- append(Output_06e_MyFinal_Stop_4bp,NA)
        MyExpectedAALength <- width(MyRefAA_HXB2[19])
        Output_06b2_MyAA_Expected_len <- append(Output_06b2_MyAA_Expected_len,MyExpectedAALength)
        Output_06f_QC3_AA_TooShort <- append(Output_06f_QC3_AA_TooShort,"Failed QC3; AA_TooShort")
        Output_07a_AA_trimmed_HXB2 <- append(Output_07a_AA_trimmed_HXB2,NA)
        Output_07b_AA_trimmed_sample <- append(Output_07b_AA_trimmed_sample,NA)
      }
    }
  } else {
    Output_01_MyWorkingRefNuc_Name <- append(Output_01_MyWorkingRefNuc_Name,"HXB2_TAT-slab_REV-slab_failed")
    Output_01a_MyWorkingRefNuc_UnAlign_Seq <- append(Output_01a_MyWorkingRefNuc_UnAlign_Seq,NA)
    Output_02_MyWorkingPreAlnNuc_Name <- append(Output_02_MyWorkingPreAlnNuc_Name,NA)
    Output_02a_QC0_TargetSeqATCG <- append(Output_02a_QC0_TargetSeqATCG,NA)
    Output_03a_MyWorkingAln_NucRef <- append(Output_03a_MyWorkingAln_NucRef,NA)
    Output_03b_MyWorkingAln_NucTarget <- append(Output_03b_MyWorkingAln_NucTarget,NA)
    Output_03a1_MyWorkingAln_NucRef_Len <- append(Output_03a1_MyWorkingAln_NucRef_Len,NA)
    Output_03a2_MyWorkingAln_NucTarget_Len <- append(Output_03a2_MyWorkingAln_NucTarget_Len,NA)
    Output_04a_MyIndex_start <- append(Output_04a_MyIndex_start,NA)
    Output_04b_MyIndex_end <- append(Output_04b_MyIndex_end,NA)
    Output_04c_MyIndex_ATCG_start <- append(Output_04c_MyIndex_ATCG_start,NA)
    Output_04d_MyIndex_ATCG_end <- append(Output_04d_MyIndex_ATCG_end,NA)
    Output_05a_QC1_lentrim1vs2 <- append(Output_05a_QC1_lentrim1vs2,NA)
    Output_05b_QC2a_PercIdentity <- append(Output_05b_QC2a_PercIdentity,NA)
    Output_05c_QC2b_PerIden_Cutoff_Nuc <- append(Output_05c_QC2b_PerIden_Cutoff_Nuc,NA)
    Output_06a_MyFinal_AA <- append(Output_06a_MyFinal_AA,NA)
    Output_06b_MyFinal_AA_len <- append(Output_06b_MyFinal_AA_len,NA)
    Output_06b2_MyAA_Expected_len <- append(Output_06b2_MyAA_Expected_len,NA)
    Output_06c_MyFinal_Nuc <- append(Output_06c_MyFinal_Nuc,NA)
    Output_06d_MyFinal_Nuc_len <- append(Output_06d_MyFinal_Nuc_len,NA)
    Output_06e_MyFinal_Stop_4bp <- append(Output_06e_MyFinal_Stop_4bp,NA)
    Output_06f_QC3_AA_TooShort <- append(Output_06f_QC3_AA_TooShort,NA)
    Output_07a_AA_trimmed_HXB2 <- append(Output_07a_AA_trimmed_HXB2,NA)
    Output_07b_AA_trimmed_sample <- append(Output_07b_AA_trimmed_sample,NA)
  }
  
  
  
  #######################################################
  
  length(Output_01_MyWorkingRefNuc_Name)
  length(Output_01a_MyWorkingRefNuc_UnAlign_Seq)
  length(Output_02_MyWorkingPreAlnNuc_Name)
  length(Output_02a_QC0_TargetSeqATCG)
  length(Output_03a_MyWorkingAln_NucRef) #36 OK :)
  length(Output_03b_MyWorkingAln_NucTarget) #36 OK :)
  length(Output_03a1_MyWorkingAln_NucRef_Len)  #36 OK :)
  length(Output_03a2_MyWorkingAln_NucTarget_Len) #36 OK :)
  length(Output_04a_MyIndex_start) #36 OK :)
  length(Output_04b_MyIndex_end) #36 OK :)
  length(Output_04c_MyIndex_ATCG_start) #36 OK :)
  length(Output_04d_MyIndex_ATCG_end) #36 OK :)
  length(Output_05a_QC1_lentrim1vs2) 
  length(Output_05b_QC2a_PercIdentity)
  length(Output_05c_QC2b_PerIden_Cutoff_Nuc)
  length(Output_06a_MyFinal_AA)
  length(Output_06b_MyFinal_AA_len)
  length(Output_06b2_MyAA_Expected_len)
  length(Output_06c_MyFinal_Nuc)
  length(Output_06d_MyFinal_Nuc_len)
  length(Output_06e_MyFinal_Stop_4bp)
  length(Output_06f_QC3_AA_TooShort)
  length(Output_07a_AA_trimmed_HXB2)
  length(Output_07b_AA_trimmed_sample)
  length(Output_08a_ATG)
  length (Output_08b_QC4_ATG)
  
  
  mydf <- data.frame(
    Output_01_MyWorkingRefNuc_Name,
    Output_01a_MyWorkingRefNuc_UnAlign_Seq,
    Output_02_MyWorkingPreAlnNuc_Name,
    Output_02a_QC0_TargetSeqATCG,
    Output_03a_MyWorkingAln_NucRef,
    Output_03b_MyWorkingAln_NucTarget,
    Output_03a1_MyWorkingAln_NucRef_Len,
    Output_03a2_MyWorkingAln_NucTarget_Len,
    Output_04a_MyIndex_start,
    Output_04b_MyIndex_end,
    Output_04c_MyIndex_ATCG_start,
    Output_04d_MyIndex_ATCG_end,
    Output_05a_QC1_lentrim1vs2,
    Output_05b_QC2a_PercIdentity,
    Output_05c_QC2b_PerIden_Cutoff_Nuc,
    Output_06a_MyFinal_AA,
    Output_06b_MyFinal_AA_len,
    Output_06b2_MyAA_Expected_len,
    Output_06c_MyFinal_Nuc,
    Output_06d_MyFinal_Nuc_len,
    Output_06e_MyFinal_Stop_4bp,
    Output_06f_QC3_AA_TooShort,
    Output_07a_AA_trimmed_HXB2,
    Output_07b_AA_trimmed_sample
  )
  
  #### QC4: Check for ATG start ###
  Output_08a_ATG <- as.character(substring(Output_07b_AA_trimmed_sample,1,1))
  #Only export QC4 for gag, vif, vpr, vpu, tat, rev, env, nef)
  Output_08b_QC4_ATG <- c()
  for (q in 1:length(Output_08a_ATG)){
    if (mydf$Output_01_MyWorkingRefNuc_Name[q]=="HXB2_gag_790-2289" | mydf$Output_01_MyWorkingRefNuc_Name[q]=="HXB2_vif_5041-5616" | mydf$Output_01_MyWorkingRefNuc_Name[q]=="HXB2_vpr_5559-5847mod" | mydf$Output_01_MyWorkingRefNuc_Name[q]=="HXB2_vpu_6062-6307" | mydf$Output_01_MyWorkingRefNuc_Name[q]=="HXB2_env_6225-8792" | mydf$Output_01_MyWorkingRefNuc_Name[q]=="HXB2_nef_8797-9414" | mydf$Output_01_MyWorkingRefNuc_Name[q]=="HXB2_TAT_slab" | mydf$Output_01_MyWorkingRefNuc_Name[q]=="HXB2_REV_slab") {
      if (is.na(Output_08a_ATG[q])){
        Output_08b_QC4_ATG <- append(Output_08b_QC4_ATG,NA)
      } else if (Output_08a_ATG[q]=="M") {
        Output_08b_QC4_ATG <- append(Output_08b_QC4_ATG,"Passed QC4; M")
      } else if (Output_08a_ATG[q]!="M") {
        Output_08b_QC4_ATG <- append(Output_08b_QC4_ATG,"Failed QC4; Not M")
      } else {
        Output_08b_QC4_ATG <- append(Output_08b_QC4_ATG,NA)
      }
    } else {
      Output_08b_QC4_ATG <- append(Output_08b_QC4_ATG,NA)
    }
  }
  mydf <- cbind(mydf
                ,Output_08a_ATG
                ,Output_08b_QC4_ATG)
  
  write.csv(mydf,file=Output_df,row.names=FALSE)
  
}



#Combine the df files in Folder "ToHIVSeqinR_DB#
AllFiles <- list.files("./ToHIVSeqinR_DB")
Output_CSV <- paste("Output_HIVSeqinR_",current_pipeline_version,"_mydf.csv",sep="")
mydf <- data.frame()
for (n in 1:length(AllFiles)){  
  #Step01. Read csv
  InputCSVFileName <- AllFiles[n]
  MyCSV <- read.csv(paste("./ToHIVSeqinR_DB/",InputCSVFileName,sep=""),header=TRUE)
  #View(MyCSV)
  #Step02.  Append to a master data.frame object
  mydf <- rbind(mydf,MyCSV)
  #View(MyAllCSV)
}
write.csv(mydf,file=Output_CSV)


#################################################################
#############Graphing############################################


# Conditions
# QC2 NA QC3 pass (TAT/REV slab) #>>>> light green "darkolivegreen4"
# QC2 NA QC3 fail (TAT/REV slab) #>>>>> red "red"
# QC2 pass QC3 pass green #>>>>> darkgreen "darkgreen "
# QC2 fail  #>>>>> grey "grey" #QC2 - likely missing this region
# QC3 fail  #>>>>> red "red"  #QC3 - likely premature stop
# QC2 pass QC3 NA (TAT/REV exon1 or 2) #>>>>> light green "darkolivegreen4"
# QC2 fail QC3 NA (TAT/REV exon1 or 2) #>>>>> red "red"
# else #>>>>> black "black"

Output_01_MyWorkingRefNuc_Name <- mydf$Output_01_MyWorkingRefNuc_Name
Output_02_MyWorkingPreAlnNuc_Name <- mydf$Output_02_MyWorkingPreAlnNuc_Name
Output_05c_QC2b_PerIden_Cutoff_Nuc <- mydf$Output_05c_QC2b_PerIden_Cutoff_Nuc
Output_06f_QC3_AA_TooShort <- mydf$Output_06f_QC3_AA_TooShort

#Step 01.  Create a vector indicating sample/genes that failed QC2 "too many mismatches" AND/or QC3 "AA too short" with color (darkgreen/red/orange)
Output_08a_QC2n3_COLOR <- c()
for (p in 1:length(Output_05c_QC2b_PerIden_Cutoff_Nuc)) {
  if (is.na(Output_05c_QC2b_PerIden_Cutoff_Nuc[p]) | is.na(Output_06f_QC3_AA_TooShort[p])){
    if (is.na(Output_05c_QC2b_PerIden_Cutoff_Nuc[p])) {
      if (substr(Output_06f_QC3_AA_TooShort[p],1,6)=="Passed"){
        Output_08a_QC2n3_COLOR <- append(Output_08a_QC2n3_COLOR,"darkolivegreen4")
        #print("darkolivegreen4")
      } else if (substr(Output_06f_QC3_AA_TooShort[p],1,6)=="Failed"){
        Output_08a_QC2n3_COLOR <- append(Output_08a_QC2n3_COLOR,"red")
        #print("red")
      } else {
        Output_08a_QC2n3_COLOR <- append(Output_08a_QC2n3_COLOR,"black")
        #print("black")
      }
    } else if (is.na(Output_06f_QC3_AA_TooShort[p])){
      if (substr(Output_05c_QC2b_PerIden_Cutoff_Nuc[p],1,6)=="Passed"){
        Output_08a_QC2n3_COLOR <- append(Output_08a_QC2n3_COLOR,"darkolivegreen4")
        #print("darkolivegreen4")
      } else if (substr(Output_05c_QC2b_PerIden_Cutoff_Nuc[p],1,6)=="Failed"){
        Output_08a_QC2n3_COLOR <- append(Output_08a_QC2n3_COLOR,"grey")
        #print("red")
      } else {
        Output_08a_QC2n3_COLOR <- append(Output_08a_QC2n3_COLOR,"black")
        #print("black")
      }
    }
  } else if (substr(Output_05c_QC2b_PerIden_Cutoff_Nuc[p],1,6)=="Passed" & substr(Output_06f_QC3_AA_TooShort[p],1,6)=="Passed") { #QC2 pass QC 3 pass #Gene is intact#
    Output_08a_QC2n3_COLOR <- append(Output_08a_QC2n3_COLOR,"darkgreen")
    #print("darkgreen")
  } else if (substr(Output_05c_QC2b_PerIden_Cutoff_Nuc[p],1,6)=="Failed" & substr(Output_06f_QC3_AA_TooShort[p],1,6)=="Passed") {  #QC2 fail QC3 passed #likely low PercIdent#
    Output_08a_QC2n3_COLOR <- append(Output_08a_QC2n3_COLOR,"darkolivegreen4")
  } else if (substr(Output_05c_QC2b_PerIden_Cutoff_Nuc[p],1,6)=="Passed" & substr(Output_06f_QC3_AA_TooShort[p],1,6)=="Failed") {  #QC2 pass QC3 fail #RED# - likely premature stop
    Output_08a_QC2n3_COLOR <- append(Output_08a_QC2n3_COLOR,"red")
    #print("red")
  } else if (substr(Output_05c_QC2b_PerIden_Cutoff_Nuc[p],1,6)=="Failed" & substr(Output_06f_QC3_AA_TooShort[p],1,6)=="Failed") {  #QC2 fail QC3 fail #RED# - likely missing the region
    Output_08a_QC2n3_COLOR <- append(Output_08a_QC2n3_COLOR,"grey")
    #print("red")
  } else {
    Output_08a_QC2n3_COLOR <- append(Output_08a_QC2n3_COLOR,"black")
  }
}
#length(Output_08a_QC2n3_COLOR)  #50 correct


start_x <- c(1,790,2085,2253,2550,3870,4230,5041,5559,5831,5970,6062,6225,6315,7110,7758,8379,8379,8797,5831,5970)
end_x <- c(9719,2289,5093,2549,3869,4229,5093,5616,5847,6045,6045,6307,8792,7757,7217,8792,8466,8650,9414,(5831+303),(5970+348))
#length(start_x) #20
#legnth(end_x) #20
mydf_plot <- data.frame(Output_01_MyWorkingRefNuc_Name,Output_02_MyWorkingPreAlnNuc_Name,Output_08a_QC2n3_COLOR)
MySampleList <- unique(Output_02_MyWorkingPreAlnNuc_Name)
# View(mydf_plot)

### Here ###

pdf(file=paste("Output_HIVSeqinR_",current_pipeline_version,".pdf",sep=""))
for (q in 1:length(MySampleList)){ #for each A1 and A2
  MySampleName <- MySampleList[q]
  mydf_plot_per_sample_2 <- subset(mydf_plot,Output_02_MyWorkingPreAlnNuc_Name==MySampleList[q])
  #mydf_plot_per_sample_2 <- subset(mydf_plot_per_sample_1,Output_01_MyWorkingRefNuc_Name!="HXB2_pol_2085-5093")
  mydf_plot_per_sample_3 <- subset(mydf_plot_per_sample_2,Output_01_MyWorkingRefNuc_Name!="HXB2_PR-RT_2253-3869")
  mydf_plot_per_sample_4 <- subset(mydf_plot_per_sample_3,Output_01_MyWorkingRefNuc_Name!="HXB2_PR-RNaseH_2253-4229")
  mydf_plot_per_sample_5 <- subset(mydf_plot_per_sample_4,Output_01_MyWorkingRefNuc_Name!="HXB2_PR-INT_2253-5093")
  mydf_plot_per_sample_6 <- subset(mydf_plot_per_sample_5,Output_01_MyWorkingRefNuc_Name!="HXB2_RT-INT_3720-5093")
  #mydf_plot_per_sample_color <- toString(mydf_plot_per_sample_6$Output_08a_QC2n3_COLOR) #class: factor
  #class(mydf_plot_per_sample_color)
  #nrow(mydf_plot_per_sample_6) #20
  #Step 01. Generate plot frame
  #plot.new()
  #pdf(file=paste(MySampleName,".pdf",sep=""),width=4000,height=3500)
  x=0
  y=0
  plot(x,y,xlim=c(0,10000),ylim=c(0,30000))
  title(MySampleName)
  #Step 02. Draw HXB2 and label; define start_y 
  start_y <- 28000
  rect(1,start_y,9719,start_y+200,density=200,angle=180,col="black") #rect(xleft,ybottom,xright,ytop)  #density is the "fill" in "lines per inch", default = NULL; angle is the angle of the shading line, default = 45
  text(5000,start_y+800,labels="HXB2_K03455_FULL_1-9717 bp",col="black",cex=0.8)
  #Step 03. Draw a legend
  #Method 1.  Only display legend for colors in the plot
  #legend_color <- unique(mydf_plot_per_sample_6$Output_08a_QC2n3_COLOR)
  #length(legend_color)
  #Method 2.  Display all possible colors in the script
  legend_color <- c("darkgreen","darkolivegreen4","grey","red","black")
  legend <- c("Passed QC2 and QC3","Passed only QC2 (TAT/REV) or only QC3","potential deletion (Failed QC2 and QC3)","premature stop codon (Passed QC2 Failed QC3)","unknown error")
  legend_start_y <- 12500
  for (c in 1:length(legend_color)){
    legend_start_y_current <- legend_start_y-(1200*(c-1))
    rect(0,legend_start_y_current,500,
         legend_start_y_current-500,density=200,angle=180,
         col=legend_color[c])
    #rect(xleft, ybottom, xright, ytop)
    text(2500,legend_start_y_current-250,
         labels=legend[c],
         col="black",
         cex=0.5)
  }
  text(2500,legend_start_y_current-1500,labels="QC2: <70% seq identity against HXB2",cex=0.5)
  text(2500,legend_start_y_current-2500,labels="QC3: amino acid length =<95% or >=120% than expected",cex=0.5)
  #Step 03. Loop through the genes
  for (g in 2:length(mydf_plot_per_sample_6$Output_08a_QC2n3_COLOR)){
    start_y_current <- (start_y-2000)-(1000*(g-1))
    #print(start_y_current)
    rect(start_x[g],start_y_current,end_x[g],
         start_y_current+200,density=200,angle=180,
         col=toString(mydf_plot_per_sample_6$Output_08a_QC2n3_COLOR[g]))
    #     print(start_x[g])
    #     print(start_y_current)
    #     print(end_x[g])
    #     print(start_y_current+200)
    #     print(toString(mydf_plot_per_sample$Output_08a_QC2n3_COLOR[g]))
    text(start_x[g],start_y_current+500,
         labels=substr(toString(mydf_plot_per_sample_6$Output_01_MyWorkingRefNuc_Name[g]),
                       6,nchar(toString(mydf_plot_per_sample_6$Output_01_MyWorkingRefNuc_Name[g]))),
         col=toString(mydf_plot_per_sample_6$Output_08a_QC2n3_COLOR[g]),
         cex=0.5)
  }
}
dev.off()
graphics.off() 
### Thanks Oyoung!


##########################################
#20171127 HIVSeqinR Pivot gag, pol, env #
##########################################

#Reload HIVSeqinR Table
MyHIVSeqinR_OutputCSV_DF_Uniq <- as.character(unique(mydf$Output_02_MyWorkingPreAlnNuc_Name))
MyGAGPassFail <- c()
MyPOLPassFail <- c()
MyPRPassFail <- c()
MyRTPassFail <- c()
MyRNaseHPassFail <- c()
MyINTPassFail <- c()
MyENVPassFail <- c()
MyGP41PassFail <- c()
MyAnyATGFailed <- c()
MyGAGLenUntrimmed <- c()
MyPOLLenUntrimmed <- c()
MyPRLenUntrimmed <- c()
MyRTLenUntrimmed <- c()
MyRNaseHLenUntrimmed <- c()
MyINTLenUntrimmed <- c()
MyENVLenUntrimmed <- c()
MyGP41LenUntrimmed <- c()
MyGAGLenHXB2 <- c()
MyPOLLenHXB2 <- c()
MyPRLenHXB2 <- c()
MyRTLenHXB2 <- c()
MyRNaseHLenHXB2 <- c()
MyINTLenHXB2 <- c()
MyENVLenHXB2 <- c()
MyGP41LenHXB2 <- c()
for (u in 1:length(MyHIVSeqinR_OutputCSV_DF_Uniq)){
  MyHIVSeqinR_CurrentDF_all <- subset(mydf,mydf$Output_02_MyWorkingPreAlnNuc_Name==MyHIVSeqinR_OutputCSV_DF_Uniq[u])
  MyHIVSeqinR_CurrentDF_GAG <- subset(MyHIVSeqinR_CurrentDF_all,MyHIVSeqinR_CurrentDF_all$Output_01_MyWorkingRefNuc_Name=="HXB2_gag_790-2289")
  MyHIVSeqinR_CurrentDF_POL <- subset(MyHIVSeqinR_CurrentDF_all,MyHIVSeqinR_CurrentDF_all$Output_01_MyWorkingRefNuc_Name=="HXB2_pol_2085-5093")
  MyHIVSeqinR_CurrentDF_PR <- subset(MyHIVSeqinR_CurrentDF_all,MyHIVSeqinR_CurrentDF_all$Output_01_MyWorkingRefNuc_Name=="HXB2_PR_2253-2549")
  MyHIVSeqinR_CurrentDF_RT <- subset(MyHIVSeqinR_CurrentDF_all,MyHIVSeqinR_CurrentDF_all$Output_01_MyWorkingRefNuc_Name=="HXB2_RT_2550-3869")
  MyHIVSeqinR_CurrentDF_RNaseH <- subset(MyHIVSeqinR_CurrentDF_all,MyHIVSeqinR_CurrentDF_all$Output_01_MyWorkingRefNuc_Name=="HXB2_RNaseH_3870-4229")
  MyHIVSeqinR_CurrentDF_INT <- subset(MyHIVSeqinR_CurrentDF_all,MyHIVSeqinR_CurrentDF_all$Output_01_MyWorkingRefNuc_Name=="HXB2_INT_4230-5093")
  MyHIVSeqinR_CurrentDF_ENV <- subset(MyHIVSeqinR_CurrentDF_all,MyHIVSeqinR_CurrentDF_all$Output_01_MyWorkingRefNuc_Name=="HXB2_env_6225-8792")
  MyHIVSeqinR_CurrentDF_GP41 <- subset(MyHIVSeqinR_CurrentDF_all,MyHIVSeqinR_CurrentDF_all$Output_01_MyWorkingRefNuc_Name=="HXB2_gp41_7758-8792")
  MyGAGPassFail <- append(MyGAGPassFail,as.character(MyHIVSeqinR_CurrentDF_GAG[1,"Output_06f_QC3_AA_TooShort"]))
  MyPOLPassFail <- append(MyPOLPassFail,as.character(MyHIVSeqinR_CurrentDF_POL[1,"Output_06f_QC3_AA_TooShort"]))
  MyPRPassFail <- append(MyPRPassFail,as.character(MyHIVSeqinR_CurrentDF_PR[1,"Output_06f_QC3_AA_TooShort"]))
  MyRTPassFail <- append(MyRTPassFail,as.character(MyHIVSeqinR_CurrentDF_RT[1,"Output_06f_QC3_AA_TooShort"]))
  MyRNaseHPassFail <- append(MyRNaseHPassFail,as.character(MyHIVSeqinR_CurrentDF_RNaseH[1,"Output_06f_QC3_AA_TooShort"]))
  MyINTPassFail <- append(MyINTPassFail,as.character(MyHIVSeqinR_CurrentDF_INT[1,"Output_06f_QC3_AA_TooShort"]))
  MyENVPassFail <- append(MyENVPassFail,as.character(MyHIVSeqinR_CurrentDF_ENV[1,"Output_06f_QC3_AA_TooShort"]))
  MyGP41PassFail <- append(MyGP41PassFail,as.character(MyHIVSeqinR_CurrentDF_GP41[1,"Output_06f_QC3_AA_TooShort"]))
  MyAnyATGFailed <- any((substring(MyHIVSeqinR_CurrentDF_all$Output_08b_QC4_ATG,1,6)=="Failed")==TRUE,na.rm=TRUE)
  MyGAGLenUntrimmed <- append(MyGAGLenUntrimmed,as.character(MyHIVSeqinR_CurrentDF_GAG[1,"Output_06b_MyFinal_AA_len"]))
  MyPOLLenUntrimmed <- append(MyPOLLenUntrimmed,as.character(MyHIVSeqinR_CurrentDF_POL[1,"Output_06b_MyFinal_AA_len"]))
  MyPRLenUntrimmed <- append(MyPRLenUntrimmed,as.character(MyHIVSeqinR_CurrentDF_PR[1,"Output_06b_MyFinal_AA_len"]))
  MyRTLenUntrimmed <- append(MyRTLenUntrimmed,as.character(MyHIVSeqinR_CurrentDF_RT[1,"Output_06b_MyFinal_AA_len"]))
  MyRNaseHLenUntrimmed <- append(MyRNaseHLenUntrimmed,as.character(MyHIVSeqinR_CurrentDF_RNaseH[1,"Output_06b_MyFinal_AA_len"]))
  MyINTLenUntrimmed <- append(MyINTLenUntrimmed,as.character(MyHIVSeqinR_CurrentDF_INT[1,"Output_06b_MyFinal_AA_len"]))
  MyENVLenUntrimmed <- append(MyENVLenUntrimmed,as.character(MyHIVSeqinR_CurrentDF_ENV[1,"Output_06b_MyFinal_AA_len"]))
  MyGP41LenUntrimmed <- append(MyGP41LenUntrimmed,as.character(MyHIVSeqinR_CurrentDF_GP41[1,"Output_06b_MyFinal_AA_len"]))
  MyGAGLenHXB2 <- append(MyGAGLenHXB2,as.character(MyHIVSeqinR_CurrentDF_GAG[1,"Output_06b2_MyAA_Expected_len"]))
  MyPOLLenHXB2 <- append(MyPOLLenHXB2,as.character(MyHIVSeqinR_CurrentDF_POL[1,"Output_06b2_MyAA_Expected_len"]))
  MyPRLenHXB2 <- append(MyPRLenHXB2,as.character(MyHIVSeqinR_CurrentDF_PR[1,"Output_06b2_MyAA_Expected_len"]))
  MyRTLenHXB2 <- append(MyRTLenHXB2,as.character(MyHIVSeqinR_CurrentDF_RT[1,"Output_06b2_MyAA_Expected_len"]))
  MyRNaseHLenHXB2 <- append(MyRNaseHLenHXB2,as.character(MyHIVSeqinR_CurrentDF_RNaseH[1,"Output_06b2_MyAA_Expected_len"]))
  MyINTLenHXB2 <- append(MyINTLenHXB2,as.character(MyHIVSeqinR_CurrentDF_INT[1,"Output_06b2_MyAA_Expected_len"]))
  MyENVLenHXB2 <- append(MyENVLenHXB2,as.character(MyHIVSeqinR_CurrentDF_ENV[1,"Output_06b2_MyAA_Expected_len"]))
  MyGP41LenHXB2 <- append(MyGP41LenHXB2,as.character(MyHIVSeqinR_CurrentDF_GP41[1,"Output_06b2_MyAA_Expected_len"]))
}

MyHIVSeqinR_OutputPivot <- cbind(
  MyHIVSeqinR_OutputCSV_DF_Uniq
  ,MyGAGPassFail
  ,MyPOLPassFail
  ,MyPRPassFail
  ,MyRTPassFail
  ,MyRNaseHPassFail
  ,MyINTPassFail
  ,MyENVPassFail
  ,MyGP41PassFail
  ,MyAnyATGFailed
  ,MyGAGLenUntrimmed
  ,MyPOLLenUntrimmed
  ,MyPRLenUntrimmed
  ,MyRTLenUntrimmed
  ,MyRNaseHLenUntrimmed
  ,MyINTLenUntrimmed
  ,MyENVLenUntrimmed
  ,MyGP41LenUntrimmed
  ,MyGAGLenHXB2
  ,MyPOLLenHXB2
  ,MyPRLenHXB2
  ,MyRTLenHXB2
  ,MyRNaseHLenHXB2
  ,MyINTLenHXB2
  ,MyENVLenHXB2
  ,MyGP41LenHXB2
)

#View(MyHIVSeqinR_OutputCSV_DF)
#View(MyHIVSeqinR_OutputPivot)
write.csv(MyHIVSeqinR_OutputPivot,file="Output_MyHIVSeqinR_Pivot.csv")

#############################
# Part II. R04_BigSummary_ver45_PSC5DEFECT.R
#############################

#Step 00. Setup R #PC@work R version 3.2.2
rm(list= ls()[!(ls() %in% c("MyBlastnDir","StartSysTime","Primer2ndF_HXB2","Primer2ndR_HXB2"))])
closeAllConnections() 
graphics.off() 
cat("\014") #send CTRL+L to R, clean console
Sys.setenv(TZ="EST")
MyWD <- getwd()
setwd(MyWD)
library(Biostrings) #PC@work Biostrings_2.38.1
#library(muscle) #PC@work muscle_3.12.0
#library(ape)
#sessionInfo()

DF_FINAL <- read.csv("Output_MyBigSummary_DF_FINAL_half.csv")[,-1]
MyHIVSeqinRPivot <- read.csv("Output_MyHIVSeqinR_Pivot.csv")
MyHIVSeqinRPivot <- MyHIVSeqinRPivot[,2:ncol(MyHIVSeqinRPivot)]
DF_FINAL <- merge(DF_FINAL,MyHIVSeqinRPivot,by.x="SEQID",by.y="MyHIVSeqinR_OutputCSV_DF_Uniq",all=TRUE)
#View(MyHIVSeqinRPivot)
#View(DF_FINAL)

##################################################################


## FINAL VERDICT ##

## R script for VERDICT; order matters ##
MyVerdict <- c()
for (m in 1:length(DF_FINAL$MyHIVTF_FINAL)){
  if ((is.na(DF_FINAL$MyHIVTF_FINAL[m])==FALSE & DF_FINAL$MyHIVTF_FINAL[m]!=TRUE)){
    MyVerdict <- append(MyVerdict,"NonHIV")
  } else if (is.na(DF_FINAL$MyLenFinal[m])==FALSE & DF_FINAL$MyLenFinal[m]<8000){
    MyVerdict <- append(MyVerdict,"LargeDeletion")
  } else if (is.na(DF_FINAL$MyStrandBlastHXB2_FINAL[m])==FALSE & DF_FINAL$MyStrandBlastHXB2_FINAL[m]=="mix"){
    MyVerdict <- append(MyVerdict,"InternalInversion")
  } else if (is.na(DF_FINAL$MyStrandBlastHXB2_FINAL[m])==FALSE & DF_FINAL$MyStrandBlastHXB2_FINAL[m]=="plusScramble") {
    MyVerdict <- append(MyVerdict,"ScramblePlus")
  } else if (is.na(DF_FINAL$MyStrandBlastHXB2_FINAL[m])==FALSE & DF_FINAL$MyStrandBlastHXB2_FINAL[m]=="minusScramble") {
    MyVerdict <- append(MyVerdict,"ScrambleMinus")
  } else if (is.na(DF_FINAL$MyStrandBlastHXB2_FINAL[m])==FALSE & DF_FINAL$MyStrandBlastHXB2_FINAL[m]=="Check") {
    MyVerdict <- append(MyVerdict,"ScrambleCheck")
  } else if (is.na(DF_FINAL$MyHypermutTF[m])==FALSE & DF_FINAL$MyHypermutTF[m]=="TRUE"){
    MyVerdict <- append(MyVerdict,"Hypermut")
    #For everything without primers: InferredXXX#
  } else if (is.na(DF_FINAL$MyGAGPassFail[m])==FALSE & (DF_FINAL$MyPrimerTF[m]=="NoPrimer") ){
    #If gag has ATG, all GAG/POL/ENV has to pass to be inferred intact
    if (is.na(DF_FINAL$MyGAGPassFail[m])==FALSE & (DF_FINAL$MyGagATG[m]==TRUE)){
      if (is.na(DF_FINAL$MyGAGPassFail[m])==FALSE & (substring(DF_FINAL$MyGAGPassFail[m],1,6)=="Passed" & substring(DF_FINAL$MyPOLPassFail[m],1,6)=="Passed" & substring(DF_FINAL$MyINTPassFail[m],1,6)=="Passed" & substring(DF_FINAL$MyENVPassFail[m],1,6)=="Passed" & substring(DF_FINAL$MyGP41PassFail[m],1,6)=="Passed")){
        MyVerdict <- append(MyVerdict,"Inferred_Intact")
      } else if (is.na(DF_FINAL$MyGAGPassFail[m])==FALSE & (substring(DF_FINAL$MyGAGPassFail[m],1,6)=="Failed" | substring(DF_FINAL$MyPOLPassFail[m],1,6)=="Failed" | substring(DF_FINAL$MyINTPassFail[m],1,6)=="Failed" | substring(DF_FINAL$MyENVPassFail[m],1,6)=="Failed" | substring(DF_FINAL$MyGP41PassFail[m],1,6)=="Failed")) {
        MyVerdict <- append(MyVerdict,"Inferred_PrematureStopORInframeDEL")
      } else {
        MyVerdict <- append(MyVerdict,"Inferred_Check#1_NoPrimerGagATGNoInferredIntactNoInferredPrematureStop")
      }
      #if gag has no ATG
    } else if (is.na(DF_FINAL$MyGAGPassFail[m])==FALSE & (DF_FINAL$MyGagATG[m]==FALSE)) {
      if (is.na(DF_FINAL$MyGAGPassFail[m])==FALSE & (substring(DF_FINAL$MyGAGPassFail[m],1,6)=="Passed" & substring(DF_FINAL$MyPOLPassFail[m],1,6)=="Passed" & substring(DF_FINAL$MyINTPassFail[m],1,6)=="Passed" & substring(DF_FINAL$MyENVPassFail[m],1,6)=="Passed" & substring(DF_FINAL$MyGP41PassFail[m],1,6)=="Passed")){
        MyVerdict <- append(MyVerdict,"Inferred_Intact_GagNoATG")
      } else if (is.na(DF_FINAL$MyGAGPassFail[m])==FALSE & (substring(DF_FINAL$MyGAGPassFail[m],1,6)=="Failed" & substring(DF_FINAL$MyPOLPassFail[m],1,6)=="Passed" & substring(DF_FINAL$MyINTPassFail[m],1,6)=="Passed" & substring(DF_FINAL$MyENVPassFail[m],1,6)=="Passed" & substring(DF_FINAL$MyGP41PassFail[m],1,6)=="Passed") & (as.numeric(as.character(DF_FINAL$MyPsiDELrelative[m])) + as.numeric(as.character(DF_FINAL$MyPsiINrelatiive[m])))>=15 ) {
        MyVerdict <- append(MyVerdict,"Inferred_Intact_NoGag")
      } else if (is.na(DF_FINAL$MyGAGPassFail[m])==FALSE & (substring(DF_FINAL$MyGAGPassFail[m],1,6)=="Failed" & (substring(DF_FINAL$MyPOLPassFail[m],1,6)=="Passed" & substring(DF_FINAL$MyINTPassFail[m],1,6)=="Passed" & substring(DF_FINAL$MyENVPassFail[m],1,6)=="Passed" & substring(DF_FINAL$MyGP41PassFail[m],1,6)=="Passed"))) {
        MyVerdict <- append(MyVerdict,"Inferred_PrematureStopORInframeDEL_GagNoATGandFailed")
      } else if (is.na(DF_FINAL$MyGAGPassFail[m])==FALSE & (substring(DF_FINAL$MyPOLPassFail[m],1,6)=="Failed" | substring(DF_FINAL$MyINTPassFail[m],1,6)=="Failed" | substring(DF_FINAL$MyENVPassFail[m],1,6)=="Failed" | substring(DF_FINAL$MyGP41PassFail[m],1,6)=="Failed")) {
        MyVerdict <- append(MyVerdict,"Inferred_PrematureStopORInframeDEL_GagNoATG")
      } else {
        MyVerdict <- append(MyVerdict,"Inferred_Check#2_NoPrimerGagNoATGNoInferredIntactNoInferredPrematureStop")
      }
    } else {
      MyVerdict <- append(MyVerdict,"Inferred_Check#3_NoPrimer")
    }
    #For everything with primers #
  } else if (is.na(DF_FINAL$MyGAGPassFail[m])==FALSE & (DF_FINAL$MyPrimerTF[m]=="YesPrimer") ){
    #if gag has ATG
    if (is.na(DF_FINAL$MyGAGPassFail[m])==FALSE & (DF_FINAL$MyGagATG[m]==TRUE)) {
      if (is.na(DF_FINAL$MyGAGPassFail[m])==FALSE & (substring(DF_FINAL$MyGAGPassFail[m],1,6)=="Failed" | substring(DF_FINAL$MyPOLPassFail[m],1,6)=="Failed" | substring(DF_FINAL$MyINTPassFail[m],1,6)=="Failed" | substring(DF_FINAL$MyENVPassFail[m],1,6)=="Failed" | substring(DF_FINAL$MyGP41PassFail[m],1,6)=="Failed")) {
        MyVerdict <- append(MyVerdict,"PrematureStop_OR_AAtooLong_OR_AAtooShort")
      } else if (is.na(DF_FINAL$MyGAGPassFail[m])==FALSE & (substring(DF_FINAL$MyGAGPassFail[m],1,6)=="Passed" & substring(DF_FINAL$MyPOLPassFail[m],1,6)=="Passed" & substring(DF_FINAL$MyINTPassFail[m],1,6)=="Passed" & substring(DF_FINAL$MyENVPassFail[m],1,6)=="Passed" & substring(DF_FINAL$MyGP41PassFail[m],1,6)=="Passed")) {
        if ((as.numeric(as.character(DF_FINAL$MyPsiDELrelative[m])) + as.numeric(as.character(DF_FINAL$MyPsiINrelatiive[m])))>=15) {
          MyVerdict <- append(MyVerdict,"5DEFECT")
        } else if ( (as.numeric(as.character(DF_FINAL$MyPsiDELrelative[m])) + as.numeric(as.character(DF_FINAL$MyPsiINrelatiive[m])))<15 ) {
          MyVerdict <- append(MyVerdict,"Intact")
        } else {
          MyVerdict <- append(MyVerdict,"Check#4_YesPrimerYesGagATGNoPrematureStopNot5DEFECTNotIntact")
        }
      } else {
        MyVerdict <- append(MyVerdict,"Check#5_YesPrimerYesGagATG")
      }
      #if gag has no ATG
    } else if (is.na(DF_FINAL$MyGAGPassFail[m])==FALSE & (DF_FINAL$MyGagATG[m]==FALSE) ) {
      if (is.na(DF_FINAL$MyGAGPassFail[m])==FALSE & (substring(DF_FINAL$MyPOLPassFail[m],1,6)=="Failed" | substring(DF_FINAL$MyINTPassFail[m],1,6)=="Failed" | substring(DF_FINAL$MyENVPassFail[m],1,6)=="Failed" | substring(DF_FINAL$MyGP41PassFail[m],1,6)=="Failed") & (as.numeric(as.character(DF_FINAL$MyPsiDELrelative[m])) + as.numeric(as.character(DF_FINAL$MyPsiINrelatiive[m])))>=15) {
        MyVerdict <- append(MyVerdict,"5DFECT_IntoGag")
      } else if (is.na(DF_FINAL$MyGAGPassFail[m])==FALSE & (substring(DF_FINAL$MyPOLPassFail[m],1,6)=="Failed" | substring(DF_FINAL$MyINTPassFail[m],1,6)=="Failed" | substring(DF_FINAL$MyENVPassFail[m],1,6)=="Failed" | substring(DF_FINAL$MyGP41PassFail[m],1,6)=="Failed")) {
        MyVerdict <- append(MyVerdict,"PrematureStop_OR_AAtooLong_OR_AAtooShort_GagNoATG")
      } else if (is.na(DF_FINAL$MyGAGPassFail[m])==FALSE & (substring(DF_FINAL$MyPOLPassFail[m],1,6)=="Passed" & substring(DF_FINAL$MyINTPassFail[m],1,6)=="Passed" & substring(DF_FINAL$MyENVPassFail[m],1,6)=="Passed" & substring(DF_FINAL$MyGP41PassFail[m],1,6)=="Passed")) {
        if ( (substring(DF_FINAL$MyGAGPassFail[m],1,6)=="Passed") & (as.numeric(as.character(DF_FINAL$MyPsiDELrelative[m])) + as.numeric(as.character(DF_FINAL$MyPsiINrelatiive[m])))>=15 ) {
          MyVerdict <- append(MyVerdict,"5DEFECT_GagNoATGGagPassed")
        } else if ( (substring(DF_FINAL$MyGAGPassFail[m],1,6)=="Passed") & (as.numeric(as.character(DF_FINAL$MyPsiDELrelative[m])) + as.numeric(as.character(DF_FINAL$MyPsiINrelatiive[m])))<15 ) {
          MyVerdict <- append(MyVerdict,"Intact_GagNoATG")
        } else if ((substring(DF_FINAL$MyGAGPassFail[m],1,6)=="Failed") & (as.numeric(as.character(DF_FINAL$MyPsiDELrelative[m])) + as.numeric(as.character(DF_FINAL$MyPsiINrelatiive[m])))>=15 ) {
          MyVerdict <- append(MyVerdict,"5DEFECT_GagNoATGGagFailed")
        } else if ((substring(DF_FINAL$MyGAGPassFail[m],1,6)=="Failed") & (as.numeric(as.character(DF_FINAL$MyPsiDELrelative[m])) + as.numeric(as.character(DF_FINAL$MyPsiINrelatiive[m])))<15 ) {
          MyVerdict <- append(MyVerdict,"Check#6a_YesPrimerGagNoATGGagFailed<15bp")
        } else {
          MyVerdict <- append(MyVerdict,"Check#6_YesPrimerGagNoATGNo5DEFECTNotIntact")
        }
      } else {
        MyVerdict <- append(MyVerdict,"Check#7_YesPrimerGagNoATG")
      }
    } else {
      MyVerdict <- append(MyVerdict,"Check#8_YesPrimer")
    }
  } else {
    MyVerdict <- append(MyVerdict,"Check#9_NoCategory")
  }
}

#Look for "Check" to find logic errors
#length(MyVerdict)
#View(as.matrix(MyVerdict))
#View(DF_FINAL[671,])

DF_FINAL <- cbind(DF_FINAL,MyVerdict)

#View(DF_FINAL)

#################################################################

### Step 05. Merge QC Output ###

if (any(list.files()=="RAW_QC") == FALSE) {
  #Skip this# 
} else {
  #merge the output from QCPipeline
  MyQCFile <- "Output_QC_ver05.csv"
  MyQCDF_RAW <- read.csv(MyQCFile)
  
  #Step05a. Create linking UniqID_Contig so that it will match with DF_FINAL
  MyUniqSampleContigID_Len <- mapply(paste,MyQCDF_RAW$Output01b_CurrentSample,"_",MyQCDF_RAW$Output02a1_CurrentContigSimple,"_",MyQCDF_RAW$Output04b_TotalContigLength,sep="") #n1583
  MyQCDF_RAW <- cbind(MyQCDF_RAW,MyUniqSampleContigID_Len)
  
  #Merge
  DF_FINAL <- merge(DF_FINAL,MyQCDF_RAW,by="MyUniqSampleContigID_Len",all=TRUE) 
}

##################################################################

## When all looks ready ##

write.csv(DF_FINAL,paste("Output_MyBigSummary_DF_FINAL.csv",sep=""))

##################################################################

## Copy files into final folders ##

MyResultFiles_Temp <- c("Output_Blastn_HXB2MEGA28_tabdelim.txt"
  ,"Output_Concatenate_RC.csv"
  ,"Output_Concatenate_RC.fasta"
  ,"Output_Concatenate.csv"
  ,"Output_Concatenate.fasta"
  ,"Output_HIV_Over8000_toHypermut_ver05fixscramble.fasta"
  ,"Output_HIV_Over8000_toHypermut_ver05fixscramble.fasta_Raligned.fasta_GLtrimmed_psigag.fasta"
  ,"Output_HIV_Over8000_toHypermut_ver05fixscramble.fasta_Raligned.fasta_GLtrimmed.fasta"
  ,"Output_HIVSeqinR_ver02.6_muscle_nuc_log.txt"
  ,"Output_HIVSeqinR_ver02.6_MyWorkingAln_AA_muscleAlign_oneb4stop.fasta"
  ,"Output_HIVSeqinR_ver02.6_MyWorkingAln_Nuc_muscleAlign_full_untrimmed.fasta"
  ,"Output_Hypermut01_REFCategory.txt"
  ,"Output_Hypermut02_MyFASTA_Reverted.fasta"
  ,"Output_Hypermut03_MyDF_Summary.csv"
  ,"Output_MyBigSummary_DF_FINAL_half.csv"
  ,"Output_MyHIVSeqinR_Pivot.csv"
)
MyResultFiles_Temp_1 <- list.files("./ToHIVSeqinR")
MyResultFiles_Temp_2 <- list.files("./ToHIVSeqinR_DB")

MyResultFiles_Final <- c("Output_HIVSeqinR_ver02.6_mydf.csv"
  ,"Output_HIVSeqinR_ver02.6_MyWorkingAln_Nuc_muscleAlign_trimmed.fasta"
  ,"Output_HIVSeqinR_ver02.6_MyWorkingAln_NUC2AA_muscleAlign_trimmed.fasta"
  ,"Output_HIVSeqinR_ver02.6.pdf"
  ,"Output_MyBigSummary_DF_FINAL.csv"
  ,"Output_HIV_Over8000_toHypermut_ver05fixscramble.fasta_Raligned.fasta"
)

#Temp files
dir.create("./Results_Intermediate")
dir.create("./Results_Intermediate/ToHIVSeqinR")
dir.create("./Results_Intermediate/ToHIVSeqinR_DB")
#Copy
file.copy(MyResultFiles_Temp,"./Results_Intermediate")
file.copy(paste("./ToHIVSeqinR/",MyResultFiles_Temp_1,sep=""),"./Results_Intermediate/ToHIVSeqinR")
file.copy(paste("./ToHIVSeqinR_DB/",MyResultFiles_Temp_2,sep=""),"./Results_Intermediate/ToHIVSeqinR_DB")
#Delete
file.remove(MyResultFiles_Temp)
unlink("./ToHIVSeqinR/",recursive=TRUE)
unlink("./ToHIVSeqinR_DB/",recursive=TRUE)

#QC file
if (any(list.files()=="Output_QC_ver05.csv")==TRUE){
  file.copy("Output_QC_ver05.csv","./Results_Intermediate")
  file.remove("Output_QC_ver05.csv")
} else {
  #do nothing
}

#Final files
dir.create("./Results_Final")
file.copy(MyResultFiles_Final,"./Results_Final")
file.remove(MyResultFiles_Final)

##################################################################

EndSysTime <- Sys.time()
print(paste("HIVSeqinR started running",StartSysTime))
print(paste("HIVSeqinR finished running",EndSysTime))

