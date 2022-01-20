#Balaji Sundararaman
#10-20-21
#Fixed Alisa's wrong math for endogenous calculation
#adopted from Alisa V's version of "endDNA_by_AK_complex_AV.R"
#Using Picard MarkDup %_Duplication as followed by Jonas and others

args = commandArgs(trailingOnly=TRUE)
if (length(args) == 0){
  cat("Syntax: Rscript shortReads_QCmetrics.R [path to seqprep files] [path to flagstats files] [path to lgdist folders]\n")
  cat("Example: Rscript shortReads_QCmetrics.R ./SeqPrep_output/ ./BWA_analyses/ ./MapDamage_output/\n")
  quit()
}
#directory paths
seqPrepPath = args[1]
bwaPath = args[2]
fragLenPaths <- list.dirs(path = args[3], full.names = TRUE, recursive = TRUE)

#files
seqPrepFiles <- list.files(path=seqPrepPath, pattern="*SeqPrep_output.txt", full.names=T, recursive=FALSE)
allReadFSFiles <- list.files(path=bwaPath, pattern="*all_reads.complexity_filtered_EquCab2.sorted.flagstats.txt", full.names=T, recursive=FALSE)

cat ("checking number of found files\n")
#cat(seqPrepFiles)
#cat(allReadFSFiles)
#cat(paste0(length(seqPrepFiles), length(allReadFSFiles), length(fragLenPaths)))

if (length(seqPrepFiles) == length(allReadFSFiles) && length(allReadFSFiles) == length(fragLenPaths)-1)
{
  cat (paste0("OK seqpreps:", length(seqPrepFiles), ", AllReadFlagStats:", length(allReadFSFiles), ", "fraglenfolders:", length(fragLenPaths)-1, "\n"))
  #create output matrix
  output <- matrix(ncol=21, nrow=length(seqPrepFiles))
  #go through each sample
  i = 0
  for (f in seqPrepFiles)
  {
    i = i + 1
    #obtain sample name
    fileName <- unlist(strsplit(f, "/"))
    sampleName <- unlist(strsplit(unlist(fileName[length(fileName)]), "_SeqPrep_output.txt"))
    cat(paste0("sample name:", sampleName, "\n"))

    #just to make sure on correct pairing of flagstats and seqprep files:
    cat(paste0("SeqPrep filename:", seqPrepFiles[grep(sampleName, seqPrepFiles)], "\n"))
    cat(paste0("All reads flagstats filename:", allReadFSFiles[grep(sampleName, allReadFSFiles)], "\n"))
    
    #look up files with sample name and read them
    seqPrep <- read.table(seqPrepFiles[grep(sampleName, seqPrepFiles)], fill=TRUE, header = FALSE, sep = ":")
    allReadFlagstat <- read.table(allReadFSFiles[grep(sampleName, allReadFSFiles)], fill = TRUE, header = FALSE, sep = "+")
    
    #convert flagstats first column to numeric since the last row has chars...
    allReadFlagstat$V1 <- as.numeric(as.character(allReadFlagstat$V1))
    
    #alisa's lgdistribution script
    fragLenDistFilename <- paste0(fragLenPaths[grep(sampleName, fragLenPaths)],"/lgdistribution.txt")
    dat <- read.table(fragLenDistFilename, header=T)
    mean = sum(dat$Length * dat$Occurences) / sum(dat$Occurences)
    SD = sqrt(sum((dat$Length - mean)**2 * dat$Occurences) / (sum(dat$Occurences)-1))

    #balaji's endogenous math corrected version of Alisa's script
    #Read numbers
    RawReadPairs <- seqPrep$V2[2]
    Discarded <- seqPrep$V2[5]
    ReadPairsUsed <- (RawReadPairs-Discarded)
    ReadPairsMerged <- seqPrep$V2[3]
    ReadsMergedPCT <-  ReadPairsMerged/ReadPairsUsed
    MergedCplexFilterReads <- seqPrep$V2[7] # Complexity filtered merged read count is appended to the SeqPrep stat.txt by the reads processing bash script
    MergCplexFilterPCT <- MergedCplexFilterReads/ReadPairsMerged
    UmergedCplexFilterReadPairs <- seqPrep$V2[8] # Complexity filtered unmerged read count is appended to the SeqPrep stat.txt by the reads processing bash script
    UnmergeCplexFilterPCT <- UmergedCplexFilterReadPairs/(ReadPairsUsed-ReadPairsMerged)
    TotalCplexFilterReads <- MergedCplexFilterReads + (2*UmergedCplexFilterReadPairs)

    #mapping numbers
    #AllMapped <- allReadFlagstat$V1[3]
    #ProperlyPaired <- allReadFlagstat$V1[7]
    #Singletons <- allReadFlagstat$V1[9]
    UnmergMapped <- allReadFlagstat$V1[7] #Properly_Paired_mapped_reads
    UnmergMapPCT <- UnmergMapped/(2*UmergedCplexFilterReadPairs)
    MergMapped <- allReadFlagstat$V1[12] # Merged mapped read count is appended to the 'all_read_flagstat.txt' by the read processing bash script.
    MergMapPCT <- MergMapped/MergedCplexFilterReads
    AllMapped <- MergMapped+UnmergMapped
    Endog <- AllMapped/TotalCplexFilterReads
    Dups <- allReadFlagstat$V1[13] #Picard's MarkDup PERCENT_DUPLICATION is appended to the 'all_read_flagstat.txt' by the read processing bash script.
    ComplexityPCT <- 1-allReadFlagstat$V1[13] #Picard's MarkDup PERCENT_DUPLICATION is appended to the 'all_read_flagstat.txt' by the read processing bash script.
    
    cat(paste0("Endog:", Endog,"\n"))
    #put results into the matrix
    output[i,] <- c(sampleName, RawReadPairs, Discarded, ReadPairsUsed, ReadPairsMerged, ReadsMergedPCT, MergedCplexFilterReads, MergCplexFilterPCT, UmergedCplexFilterReadPairs, UnmergeCplexFilterPCT, TotalCplexFilterReads, MergMapped, MergMapPCT, UnmergMapped, UnmergMapPCT, AllMapped, Endog, Dups, ComplexityPCT, mean, SD)
  }
} else {
  cat ("Error, different number of files/folders found\n")
  #cat (paste0("seqpreps:", length(filessp), ", flagstats:", length(filesfs),", lgdistfolders:", length(lgdistpath)-1, "\n"))
  quit()
}

#transform matrix to dataframe
output <- data.frame(output)

#label columns
colnames(output) <- c("Sample", "Raw_ReadPairs", "ReadPairs_Discarded", "ReadPairs_Used", "ReadPairs_Merged", "ReadPairs_Merged_PCT", "Merged_Complexity_Filtered_Reads", "Merged_Complexity_Filtered_Reads_PCT", "Unmerged_Complexity_Filtered_ReadPairs", "Unmerged_Complexity_Filtered_ReadPairs_PCT", "Total_Complexity_Filtered_Reads", "Merged_Reads_Mapped", "Merged_Reads_Mapped_PCT", "Unmerged_Reads_Mapped", "Unmerged_Reads_Mapped_PCT", "Mapped_Reads_Total", "Endogeneous_PCT", "Picard_Duplication_PCT", "Complexity_PCT", "Frag_Length_Mean", "Frag_Length_SD")

#wtite output file
wd <- basename(getwd())
write.csv(output, file=paste(wd, "read_mapping_metrics.csv", sep="."))

cat("done\n")
