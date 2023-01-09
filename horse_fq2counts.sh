#!/usr/local/bin
#Balaji Sundararaman
#03-30-2022

#usage: bash horse_fq2counts.sh sample_name probe_name raw_R1.fq raw_R2.fq

mkdir $1
cd $1
mkdir aln
mkdir seq
mkdir captureCounts

SEQPREP_LOCATION=/path/to/SeqPrep2/      # path for SeqPrep
REFERENCE=/path/to/equcab2_reference/whole_genome_horse.fa # Remember to include the full path for EquCab2 reference

if [ $2 == "CNER" ]; then
    PROBEbed=/path/to/ecab2_23999probes_23771SNPs_CNERs.bed
    PROBElist=/path/to/ecab2_23999probes_23771SNPs_CNERs_interval_list
elif [ $2 == "Arbor" ]; then
    PROBEbed=/path/to/ecab2_Arbor_59292baits_22617SNPs_myBaits.bed
    PROBElist=/path/to/ecab2_Arbor_59292baits_22617SNPs_myBaits_interval_list
fi

SNPbed=/path/to/ecab2_common_22619SNPs.bed
SNPlist=/path/to/ecab2_common_22619SNPs_interval_list
CNERsnps=/path/to/ecab2_23999probes_23771SNPs.bed
UpdwBED=/path/to/ecab2_common_22619SNPs_100updown.bed

if [ ! -f ./seq/$1.mrgReads.fq.gz ]
then
    echo "Running SeqPrep for"$1
    ### SeqPrep -- removing adapters and merging reads. Overlap (-o) is set to 15, minimum length (-L) is set to 25, minimum quality (-q) is set to 15, mismatch fraction set to 0.05.
    $SEQPREP_LOCATION/SeqPrep2 -f $3 -r $4 -1 ./seq/$1.R1.fq.gz -2 ./seq/$1.R2.fq.gz -s ./seq/$1.mrgReads.fq.gz -q 15 -L 30 -o 15 -m 0.05 \
			       -A AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -B AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -C ATCTCGTATGCCGTCTTCTGCTTG -D GATCTCGGTGGTCGCCGTATCATT >& ./seq/$1.SeqPrep_output.txt
    wait
else
    echo "SeqPrep already done for "$1
fi

### BWA -- aligning reads to a reference sequence
if [ ! -f ./aln/$1.all_reads.mrkd.s.bam.bai ]
then
    echo "BWA aligning"$1
    bwa aln -l 16500 -n 0.01 -o 2 -t 16 $REFERENCE ./seq/$1.mrgReads.fq.gz > ./aln/$1.mrgReads.sai
    wait
    bwa samse $REFERENCE ./aln/$1.mrgReads.sai ./seq/$1.mrgReads.fq.gz | samtools view -Sbh - | samtools sort - -@ 8 -o ./aln/$1.mrgReads.s.bam
    wait
    bwa aln -l 16500 -n 0.01 -o 2 -t 16 $REFERENCE ./seq/$1.R1.fq.gz > ./aln/$1.R1.sai
    wait
    bwa aln -l 16500 -n 0.01 -o 2 -t 16 $REFERENCE ./seq/$1.R2.fq.gz > ./aln/$1.R2.sai
    wait
    bwa sampe $REFERENCE ./aln/$1.R1.sai ./aln/$1.R2.sai ./seq/$1.R1.fq.gz ./seq/$1.R2.fq.gz | samtools view -Sbh - | samtools sort - -@ 8 -o ./aln/$1.unmrgReads.s.bam
    wait
    samtools index ./aln/$1.mrgReads.s.bam
    wait
    samtools index ./aln/$1.unmrgReads.s.bam
    wait
    samtools merge ./aln/$1.all_reads.bam ./aln/$1.mrgReads.s.bam ./aln/$1.unmrgReads.s.bam
    wait
    samtools sort -@ 8 ./aln/$1.all_reads.bam -o ./aln/$1.all_reads.s.bam
    wait
    samtools index ./aln/$1.all_reads.s.bam
    wait
    java -jar /path/to/picard.jar MarkDuplicates I=./aln/$1.all_reads.s.bam O=./aln/$1.all_reads.mrkd.bam M=./aln/$1.markdup_metrics.txt REMOVE_DUPLICATES=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=800 VALIDATION_STRINGENCY=LENIENT
    wait
    samtools sort -@ 8 ./aln/$1.all_reads.mrkd.bam -o ./aln/$1.all_reads.mrkd.s.bam
    wait
    samtools index ./aln/$1.all_reads.mrkd.s.bam
    wait
else
    echo "All alignment files found for "$1
fi
wait

###Capture counting
if [ ! -f ./captureCounts/$1.$2.captureCounts.csv ]; then
    echo "Counting $1 for $2 capture"
    echo "Sample,"$1"_"$2 > ./captureCounts/$1.$2.captureCounts.csv
    zcat $3 | wc -l | awk '{print"Raw_Reads,"$1/2}' >> ./captureCounts/$1.$2.captureCounts.csv
    samtools view -c ./aln/$1.mrgReads.s.bam | awk '{print"MergedReads_fromBam,"$1}' >> ./captureCounts/$1.$2.captureCounts.csv
    samtools view -c -F0x4 ./aln/$1.mrgReads.s.bam | awk '{print"MergedReads_totalMapped,"$1}' >> ./captureCounts/$1.$2.captureCounts.csv
    samtools view -c ./aln/$1.unmrgReads.s.bam | awk '{print"UnmergedReads_fromBam,"$1}' >> ./captureCounts/$1.$2.captureCounts.csv
    samtools view -c -F0x4 ./aln/$1.unmrgReads.s.bam | awk '{print"UnmergedReads_totalMapped,"$1}' >> ./captureCounts/$1.$2.captureCounts.csv
    samtools view -c -F0x4 -f3 ./aln/$1.unmrgReads.s.bam | awk '{print"UnmergedReads_properlyMapped,"$1}' >> ./captureCounts/$1.$2.captureCounts.csv
    samtools view -c ./aln/$1.all_reads.s.bam | awk '{print"AllReads_fromBam,"$1}' >> ./captureCounts/$1.$2.captureCounts.csv
    samtools view -c -F0x4 ./aln/$1.all_reads.s.bam | awk '{print"AllReads_totalMapped,"$1}' >> ./captureCounts/$1.$2.captureCounts.csv
    samtools view -c ./aln/$1.all_reads.mrkd.s.bam | awk '{print"UniqueReads_fromBam,"$1}' >> ./captureCounts/$1.$2.captureCounts.csv
    samtools view -c -F0x4 ./aln/$1.all_reads.mrkd.s.bam | awk '{print"AllReads_uniqueMapped,"$1}' >> ./captureCounts/$1.$2.captureCounts.csv
    grep '^LIB' -A1 ./aln/$1.markdup_metrics.txt | grep -v '^LIB' | awk '{print "Picard_MarkDup_PCT,"$10}' >> ./captureCounts/$1.$2.captureCounts.csv
    wait
    
    bedtools multicov -bams ./aln/$1.mrgReads.s.bam -bed $SNPbed | awk '{sum+=$7} END{print "MergedReads_SNPcounts,"sum}' >> ./captureCounts/$1.$2.captureCounts.csv
    bedtools multicov -bams ./aln/$1.unmrgReads.s.bam -bed $SNPbed | awk '{sum+=$7} END{print "UnMergedReads_SNPcounts,"sum}' >> ./captureCounts/$1.$2.captureCounts.csv
    
    bedtools multicov -bams ./aln/$1.mrgReads.s.bam -bed $PROBEbed | awk '{sum+=$7} END{print "MergedReads_ProbeCounts,"sum}' >> ./captureCounts/$1.$2.captureCounts.csv
    bedtools multicov -bams ./aln/$1.unmrgReads.s.bam -bed $PROBEbed | awk '{sum+=$7} END{print "UnMergedReads_ProbeCounts,"sum}' >> ./captureCounts/$1.$2.captureCounts.csv

    bedtools multicov -bams ./aln/$1.all_reads.s.bam -bed $PROBEbed | awk '{sum+=$7} END{print "AllReads_ProbeCounts,"sum}' >> ./captureCounts/$1.$2.captureCounts.csv
    bedtools multicov -bams ./aln/$1.all_reads.mrkd.s.bam -bed $PROBEbed | awk '{sum+=$7} END{print "UniqueReads_ProbeCounts,"sum}' >> ./captureCounts/$1.$2.captureCounts.csv
    wait
    
    bedtools multicov -bams ./aln/$1.all_reads.s.bam -bed $SNPbed > ./$1.$2.AllReads_22619SNPs_depth.tsv
    bedtools multicov -bams ./aln/$1.all_reads.mrkd.s.bam -bed $SNPbed > ./$1.$2.UniqueReads_22619SNPs_depth.tsv
    wait
    
    awk '{sum+=$7} END{print "AllReads_SNPcounts,"sum}' ./$1.$2.AllReads_22619SNPs_depth.tsv >> ./captureCounts/$1.$2.captureCounts.csv
    awk '{sum+=$7} END{print "UniqueReads_SNPcounts,"sum}' ./$1.$2.UniqueReads_22619SNPs_depth.tsv  >> ./captureCounts/$1.$2.captureCounts.csv
    awk ' BEGIN {count=0;}  { if ($7 == 0) count+=1} END {print "AllReads_zeroCovSNPs,"count}' ./$1.$2.AllReads_22619SNPs_depth.tsv >> ./captureCounts/$1.$2.captureCounts.csv
    awk ' BEGIN {count=0;}  { if ($7 == 0) count+=1} END {print "UniqueReads_zeroCovSNPs,"count}' ./$1.$2.UniqueReads_22619SNPs_depth.tsv  >> ./captureCounts/$1.$2.captureCounts.csv
    wait
    
    for p in {1..5..1}
    do
	awk -v d=$p ' BEGIN {count=0;}  { if ($7 >= d) count+=1} END {print "AllReads_SNPs_" d "X,"count}' ./$1.$2.AllReads_22619SNPs_depth.tsv >> ./captureCounts/$1.$2.captureCounts.csv
	awk -v d=$p ' BEGIN {count=0;}  { if ($7 >= d) count+=1} END {print "UniqueReads_SNPs_" d "X,"count}' ./$1.$2.UniqueReads_22619SNPs_depth.tsv  >> ./captureCounts/$1.$2.captureCounts.csv
    done
    wait

    for q in {10..50..10}
    do
	awk -v t=$q ' BEGIN {count=0;}  { if ($7 >= t) count+=1} END {print "AllReads_SNPs_" t "X,"count}' ./$1.$2.AllReads_22619SNPs_depth.tsv >> ./captureCounts/$1.$2.captureCounts.csv
	awk -v t=$q ' BEGIN {count=0;}  { if ($7 >= t) count+=1} END {print "UniqueReads_SNPs_" t "X,"count}' ./$1.$2.UniqueReads_22619SNPs_depth.tsv  >> ./captureCounts/$1.$2.captureCounts.csv
    done
    wait

    awk ' BEGIN {count=0;} {count+=$7} END {print "AllReads_SNPmeanCov,"count/NR}' ./$1.$2.AllReads_22619SNPs_depth.tsv >> ./captureCounts/$1.$2.captureCounts.csv
    awk ' BEGIN {count=0;} {count+=$7} END {print "UniqueReads_SNPmeanCov,"count/NR}' ./$1.$2.UniqueReads_22619SNPs_depth.tsv  >> ./captureCounts/$1.$2.captureCounts.csv
    wait
    
    cut -f7 ./$1.$2.AllReads_22619SNPs_depth.tsv | sort -n | awk '{a[NR]=$1} END {if(NR%2==1) print "AllReads_SNPmedianCov,"a[(NR+1)/2]; else print "AllReads_SNPmedianCov,"(a[NR/2]+a[NR/2+1])/2}' >> ./captureCounts/$1.$2.captureCounts.csv
    cut -f7 ./$1.$2.UniqueReads_22619SNPs_depth.tsv | sort -n | awk '{a[NR]=$1} END {if(NR%2==1) print "UnqiueReads_SNPmedianCov,"a[(NR+1)/2]; else print "UnqiueReads_SNPmedianCov,"(a[NR/2]+a[NR/2+1])/2}' >> ./captureCounts/$1.$2.captureCounts.csv
    wait
    
    java -jar /path/to/picard.jar CollectHsMetrics I=./aln/$1.all_reads.mrkd.s.bam O=$1.temp1.tsv R=$REFERENCE BAIT_INTERVALS=$PROBElist TARGET_INTERVALS=$PROBElist PER_TARGET_COVERAGE=$1.$2.perProbe_coverage.tsv VALIDATION_STRINGENCY=SILENT
    wait
    java -jar /path/to/picard.jar CollectHsMetrics I=./aln/$1.all_reads.mrkd.s.bam O=$1.temp2.tsv R=$REFERENCE BAIT_INTERVALS=$PROBElist TARGET_INTERVALS=$SNPlist VALIDATION_STRINGENCY=SILENT
    wait
    
    grep -A 1 '^BAIT' $1.temp2.tsv | sed 's:\t:,:g' | sed "s:$:$1:g" > $1.$2.picard.csv
    wait
fi
wait

echo "Plotting SNP coverages"
if [ $2 == "CNER" ]; then
    bedtools multicov -bams ./aln/$1.all_reads.mrkd.s.bam -bed $CNERsnps > ./$1.$2.UniqueReads_23771SNPs_depth.tsv
    wait
    python3 /path/to/cner_len_SNPdepth.py ./$1.$2.UniqueReads_23771SNPs_depth.tsv $1"_"$2
fi
wait
bedtools coverage -a $UpdwBED -b ./aln/$1.all_reads.mrkd.s.bam -d > ./$1.$2.22619SNPs_100updown.tsv
wait
python3 /path/to/gc_cover_plot.py $1.$2.perProbe_coverage.tsv $1"_"$2
python3 /path/to/plot_SNP_100updown_coverage.py ./$1.$2.22619SNPs_100updown.tsv $1"_"$2
python3 /path/to/targetSNP_covDepth_PCTplot.py ./$1.$2.UniqueReads_22619SNPs_depth.tsv $1"_"$2
wait
echo "Done analysing $1 for $2 captures"
mkdir figures
mv *png ./figures
mv *tsv ./captureCounts
wait
echo "Removing temp and intermediate files"
rm ./aln/$1.*sai
rm ./aln/$1.*all_reads.bam
rm ./aln/$1.*all_reads.mrkd.bam
rm $1.temp1.tsv && $1.temp2.tsv
wait
echo "All done!"
