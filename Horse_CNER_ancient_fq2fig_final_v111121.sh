#!/usr/local/bin
#Balaji Sundararaman
#11-11-2021
#Usage bash Horse_CNER_ancient_fq2fig_final_v111121.sh <directory_to_ouput> <sample_text_file>

cd $1
mkdir BWA_analyses
mkdir SeqPrep_output
mkdir MapDamage_output

#SEQPREP2 envelopes
SEQPREP_MIN_LENGTH=27        # Removes unmappably short reads.
SEQPREP_OVERLAP=15           # Allows for confident merging of reads. Can be reduced to 10 if needed.
SEQPREP_LOCATION=/path/to/SeqPrep2-master    # Refer: https://github.com/jeizenga/SeqPrep2

# BWA and SAMtools envelopes
REFERENCE_SEQUENCE=/path/to/ecab2_reference/whole_genome_horse.fa 		# Remember to include the full path

# Prinseq envelopes
PRINSEQ_LITE=/path/to/prinseq-lite.pl 	# Refer: https://github.com/uwb-linux/prinseq
COMPLEXITY_METHOD=dust			# dust is the standard used by BLAST. The entropy method is an alternative, but is not widely used.
COMPLEXITY_THRESHOLD=7			# Pretty low, but follows the PRINSEQ_LITE manual. Recommends 70 if using entropy.
COMBINE_PAIRED_END_READS=/path/to/combinePairedEndReads.pl	# Refer: https://github.com/uwb-linux/prinseq
SPLIT_PAIRED_END_READS=/path/to/splitPairedEndReads.pl		# Refer: https://github.com/uwb-linux/prinseq

fn=$2
while read ln; do
    echo "$ln"
    RAW_R1=`echo $ln | cut -d ',' -f 1`
    RAW_R2=`echo $ln | cut -d ',' -f 2`
    SAMPLE=`echo $ln | cut -d ',' -f 3`

    if [ ! -f ./SeqPrep_output/${SAMPLE}_SeqPrep_output.txt ]
    then

### SeqPrep -- removing adapters and merging reads. Overlap (-o) is set to 15, minimum length (-l) is set to 27.
	${SEQPREP_LOCATION}/SeqPrep2 -f $RAW_R1 -r $RAW_R2 -1 ./SeqPrep_output/${SAMPLE}_R1_unmerged.fastq.gz -2 ./SeqPrep_output/${SAMPLE}_R2_unmerged.fastq.gz -q 15 -L ${SEQPREP_MIN_LENGTH} -A AGATCGGAAGAGCACACGTC -B AGATCGGAAGAGCGTCGTGT -s ./SeqPrep_output/${SAMPLE}_merged.fastq.gz -o ${SEQPREP_OVERLAP} -d 1 -C ATCTCGTATGCCGTCTTCTGCTTG -D GATCTCGGTGGTCGCCGTATCATT >& ./SeqPrep_output/${SAMPLE}_SeqPrep_output.txt
	gunzip ./SeqPrep_output/${SAMPLE}*merged.fastq.gz

### Filtering reads for low complexity sequences
	# Remove low complexity reads from merged files
	perl ${PRINSEQ_LITE} -fastq ./SeqPrep_output/${SAMPLE}_merged.fastq -out_good ./SeqPrep_output/${SAMPLE}_merged.complexity_filtered -out_bad null -lc_method ${COMPLEXITY_METHOD} -lc_threshold ${COMPLEXITY_THRESHOLD} -line_width 0

	# Remove low complexity reads from unmerged files
	perl ${COMBINE_PAIRED_END_READS} ./SeqPrep_output/${SAMPLE}_R1_unmerged.fastq ./SeqPrep_output/${SAMPLE}_R2_unmerged.fastq ./SeqPrep_output/${SAMPLE}_unmerged_combined.fastq
	wait
	perl ${PRINSEQ_LITE} -fastq ./SeqPrep_output/${SAMPLE}_unmerged_combined.fastq -out_good ./SeqPrep_output/${SAMPLE}_unmerged_combined.complexity_filtered -out_bad null -lc_method ${COMPLEXITY_METHOD} -lc_threshold ${COMPLEXITY_THRESHOLD} -line_width 0
	wait
	perl ${SPLIT_PAIRED_END_READS} ./SeqPrep_output/${SAMPLE}_unmerged_combined.complexity_filtered.fastq
	wait
	mv ./SeqPrep_output/${SAMPLE}_unmerged_combined.complexity_filtered.fastq_1 ./SeqPrep_output/${SAMPLE}_R1_unmerged.complexity_filtered.fastq
	wait
	mv ./SeqPrep_output/${SAMPLE}_unmerged_combined.complexity_filtered.fastq_2 ./SeqPrep_output/${SAMPLE}_R2_unmerged.complexity_filtered.fastq

	wc -l ./SeqPrep_output/${SAMPLE}_merged.complexity_filtered.fastq | awk '{print "PrinSeq_merged_complexity_filtered_reads:\t" $1/4}' >>./SeqPrep_output/${SAMPLE}_SeqPrep_output.txt
	wc -l ./SeqPrep_output/${SAMPLE}_R1_unmerged.complexity_filtered.fastq | awk '{print "PrinSeq_R1_unmerged_complexity_filtered_reads:\t" $1/4}' >>./SeqPrep_output/${SAMPLE}_SeqPrep_output.txt
    fi

    ### BWA -- aligning reads to a reference sequence
    if [ ! -f ./BWA_analyses/${SAMPLE}_all_reads.complexity_filtered_EquCab2.sorted.flagstats.txt ]
    then
	/usr/bin/bwa aln -l 16500 -n 0.01 -o 2 -t16 ${REFERENCE_SEQUENCE} ./SeqPrep_output/${SAMPLE}_merged.complexity_filtered.fastq > ./BWA_analyses/${SAMPLE}_merged.complexity_filtered_EquCab2.sai
	/usr/bin/bwa samse ${REFERENCE_SEQUENCE} ./BWA_analyses/${SAMPLE}_merged.complexity_filtered_EquCab2.sai ./SeqPrep_output/${SAMPLE}_merged.complexity_filtered.fastq | /usr/bin/samtools view -Sbh - | /usr/bin/samtools sort - ./BWA_analyses/${SAMPLE}_merged.complexity_filtered_EquCab2.sorted

	/usr/bin/bwa aln -l 16500 -n 0.01 -o 2 -t16 ${REFERENCE_SEQUENCE} ./SeqPrep_output/${SAMPLE}_R1_unmerged.complexity_filtered.fastq > ./BWA_analyses/${SAMPLE}_R1_unmerged.complexity_filtered_EquCab2.sai
	/usr/bin/bwa aln -l 16500 -n 0.01 -o 2 -t16 ${REFERENCE_SEQUENCE} ./SeqPrep_output/${SAMPLE}_R2_unmerged.complexity_filtered.fastq > ./BWA_analyses/${SAMPLE}_R2_unmerged.complexity_filtered_EquCab2.sai
	/usr/bin/bwa sampe ${REFERENCE_SEQUENCE} ./BWA_analyses/${SAMPLE}_R1_unmerged.complexity_filtered_EquCab2.sai ./BWA_analyses/${SAMPLE}_R2_unmerged.complexity_filtered_EquCab2.sai ./SeqPrep_output/${SAMPLE}_R1_unmerged.complexity_filtered.fastq ./SeqPrep_output/${SAMPLE}_R2_unmerged.complexity_filtered.fastq | /usr/bin/samtools view -Sbh - | /usr/bin/samtools sort - ./BWA_analyses/${SAMPLE}_unmerged.complexity_filtered_pe_EquCab2.sorted

	/usr/bin/samtools merge ./BWA_analyses/${SAMPLE}_all_reads.complexity_filtered_EquCab2.bam ./BWA_analyses/${SAMPLE}_merged.complexity_filtered_EquCab2.sorted.bam ./BWA_analyses/${SAMPLE}_unmerged.complexity_filtered_pe_EquCab2.sorted.bam

	/usr/bin/samtools sort ./BWA_analyses/${SAMPLE}_all_reads.complexity_filtered_EquCab2.bam ./BWA_analyses/${SAMPLE}_all_reads.complexity_filtered_EquCab2.sorted
	
	java -jar ~/bin/picard.jar MarkDuplicates -I ./BWA_analyses/${SAMPLE}_all_reads.complexity_filtered_EquCab2.sorted.bam -O ./BWA_analyses/${SAMPLE}_all_reads.complexity_filtered_EquCab2.sorted.rmdup.bam -METRICS_FILE ./BWA_analyses/${SAMPLE}_markdup_metrics.txt -REMOVE_DUPLICATES true -MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 800 -VALIDATION_STRINGENCY LENIENT

	/usr/bin/samtools index ./BWA_analyses/${SAMPLE}_all_reads.complexity_filtered_EquCab2.sorted.rmdup.bam

# Generate statistics like number mapped, duplicates, and number aligned to each chromosome
	/usr/bin/samtools flagstat ./BWA_analyses/${SAMPLE}_all_reads.complexity_filtered_EquCab2.sorted.bam > ./BWA_analyses/${SAMPLE}_all_reads.complexity_filtered_EquCab2.sorted.flagstats.txt

	/usr/bin/samtools view -F4 -c ./BWA_analyses/${SAMPLE}_merged.complexity_filtered_EquCab2.sorted.bam | awk '{print $1 " + Merged_mapped_reads"}' >> ./BWA_analyses/${SAMPLE}_all_reads.complexity_filtered_EquCab2.sorted.flagstats.txt

	grep '^LIB' -A1 ./BWA_analyses/${SAMPLE}_markdup_metrics.txt | grep -v '^LIB' | awk '{print $10 " + Picard_MarkDup_PCT_Duplication"}' >> ./BWA_analyses/${SAMPLE}_all_reads.complexity_filtered_EquCab2.sorted.flagstats.txt

### mapDamage to assess damage rates from all aligned and duplicate-removed data and to draw fragment length distributions of aligned data
	/soe/pheintzman/bin/mapDamage -i ./BWA_analyses/${SAMPLE}_all_reads.complexity_filtered_EquCab2.sorted.bam -r ${REFERENCE_SEQUENCE} --merge-reference-sequences -l 150 -d ./MapDamage_output/mapDamage_${SAMPLE} -y 0.5 -m 25 -t ${SAMPLE}
    fi

###CNER capture metrics using Picard tool
    if [ ! -f ${SAMPLE}.picard_temp.csv ]
    then
	bedtools multicov -q 20 -bams ./BWA_analyses/${SAMPLE}_all_reads.complexity_filtered_EquCab2.sorted.rmdup.bam -bed /projects/redser3-notbackedup/projects/bsundara/probes/ecab2/ecab2_23999baits_23771SNPs_targets.bed > ${SAMPLE}.all23771SNPs_depth.tsv

	bedtools coverage -a /projects/redser2/projects/bsundara/probes/ecab2/ecab2_all23771SNPs_100updown.bed -b ./BWA_analyses/${SAMPLE}_all_reads.complexity_filtered_EquCab2.sorted.rmdup.bam -d > ${SAMPLE}.all23771SNPs_100updown.tsv

	java -jar ~/bin/picard.jar CollectHsMetrics -I ./BWA_analyses/${SAMPLE}_all_reads.complexity_filtered_EquCab2.sorted.rmdup.bam -O ${SAMPLE}.temp1.tsv -R ${REFERENCE_SEQUENCE} -BAIT_INTERVALS /projects/redser3-notbackedup/projects/bsundara/probes/ecab2/ecab2_23999baits_23771SNPs_baits_interval_list -TARGET_INTERVALS /projects/redser3-notbackedup/projects/bsundara/probes/ecab2/ecab2_23999baits_23771SNPs_baits_interval_list -PER_TARGET_COVERAGE ${SAMPLE}.ecab2_perCNER.tsv -VALIDATION_STRINGENCY SILENT

	java -jar ~/bin/picard.jar CollectHsMetrics -I ./BWA_analyses/${SAMPLE}_all_reads.complexity_filtered_EquCab2.sorted.rmdup.bam -O ${SAMPLE}.temp2.tsv -R ${REFERENCE_SEQUENCE} -BAIT_INTERVALS /projects/redser3-notbackedup/projects/bsundara/probes/ecab2/ecab2_23999baits_23771SNPs_baits_interval_list -TARGET_INTERVALS /projects/redser3-notbackedup/projects/bsundara/probes/ecab2/ecab2_23999baits_23771SNPs_targets_interval_list -VALIDATION_STRINGENCY SILENT

	grep -A 1 '^BAIT' ${SAMPLE}.temp2.tsv | sed 's:\t:,:g' | sed "s:$:$SAMPLE:g" > ${SAMPLE}.picard_temp.csv #>> Picard_metrics.csv
	
	###Plot coverage metrics using custom py scripts
	python3 /path/to/gc_cover_plot_final.py ${SAMPLE}.ecab2_perCNER.tsv ${SAMPLE} >> ${SAMPLE}.pipe_log.txt
	python3 /path/to/cner_len_SNPdepth_final.py ${SAMPLE}.all23771SNPs_depth.tsv ${SAMPLE} >> ${SAMPLE}.pipe_log.txt
	python3 /path/to/plot_SNP_100updown_coverage_final.py ${SAMPLE}.all23771SNPs_100updown.tsv ${SAMPLE}
	python3 /path/to/targetSNP_covDepth_PCTplot_final.py ${SAMPLE}.all23771SNPs_depth.tsv ${SAMPLE} >> ${SAMPLE}.pipe_log.txt
    fi

done < $fn

rm -f ./SeqPrep_output/*_unmerged_combined*.fastq
gzip ./SeqPrep_output/*.fastq

rm -rf ./BWA_analyses/*.sai
rm -rf ./BWA_analyses/*_all_reads.complexity_filtered_EquCab2.bam

cat *picard_temp.csv | sed '3~2d' > Picard_metrics.csv
rm *temp*.tsv

mkdir ./CNER_CaptureAnalyses
mv *.*sv ./CNER_CaptureAnalyses
mv *.png ./CNER_CaptureAnalyses

#For final mapping metrics of all samples consolidation:
Rscript /path/to/shortReads_QCmetrics_v2.R ./SeqPrep_output/ ./BWA_analyses/ ./MapDamage_output/
