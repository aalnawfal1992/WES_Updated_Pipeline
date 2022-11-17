##!/bin/bash


Batch_Id='/home/bioinformatics2/Desktop/ClinicalExome_Batch1_01-11-2022'
###Change it to active batch using 'pwd'
ref='/home/bioinformatics2/Genomics_Reference_Files/Human_Reference_Genome/UCSC/UCSC_Homo_Sapiens_GRCH37_hg19.fasta'
bed='/home/bioinformatics2/Genomics_Reference_Files/Human_Reference_Genome/UCSC/bed_files/S04380110_Padded.bed'
bait='/home/bioinformatics2/Genomics_Reference_Files/Human_Reference_Genome/UCSC/bed_files/S04380110_baits.interval_list' #Generated from S04380110_Covered.bed
target='/home/bioinformatics2/Genomics_Reference_Files/Human_Reference_Genome/UCSC/bed_files/S04380110_targets.interval_list' #Generated from S04380110_Regions.bed

	for sample in $Batch_Id;
    do
            echo "Task Started:"  `date '+%T'`
            find $Batch_Id/*/*.fq.gz -type f | parallel -j 128 -k md5sum > $Batch_Id/md5_After.txt
            echo "MD5 Generated for all Fastq Files and Stored in "${PWD##*/}"/md5_After.txt" `date '+%T'`
                
    done 

for sample in "$Batch_Id"/*;
do
	if [ -d "$sample" ]
	then
		patirnt_MRN=$(basename "$sample")
		
	echo '### Cleaning Fastq Files started in' `date '+%T'` 'For' $patirnt_MRN '###'
	/home/bioinformatics2/fastp -w 16 --detect_adapter_for_pe --html "$sample/$patirnt_MRN" -y -i "$sample"/*_1.fq.gz -o "$sample/$patirnt_MRN"_1_out.fq.gz -I "$sample"/*_2.fq.gz -O "$sample/$patirnt_MRN"_2_out.fq.gz
	echo '### Cleaning Fastq Files Finished in' `date '+%T'` 'For' $patirnt_MRN '###'


	echo '### Mapping/Post-Mapping process Started in' `date '+%T'` 'For' $patirnt_MRN '###'
	/home/bioinformatics2/bwa/bwa mem -t 128 -R "@RG\tID:Novaseq\tSM:"$sample/$patirnt_MRN"\tPL:ILLUMINA\tPI:330" $ref "$sample/$patirnt_MRN"_1_out.fq.gz "$sample/$patirnt_MRN"_2_out.fq.gz | samtools sort -T "$sample" -@ 128 - > "$sample/$patirnt_MRN".bam
	echo '### Mapping/Post-Mapping process Finished in' `date '+%T'` 'For' $patirnt_MRN '###'


	echo "Runing FastQC for All Fastq Files" `date '+%T'`
    mkdir "${PWD##*/}"_FastQC
    perl /home/bioinformatics2/FastQC/fastqc -t 256 $Batch_Id/*/*.fq.gz -o "${PWD##*/}"_FastQC/
    echo "FastQC Result Stored in "${PWD##*/}"_FastQC" `date '+%T'`

	echo '### Deduplication start in' `date '+%T'` 'For' $patirnt_MRN '###'
	/home/bioinformatics2/sambamba-0.8.2 markdup -t 128 -r "$sample/$patirnt_MRN".bam "$sample/$patirnt_MRN"_PCRFree.bam
	echo '### Deduplication Finished in' `date '+%T'` 'For' $patirnt_MRN '###'

	samtools index -@ 128 "$sample/$patirnt_MRN"_PCRFree.bam
	rm "$sample/$patirnt_MRN".bam
    rm "$sample/$patirnt_MRN"_1_out.fq.gz "$sample/$patirnt_MRN"_2_out.fq.gz
	fi
done


	echo '### Post Aligment QC start in' `date '+%T'` '###'
	ls "$Batch_Id"/*/*_PCRFree.bam &>> ListOfBAMs.txt
	cat ListOfBAMs.txt | parallel -j 10 java -jar /home/bioinformatics2/picard.jar CollectHsMetrics -I {} -O {.}_metrics.txt -R $ref --BAIT_INTERVALS $bait --TARGET_INTERVALS $target
	cat ListOfBAMs.txt | parallel -j 100% /home/bioinformatics2/sambamba-0.8.2 flagstat -t 128 {} {.}.txt
	echo '### Post Aligment QC finished in' `date '+%T'` '###'


	echo '### Variant Calling Started in' `date '+%T'` '###'
	cat ListOfBAMs.txt | parallel -j 20 java -jar /home/bioinformatics2/gatk-4.2.6.1/gatk-package-4.2.6.1-local.jar HaplotypeCaller --native-pair-hmm-threads 128 -I {} -O {.}.vcf -R $ref -L $bed
	echo '### Processing finished in' `date '+%T'` '###'

exit 0 	


	

	