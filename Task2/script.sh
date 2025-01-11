###################################################### VARIANT CALLING STEPS ####################################################################


#Directories
# - project_dir(/mnt/c/project/Pupil_Bio/task_2)
#    - ref
#    - reads
#    - results
#    - aligned
#    - tmp


# Set directories path
ref=/mnt/c/project/Pupil_Bio/task_2/ref/hg38.fa  
project_dir=/mnt/c/project/Pupil_Bio/task_2
aligned_reads=$project_dir/aligned
reads=$project_dir/reads
results=$project_dir/results



# -------------------
# STEP 1: QC - Run fastqc 
# -------------------

echo "STEP 1: QC - Run fastqc"

#for i in *.gz; do fastqc $i -o fastqc/; done 

fastqc ${reads}/PA220KH-lib09-P19-Tumor_S2_L001_R1_001.fastq -o ${reads}/
fastqc ${reads}/PA220KH-lib09-P19-Tumor_S2_L001_R2_001.fastq -o ${reads}/
fastqc ${reads}/PA221MH-lib09-P19-Norm_S1_L001_R1_001.fastq -o ${reads}/
fastqc ${reads}/PA221MH-lib09-P19-Norm_S1_L001_R2_001.fastq -o ${reads}/

# ------------------------------
#STEP 1.a: MUTIQC
# -------------------------------
echo "summarizing fastqc results using multiqc"

mutiqc ${reads}/

# --------------------------------
# STEP 1.b: Trimming(Trimmomatic)
# --------------------------------
echo "trimming the adapters"

trimmomatic PE -threads 4 ${reads}/PA220KH-lib09-P19-Tumor_S2_L001_R1_001.fastq ${reads}/PA220KH-lib09-P19-Tumor_S2_L001_R2_001.fastq \
${reads}/Tumor_paired_forward.fastq /dev/null \
${reads}/Tumor_paired_reverse.fastq /dev/null \
ILLUMINACLIP:adapters.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:36

trimmomatic PE -threads 4 ${reads}/PA221MH-lib09-P19-Norm_S1_L001_R1_001.fastq  ${reads}/PA221MH-lib09-P19-Norm_S1_L001_R2_001.fastq \
${reads}/Normal_paired_forward.fastq /dev/null \
${reads}/Normal_paired_reverse.fastq /dev/null \
ILLUMINACLIP:adapters.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:36

# /dev/null discards unpaired reads instead of saving them to a file.

# --------------------------------------
# STEP 2: Map to reference using BWA-MEM
# --------------------------------------

echo "STEP 2: Map to reference using BWA-MEM"

# BWA index reference 
bwa index ${ref}

# BWA alignment
bwa mem -t 4 -R "@RG\tID:Normal\tPL:ILLUMINA\tSM:Normal" ${ref}${reads}/Normal_paired_forward.fastq ${reads}/Normal_paired_reverse.fastq > ${aligned_reads}/Normal_paired.sam
bwa mem -t 4 -R "@RG\tID:Tumor\tPL:ILLUMINA\tSM:Tumor" ${ref} ${reads}/Tumor_paired_forward.fastq ${reads}/Tumor_paired_reverse.fastq > ${aligned_reads}/Tumor_paired.sam

#---------------------------------
# Activate the Conda environment
#---------------------------------
echo "Activating Conda environment..."

source /mnt/c/project/conda/Anaconda3-2023.09-0-Linux-x86_64.sh  # Adjust this path based on your Conda installation
conda activate gatk4_env

echo "gatk4 activated" 

# -----------------------------------------
# STEP 3: Mark Duplicates and Sort - GATK4
# -----------------------------------------

echo "STEP 3: Mark Duplicates and Sort - GATK4"


gatk MarkDuplicatesSpark -I ${aligned_reads}/Normal_paired.sam -O ${aligned_reads}/Tumor_sorted_dedup_reads.bam
gatk MarkDuplicatesSpark -I ${aligned_reads}/Tumor_paired.sam -O ${aligned_reads}/Normal_sorted_dedup_reads.bam




# -----------------------------------------------
# STEP 4: Collect Alignment & Insert Size Metrics
# -----------------------------------------------


echo "STEP 5: Collect Alignment & Insert Size Metrics"

gatk CollectAlignmentSummaryMetrics R=${ref} I= ${aligned_reads}/Normal_sorted_dedup_reads.bam  O=${aligned_reads}/normal_alignment_metrics.txt
gatk CollectInsertSizeMetrics INPUT= ${aligned_reads}/Normal_sorted_dedup_reads.bam OUTPUT=${aligned_reads}/Normal_insert_size_metrics.txt HISTOGRAM_FILE=${aligned_reads}/Normal_insert_size_histogram.pdf


gatk CollectAlignmentSummaryMetrics R=${ref} I=${aligned_reads}/Tumor_sorted_dedup_reads.bam  O=${aligned_reads}/Tumor_alignment_metrics.txt
gatk CollectInsertSizeMetrics INPUT=${aligned_reads}/Tumor_sorted_dedup_reads.bam OUTPUT=${aligned_reads}/Tumor_insert_size_metrics.txt HISTOGRAM_FILE=${aligned_reads}/Tumor_insert_size_histogram.pdf


#--------------------------------------------------------------------------------------------------------------------------------------
# VARIANT CALLING USING MUTECT2
# -------------------------------------------------------------------------------------------------------------------------------------












