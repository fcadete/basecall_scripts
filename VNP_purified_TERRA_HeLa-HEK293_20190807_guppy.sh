#!/bin/bash
#SBATCH --job-name=guppy_VNP_purified_TERRA_HeLa-HEK293_20190807
#SBATCH --time=48:00:00
#SBATCH --output=guppy_VNP_purified_TERRA_HeLa-HEK293_20190807.output
#SBATCH --mem-per-cpu=6G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --image=docker:fcadete/guppy_original_barcodes:latest
#SBATCH --workdir=/mnt/nfs/lobo/SALMEIDA-NFS/fcadete/scratch/minion/

#srun -n1 -N1 shifter /ont-guppy-cpu/bin/guppy_basecaller \
#                                     -c /ont-guppy-cpu/data/dna_r9.4.1_450bps.cfg \
#                                     --num_callers 40 \
#                                     --qscore_filtering \
#                                     --recursive \
#                                     -i 190807_Seventh_run/190807_VNP-TERRA_HeLa_HEK293/HeLA_HEK_barcoded/20190807_1527_MN29796_FAL00837_baf29650/ \
#                                     -s 190807_Seventh_run/190807_VNP-TERRA_HeLa_HEK293/HeLA_HEK_barcoded/20190807_1527_MN29796_FAL00837_baf29650/fastq_guppy
#
#srun -n1 -N1 shifter Rscript /MinIONQC.R -i 190807_Seventh_run/190807_VNP-TERRA_HeLa_HEK293/HeLA_HEK_barcoded/20190807_1527_MN29796_FAL00837_baf29650/fastq_guppy/sequencing_summary.txt \
#                                                     -o 190807_Seventh_run/190807_VNP-TERRA_HeLa_HEK293/HeLA_HEK_barcoded/20190807_1527_MN29796_FAL00837_baf29650/fastq_guppy/qc_plots
#
#
#srun -n1 -N1 shifter /Porechop/porechop-runner.py -i 190807_Seventh_run/190807_VNP-TERRA_HeLa_HEK293/HeLA_HEK_barcoded/20190807_1527_MN29796_FAL00837_baf29650/fastq_guppy/pass/ \
#                                          -b 190807_Seventh_run/190807_VNP-TERRA_HeLa_HEK293/HeLA_HEK_barcoded/20190807_1527_MN29796_FAL00837_baf29650/fastq_guppy/barcode_separated \
#                                          --check_reads 50000 -t 40 --verbosity 3 > guppy_VNP_purified_TERRA_HeLa-HEK293_20190807.porechop_output
#
#srun -n1 -N1 shifter python3 basecall_scripts/parse_porechop_output.py guppy_VNP_purified_TERRA_HeLa-HEK293_20190807.porechop_output
#
#srun -n1 -N1 shifter Rscript -e "library(dplyr); library('readr'); library('stringr'); latest_seq_barcodes <- read_tsv('guppy_VNP_purified_TERRA_HeLa-HEK293_20190807.porechop_finalcall_table'); barcode_info <- read_tsv('190807_Seventh_run/Barcodes_clean.txt', col_names = c('cell_type', 'date', 'barcode')) %>% mutate(barcode = str_replace(barcode, 'NB', 'BC')); full_join(barcode_info, latest_seq_barcodes %>% count(final_barcode_call), by = c('barcode' = 'final_barcode_call')) %>% write_tsv(path = 'guppy_VNP_purified_TERRA_HeLa-HEK293_20190807.barcode_counts_output') "
#

barcoded=( "BC05" "BC06" "BC07" "BC08"
           "BC09" "BC10" "BC11" "BC12"
           "none" )

mkdir HeLa-HEK_barcoded/terra_primer_separated

for barcode in "${barcoded[@]}"; do
   
   srun -n1 -N1 shifter /FastQC/fastqc -k 6 --nano --threads 40 190807_Seventh_run/190807_VNP-TERRA_HeLa_HEK293/HeLA_HEK_barcoded/20190807_1527_MN29796_FAL00837_baf29650/fastq_guppy/barcode_separated/$barcode.fastq
   
   srun -n1 -N1 shifter --image=fcadete/guppy /Porechop/porechop-runner.py -i 190807_Seventh_run/190807_VNP-TERRA_HeLa_HEK293/HeLA_HEK_barcoded/20190807_1527_MN29796_FAL00837_baf29650/fastq_guppy/barcode_separated/$barcode.fastq \
                                             -b HeLa-HEK_barcoded/terra_primer_separated/$barcode \
                                             --barcode_diff 1 -t 40 --verbosity 3 > HeLa-HEK_barcoded/terra_primer_separated/$barcode.porechop_output
   
   srun -n1 -N1 shifter python3 basecall_scripts/parse_porechop_output.py HeLa-HEK_barcoded/terra_primer_separated/$barcode.porechop_output
   
   
   for primered_file in HeLa-HEK_barcoded/terra_primer_separated/$barcode/*fastq; do
      
      srun -n1 -N1 shifter bash -c "/minimap2/minimap2 -a -x map-ont -t 40 \
                                           references/ConcatenatedFASTAAassemblies_hTel.txt \
                                           $primered_file | \
                                           samtools view -b - > ${primered_file/.fastq/_on_rhietman_mapont.bam}"
      
      
      srun -n1 -N1 shifter bash -c "samtools sort -@ 40 ${primered_file/.fastq/_on_rhietman_mapont.bam} | samtools view -bq 30 - | samtools depth - > ${primered_file/.fastq/_on_rhietman_mapont_q30.depth}"
      
      srun -n1 -N1 shifter bash -c "samtools sort -@ 40 ${primered_file/.fastq/_on_rhietman_mapont.bam} | samtools view -bF 256 - | samtools depth - > ${primered_file/.fastq/_on_rhietman_mapont_primary.depth}"
      
      srun -n1 -N1 shifter bash -c "samtools sort -@ 40 ${primered_file/.fastq/_on_rhietman_mapont.bam} | samtools view -bq 30 - > ${primered_file/.fastq/_on_rhietman_mapont_q30.bam}; samtools index ${primered_file/.fastq/_on_rhietman_mapont_q30.bam}"
      
      srun -n1 -N1 shifter bash -c "samtools sort -@ 40 ${primered_file/.fastq/_on_rhietman_mapont.bam} | samtools view -bF 256 - > ${primered_file/.fastq/_on_rhietman_mapont_primary.bam}; samtools index ${primered_file/.fastq/_on_rhietman_mapont_primary.bam}"
      
   done
   
done


#
## Map all reads, without using porechop to identify the adapters
## The rationale is that porechop is removing the adapters and, with them, useful telomeric sequences that we want to use.
#srun -n1 -N1 shifter bash -c "cat 190807_Seventh_run/190807_VNP-TERRA_HeLa_HEK293/HeLA_HEK_barcoded/20190807_1527_MN29796_FAL00837_baf29650/fastq_guppy/pass/*fastq |
#                              /minimap2/minimap2 -a -x map-ont -t 40 \
#                                                references/ConcatenatedFASTAAassemblies_hTel.txt - |
#                              samtools view -b - > HeLa-HEK_barcoded/alignment_outputs_no_barcoding/VNP_purified_TERRA_HeLa-HEK293_20190807_on_rhietman_mapont.bam"
#
#srun -n1 -N1 shifter bash -c "samtools sort -@ 40 HeLa-HEK_barcoded/alignment_outputs_no_barcoding/VNP_purified_TERRA_HeLa-HEK293_20190807_on_rhietman_mapont.bam | samtools view -bq 30 - | samtools depth - > HeLa-HEK_barcoded/alignment_outputs_no_barcoding/VNP_purified_TERRA_HeLa-HEK293_20190807_on_rhietman_mapont_q30.depth"
#
#srun -n1 -N1 shifter bash -c "samtools sort -@ 40 HeLa-HEK_barcoded/alignment_outputs_no_barcoding/VNP_purified_TERRA_HeLa-HEK293_20190807_on_rhietman_mapont.bam | samtools view -bF 256 - | samtools depth - > HeLa-HEK_barcoded/alignment_outputs_no_barcoding/VNP_purified_TERRA_HeLa-HEK293_20190807_on_rhietman_mapont_primary.depth"
#
#srun -n1 -N1 shifter bash -c "samtools sort -@ 40 HeLa-HEK_barcoded/alignment_outputs_no_barcoding/VNP_purified_TERRA_HeLa-HEK293_20190807_on_rhietman_mapont.bam | samtools view -bq 30 - > HeLa-HEK_barcoded/alignment_outputs_no_barcoding/VNP_purified_TERRA_HeLa-HEK293_20190807_on_rhietman_mapont_q30.bam; samtools index HeLa-HEK_barcoded/alignment_outputs_no_barcoding/VNP_purified_TERRA_HeLa-HEK293_20190807_on_rhietman_mapont_q30.bam"
#
#srun -n1 -N1 shifter bash -c "samtools sort -@ 40 HeLa-HEK_barcoded/alignment_outputs_no_barcoding/VNP_purified_TERRA_HeLa-HEK293_20190807_on_rhietman_mapont.bam | samtools view -bF 256 - > HeLa-HEK_barcoded/alignment_outputs_no_barcoding/VNP_purified_TERRA_HeLa-HEK293_20190807_on_rhietman_mapont_primary.bam; samtools index HeLa-HEK_barcoded/alignment_outputs_no_barcoding/VNP_purified_TERRA_HeLa-HEK293_20190807_on_rhietman_mapont_primary.bam"
#

echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,NodeList,NNodes,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID

