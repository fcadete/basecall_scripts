#!/bin/bash
#SBATCH --job-name=guppy_VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded
#SBATCH --time=48:00:00
#SBATCH --output=guppy_VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded.output
#SBATCH --mem-per-cpu=6G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --image=docker:fcadete/guppy_original_barcodes:latest
#SBATCH --workdir=/mnt/nfs/lobo/SALMEIDA-NFS/fcadete/scratch/minion/

#srun -n1 -N1 shifter /ont-guppy-cpu/bin/guppy_basecaller \
#                                     -c /ont-guppy-cpu/data/dna_r9.4.1_450bps.cfg \
#                                     --num_callers 32 \
#                                     --qscore_filtering \
#                                     --recursive \
#                                     -i 190723_Sixth_run/190723_VNP-TERRA_USOS_SID4_1/U2OS_SID4_PlusMinusDOX_Barcoded/20190723_2008_MN29796_FAK86261_b4c588a4/ \
#                                     -s 190723_Sixth_run/190723_VNP-TERRA_USOS_SID4_1/U2OS_SID4_PlusMinusDOX_Barcoded/20190723_2008_MN29796_FAK86261_b4c588a4/fastq_guppy
#
#srun -n1 -N1 shifter Rscript /MinIONQC.R -i 190723_Sixth_run/190723_VNP-TERRA_USOS_SID4_1/U2OS_SID4_PlusMinusDOX_Barcoded/20190723_2008_MN29796_FAK86261_b4c588a4/fastq_guppy/sequencing_summary.txt \
#                                                     -o 190723_Sixth_run/190723_VNP-TERRA_USOS_SID4_1/U2OS_SID4_PlusMinusDOX_Barcoded/20190723_2008_MN29796_FAK86261_b4c588a4/fastq_guppy/qc_plots
#
#
#srun -n1 -N1 shifter /Porechop/porechop-runner.py -i 190723_Sixth_run/190723_VNP-TERRA_USOS_SID4_1/U2OS_SID4_PlusMinusDOX_Barcoded/20190723_2008_MN29796_FAK86261_b4c588a4/fastq_guppy/pass/ \
#                                          -b 190723_Sixth_run/190723_VNP-TERRA_USOS_SID4_1/U2OS_SID4_PlusMinusDOX_Barcoded/20190723_2008_MN29796_FAK86261_b4c588a4/fastq_guppy/barcode_separated \
#                                          --check_reads 50000 -t 32 --verbosity 3 > guppy_VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded.porechop_output
#
#srun -n1 -N1 shifter python3 basecall_scripts/parse_porechop_output.py guppy_VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded.porechop_output
#

barcoded=( "BC01" "BC02" "BC03" "BC04" "none" )


for barcode in "${barcoded[@]}"; do

srun -n1 -N1 shifter /FastQC/fastqc -k 6 --nano --threads 32 190723_Sixth_run/190723_VNP-TERRA_USOS_SID4_1/U2OS_SID4_PlusMinusDOX_Barcoded/20190723_2008_MN29796_FAK86261_b4c588a4/fastq_guppy/barcode_separated/$barcode.fastq

srun -n1 -N1 shifter bash -c "/minimap2/minimap2 -a -x map-ont -t 32 \
                                     references/ConcatenatedFASTAAassemblies_hTel.txt \
                                     190723_Sixth_run/190723_VNP-TERRA_USOS_SID4_1/U2OS_SID4_PlusMinusDOX_Barcoded/20190723_2008_MN29796_FAK86261_b4c588a4/fastq_guppy/barcode_separated/$barcode.fastq | \
                                     samtools view -b - > TALEs/alignment_outputs/VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_${barcode}_on_rhietman_mapont.bam"


srun -n1 -N1 shifter bash -c "samtools sort -@ 32 TALEs/alignment_outputs/VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_${barcode}_on_rhietman_mapont.bam | samtools view -bq 30 - | samtools depth - > TALEs/alignment_outputs/VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_${barcode}_on_rhietman_mapont_q30.depth"

srun -n1 -N1 shifter bash -c "samtools sort -@ 32 TALEs/alignment_outputs/VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_${barcode}_on_rhietman_mapont.bam | samtools view -bF 256 - | samtools depth - > TALEs/alignment_outputs/VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_${barcode}_on_rhietman_mapont_primary.depth"

srun -n1 -N1 shifter bash -c "samtools sort -@ 32 TALEs/alignment_outputs/VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_${barcode}_on_rhietman_mapont.bam | samtools view -bq 30 - > TALEs/alignment_outputs/VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_${barcode}_on_rhietman_mapont_q30.bam; samtools index alignment_outputs/VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_${barcode}_on_rhietman_mapont_q30.bam"

srun -n1 -N1 shifter bash -c "samtools sort -@ 32 TALEs/alignment_outputs/VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_${barcode}_on_rhietman_mapont.bam | samtools view -bF 256 - > TALEs/alignment_outputs/VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_${barcode}_on_rhietman_mapont_primary.bam; samtools index alignment_outputs/VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_${barcode}_on_rhietman_mapont_primary.bam"

done


#
## Map all reads, without using porechop to identify the adapters
## The rationale is that porechop is removing the adapters and, with them, useful telomeric sequences that we want to use.
#srun -n1 -N1 shifter bash -c "cat 190723_Sixth_run/190723_VNP-TERRA_USOS_SID4_1/U2OS_SID4_PlusMinusDOX_Barcoded/20190723_2008_MN29796_FAK86261_b4c588a4/fastq_guppy/pass/*fastq |
#                              /minimap2/minimap2 -a -x map-ont -t 32 \
#                                                references/ConcatenatedFASTAAassemblies_hTel.txt - |
#                              samtools view -b - > TALEs/alignment_outputs_no_barcoding/VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_on_rhietman_mapont.bam"
#
#srun -n1 -N1 shifter bash -c "samtools sort -@ 32 TALEs/alignment_outputs_no_barcoding/VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_on_rhietman_mapont.bam | samtools view -bq 30 - | samtools depth - > TALEs/alignment_outputs_no_barcoding/VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_on_rhietman_mapont_q30.depth"
#
#srun -n1 -N1 shifter bash -c "samtools sort -@ 32 TALEs/alignment_outputs_no_barcoding/VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_on_rhietman_mapont.bam | samtools view -bF 256 - | samtools depth - > TALEs/alignment_outputs_no_barcoding/VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_on_rhietman_mapont_primary.depth"
#
#srun -n1 -N1 shifter bash -c "samtools sort -@ 32 TALEs/alignment_outputs_no_barcoding/VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_on_rhietman_mapont.bam | samtools view -bq 30 - > TALEs/alignment_outputs_no_barcoding/VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_on_rhietman_mapont_q30.bam; samtools index alignment_outputs_no_barcoding/VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_on_rhietman_mapont_q30.bam"
#
#srun -n1 -N1 shifter bash -c "samtools sort -@ 32 TALEs/alignment_outputs_no_barcoding/VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_on_rhietman_mapont.bam | samtools view -bF 256 - > TALEs/alignment_outputs_no_barcoding/VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_on_rhietman_mapont_primary.bam; samtools index alignment_outputs_no_barcoding/VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_on_rhietman_mapont_primary.bam"
#

echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,NodeList,NNodes,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID

