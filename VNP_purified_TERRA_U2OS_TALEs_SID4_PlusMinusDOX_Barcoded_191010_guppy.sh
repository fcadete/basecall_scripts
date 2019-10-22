#!/bin/bash
#SBATCH --job-name=guppy_VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_191010
#SBATCH --time=48:00:00
#SBATCH --output=guppy_VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_191010.output
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
#                                     -i 191010_Eighth_run/191010_TALES/NLS3_SID4/20191010_1549_MN29796_FAL24491_be365db2/ \
#                                     -s 191010_Eighth_run/191010_TALES/NLS3_SID4/20191010_1549_MN29796_FAL24491_be365db2/fastq_guppy
#
#srun -n1 -N1 shifter Rscript /MinIONQC.R -i 191010_Eighth_run/191010_TALES/NLS3_SID4/20191010_1549_MN29796_FAL24491_be365db2/fastq_guppy/sequencing_summary.txt \
#                                                     -o 191010_Eighth_run/191010_TALES/NLS3_SID4/20191010_1549_MN29796_FAL24491_be365db2/fastq_guppy/qc_plots
#
#
#srun -n1 -N1 shifter /Porechop/porechop-runner.py -i 191010_Eighth_run/191010_TALES/NLS3_SID4/20191010_1549_MN29796_FAL24491_be365db2/fastq_guppy/pass/ \
#                                          -b 191010_Eighth_run/191010_TALES/NLS3_SID4/20191010_1549_MN29796_FAL24491_be365db2/fastq_guppy/barcode_separated \
#                                          --check_reads 50000 -t 40 --verbosity 3 > guppy_VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_191010.porechop_output
#
#srun -n1 -N1 shifter python3 basecall_scripts/parse_porechop_output.py guppy_VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_191010.porechop_output
#

barcoded=( "BC01" "BC02" "BC03" "BC04" "BC05" "BC06" "BC07" "BC08" "BC09" "BC10" "BC11" "BC12" "none" )


for barcode in "${barcoded[@]}"; do

#srun -n1 -N1 shifter /FastQC/fastqc -k 6 --nano --threads 40 191010_Eighth_run/191010_TALES/NLS3_SID4/20191010_1549_MN29796_FAL24491_be365db2/fastq_guppy/barcode_separated/$barcode.fastq
#
#srun -n1 -N1 shifter bash -c "/minimap2/minimap2 -a -x map-ont -t 40 \
#                                     references/ConcatenatedFASTAAassemblies_hTel.txt \
#                                     191010_Eighth_run/191010_TALES/NLS3_SID4/20191010_1549_MN29796_FAL24491_be365db2/fastq_guppy/barcode_separated/$barcode.fastq | \
#                                     samtools view -b - > TALEs/alignment_outputs/VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_191010_${barcode}_on_rhietman_mapont.bam"
#
#
#srun -n1 -N1 shifter bash -c "samtools sort -@ 40 TALEs/alignment_outputs/VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_191010_${barcode}_on_rhietman_mapont.bam | samtools view -bq 30 - | samtools depth - > TALEs/alignment_outputs/VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_191010_${barcode}_on_rhietman_mapont_q30.depth"
#
#srun -n1 -N1 shifter bash -c "samtools sort -@ 40 TALEs/alignment_outputs/VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_191010_${barcode}_on_rhietman_mapont.bam | samtools view -bF 256 - | samtools depth - > TALEs/alignment_outputs/VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_191010_${barcode}_on_rhietman_mapont_primary.depth"
#
#srun -n1 -N1 shifter bash -c "samtools sort -@ 40 TALEs/alignment_outputs/VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_191010_${barcode}_on_rhietman_mapont.bam | samtools view -bq 30 - > TALEs/alignment_outputs/VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_191010_${barcode}_on_rhietman_mapont_q30.bam; samtools index TALEs/alignment_outputs/VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_191010_${barcode}_on_rhietman_mapont_q30.bam"
#
#srun -n1 -N1 shifter bash -c "samtools sort -@ 40 TALEs/alignment_outputs/VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_191010_${barcode}_on_rhietman_mapont.bam | samtools view -bF 256 - > TALEs/alignment_outputs/VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_191010_${barcode}_on_rhietman_mapont_primary.bam; samtools index TALEs/alignment_outputs/VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_191010_${barcode}_on_rhietman_mapont_primary.bam"
#
#
#srun -n1 -N1 shifter --image=docker:fcadete/guppy:latest /Porechop/porechop-runner.py -i 191010_Eighth_run/191010_TALES/NLS3_SID4/20191010_1549_MN29796_FAL24491_be365db2/fastq_guppy/barcode_separated/$barcode.fastq \
#                                          -b 191010_Eighth_run/191010_TALES/NLS3_SID4/20191010_1549_MN29796_FAL24491_be365db2/fastq_guppy/barcode_separated/$barcode \
#                                          --barcode_diff 1 --check_reads 50000 -t 40 --verbosity 3 > TALEs/primer_separation/guppy_VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_191010_$barcode.porechop_output
#
#srun -n1 -N1 shifter python3 basecall_scripts/parse_porechop_output.py TALEs/primer_separation/guppy_VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_191010_$barcode.porechop_output
#

mkdir TALEs/alignment_outputs_primered
   
   primers=("Joana_TERRA_VNP_no_PCR" "Joana_polyA_VNP_no_PCR" "none")
   
   for primer in "${primers[@]}"; do
   
   srun -n1 -N1 shifter bash -c "/minimap2/minimap2 -a -x map-ont -t 40 \
                                        references/ConcatenatedFASTAAassemblies_hTel.txt \
                                        191010_Eighth_run/191010_TALES/NLS3_SID4/20191010_1549_MN29796_FAL24491_be365db2/fastq_guppy/barcode_separated/$barcode/$primer.fastq | \
                                        samtools view -b - > TALEs/alignment_outputs_primered/VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_191010_${barcode}_${primer}_on_rhietman_mapont.bam"
   
   
   srun -n1 -N1 shifter bash -c "samtools sort -@ 40 TALEs/alignment_outputs_primered/VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_191010_${barcode}_${primer}_on_rhietman_mapont.bam | samtools view -bq 30 - | samtools depth - > TALEs/alignment_outputs_primered/VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_191010_${barcode}_${primer}_on_rhietman_mapont_q30.depth"
   
   srun -n1 -N1 shifter bash -c "samtools sort -@ 40 TALEs/alignment_outputs_primered/VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_191010_${barcode}_${primer}_on_rhietman_mapont.bam | samtools view -bF 256 - | samtools depth - > TALEs/alignment_outputs_primered/VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_191010_${barcode}_${primer}_on_rhietman_mapont_primary.depth"
   
   srun -n1 -N1 shifter bash -c "samtools sort -@ 40 TALEs/alignment_outputs_primered/VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_191010_${barcode}_${primer}_on_rhietman_mapont.bam | samtools view -bq 30 - > TALEs/alignment_outputs_primered/VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_191010_${barcode}_${primer}_on_rhietman_mapont_q30.bam; samtools index TALEs/alignment_outputs_primered/VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_191010_${barcode}_${primer}_on_rhietman_mapont_q30.bam"
   
   srun -n1 -N1 shifter bash -c "samtools sort -@ 40 TALEs/alignment_outputs_primered/VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_191010_${barcode}_${primer}_on_rhietman_mapont.bam | samtools view -bF 256 - > TALEs/alignment_outputs_primered/VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_191010_${barcode}_${primer}_on_rhietman_mapont_primary.bam; samtools index TALEs/alignment_outputs_primered/VNP_purified_TERRA_U2OS_TALEs_SID4_PlusMinusDOX_Barcoded_191010_${barcode}_${primer}_on_rhietman_mapont_primary.bam"
   
   
   done
   
done



echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,NodeList,NNodes,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID

