#!/bin/bash
#SBATCH --job-name=guppy_VNP_pA-TERRA_20190327
#SBATCH --time=4:00:00
#SBATCH --output=guppy_VNP_pA-TERRA_20190327.output
#SBATCH --mem-per-cpu=6G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --image=docker:fcadete/guppy:latest
#SBATCH --workdir=/mnt/nfs/lobo/SALMEIDA-NFS/fcadete/scratch/minion/

#srun -n1 -N1 --exclusive shifter /ont-guppy-cpu/bin/guppy_basecaller \
#                                     -c /ont-guppy-cpu/data/dna_r9.4.1_450bps.cfg \
#                                     --num_callers 32 \
#                                     --qscore_filtering \
#                                     --recursive \
#                                     -i 20190327_first_run/VNP_polyA-TERRA/VNP_pA-TERRA_20190327/20190327_1801_MN29796_FAK43621_7641ca36/ \
#                                     -s 20190327_first_run/VNP_polyA-TERRA/VNP_pA-TERRA_20190327/20190327_1801_MN29796_FAK43621_7641ca36/fastq_guppy
#
#srun -n1 -N1 --exclusive shifter Rscript /MinIONQC.R -i 20190327_first_run/VNP_polyA-TERRA/VNP_pA-TERRA_20190327/20190327_1801_MN29796_FAK43621_7641ca36/fastq_guppy/sequencing_summary.txt \
#                                                     -o 20190327_first_run/VNP_polyA-TERRA/VNP_pA-TERRA_20190327/20190327_1801_MN29796_FAK43621_7641ca36/fastq_guppy/qc_plots
#
#
#srun -n1 -N1 --exclusive shifter /Porechop/porechop-runner.py -i 20190327_first_run/VNP_polyA-TERRA/VNP_pA-TERRA_20190327/20190327_1801_MN29796_FAK43621_7641ca36/fastq_guppy/pass/ \
#                                          -b 20190327_first_run/VNP_polyA-TERRA/VNP_pA-TERRA_20190327/20190327_1801_MN29796_FAK43621_7641ca36/fastq_guppy/barcode_separated \
#                                          --barcode_diff 1 --adapter_threshold 50.0 \
#                                          -t 32 --verbosity 3 > guppy_VNP_pA-TERRA_20190327.porechop_output
#
#srun -n1 -N1 --exclusive shifter python3 basecall_scripts/parse_porechop_output.py guppy_VNP_pA-TERRA_20190327.porechop_output
#
#
#barcoded=( "Joana_polyA_VNP" "Joana_polyA_VNP_no_PCR" "Joana_TERRA_VNP" "Joana_TERRA_VNP_no_PCR" "none" )
#
#
#for barcode in "${barcoded[@]}"; do
#
#srun -n1 -N1 --exclusive shifter /FastQC/fastqc -k 6 --nano --threads 32 20190327_first_run/VNP_polyA-TERRA/VNP_pA-TERRA_20190327/20190327_1801_MN29796_FAK43621_7641ca36/fastq_guppy/barcode_separated/$barcode.fastq
#
#srun -n1 -N1 --exclusive shifter /minimap2/minimap2 -a -x map-ont -t 32 \
#                                     references/ConcatenatedFASTAAassemblies_hTel.txt \
#                                     20190327_first_run/VNP_polyA-TERRA/VNP_pA-TERRA_20190327/20190327_1801_MN29796_FAK43621_7641ca36/fastq_guppy/barcode_separated/$barcode.fastq |
#                                     samtools view -b - > alignment_outputs/VNP_pA-TERRA_20190327_${barcode}_on_rhietman_mapont.bam
#
#srun -n1 -N1 --exclusive shifter /minimap2/minimap2 -a -x splice -t 32 \
#                                     references/ConcatenatedFASTAAassemblies_hTel.txt \
#                                     20190327_first_run/VNP_polyA-TERRA/VNP_pA-TERRA_20190327/20190327_1801_MN29796_FAK43621_7641ca36/fastq_guppy/barcode_separated/$barcode.fastq |
#                                     samtools view -b - > alignment_outputs/VNP_pA-TERRA_20190327_${barcode}_on_rhietman_splice.bam
#
#srun -n1 -N1 --exclusive shifter /minimap2/minimap2 -a -x map-ont -t 32 \
#                                     references/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
#                                     20190327_first_run/VNP_polyA-TERRA/VNP_pA-TERRA_20190327/20190327_1801_MN29796_FAK43621_7641ca36/fastq_guppy/barcode_separated/$barcode.fastq |
#                                     samtools view -b - > alignment_outputs/VNP_pA-TERRA_20190327_${barcode}_on_GRCh38_mapont.bam
#
#srun -n1 -N1 --exclusive shifter /minimap2/minimap2 -a -x splice -t 32 \
#                                     references/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
#                                     20190327_first_run/VNP_polyA-TERRA/VNP_pA-TERRA_20190327/20190327_1801_MN29796_FAK43621_7641ca36/fastq_guppy/barcode_separated/$barcode.fastq |
#                                     samtools view -b - > alignment_outputs/VNP_pA-TERRA_20190327_${barcode}_on_GRCh38_splice.bam
#
#
#
#srun -n1 -N1 --exclusive shifter bash -c "samtools sort -@ 32 alignment_outputs/VNP_pA-TERRA_20190327_${barcode}_on_rhietman_mapont.bam | samtools view -bq 30 - | samtools depth - > alignment_outputs/VNP_pA-TERRA_20190327_${barcode}_on_rhietman_mapont_q30.depth"
#srun -n1 -N1 --exclusive shifter bash -c "samtools sort -@ 32 alignment_outputs/VNP_pA-TERRA_20190327_${barcode}_on_rhietman_splice.bam | samtools view -bq 30 - | samtools depth - > alignment_outputs/VNP_pA-TERRA_20190327_${barcode}_on_rhietman_splice_q30.depth"
#srun -n1 -N1 --exclusive shifter bash -c "samtools sort -@ 32 alignment_outputs/VNP_pA-TERRA_20190327_${barcode}_on_GRCh38_mapont.bam | samtools view -bq 30 - | samtools depth - > alignment_outputs/VNP_pA-TERRA_20190327_${barcode}_on_GRCh38_mapont_q30.depth"
#srun -n1 -N1 --exclusive shifter bash -c "samtools sort -@ 32 alignment_outputs/VNP_pA-TERRA_20190327_${barcode}_on_GRCh38_splice.bam | samtools view -bq 30 - | samtools depth - > alignment_outputs/VNP_pA-TERRA_20190327_${barcode}_on_GRCh38_splice_q30.depth"
#
#srun -n1 -N1 --exclusive shifter bash -c "samtools sort -@ 32 alignment_outputs/VNP_pA-TERRA_20190327_${barcode}_on_rhietman_mapont.bam | samtools view -bF 256 - | samtools depth - > alignment_outputs/VNP_pA-TERRA_20190327_${barcode}_on_rhietman_mapont_primary.depth"
#srun -n1 -N1 --exclusive shifter bash -c "samtools sort -@ 32 alignment_outputs/VNP_pA-TERRA_20190327_${barcode}_on_rhietman_splice.bam | samtools view -bF 256 - | samtools depth - > alignment_outputs/VNP_pA-TERRA_20190327_${barcode}_on_rhietman_splice_primary.depth"
#srun -n1 -N1 --exclusive shifter bash -c "samtools sort -@ 32 alignment_outputs/VNP_pA-TERRA_20190327_${barcode}_on_GRCh38_mapont.bam | samtools view -bF 256 - | samtools depth - > alignment_outputs/VNP_pA-TERRA_20190327_${barcode}_on_GRCh38_mapont_primary.depth"
#srun -n1 -N1 --exclusive shifter bash -c "samtools sort -@ 32 alignment_outputs/VNP_pA-TERRA_20190327_${barcode}_on_GRCh38_splice.bam | samtools view -bF 256 - | samtools depth - > alignment_outputs/VNP_pA-TERRA_20190327_${barcode}_on_GRCh38_splice_primary.depth"
#
#srun -n1 -N1 --exclusive shifter bash -c "samtools sort -@ 32 alignment_outputs/VNP_pA-TERRA_20190327_${barcode}_on_rhietman_mapont.bam | samtools view -bq 30 - > alignment_outputs/VNP_pA-TERRA_20190327_${barcode}_on_rhietman_mapont_q30.bam; samtools index alignment_outputs/VNP_pA-TERRA_20190327_${barcode}_on_rhietman_mapont_q30.bam"
#srun -n1 -N1 --exclusive shifter bash -c "samtools sort -@ 32 alignment_outputs/VNP_pA-TERRA_20190327_${barcode}_on_rhietman_splice.bam | samtools view -bq 30 - > alignment_outputs/VNP_pA-TERRA_20190327_${barcode}_on_rhietman_splice_q30.bam; samtools index alignment_outputs/VNP_pA-TERRA_20190327_${barcode}_on_rhietman_splice_q30.bam"
#srun -n1 -N1 --exclusive shifter bash -c "samtools sort -@ 32 alignment_outputs/VNP_pA-TERRA_20190327_${barcode}_on_GRCh38_mapont.bam | samtools view -bq 30 - > alignment_outputs/VNP_pA-TERRA_20190327_${barcode}_on_GRCh38_mapont_q30.bam; samtools index alignment_outputs/VNP_pA-TERRA_20190327_${barcode}_on_GRCh38_mapont_q30.bam"
#srun -n1 -N1 --exclusive shifter bash -c "samtools sort -@ 32 alignment_outputs/VNP_pA-TERRA_20190327_${barcode}_on_GRCh38_splice.bam | samtools view -bq 30 - > alignment_outputs/VNP_pA-TERRA_20190327_${barcode}_on_GRCh38_splice_q30.bam; samtools index alignment_outputs/VNP_pA-TERRA_20190327_${barcode}_on_GRCh38_splice_q30.bam"
#
#
#srun -n1 -N1 --exclusive shifter bash -c "samtools sort -@ 32 alignment_outputs/VNP_pA-TERRA_20190327_${barcode}_on_rhietman_mapont.bam | samtools view -bF 256 - > alignment_outputs/VNP_pA-TERRA_20190327_${barcode}_on_rhietman_mapont_primary.bam; samtools index alignment_outputs/VNP_pA-TERRA_20190327_${barcode}_on_rhietman_mapont_primary.bam"
#srun -n1 -N1 --exclusive shifter bash -c "samtools sort -@ 32 alignment_outputs/VNP_pA-TERRA_20190327_${barcode}_on_rhietman_splice.bam | samtools view -bF 256 - > alignment_outputs/VNP_pA-TERRA_20190327_${barcode}_on_rhietman_splice_primary.bam; samtools index alignment_outputs/VNP_pA-TERRA_20190327_${barcode}_on_rhietman_splice_primary.bam"
#srun -n1 -N1 --exclusive shifter bash -c "samtools sort -@ 32 alignment_outputs/VNP_pA-TERRA_20190327_${barcode}_on_GRCh38_mapont.bam | samtools view -bF 256 - > alignment_outputs/VNP_pA-TERRA_20190327_${barcode}_on_GRCh38_mapont_primary.bam; samtools index alignment_outputs/VNP_pA-TERRA_20190327_${barcode}_on_GRCh38_mapont_primary.bam"
#srun -n1 -N1 --exclusive shifter bash -c "samtools sort -@ 32 alignment_outputs/VNP_pA-TERRA_20190327_${barcode}_on_GRCh38_splice.bam | samtools view -bF 256 - > alignment_outputs/VNP_pA-TERRA_20190327_${barcode}_on_GRCh38_splice_primary.bam; samtools index alignment_outputs/VNP_pA-TERRA_20190327_${barcode}_on_GRCh38_splice_primary.bam"
#
#
#
#done
#
#srun -n1 -N1 --exclusive shifter /bbmap/kmercountexact.sh in=20190327_first_run/VNP_polyA-TERRA/VNP_pA-TERRA_20190327/20190327_1801_MN29796_FAK43621_7641ca36/fastq_guppy/pass_trimmed.fastq \
#                                                          out=20190327_first_run/VNP_polyA-TERRA/VNP_pA-TERRA_20190327/20190327_1801_MN29796_FAK43621_7641ca36/fastq_guppy/pass_trimmed.kmer_count \
#                                                          fastadump=f k=6 threads=32 qin=33 rcomp=f overwrite=t
#
#srun -n1 -N1 --exclusive shifter /bbmap/commonkmers.sh in=20190327_first_run/VNP_polyA-TERRA/VNP_pA-TERRA_20190327/20190327_1801_MN29796_FAK43621_7641ca36/fastq_guppy/pass_trimmed.fastq \
#                                                          out=20190327_first_run/VNP_polyA-TERRA/VNP_pA-TERRA_20190327/20190327_1801_MN29796_FAK43621_7641ca36/fastq_guppy/pass_trimmed.kmer_common \
#                                                          fastadump=f k=6 threads=32 qin=33 display=100 count=t overwrite=t
#
#srun -n1 -N1 --exclusive shifter /bbmap/bbduk.sh in=20190327_first_run/VNP_polyA-TERRA/VNP_pA-TERRA_20190327/20190327_1801_MN29796_FAK43621_7641ca36/fastq_guppy/pass_trimmed.fastq \
#                                                 literal=TTAGGGTTAGGGTTAGGG k=18 mm=f hdist=2 rcomp=f qin=33 \
#                                                 outm=20190327_first_run/VNP_polyA-TERRA/VNP_pA-TERRA_20190327/20190327_1801_MN29796_FAK43621_7641ca36/fastq_guppy/pass_trimmed_3joint_repeats.fastq
#
#srun -n1 -N1 --exclusive grep -o -n 'TT[AG]GGG' 20190327_first_run/VNP_polyA-TERRA/VNP_pA-TERRA_20190327/20190327_1801_MN29796_FAK43621_7641ca36/fastq_guppy/pass_trimmed_3joint_repeats.fastq | cut -d : -f 1 | uniq -c > 20190327_first_run/VNP_polyA-TERRA/VNP_pA-TERRA_20190327/20190327_1801_MN29796_FAK43621_7641ca36/fastq_guppy/pass_trimmed_3joint_repeats.TTAGGG_counts


# Map reads to GRCh38, then use the ones that do NOT map to map agains the Riethman lab assembly
srun -n1 -N1 --exclusive shifter bash -c "cat 20190327_first_run/VNP_polyA-TERRA/VNP_pA-TERRA_20190327/20190327_1801_MN29796_FAK43621_7641ca36/fastq_guppy/barcode_separated/*fastq |
                                          /minimap2/minimap2 -a -x map-ont -t 32 \
                                              references/Homo_sapiens.GRCh38.dna.primary_assembly.fa - |
                                          samtools view -f 4 -b - | samtools sort -n - | samtools fastq - |
                                          /minimap2/minimap2 -a -x map-ont -t 32 \
                                              references/ConcatenatedFASTAAassemblies_hTel.txt - | samtools view -b - > alignment_outputs/VNP_pA-TERRA_20190327_maps_only_to_riethman_assembly.bam &&
                                          samtools sort -@ 32 alignment_outputs/VNP_pA-TERRA_20190327_maps_only_to_riethman_assembly.bam | samtools view -bq 30 - | samtools depth - > alignment_outputs/VNP_pA-TERRA_20190327_maps_only_to_riethman_assembly_q30.depth"


echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,NodeList,NNodes,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID

