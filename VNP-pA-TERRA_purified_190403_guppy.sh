#!/bin/bash
#SBATCH --job-name=guppy_VNP-pA-TERRA_purified_190403
#SBATCH --time=4:00:00
#SBATCH --output=guppy_VNP-pA-TERRA_purified_190403.output
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
#                                     -i 20190403_second_run/VNP-pA-TERRA_purified_190403/VNP-pA-TERRA_purified_190403/20190403_1419_MN29796_FAK43621_aaa66819/ \
#                                     -s 20190403_second_run/VNP-pA-TERRA_purified_190403/VNP-pA-TERRA_purified_190403/20190403_1419_MN29796_FAK43621_aaa66819/fastq_guppy
#
#srun -n1 -N1 --exclusive shifter Rscript /MinIONQC.R -i 20190403_second_run/VNP-pA-TERRA_purified_190403/VNP-pA-TERRA_purified_190403/20190403_1419_MN29796_FAK43621_aaa66819/fastq_guppy/sequencing_summary.txt \
#                                                     -o 20190403_second_run/VNP-pA-TERRA_purified_190403/VNP-pA-TERRA_purified_190403/20190403_1419_MN29796_FAK43621_aaa66819/fastq_guppy/qc_plots
#
#
srun -n1 -N1 --exclusive shifter /Porechop/porechop-runner.py -i 20190403_second_run/VNP-pA-TERRA_purified_190403/VNP-pA-TERRA_purified_190403/20190403_1419_MN29796_FAK43621_aaa66819/fastq_guppy/pass/ \
                                          -o 20190403_second_run/VNP-pA-TERRA_purified_190403/VNP-pA-TERRA_purified_190403/20190403_1419_MN29796_FAK43621_aaa66819/fastq_guppy/barcode_separated \
                                          --barcode_diff 1 \
                                          -t 32 --verbosity 3 > guppy_VNP-pA-TERRA_purified_190403.porechop_output
#
#srun -n1 -N1 --exclusive shifter /FastQC/fastqc -k 6 --nano --threads 32 20190403_second_run/VNP-pA-TERRA_purified_190403/VNP-pA-TERRA_purified_190403/20190403_1419_MN29796_FAK43621_aaa66819/fastq_guppy/pass_trimmed.fastq
#
#srun -n1 -N1 --exclusive shifter /minimap2/minimap2 -x map-ont -t 32 \
#                                     references/ConcatenatedFASTAAassemblies_hTel.txt \
#                                     20190403_second_run/VNP-pA-TERRA_purified_190403/VNP-pA-TERRA_purified_190403/20190403_1419_MN29796_FAK43621_aaa66819/fastq_guppy/pass_trimmed.fastq > \
#                                     alignment_outputs/VNP-pA-TERRA_purified_190403_on_rhietman_mapont.sam
#
#srun -n1 -N1 --exclusive shifter /minimap2/minimap2 -x splice -t 32 \
#                                     references/ConcatenatedFASTAAassemblies_hTel.txt \
#                                     20190403_second_run/VNP-pA-TERRA_purified_190403/VNP-pA-TERRA_purified_190403/20190403_1419_MN29796_FAK43621_aaa66819/fastq_guppy/pass_trimmed.fastq > \
#                                     alignment_outputs/VNP-pA-TERRA_purified_190403_on_rhietman_splice.sam
#
#srun -n1 -N1 --exclusive shifter /minimap2/minimap2 -x map-ont -t 32 \
#                                     references/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
#                                     20190403_second_run/VNP-pA-TERRA_purified_190403/VNP-pA-TERRA_purified_190403/20190403_1419_MN29796_FAK43621_aaa66819/fastq_guppy/pass_trimmed.fastq > \
#                                     alignment_outputs/VNP-pA-TERRA_purified_190403_on_GRCh38_mapont.sam
#
#srun -n1 -N1 --exclusive shifter /minimap2/minimap2 -x splice -t 32 \
#                                     references/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
#                                     20190403_second_run/VNP-pA-TERRA_purified_190403/VNP-pA-TERRA_purified_190403/20190403_1419_MN29796_FAK43621_aaa66819/fastq_guppy/pass_trimmed.fastq > \
#                                     alignment_outputs/VNP-pA-TERRA_purified_190403_on_GRCh38_splice.sam
#
#
#srun -n1 -N1 --exclusive shifter bash -c "samtools view -H -T references/ConcatenatedFASTAAassemblies_hTel.txt alignment_outputs/VNP-pA-TERRA_purified_190403_on_rhietman_mapont.sam > alignment_outputs/VNP-pA-TERRA_purified_190403_on_rhietman_mapont.header \
#                                          && cat alignment_outputs/VNP-pA-TERRA_purified_190403_on_rhietman_mapont.header alignment_outputs/VNP-pA-TERRA_purified_190403_on_rhietman_mapont.sam | samtools depth - > alignment_outputs/VNP-pA-TERRA_purified_190403_on_rhietman_mapont.depth"
#srun -n1 -N1 --exclusive shifter bash -c "samtools view -H -T references/ConcatenatedFASTAAassemblies_hTel.txt alignment_outputs/VNP-pA-TERRA_purified_190403_on_rhietman_splice.sam > alignment_outputs/VNP-pA-TERRA_purified_190403_on_rhietman_splice.header \
#                                          && cat alignment_outputs/VNP-pA-TERRA_purified_190403_on_rhietman_splice.header alignment_outputs/VNP-pA-TERRA_purified_190403_on_rhietman_splice.sam | samtools depth - > alignment_outputs/VNP-pA-TERRA_purified_190403_on_rhietman_splice.depth"
#srun -n1 -N1 --exclusive shifter bash -c "samtools view -H -T references/Homo_sapiens.GRCh38.dna.primary_assembly.fa alignment_outputs/VNP-pA-TERRA_purified_190403_on_GRCh38_mapont.sam > alignment_outputs/VNP-pA-TERRA_purified_190403_on_GRCh38_mapont.header \
#                                          && cat alignment_outputs/VNP-pA-TERRA_purified_190403_on_GRCh38_mapont.header alignment_outputs/VNP-pA-TERRA_purified_190403_on_GRCh38_mapont.sam | samtools depth - > alignment_outputs/VNP-pA-TERRA_purified_190403_on_GRCh38_mapont.depth"
#srun -n1 -N1 --exclusive shifter bash -c "samtools view -H -T references/Homo_sapiens.GRCh38.dna.primary_assembly.fa alignment_outputs/VNP-pA-TERRA_purified_190403_on_GRCh38_splice.sam > alignment_outputs/VNP-pA-TERRA_purified_190403_on_GRCh38_splice.header \
#                                          && cat alignment_outputs/VNP-pA-TERRA_purified_190403_on_GRCh38_splice.header alignment_outputs/VNP-pA-TERRA_purified_190403_on_GRCh38_splice.sam | samtools depth - > alignment_outputs/VNP-pA-TERRA_purified_190403_on_GRCh38_splice.depth"
#
#srun -n1 -N1 --exclusive shifter /bbmap/kmercountexact.sh in=20190403_second_run/VNP-pA-TERRA_purified_190403/VNP-pA-TERRA_purified_190403/20190403_1419_MN29796_FAK43621_aaa66819/fastq_guppy/pass_trimmed.fastq \
#                                                          out=20190403_second_run/VNP-pA-TERRA_purified_190403/VNP-pA-TERRA_purified_190403/20190403_1419_MN29796_FAK43621_aaa66819/fastq_guppy/pass_trimmed.kmer_count \
#                                                          fastadump=f k=6 threads=32 qin=33 overwrite=t rcomp=f
#
#srun -n1 -N1 --exclusive shifter /bbmap/commonkmers.sh in=20190403_second_run/VNP-pA-TERRA_purified_190403/VNP-pA-TERRA_purified_190403/20190403_1419_MN29796_FAK43621_aaa66819/fastq_guppy/pass_trimmed.fastq \
#                                                          out=20190403_second_run/VNP-pA-TERRA_purified_190403/VNP-pA-TERRA_purified_190403/20190403_1419_MN29796_FAK43621_aaa66819/fastq_guppy/pass_trimmed.kmer_common \
#                                                          fastadump=f k=6 threads=32 qin=33 display=100 count=t overwrite=t
#
#srun -n1 -N1 --exclusive shifter /bbmap/bbduk.sh in=20190403_second_run/VNP-pA-TERRA_purified_190403/VNP-pA-TERRA_purified_190403/20190403_1419_MN29796_FAK43621_aaa66819/fastq_guppy/pass_trimmed.fastq \
#                                                 literal=TTAGGGTTAGGGTTAGGG k=18 mm=f hdist=2 rcomp=f qin=33 \
#                                                 outm=20190403_second_run/VNP-pA-TERRA_purified_190403/VNP-pA-TERRA_purified_190403/20190403_1419_MN29796_FAK43621_aaa66819/fastq_guppy/pass_trimmed_3joint_repeats.fastq


srun -n1 -N1 --exclusive grep -o -n 'TT[AG]GGG' 20190403_second_run/VNP-pA-TERRA_purified_190403/VNP-pA-TERRA_purified_190403/20190403_1419_MN29796_FAK43621_aaa66819/fastq_guppy/pass_trimmed_3joint_repeats.fastq | cut -d : -f 1 | uniq -c > 20190403_second_run/VNP-pA-TERRA_purified_190403/VNP-pA-TERRA_purified_190403/20190403_1419_MN29796_FAK43621_aaa66819/fastq_guppy/pass_trimmed_3joint_repeats.TTAGGG_counts

echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,NodeList,NNodes,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID

