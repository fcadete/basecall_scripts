#!/bin/bash
#SBATCH --job-name=extract_soft_clippings_from_subtel_ends_TALEs_primered
#SBATCH --time=0:30:00
#SBATCH --output=extract_soft_clippings_from_subtel_ends_TALEs_primered.output
#SBATCH --mem-per-cpu=8G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --image=docker:fcadete/guppy:latest
#SBATCH --workdir=/mnt/nfs/lobo/SALMEIDA-NFS/fcadete/scratch/minion/

mkdir TALEs/soft_clipped_ends_primered

for f in TALEs/alignment_outputs_primered/*on_rhietman_mapont_primary.bam
do
  filename=$(basename $f)
  bam_output_filename=${filename/primary.bam/primary_subtelomere_start.bam}
  fq_output_filename=${filename/primary.bam/primary_subtelomere_start_softclipped.fq.gz}
  srun -n1 -N1 shifter samtools view -L subtelomere_starts_1_5000.bed -b $f > TALEs/soft_clipped_ends_primered/$bam_output_filename
  srun -n1 -N1 ~/SE-MEI/extractSoftclipped TALEs/soft_clipped_ends_primered/$bam_output_filename > TALEs/soft_clipped_ends_primered/$fq_output_filename
done

