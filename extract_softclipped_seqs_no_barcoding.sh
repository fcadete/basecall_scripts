#!/bin/bash
#SBATCH --job-name=extract_soft_clippings_from_subtel_ends
#SBATCH --time=0:08:00
#SBATCH --output=extract_soft_clippings_from_subtel_ends.output
#SBATCH --mem-per-cpu=8G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --image=docker:fcadete/guppy:latest
#SBATCH --workdir=/mnt/nfs/lobo/SALMEIDA-NFS/fcadete/scratch/minion/

mkdir soft_clipped_ends_no_barcoding

for f in alignment_outputs_no_barcoding/*on_rhietman_mapont_primary.bam
do
  filename=$(basename $f)
  bam_output_filename=${filename/primary.bam/primary_subtelomere_start.bam}
  fq_output_filename=${filename/primary.bam/primary_subtelomere_start_softclipped.fq.gz}
  srun -n1 -N1 shifter samtools view -L subtelomere_starts_1_5000.bed -b $f > soft_clipped_ends_no_barcoding/$bam_output_filename
  srun -n1 -N1 ~/SE-MEI/extractSoftclipped soft_clipped_ends_no_barcoding/$bam_output_filename > soft_clipped_ends_no_barcoding/$fq_output_filename
done

