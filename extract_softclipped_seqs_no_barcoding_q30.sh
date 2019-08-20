#!/bin/bash
#SBATCH --job-name=extract_soft_clippings_from_subtel_ends
#SBATCH --time=0:08:00
#SBATCH --output=extract_soft_clippings_from_subtel_ends.output
#SBATCH --mem-per-cpu=8G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --image=docker:fcadete/guppy:latest
#SBATCH --workdir=/mnt/nfs/lobo/SALMEIDA-NFS/fcadete/scratch/minion/

mkdir soft_clipped_ends_no_barcoding_q30
mkdir soft_clipped_ends_no_barcoding_q30/by_subtelomere

for f in alignment_outputs_no_barcoding/*on_rhietman_mapont_q30.bam
do
  filename=$(basename $f)
  bam_output_filename=${filename/q30.bam/q30_subtelomere_start.bam}
  fq_output_filename=${filename/q30.bam/q30_subtelomere_start_softclipped.fq.gz}
#  srun -n1 -N1 shifter samtools view -L subtelomere_starts_1_5000.bed -b $f > soft_clipped_ends_no_barcoding_q30/$bam_output_filename
#  srun -n1 -N1 shifter samtools index soft_clipped_ends_no_barcoding_q30/$bam_output_filename
#  srun -n1 -N1 ~/SE-MEI/extractSoftclipped soft_clipped_ends_no_barcoding_q30/$bam_output_filename > soft_clipped_ends_no_barcoding_q30/$fq_output_filename

  while read line; do
    array=($line)
    subtelomere=${array[0]}
    echo $array
    echo ${subtelomere}
    echo $filename
 #   srun -n1 -N1 shifter samtools view -b soft_clipped_ends_no_barcoding_q30/$bam_output_filename $subtelomere > soft_clipped_ends_no_barcoding_q30/by_subtelomere/${subtelomere}_${bam_output_filename} </dev/null
 #   srun -n1 -N1 ~/SE-MEI/extractSoftclipped soft_clipped_ends_no_barcoding_q30/by_subtelomere/${subtelomere}_${bam_output_filename} > soft_clipped_ends_no_barcoding_q30/by_subtelomere/${subtelomere}_${fq_output_filename} </dev/null
 #   srun -n1 -N1 shifter /seqtk/seqtk seq -A soft_clipped_ends_no_barcoding_q30/by_subtelomere/${subtelomere}_${fq_output_filename} > soft_clipped_ends_no_barcoding_q30/by_subtelomere/${subtelomere}_${fq_output_filename/fq.gz/fasta} </dev/null
    srun -n1 -N1 ~/muscle3.8.31_i86linux64 -in soft_clipped_ends_no_barcoding_q30/by_subtelomere/${subtelomere}_${fq_output_filename/fq.gz/fasta} -out soft_clipped_ends_no_barcoding_q30/by_subtelomere/${subtelomere}_${fq_output_filename/fq.gz/aln} </dev/null
  done < subtelomere_starts_1_5000.bed

done

