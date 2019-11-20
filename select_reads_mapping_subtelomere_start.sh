#!/bin/bash
#SBATCH --job-name=select_reads_mapping_to_very_end
#SBATCH --time=0:08:00
#SBATCH --output=select_reads_mapping_to_very_end.output
#SBATCH --mem-per-cpu=8G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --image=docker:fcadete/guppy:latest
#SBATCH --workdir=/mnt/nfs/lobo/SALMEIDA-NFS/fcadete/scratch/minion/

#mkdir reads_mapping_to_very_end
#
#for f in alignment_outputs/*on_rhietman_mapont_q30.bam
#do
#  filename=$(basename $f)
#  bam_output_filename=${filename/q30.bam/q30_very_end.bam}
#  depth_output_filename=${filename/q30.bam/q30_very_end.depth}
#  srun -n1 -N1 shifter samtools view -L subtelomere_starts_1_10.bed -b $f > reads_mapping_to_very_end/$bam_output_filename
#  srun -n1 -N1 shifter samtools depth reads_mapping_to_very_end/$bam_output_filename > reads_mapping_to_very_end/$depth_output_filename
#done
#
#for f in alignment_outputs/*on_rhietman_mapont_primary.bam
#do
#  filename=$(basename $f)
#  bam_output_filename=${filename/primary.bam/primary_very_end.bam}
#  depth_output_filename=${filename/primary.bam/primary_very_end.depth}
#  srun -n1 -N1 shifter samtools view -L subtelomere_starts_1_10.bed -b $f > reads_mapping_to_very_end/$bam_output_filename
#  srun -n1 -N1 shifter samtools depth reads_mapping_to_very_end/$bam_output_filename > reads_mapping_to_very_end/$depth_output_filename
#done
#
#
#mkdir reads_mapping_to_very_end/HeLa-HEK_barcoded/
#
#barcode_pat='BC([0-9]+)'
#for f in HeLa-HEK_barcoded/terra_primer_separated/BC*/*on_rhietman_mapont_primary.bam
#do
#  [[ $f =~ $barcode_pat ]]
#  barcode=${BASH_REMATCH[0]}
#  mkdir reads_mapping_to_very_end/HeLa-HEK_barcoded/$barcode
#  filename=$(basename $f)
#  bam_output_filename=${filename/primary.bam/primary_very_end.bam}
#  depth_output_filename=${filename/primary.bam/primary_very_end.depth}
#  srun -n1 -N1 shifter samtools view -L subtelomere_starts_1_10.bed -b $f > reads_mapping_to_very_end/HeLa-HEK_barcoded/$barcode/$bam_output_filename
#  srun -n1 -N1 shifter samtools depth reads_mapping_to_very_end/HeLa-HEK_barcoded/$barcode/$bam_output_filename > reads_mapping_to_very_end/HeLa-HEK_barcoded/$barcode/$depth_output_filename
#done
#
#barcode_pat='BC([0-9]+)'
#for f in HeLa-HEK_barcoded/terra_primer_separated/BC*/*on_rhietman_mapont_q30.bam
#do
#  [[ $f =~ $barcode_pat ]]
#  barcode=${BASH_REMATCH[0]}
#  mkdir reads_mapping_to_very_end/HeLa-HEK_barcoded/$barcode
#  filename=$(basename $f)
#  bam_output_filename=${filename/q30.bam/q30_very_end.bam}
#  depth_output_filename=${filename/q30.bam/q30_very_end.depth}
#  srun -n1 -N1 shifter samtools view -L subtelomere_starts_1_10.bed -b $f > reads_mapping_to_very_end/HeLa-HEK_barcoded/$barcode/$bam_output_filename
#  srun -n1 -N1 shifter samtools depth reads_mapping_to_very_end/HeLa-HEK_barcoded/$barcode/$bam_output_filename > reads_mapping_to_very_end/HeLa-HEK_barcoded/$barcode/$depth_output_filename
#done
#

mkdir reads_mapping_to_very_end/TALEs_primered/

barcode_pat='BC([0-9]+)'
for f in TALEs/alignment_outputs_primered/*191010*primary.bam
do
  [[ $f =~ $barcode_pat|none ]]
  barcode=${BASH_REMATCH[0]}
  mkdir reads_mapping_to_very_end/TALEs_primered/$barcode
  filename=$(basename $f)
  bam_output_filename=${filename/primary.bam/primary_very_end.bam}
  depth_output_filename=${filename/primary.bam/primary_very_end.depth}
  srun -n1 -N1 shifter samtools view -L subtelomere_starts_1_10.bed -b $f > reads_mapping_to_very_end/TALEs_primered/$barcode/$bam_output_filename
  srun -n1 -N1 shifter samtools depth reads_mapping_to_very_end/TALEs_primered/$barcode/$bam_output_filename > reads_mapping_to_very_end/TALEs_primered/$barcode/$depth_output_filename
done

barcode_pat='BC([0-9]+)'
for f in TALEs/alignment_outputs_primered/*191010*q30.bam
do
  [[ $f =~ $barcode_pat|none ]]
  barcode=${BASH_REMATCH[0]}
  mkdir reads_mapping_to_very_end/TALEs_primered/$barcode
  filename=$(basename $f)
  bam_output_filename=${filename/q30.bam/q30_very_end.bam}
  depth_output_filename=${filename/q30.bam/q30_very_end.depth}
  srun -n1 -N1 shifter samtools view -L subtelomere_starts_1_10.bed -b $f > reads_mapping_to_very_end/TALEs_primered/$barcode/$bam_output_filename
  srun -n1 -N1 shifter samtools depth reads_mapping_to_very_end/TALEs_primered/$barcode/$bam_output_filename > reads_mapping_to_very_end/TALEs_primered/$barcode/$depth_output_filename
done



barcode_pat='BC([0-9]+)'
for f in TALEs/alignment_outputs/*191010*primary.bam
do
  [[ $f =~ $barcode_pat|none ]]
  barcode=${BASH_REMATCH[0]}
  mkdir reads_mapping_to_very_end/TALEs/$barcode
  filename=$(basename $f)
  bam_output_filename=${filename/primary.bam/primary_very_end.bam}
  depth_output_filename=${filename/primary.bam/primary_very_end.depth}
  srun -n1 -N1 shifter samtools view -L subtelomere_starts_1_10.bed -b $f > reads_mapping_to_very_end/TALEs/$barcode/$bam_output_filename
  srun -n1 -N1 shifter samtools depth reads_mapping_to_very_end/TALEs/$barcode/$bam_output_filename > reads_mapping_to_very_end/TALEs/$barcode/$depth_output_filename
done

barcode_pat='BC([0-9]+)'
for f in TALEs/alignment_outputs/*191010*q30.bam
do
  [[ $f =~ $barcode_pat|none ]]
  barcode=${BASH_REMATCH[0]}
  mkdir reads_mapping_to_very_end/TALEs/$barcode
  filename=$(basename $f)
  bam_output_filename=${filename/q30.bam/q30_very_end.bam}
  depth_output_filename=${filename/q30.bam/q30_very_end.depth}
  srun -n1 -N1 shifter samtools view -L subtelomere_starts_1_10.bed -b $f > reads_mapping_to_very_end/TALEs/$barcode/$bam_output_filename
  srun -n1 -N1 shifter samtools depth reads_mapping_to_very_end/TALEs/$barcode/$bam_output_filename > reads_mapping_to_very_end/TALEs/$barcode/$depth_output_filename
done




