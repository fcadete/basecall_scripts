#!/bin/bash
#SBATCH --job-name=telomerehunter
#SBATCH --time=48:00:00
#SBATCH --output=polyA_191023_telomerehunter.output
#SBATCH --mem-per-cpu=6G
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=6
#SBATCH --image=docker:fcadete/guppy:latest
#SBATCH --workdir=/mnt/nfs/lobo/SALMEIDA-NFS/fcadete/scratch/minion/

mkdir PCR_amplified/telomere_hunter_results
mkdir PCR_amplified/telomere_hunter_results/sorted_bams

for f in PCR_amplified/alignment_outputs/*mapont.bam; do

#   echo $f
   filename=$(basename $f)
#   echo $filename
   sorted_bam=${filename//_on_rhietman_mapont.bam/_on_rhietman_mapont_sorted.bam}
#   echo $sorted_bam
   samplename=${filename//_on_rhietman_mapont.bam/}
#   echo $samplename

   rm -r PCR_amplified/telomere_hunter_results/$samplename
   
   srun -n1 -N1 shifter --image fcadete/guppy bash -c "samtools sort -@ 6 $f > PCR_amplified/telomere_hunter_results/sorted_bams/$sorted_bam
                                                       telomerehunter -ibt PCR_amplified/telomere_hunter_results/sorted_bams/$sorted_bam \
                                                                      -o PCR_amplified/telomere_hunter_results \
                                                                      -p $samplename \
                                                                      -b subtelomere_arbitrary_bands.txt \
                                                                      --parallel \
                                                                      -rt 3" &

done
wait

mkdir PCR_amplified/telomere_hunter_results/depths_from_telomeric_reads

for f in PCR_amplified/telomere_hunter_results/*/*/*_intratelomeric.bam; do

   filename=$(basename $f)
   samplename=${filename//_filtered_intratelomeric.bam/}
   depthFile=${filename//_filtered_intratelomeric.bam/_filtered_intratelomeric_primary.depth}

   srun -n1 -N1 shifter --image fcadete/guppy bash -c "samtools view -bF 256 $f | samtools sort - | samtools depth - > PCR_amplified/telomere_hunter_results/depths_from_telomeric_reads/$depthFile" &

done
wait


echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,NodeList,NNodes,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID


