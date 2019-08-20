#!/bin/bash
#SBATCH --time=8:00:00
#SBATCH --mem-per-cpu=6G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --image=docker:fcadete/guppy:latest
#SBATCH --workdir=/mnt/nfs/lobo/SALMEIDA-NFS/fcadete/scratch/minion/


for i in $(seq 5 5 95); do

   mkdir subsample_analysis/barcoded/${i}
   mkdir subsample_analysis/barcoded/${i}/${SLURM_ARRAY_TASK_ID}
   mkdir subsample_analysis/depths/${i}
   mkdir subsample_analysis/depths/${i}/${SLURM_ARRAY_TASK_ID}
   mkdir subsample_analysis/alignments/${i}
   mkdir subsample_analysis/alignments/${i}/${SLURM_ARRAY_TASK_ID}
   mkdir subsample_analysis/reads/${i}
   mkdir subsample_analysis/porechop/${i}

   mkdir ~/minion/subsample_analysis/barcoded/${i}
   mkdir ~/minion/subsample_analysis/barcoded/${i}/${SLURM_ARRAY_TASK_ID}
   mkdir ~/minion/subsample_analysis/depths/${i}
   mkdir ~/minion/subsample_analysis/depths/${i}/${SLURM_ARRAY_TASK_ID}
   mkdir ~/minion/subsample_analysis/alignments/${i}
   mkdir ~/minion/subsample_analysis/alignments/${i}/${SLURM_ARRAY_TASK_ID}
   mkdir ~/minion/subsample_analysis/reads/${i}
   mkdir ~/minion/subsample_analysis/porechop/${i}


   srun -n1 -N1 shifter /seqtk/seqtk sample -s $SLURM_ARRAY_TASK_ID large_U2OS_sample_reads.fastq.gz $(expr 1638648 / 4 / 100 \* $i) | gzip > subsample_analysis/reads/${i}/${SLURM_ARRAY_TASK_ID}_reads.fastq.gz

   srun -n1 -N1 shifter /Porechop/porechop-runner.py -i subsample_analysis/reads/${i}/${SLURM_ARRAY_TASK_ID}_reads.fastq.gz \
                                             -b subsample_analysis/barcoded/${i}/${SLURM_ARRAY_TASK_ID}/ \
                                             --barcode_diff 1 \
                                             -t 8 --verbosity 3 > subsample_analysis/porechop/${i}/${SLURM_ARRAY_TASK_ID}.porechop_output

   srun -n1 -N1 shifter python3 basecall_scripts/parse_porechop_output.py subsample_analysis/porechop/${i}/${SLURM_ARRAY_TASK_ID}.porechop_output

   srun -n1 -N1 gzip subsample_analysis/porechop/${i}/${SLURM_ARRAY_TASK_ID}.porechop_*



   barcoded=( "Joana_polyA_VNP" "Joana_polyA_VNP_no_PCR" "Joana_TERRA_VNP" "Joana_TERRA_VNP_no_PCR" "none" )

   for barcode in "${barcoded[@]}"; do

      srun -n1 -N1 shifter bash -c "/minimap2/minimap2 -a -x map-ont -t 8 \
                                                  references/ConcatenatedFASTAAassemblies_hTel.txt \
                                                  subsample_analysis/barcoded/${i}/${SLURM_ARRAY_TASK_ID}/${barcode}.fastq.gz |
                                                samtools sort -@ 8 - | samtools view -bF 256 - > subsample_analysis/alignments/${i}/${SLURM_ARRAY_TASK_ID}/${barcode}_on_rhietman_mapont_primary.bam &&
                                                samtools depth subsample_analysis/alignments/${i}/${SLURM_ARRAY_TASK_ID}/${barcode}_on_rhietman_mapont_primary.bam | gzip > subsample_analysis/depths/${i}/${SLURM_ARRAY_TASK_ID}/${barcode}_on_rhietman_mapont_primary.depth.gz"

    done

   srun -n1 -N1 bash -c "mv subsample_analysis/barcoded/${i}/${SLURM_ARRAY_TASK_ID} ~/minion/subsample_analysis/barcoded/${i}/ ; \
                                     mv subsample_analysis/alignments/${i}/${SLURM_ARRAY_TASK_ID} ~/minion/subsample_analysis/alignments/${i}/ ; \
                                     mv subsample_analysis/depths/${i}/${SLURM_ARRAY_TASK_ID} ~/minion/subsample_analysis/depths/${i}/ ; \
                                     mv subsample_analysis/reads/${i}/${SLURM_ARRAY_TASK_ID}_reads.fastq.gz ~/minion/subsample_analysis/reads/${i}/ ; \
                                     mv subsample_analysis/porechop/${i}/${SLURM_ARRAY_TASK_ID}.porechop_* ~/minion/subsample_analysis/porechop/${i}/"

done

echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,NodeList,NNodes,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID

