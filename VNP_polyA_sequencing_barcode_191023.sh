#!/bin/bash
#SBATCH --job-name=guppy_VNP_polyA_sequencing_barcode_191023
#SBATCH --time=48:00:00
#SBATCH --output=guppy_VNP_polyA_sequencing_barcode_191023.output
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
#                                     -i 191023_Ninth_run/191023_polyA_sequencing_PCR_barcode/ \
#                                     -s 191023_Ninth_run/191023_polyA_sequencing_PCR_barcode/fastq_guppy
#
#
#srun -n1 -N1 shifter Rscript /MinIONQC.R -i 191023_Ninth_run/191023_polyA_sequencing_PCR_barcode/fastq_guppy/sequencing_summary.txt \
#                                         -o 191023_Ninth_run/191023_polyA_sequencing_PCR_barcode/fastq_guppy/qc_plots
#
#
#srun -n1 -N1 shifter /Porechop/porechop-runner.py -i 191023_Ninth_run/191023_polyA_sequencing_PCR_barcode/fastq_guppy/pass/ \
#                                          -b 191023_Ninth_run/191023_polyA_sequencing_PCR_barcode/fastq_guppy/barcode_separated \
#                                          --check_reads 50000 -t 40 --verbosity 3 > guppy_VNP_polyA_sequencing_barcode_191023.porechop_output
#
#srun -n1 -N1 shifter python3 basecall_scripts/parse_porechop_output.py guppy_VNP_polyA_sequencing_barcode_191023.porechop_output
#
#
#barcoded=( "BC01" "BC02" "BC03" "BC04" "BC05" "BC06" "BC07" "BC08" "BC09" "BC10" "BC11" "BC12" "none" )
#
#
#for barcode in "${barcoded[@]}"; do
#
#srun -n1 -N1 shifter /FastQC/fastqc -k 6 --nano --threads 40 191023_Ninth_run/191023_polyA_sequencing_PCR_barcode/fastq_guppy/barcode_separated/$barcode.fastq
#
#srun -n1 -N1 shifter bash -c "/minimap2/minimap2 -a -x map-ont -t 40 \
#                                     references/ConcatenatedFASTAAassemblies_hTel.txt \
#                                     191023_Ninth_run/191023_polyA_sequencing_PCR_barcode/fastq_guppy/barcode_separated/$barcode.fastq | \
#                                     samtools view -b - > PCR_amplified/alignment_outputs/VNP_polyA_sequencing_barcode_191023_${barcode}_on_rhietman_mapont.bam"
#
#
#srun -n1 -N1 shifter bash -c "samtools sort -@ 40 PCR_amplified/alignment_outputs/VNP_polyA_sequencing_barcode_191023_${barcode}_on_rhietman_mapont.bam | samtools view -bq 30 - | samtools depth - > PCR_amplified/alignment_outputs/VNP_polyA_sequencing_barcode_191023_${barcode}_on_rhietman_mapont_q30.depth"
#
#srun -n1 -N1 shifter bash -c "samtools sort -@ 40 PCR_amplified/alignment_outputs/VNP_polyA_sequencing_barcode_191023_${barcode}_on_rhietman_mapont.bam | samtools view -bF 256 - | samtools depth - > PCR_amplified/alignment_outputs/VNP_polyA_sequencing_barcode_191023_${barcode}_on_rhietman_mapont_primary.depth"
#
#srun -n1 -N1 shifter bash -c "samtools sort -@ 40 PCR_amplified/alignment_outputs/VNP_polyA_sequencing_barcode_191023_${barcode}_on_rhietman_mapont.bam | samtools view -bq 30 - > PCR_amplified/alignment_outputs/VNP_polyA_sequencing_barcode_191023_${barcode}_on_rhietman_mapont_q30.bam; samtools index PCR_amplified/alignment_outputs/VNP_polyA_sequencing_barcode_191023_${barcode}_on_rhietman_mapont_q30.bam"
#
#srun -n1 -N1 shifter bash -c "samtools sort -@ 40 PCR_amplified/alignment_outputs/VNP_polyA_sequencing_barcode_191023_${barcode}_on_rhietman_mapont.bam | samtools view -bF 256 - > PCR_amplified/alignment_outputs/VNP_polyA_sequencing_barcode_191023_${barcode}_on_rhietman_mapont_primary.bam; samtools index PCR_amplified/alignment_outputs/VNP_polyA_sequencing_barcode_191023_${barcode}_on_rhietman_mapont_primary.bam"
#
#
#
#mkdir 191023_Ninth_run/191023_polyA_sequencing_PCR_barcode/fastq_guppy/barcode_separated/jo_primers_verbatim
#mkdir PCR_amplified/primer_separation/jo_primers_verbatim
#
#srun -n1 -N1 shifter --image=docker:fcadete/guppy_dockerfile_pcr_polya_verbatim:latest /Porechop/porechop-runner.py -i 191023_Ninth_run/191023_polyA_sequencing_PCR_barcode/fastq_guppy/barcode_separated/$barcode.fastq \
#                                          -b 191023_Ninth_run/191023_polyA_sequencing_PCR_barcode/fastq_guppy/barcode_separated/jo_primers_verbatim/$barcode \
#                                          --barcode_diff 1 --check_reads 50000 -t 40 --verbosity 3 > PCR_amplified/primer_separation/jo_primers_verbatim/guppy_VNP_polyA_sequencing_barcode_191023_$barcode.porechop_output
#
#srun -n1 -N1 shifter python3 basecall_scripts/parse_porechop_output.py PCR_amplified/primer_separation/jo_primers_verbatim/guppy_VNP_polyA_sequencing_barcode_191023_$barcode.porechop_output
#
#mkdir 191023_Ninth_run/191023_polyA_sequencing_PCR_barcode/fastq_guppy/barcode_separated/jo_primers_revcomp
#mkdir PCR_amplified/primer_separation/jo_primers_revcomp
#
#srun -n1 -N1 shifter --image=docker:fcadete/guppy_dockerfile_pcr_polya_revcomp:latest /Porechop/porechop-runner.py -i 191023_Ninth_run/191023_polyA_sequencing_PCR_barcode/fastq_guppy/barcode_separated/$barcode.fastq \
#                                          -b 191023_Ninth_run/191023_polyA_sequencing_PCR_barcode/fastq_guppy/barcode_separated/jo_primers_revcomp/$barcode \
#                                          --barcode_diff 1 --check_reads 50000 -t 40 --verbosity 3 > PCR_amplified/primer_separation/jo_primers_revcomp/guppy_VNP_polyA_sequencing_barcode_191023_$barcode.porechop_output
#
#srun -n1 -N1 shifter python3 basecall_scripts/parse_porechop_output.py PCR_amplified/primer_separation/jo_primers_revcomp/guppy_VNP_polyA_sequencing_barcode_191023_$barcode.porechop_output
#
#mkdir 191023_Ninth_run/191023_polyA_sequencing_PCR_barcode/fastq_guppy/barcode_separated/jo_primers_forward_verbatim_reverse_complement
#mkdir PCR_amplified/primer_separation/jo_primers_forward_verbatim_reverse_complement
#
#srun -n1 -N1 shifter --image=docker:fcadete/guppy_dockerfile_pcr_polya_forward_verbatim_reverse_complement /Porechop/porechop-runner.py -i 191023_Ninth_run/191023_polyA_sequencing_PCR_barcode/fastq_guppy/barcode_separated/$barcode.fastq \
#                                          -b 191023_Ninth_run/191023_polyA_sequencing_PCR_barcode/fastq_guppy/barcode_separated/jo_primers_forward_verbatim_reverse_complement/$barcode \
#                                          --adapter_threshold 80.0 --untrimmed --barcode_diff 1 --check_reads 50000 -t 40 --verbosity 3 > PCR_amplified/primer_separation/jo_primers_forward_verbatim_reverse_complement/guppy_VNP_polyA_sequencing_barcode_191023_$barcode.porechop_output &
#
#done
#
#wait
#
#mkdir PCR_amplified/polyA_removed_reads
#mkdir PCR_amplified/polyA_removed_reads/fastqs
#
#srun -n1 -N1 shifter --image pcerqueira/bbtools:38.44 /NGStools/bbmap/bbduk.sh in=191023_Ninth_run/191023_polyA_sequencing_PCR_barcode/fastq_guppy/barcode_separated/jo_primers_forward_verbatim_reverse_complement/BC01/SSP_1q.fastq \
#                                                                               out=PCR_amplified/polyA_removed_reads/fastqs/BC01_SSP_1q_polyA_removed.fastq \
#                                                                               literal=AAAAAAAAAAA k=7 qin=33 ktrim=r
#
#srun -n1 -N1 shifter --image pcerqueira/bbtools:38.44 /NGStools/bbmap/bbduk.sh in=191023_Ninth_run/191023_polyA_sequencing_PCR_barcode/fastq_guppy/barcode_separated/jo_primers_forward_verbatim_reverse_complement/BC02/SSP_2p.fastq \
#                                                                               out=PCR_amplified/polyA_removed_reads/fastqs/BC02_SSP_2p_polyA_removed.fastq \
#                                                                               literal=AAAAAAAAAAA k=7 qin=33 ktrim=r
#
#srun -n1 -N1 shifter --image pcerqueira/bbtools:38.44 /NGStools/bbmap/bbduk.sh in=191023_Ninth_run/191023_polyA_sequencing_PCR_barcode/fastq_guppy/barcode_separated/jo_primers_forward_verbatim_reverse_complement/BC03/SSP_7p.fastq \
#                                                                               out=PCR_amplified/polyA_removed_reads/fastqs/BC03_SSP_7p_polyA_removed.fastq \
#                                                                               literal=AAAAAAAAAAA k=7 qin=33 ktrim=r
#
#srun -n1 -N1 shifter --image pcerqueira/bbtools:38.44 /NGStools/bbmap/bbduk.sh in=191023_Ninth_run/191023_polyA_sequencing_PCR_barcode/fastq_guppy/barcode_separated/jo_primers_forward_verbatim_reverse_complement/BC04/SSP_9p.fastq \
#                                                                               out=PCR_amplified/polyA_removed_reads/fastqs/BC04_SSP_9p_polyA_removed.fastq \
#                                                                               literal=AAAAAAAAAAA k=7 qin=33 ktrim=r
#
#srun -n1 -N1 shifter --image pcerqueira/bbtools:38.44 /NGStools/bbmap/bbduk.sh in=191023_Ninth_run/191023_polyA_sequencing_PCR_barcode/fastq_guppy/barcode_separated/jo_primers_forward_verbatim_reverse_complement/BC05/SSP_10q_16q.fastq \
#                                                                               out=PCR_amplified/polyA_removed_reads/fastqs/BC05_SSP_10q_16q_polyA_removed.fastq \
#                                                                               literal=AAAAAAAAAAA k=7 qin=33 ktrim=r
#
#srun -n1 -N1 shifter --image pcerqueira/bbtools:38.44 /NGStools/bbmap/bbduk.sh in=191023_Ninth_run/191023_polyA_sequencing_PCR_barcode/fastq_guppy/barcode_separated/jo_primers_forward_verbatim_reverse_complement/BC06/SSP_12q.fastq \
#                                                                               out=PCR_amplified/polyA_removed_reads/fastqs/BC06_SSP_12q_polyA_removed.fastq \
#                                                                               literal=AAAAAAAAAAA k=7 qin=33 ktrim=r
#
#srun -n1 -N1 shifter --image pcerqueira/bbtools:38.44 /NGStools/bbmap/bbduk.sh in=191023_Ninth_run/191023_polyA_sequencing_PCR_barcode/fastq_guppy/barcode_separated/jo_primers_forward_verbatim_reverse_complement/BC07/SSP_13q.fastq \
#                                                                               out=PCR_amplified/polyA_removed_reads/fastqs/BC07_SSP_13q_polyA_removed.fastq \
#                                                                               literal=AAAAAAAAAAA k=7 qin=33 ktrim=r
#
#srun -n1 -N1 shifter --image pcerqueira/bbtools:38.44 /NGStools/bbmap/bbduk.sh in=191023_Ninth_run/191023_polyA_sequencing_PCR_barcode/fastq_guppy/barcode_separated/jo_primers_forward_verbatim_reverse_complement/BC08/SSP_15q.fastq \
#                                                                               out=PCR_amplified/polyA_removed_reads/fastqs/BC08_SSP_15q_polyA_removed.fastq \
#                                                                               literal=AAAAAAAAAAA k=7 qin=33 ktrim=r
#
#srun -n1 -N1 shifter --image pcerqueira/bbtools:38.44 /NGStools/bbmap/bbduk.sh in=191023_Ninth_run/191023_polyA_sequencing_PCR_barcode/fastq_guppy/barcode_separated/jo_primers_forward_verbatim_reverse_complement/BC09/SSP_16p.fastq \
#                                                                               out=PCR_amplified/polyA_removed_reads/fastqs/BC09_SSP_16p_polyA_removed.fastq \
#                                                                               literal=AAAAAAAAAAA k=7 qin=33 ktrim=r
#
#srun -n1 -N1 shifter --image pcerqueira/bbtools:38.44 /NGStools/bbmap/bbduk.sh in=191023_Ninth_run/191023_polyA_sequencing_PCR_barcode/fastq_guppy/barcode_separated/jo_primers_forward_verbatim_reverse_complement/BC10/SSP_16q_10q.fastq \
#                                                                               out=PCR_amplified/polyA_removed_reads/fastqs/BC10_SSP_16q_10q_polyA_removed.fastq \
#                                                                               literal=AAAAAAAAAAA k=7 qin=33 ktrim=r
#
#srun -n1 -N1 shifter --image pcerqueira/bbtools:38.44 /NGStools/bbmap/bbduk.sh in=191023_Ninth_run/191023_polyA_sequencing_PCR_barcode/fastq_guppy/barcode_separated/jo_primers_forward_verbatim_reverse_complement/BC11/SSP_18p.fastq \
#                                                                               out=PCR_amplified/polyA_removed_reads/fastqs/BC11_SSP_18p_polyA_removed.fastq \
#                                                                               literal=AAAAAAAAAAA k=7 qin=33 ktrim=r
#
#srun -n1 -N1 shifter --image pcerqueira/bbtools:38.44 /NGStools/bbmap/bbduk.sh in=191023_Ninth_run/191023_polyA_sequencing_PCR_barcode/fastq_guppy/barcode_separated/jo_primers_forward_verbatim_reverse_complement/BC12/SSP_20q.fastq \
#                                                                               out=PCR_amplified/polyA_removed_reads/fastqs/BC12_SSP_20q_polyA_removed.fastq \
#                                                                               literal=AAAAAAAAAAA k=7 qin=33 ktrim=r
#
mkdir PCR_amplified/polyA_removed_reads/bams
mkdir PCR_amplified/polyA_removed_reads/depths
mkdir PCR_amplified/polyA_removed_reads/bams_to_human

for f in PCR_amplified/polyA_removed_reads/fastqs/*fastq; do
   
#   srun -n1 -N1 shifter bash -c "/minimap2/minimap2 -a -x map-ont -t 40 \
#                                        references/ConcatenatedFASTAAassemblies_hTel.txt \
#                                        $f | \
#                                        samtools view -b - > ${f//fastq/bam}"
#   
   filename="$(basename $f)"
#   
#   srun -n1 -N1 shifter bash -c "samtools sort -@ 40 ${f//fastq/bam} | samtools view -bq 30 - | samtools depth - > PCR_amplified/polyA_removed_reads/depths/${filename/.fastq/_q30.depth}"
#   
#   srun -n1 -N1 shifter bash -c "samtools sort -@ 40 ${f//fastq/bam} | samtools view -bF 256 - | samtools depth - > PCR_amplified/polyA_removed_reads/depths/${filename/.fastq/_primary.depth}"
#
   srun -n1 -N1 shifter bash -c "/minimap2/minimap2 -a -x map-ont -t 40 \
                                        references/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
                                        $f | \
                                        samtools view -b - > PCR_amplified/polyA_removed_reads/bams_to_human/${filename/.fastq/_on_human_mapont.bam}"

    srun -n1 -N1 shifter bash -c "/minimap2/minimap2 -a -x splice -t 40 \
                                        references/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
                                        $f | \
                                        samtools view -b - > PCR_amplified/polyA_removed_reads/bams_to_human/${filename/.fastq/_on_human_splice.bam}"
   
done

echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,NodeList,NNodes,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID

