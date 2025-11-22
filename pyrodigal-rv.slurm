#!/bin/bash -l
#SBATCH --nodes="1"
#SBATCH -c "10"
#SBATCH --cluster=wice
#SBATCH --time="2:00:00"
#SBATCH -A lp_jm_virome_group
#SBATCH -p batch_sapphirerapids

cd $VSC_SCRATCH/pyrodigal_gv_models/refseq_genomes

conda activate gb_submitter

pyrodigal-rv -i genomic_no_N.fna -o predictions_closed.gff -j 10 -c --pool process
pyrodigal -i genomic_no_N.fna -o vanilla_predictions_closed.gff -j 10 -g 1 -c --pool process
pyrodigal-gv -i genomic_no_N.fna -o gv_predictions_closed.gff -j 10 -c --pool process

hyperfine --export-markdown hyperfine.txt \
    --export-csv hyperfine.csv \
    --export-json hyperfine.json \
    -L cmd pyrodigal,pyrodigal-rv,pyrodigal-gv \
    --command-name pyrodigal pyrodigal-rv pyrodigal-gv \
    '{cmd} -i genomic_no_N.fna -o gff -j 10 -c --pool process'