#!/bin/bash

# DIA proteomics pipeline script
# modified from C Karlsson's script, 26.10.2018
# mhe 27.10.2018
###############################################

# log
# touch current_log.txt
# echo swath analysis log > current_log.txt

# Prep: Get Software Containers
###############################
# Requirement: Working docker installation
# e.g. https://www.digitalocean.com/community/tutorials/how-to-install-and-use-docker-on-ubuntu-16-04#step-1-%E2%80%94-installing-docker
# Test if it is working correctly:
docker run hello-world

# get the containerized tools
docker pull biocontainers/diau-umpire # also contains comet, msgf+ and tpp
docker pull openswath/openswath:0.1.2

# remove containers that may exist in these names
# watch out, modifications made inside these containers will be lost!
docker rm ddalibcreate
docker rm openswath

# spawn containers on host machine
docker run -u 0 -dit --name ddalibcreate -v /media/sf_D_DRIVE/DataAnalysis/SwathPipeline/:/data biocontainers/dia-umpire
docker run -u 0 -dit --name openswath -v /media/sf_D_DRIVE/DataAnalysis/SwathPipeline/:/data openswath/openswath:0.1.2

# say hello!
docker exec ddalibcreate echo hi there, ddalibcreate container is happy and alive
docker exec openswath echo hi there, openswath container is happy and alive

# STEP 2 & 3 of the OpenSwath pipeline

###########################################################
# STEP 2: OPENSWATH:
# Query library peptides in SWATH data by OpenSwathWorkflow
###########################################################
# Note: the following commands have to be run attached to the container
# While OpenswathWorkflow runs with docker exec openswath,
# Pyprophet has a problem with that TODO contact GR
###########################################################
docker attach openswath
cd data/

# OpenSwathWorkflow
mkdir /data/results/openswath
for file in data_dia/*ML; do \
bname=$(echo ${file##*/} | cut -f 1 -d '.'); \
OpenSwathWorkflow \
-in /data/data_dia/$bname.mzXML \
-tr /data/results/library/SSLibrary_target_decoy.pqp \
-tr_irt /data/data_library/cirtkit.TraML \
-min_upper_edge_dist 1 \
-batchSize 1000 \
-out_osw /data/results/openswath/$bname.osw \
-Scoring:stop_report_after_feature 5 \
-rt_extraction_window 600 \
-mz_extraction_window 30 \
-ppm \
-threads 4 \
-use_ms1_traces \
-Scoring:Scores:use_ms1_mi \
-Scoring:Scores:use_mi_score ; done

##########################################################
# STEP3: PYPROPHET:
# Score and filter OpenSwath results
##########################################################
# Note: This is the "new" pipeline training a global model
# that is then applied to each individual run (stabilized srocing)
# Error models are learned and applied on precursor, peptide and protein level
##############################################################################
# create result folders
mkdir results/pyprophet
mkdir results/pyprophet/jumbomodel

# Train Model: Subsample over all runs and create jumbo model
################################################
pyprophet merge --out=/data/results/pyprophet/jumbomodel/subsampled.osw \
--subsample_ratio=0.5 /data/results/openswath/*.osw
pyprophet score --threads 4 --in=/data/results/pyprophet/jumbomodel/subsampled.osw \
--level=ms1ms2
# for inspection by humans:
pyprophet export --in=/data/results/pyprophet/jumbomodel/subsampled.osw \
--out=/data/results/pyprophet/jumbomodel/subsampled.tsv --format=legacy_merged \
--no-ipf

# Apply: Score all runs using the jumbo model
################################################
pyprophet merge --out=/data/results/pyprophet/allruns.osw \
--subsample_ratio=1 /data/results/openswath/collinsb*.osw
pyprophet score --threads 4 --in=/data/results/pyprophet/allruns.osw \
--level=ms1ms2

pyprophet export --in=/data/results/pyprophet/allruns.osw \
--out=/data/results/pyprophet/allruns.tsv --format=legacy_merged \
--no-ipf

# Create statistical models on peptide and protein level and over
# different contexts (Within each run or over the full dataset)
pyprophet \
peptide --in=/data/results/pyprophet/allruns.osw --context=run-specific \
peptide --in=/data/results/pyprophet/allruns.osw --context=experiment-wide \
peptide --in=/data/results/pyprophet/allruns.osw --context=global \

pyprophet export --in=/data/results/pyprophet/allruns.osw \
--out=/data/results/pyprophet/allruns_pepQvals.tsv --format=legacy_merged \
--no-ipf

pyprophet \
protein --in=/data/results/pyprophet/allruns.osw --context=run-specific \
protein --in=/data/results/pyprophet/allruns.osw --context=experiment-wide \
protein --in=/data/results/pyprophet/allruns.osw --context=global \

pyprophet export --in=/data/results/pyprophet/allruns.osw \
--out=/data/results/pyprophet/allruns_pepQvals_protQvals.tsv --format=legacy_merged \
--no-ipf

pyprophet export --in=/data/results/pyprophet/allruns.osw \
--out=/data/results/pyprophet/allruns_pepQvals_protQvals_matrix.tsv --format=matrix \
--no-ipf

# export scoring pdf
pyprophet export \
--in=/data/results/pyprophet/allruns.osw \
--format=score_plots

# TRIC
# Feature Alignment based on nonlinear RT alignment
feature_alignment.py \
--in /data/results/pyprophet/allruns_pepQvals_protQvals.tsv \
--out /data/results/pyprophet/feature_alignment.tsv \
--out_matrix /data/results/pyprophet/feature_alignment_matrix.tsv \
--method LocalMST \
--realign_method lowess \
--max_rt_diff 60 \
--mst:useRTCorrection True \
--mst:Stdev_multiplier 3.0 \
--target_fdr 0.01 \
--fdr_cutoff 0.01 \
--max_fdr_quality 0.05

exit

# Congrats, you ran the SwathPipeline!
# Check output in results folder