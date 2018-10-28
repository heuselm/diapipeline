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

# STEP 1: DDALIBCREATE:
# DDA search to create sample-specific library
##############################################
mkdir results
mkdir results/library
for file in data_dda/*.mzXML; do \
docker exec ddalibcreate comet -Pparams/comet.params $file ;done

# copy result files to result folder
for file in data_dda/*pep.xml; do \
bname=$(echo ${file##*/} | cut -f 1 -d '.') && \
mv data_dda/$bname.pep.xml results/library && \
mv data_dda/$bname.txt results/library ; done

# run tpp peptideprophet per each DDA MSrun
for file in results/library/*pep.xml; do \
fname=$(echo ${file##*/}); \
docker exec ddalibcreate xinteract -dDECOY_ -OARPld -Nresults/library/"interact-"$fname $file; done

# OPTIONAL, when desired run proteinProphet and conversion to mzIdentML for
# full submission to PRIDE 
# docker exec ddalibcreate idconvert --mzIdentML results/library/interact*.pep.xml
# .mzid files will be in top folder

# combine by InterProphetParser (iprophet)
docker exec ddalibcreate InterProphetParser DECOY=DECOY_ \
results/library/interact*.pep.xml results/library/iprophet.pep.xml

docker exec ddalibcreate ProteinProphet results/library/iprophet.pep.xml \
results/library/iprophet.prot.xml IPROPHET

# Perform Mayu FDR estimation
docker exec ddalibcreate Mayu.pl -A results/library/iprophet.pep.xml \
-C data_library/library_fwd_with_decoys.fasta -E DECOY_ -I 2 -G 0.01 -H 100 \
-M Mayu -verbose -P pepFDR=0.01:td && mv Mayu* results/library
# prints out PSM results at designated FDR that can then be used to obtain the 
# IP (i-probability) necessary to filter the library

# Retrieve maximal iprobability from mayu output table used for filtering
ip_cutoff=$(cat results/library/Mayu_psm_pepFDR0.01_td_1.07.csv \
	| cut -d ',' -f5 | sort | head -n 1)
$ip_cutoff
echo $ip_cutoff > ip_cutoff.txt

# Spectrast Library generation
##############################
# extract spectra from mzXML and filter based on cutoff and remove decoys
docker exec ddalibcreate spectrast -cNresults/library/SpectrastStep1_all \
-cICID-QTOF \
-cf "Protein! ~ DECOY_" \
-cP$ip_cutoff \
-c_IRTdata_library/cirtkit.txt \
-c_IRR results/library/iprophet.pep.xml

# Build Consensus spectra
docker exec ddalibcreate spectrast -cNresults/library/SpectrastStep2_consensus \
-cICID-QTOF -cf'Protein!~DECOY' \
-cAC -cM results/library/SpectrastStep1_all.splib

# Legacy Script: Spectrast2Tsv:
# convert to tsv, picking the 4-6 highest-intense fragment ions and 
# feeding in the swath window information
# broken
# docker exec openswath spectrast2tsv.py -l 400,2000 -k openswath \
# -s b,y -w/data/data_library/swathwindows_64vw_human.txt -o 4 -n 6  \
# -a /data/results/library/sslib.tsv /data/results/library/SpectrastStep2.sptxt

# Convert to TraML, generate and append decoys and convert to .pqp
docker exec openswath TargetedFileConverter \
-in /data/results/library/SpectrastStep2_consensus.mrm \
-out /data/results/library/SSLibrary_transitionlist_all.TraML

# -in /data/results/library/SpectrastStep3.mrm \ #
docker exec openswath OpenSwathAssayGenerator \
-in /data/results/library/SpectrastStep2_consensus.mrm \
-swath_windows_file /data/data_library/swathwindows_64vw_human_wheader.txt \
-out /data/results/library/SSLibrary_target.TraML

docker exec openswath OpenSwathDecoyGenerator \
-in /data/results/library/SSLibrary_target.TraML \
-out /data/results/library/SSLibrary_target_decoy.TraML \
-method shuffle

# Convert to pqp for osw analysis
docker exec openswath TargetedFileConverter \
-in /data/results/library/SSLibrary_target_decoy.TraML \
-out /data/results/library/SSLibrary_target_decoy.pqp

# And for easy vieweing and processing in other tools to tsv
docker exec openswath TargetedFileConverter \
-in /data/results/library/SSLibrary_target_decoy.TraML \
-out /data/results/library/SSLibrary_target_decoy.tsv

#not running robustly yet..

###########################################################
# STEP 2: OPENSWATH:
# Query library peptides in SWATH data by OpenSwathWorkflow
###########################################################
# OpenSwathWorkflow
mkdir results/openswath
for file in data_dia/*ML; do \
bname=$(echo ${file##*/} | cut -f 1 -d '.'); \
docker exec openswath OpenSwathWorkflow \
-in /data/data_dia/$bname.mzXML \
-tr /data/results/library/SSLibrary_target_decoy.pqp \
-tr_irt /data/data_library/cirtkit.TraML \
-min_upper_edge_dist 1 \
-batchSize 1000 \
-out_osw /data/results/openswath/$bname.osw \
-out_tsv /data/results/openswath/$bname.tsv \
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

# train global model by subsampling from all raw files
mkdir results/pyprophet
mkdir results/pyprophet/jumbomodel
# docker exec openswath pyprophet merge --out=/data/results/pyprophet/jumbomodel/subsampled.osw \
# --subsample_ratio=0.5 /data/results/openswath/*.osw && \
# pyprophet score --threads 4 --in=/data/results/pyprophet/jumbomodel/subsampled.osw \
# --level=ms1ms2

#Error: Invalid value for "infiles": Path "/data/results/openswath/*" does not exist.
# But: Works when entering the container shell and executing the command there
# -> Maybe related to python virtualenvs?

# attach and exit, restart for next steps
# doesn#t run automatically, need to execute manually inside docker..
docker attach openswath

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