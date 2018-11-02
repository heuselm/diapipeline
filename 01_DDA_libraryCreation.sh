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
docker run -u 0 -dit --name ddalibcreate -v $PWD/:/data biocontainers/dia-umpire
docker run -u 0 -dit --name openswath -v $PWD/:/data openswath/openswath:0.1.2

# say hello!
docker exec ddalibcreate echo hi there, ddalibcreate container is happy and alive
docker exec openswath echo hi there, openswath container is happy and alive

# Local Fix 1:
docker exec ddalibcreate gunzip /data/data_dda/*.gz
# Local Fix 2: pre-run sed -i 's/DECOY/DECOY_/g' data_library/library_fwd_with_decoys.fasta
# Local Fix 3: Open permissions inside diau/libcreate container
docker exec ddalibcreate chmod 777 /usr/local/tpp/bin/*

# STEP 1: DDALIBCREATE:
# DDA search to create sample-specific library
##############################################
mkdir results
mkdir results/library
for file in data_dda/*; do \
docker exec ddalibcreate comet -Pparams/comet.params $file ;done

# copy result files to result folder
mv data_dda/*.pep.xml results/library
mv data_dda/*.txt results/library

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

# Make results open access
docker exec openswath chmod -R 777 /data/results

# done
echo "done"
