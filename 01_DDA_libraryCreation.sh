#!/bin/bash

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
docker stop ddalibcreate && \
docker stop openswath && \
docker rm ddalibcreate && \
docker rm openswath

# spawn containers on host machine
docker run -u 0 -dit --name ddalibcreate -v $PWD/:/data biocontainers/dia-umpire
docker run -u 0 -dit --name openswath -v $PWD/:/data openswath/openswath:0.1.2

# say hello!
docker exec ddalibcreate echo hi there, ddalibcreate container is happy and alive
docker exec openswath echo hi there, openswath container is happy and alive

# STEP 1: DDALIBCREATE:
# DDA search to create sample-specific library
##############################################
# Make decoy fasta database
docker exec ddalibcreate subsetdb -R -DDECOY_ data_library/library_fwd.fasta && \
cat data_library/library_fwd.fasta library_fwd.fasta.new > data_library/library_fwd_with_decoys.fasta &&\
chmod 777 data_library/library_fwd_with_decoys.fasta

# Search via Comet
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

# Generate different flavors of the Library.. (optimized for DIANN)
docker exec openswath OpenSwathAssayGenerator \
-in /data/results/library/SpectrastStep2_consensus.mrm \
-swath_windows_file /data/data_library/swathwindows.txt \
-precursor_lower_mz_limit 350 \
-precursor_upper_mz_limit 1650 \
-product_lower_mz_limit 150 \
-out /data/results/library/SSLibrary_target.TraML

docker exec openswath OpenSwathAssayGenerator \
-in /data/results/library/SpectrastStep2_consensus.mrm \
-swath_windows_file /data/data_library/swathwindows.txt \
-precursor_lower_mz_limit 350 \
-precursor_upper_mz_limit 1650 \
-product_lower_mz_limit 150 \
-out /data/results/library/SSLibrary_target_6-6t_std.tsv

docker exec openswath OpenSwathAssayGenerator \
-in /data/results/library/SpectrastStep2_consensus.mrm \
-swath_windows_file /data/data_library/swathwindows.txt \
-precursor_lower_mz_limit 350 \
-precursor_upper_mz_limit 1650 \
-product_lower_mz_limit 150 \
-enable_detection_unspecific_losses \
-out /data/results/library/SSLibrary_target_6-6t_losses.tsv

docker exec openswath OpenSwathAssayGenerator \
-in /data/results/library/SpectrastStep2_consensus.mrm \
-min_transitions 4 \
-max_transitions 6 \
-allowed_fragment_types b,y \
-enable_detection_unspecific_losses \
-precursor_lower_mz_limit 350 \
-precursor_upper_mz_limit 1650 \
-product_lower_mz_limit 150 \
-swath_windows_file /data/data_library/swathwindows.txt \
-out /data/results/library/SSLibrary_target_4-6t_losses.tsv

docker exec openswath OpenSwathAssayGenerator \
-in /data/results/library/SpectrastStep2_consensus.mrm \
-min_transitions 4 \
-max_transitions 24 \
-allowed_fragment_types b,y \
-enable_detection_unspecific_losses \
-precursor_lower_mz_limit 350 \
-precursor_upper_mz_limit 1650 \
-product_lower_mz_limit 150 \
-swath_windows_file /data/data_library/swathwindows.txt \
-out /data/results/library/SSLibrary_target_4-24t_losses.tsv

docker exec openswath OpenSwathAssayGenerator \
-in /data/results/library/SpectrastStep2_consensus.mrm \
-min_transitions 4 \
-max_transitions 24 \
-allowed_fragment_types b,y,x,z \
-enable_detection_unspecific_losses \
-precursor_lower_mz_limit 350 \
-precursor_upper_mz_limit 1650 \
-product_lower_mz_limit 150 \
-swath_windows_file /data/data_library/swathwindows.txt \
-out /data/results/library/SSLibrary_target_4-24t_losses_byxz.tsv

## Classical pathway (5600, 400-1200 mz, b,y, 6-6t)
# Convert to TraML, generate and append decoys and convert to .pqp
docker exec openswath TargetedFileConverter \
-in /data/results/library/SpectrastStep2_consensus.mrm \
-out /data/results/library/SSLibrary_transitionlist_all.TraML

# Decoy generation and conversion (only for the standard parameter library..)
docker exec openswath OpenSwathDecoyGenerator \
-in /data/results/library/SSLibrary_target.TraML \
-out /data/results/library/SSLibrary_target_decoy.TraML \
-method shuffle

# Convert to pqp for osw analysis
docker exec openswath TargetedFileConverter \
-in /data/results/library/SSLibrary_target_decoy.TraML \
-out /data/results/library/SSLibrary_target_decoy.pqp

# And for easy viewing and processing in other tools to tsv
docker exec openswath TargetedFileConverter \
-in /data/results/library/SSLibrary_target_decoy.TraML \
-out /data/results/library/SSLibrary_target_decoy.tsv

# done
echo "done"