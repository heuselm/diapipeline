# diapipeline
Targeted analysis of DIA mass spectrometry data

# Prep: Get docker and software in containers
- Requirement: Working docker installation
- e.g. https://www.digitalocean.com/community/tutorials/how-to-install-and-use-docker-on-ubuntu-16-04#step-1-%E2%80%94-installing-docker

check if it's working..
´´´docker run hello-world´´´

get the containerized tools
´´´docker pull biocontainers/diau-umpire´´´
´´´docker pull openswath/openswath:0.1.2´´´

# Step 01: Build peptide query parameter library from DDA-MS files (mzXML), script 01_DDA_libraryCreation.sh

1.2 cd diapipeline

1.3 mkdir data_dda

1.4 cp (MYDDADATA*.mzXML) data_dda/

1.5.

  a Wanna use standard settings and human swissprot db? -> ./01_DDA_libraryCreation.sh
  
  b Wanna use a different fasta protein sequence database?
  
    0 your db must contain DECOY_ prefixed decoy sequences
  
    1 cp /myfasta/xyz.fasta data_library
  
    2 mv data_library/library_fwd_with_decoys.fasta data_library/library_fwd_with_decoys_backup.fasta
  
    3 cp data_library/xyz.fasta data_library/library_fwd_with_decoys.fasta
  
    (Alternatively you can also edit the params/comet.params file to point to your fasta, but it must be in the diapipeline folder for the dockers to see it)
  
    4 -> ./01_DDA_libraryCreation.sh
  c Wanna use different parameters for the comet searches? Check out the params folder & edit params/comet.params
  
loading spectra -> good. Anything else -> not good.

1.6 Check output in results/library!

# Step 02: Query peptide abundance in DIA-MS files (mzXML), script 02_OpenSwathJumboProphetTric.sh
## Note: Pyprophet cannot be run with docker exec; it is required to attach to the docker container and issue the commands within the docker container's shell..

2.1 mkdir data_dia

2.2 cp (MYDIADATA*.mzXML) data_dia/

2.3 (Manually) Go through the steps in script 02_, attach to docker and issue the commands from within..

## Note: This workflow is based on learning one global scoring function for a given dataset. That's good to obtain stabilized scoring across heterogeneous sample sets (e.g. AP-MS baits vs. controls, fractionated samples such as SEC-SWATH-MS)

For details on parameterization, check out the original papers and docs..
http://openswath.org/

Congrats, you ran the SwathPipeline!
Check output in results folder

Inspired by & thanks to http://openswath.org/en/latest/docs/docker.html
& https://github.com/OpenSWATH/workflows
