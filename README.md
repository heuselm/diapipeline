# diapipeline
Targeted analysis of DIA mass spectrometry data

# Step 01: Build peptide query parameter library from DDA-MS files (mzXML)
1.1 pull this repository git clone https://github.com/heuselm/diapipeline.git

1.2 cd diapipeline

1.3 mkdir data_dda

1.4 cp (MYDATA*.mzXML) data_dda/

1.5.

  a Wanna use standard settings and human swissprot db? -> ./01_01_DDA_libraryCreation.sh
  
  b Wanna use a different fasta protein sequence database?
  
  0 your db must contain DECOY_ prefixed decoy sequences
  
  1 cp /myfasta/xyz.fasta data_library
  
  2 mv data_library/library_fwd_with_decoys.fasta data_library/library_fwd_with_decoys_backup.fasta
  
  3 cp data_library/xyz.fasta data_library/library_fwd_with_decoys.fasta
  
  (Alternatively you can also edit the params/comet.params file to point to your fasta, but it must be in the diapipeline folder for the dockers to see it)
  
  4 -> ./01_DDA_libraryCreation.sh
  

loading spectra -> good. Anything else -> not good.

1.6 Check output in results/library!

Inspired by & thanks to http://openswath.org/en/latest/docs/docker.html
& https://github.com/OpenSWATH/workflows
