# diapipeline

A scaffold for **containerized, library-based DIA/SWATH-MS analysis**: build a spectral
(peptide query parameter) library from DDA data, then query that library in DIA data with
the [OpenSWATH](http://openswath.org/) toolchain. All tools run in Docker, so the only
host requirement is a working Docker installation.

```
DDA mzXML ──► Comet ──► PeptideProphet / iProphet ──► Mayu (FDR) ──► SpectraST ──► OpenSwathAssayGenerator
                                                                                    + DecoyGenerator
                                                                                          │
                                                                                   assay library (.pqp)
                                                                                          │
DIA mzXML ─────────────────────────────────────────► OpenSwathWorkflow ──► PyProphet ──► TRIC alignment
                                                                                          │
                                                                              quantitative peptide/protein matrix
```

## Repository layout

| Path | Content |
|---|---|
| `01_DDA_libraryCreation.sh` | Step 1: DDA search, statistical validation, library generation |
| `02_OpenSwathJumboProphetTric.sh` | Step 2: OpenSWATH extraction, global PyProphet scoring, TRIC |
| `dockers/` | Dockerfiles for the OpenSWATH and DIA-Umpire images |
| `params/` | Comet parameter sets (QE / TripleTOF / QTOF, tryptic / semitryptic, TMT) |
| `data_library/` | Default human SwissProt FASTA (+decoys), iRT / CiRT kits, SWATH window schemes |

## Quick start

```bash
# 0. prerequisites
docker run hello-world
docker pull openswath/openswath:0.1.2
docker pull biocontainers/dia-umpire

git clone https://github.com/heuselm/diapipeline && cd diapipeline

# 1. library from DDA data
mkdir data_dda && cp /path/to/*.mzXML data_dda/
./01_DDA_libraryCreation.sh          # output: results/library/

# 2. query the library in DIA data
mkdir data_dia && cp /path/to/*.mzXML data_dia/
# step through 02_OpenSwathJumboProphetTric.sh (see note below)
```

### Customizing

- **Other FASTA:** must contain `DECOY_`-prefixed decoys and live inside the repo folder
  (Docker mounts it). Replace `data_library/library_fwd_with_decoys.fasta`, or point
  `params/comet.params` at your file.
- **Search settings:** pick or edit a parameter set in `params/`.
- **SWATH windows:** choose the matching scheme in `data_library/swathwindows_*.txt`.

### Note on step 2

PyProphet is not run via `docker exec`; attach to the container and issue the commands
in script 02 from its shell.

The workflow learns **one global scoring model per dataset** ("jumbo" PyProphet), which
stabilizes scoring across heterogeneous sample sets such as AP-MS baits vs. controls or
fractionated samples (e.g. SEC-SWATH-MS).

## Status

A scaffold from 2017–2020 OpenSWATH work, kept for reference and reuse. Tool versions are
pinned to the containers above; for current versions see
[OpenSWATH](http://openswath.org/) and the
[OpenSWATH Snakemake workflows](https://github.com/OpenSWATH/workflows).

## References

- Röst et al. (2014) OpenSWATH. *Nat Biotechnol* 32:219–223.
- Rosenberger et al. (2017) Statistical control of peptide and protein error rates in
  large-scale targeted DIA analyses. *Nat Methods* 14:921–927.
- Röst et al. (2016) TRIC: alignment of targeted DIA-MS data. *Nat Methods* 13:777–783.

---
<sub>README optimized with Claude Opus 5.5.</sub>
