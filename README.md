# Sepsis Atlas

Publication-oriented analysis code for construction of the
Sepsis Atlas from publicly available single-cell RNA-sequencing
datasets.

## Contents

- `code/R/`: reusable preprocessing and quality-control functions.
- `code/loaders/`: study-specific loading procedures.
- `code/canonical_metadata_freeze_v2/`: metadata harmonization.
- `code/gene_feature_mapping_v1/`: common feature-space construction.
- `code/full_atlas_primary_integration_v1/`: primary Atlas integration.
- `code/full_atlas_annotation_v1/`: annotation, validation,
  reconciliation, global QC, and final annotation freezing.
- `EXECUTION_ORDER.tsv`: analysis workflow order.
- `environment/`: software-environment information.
- `config/paths.example.env`: portable path configuration.

## Input data

Raw sequencing reads are not redistributed with this code package.
They should be obtained from the original GEO/SRA repositories.

Cell Ranger output locations can be configured with:

    export ATLAS_CELLRANGER_OUTPUT_ROOT=/path/to/cellranger_count/output

## Processed expression data

The public expression matrices contain final-QC,
SoupX-corrected integer RNA counts. They are processed counts,
not raw sequencing counts. No normalized-expression matrix is
distributed as the primary public MEX payload.

## Final Atlas

The frozen Atlas contains 665,816 cells from 158 included
libraries, represented in a common 38,606-feature RNA space
across nine source studies.

Libraries with fewer than 200 cells remaining after final QC were
excluded from Atlas integration.

## Reproducibility

See `EXECUTION_ORDER.tsv`, `environment/`, `MANIFEST.tsv`,
`SYNTAX_AUDIT.tsv`, `PUBLIC_PATH_AUDIT.tsv`, and
`SHA256SUMS_RELEASE.txt`.

## Downstream reuse resources

A companion reuse package provides compartment-wise batch-corrected reference
matrices for annotation/mapping, rich metadata, and a study-balanced
deconvolution reference. See `docs/P05E_REUSE_DATA_USE_GUIDE.txt` and
`docs/RELEASE_ARCHITECTURE.md`.

Corrected expression is not intended for differential-expression inference;
native final-QC SoupX-corrected integer counts remain authoritative for DE
and pseudobulk analyses.

## License

The analysis code in this repository is released under the MIT License.

The harmonized and derived Sepsis Atlas data release is distributed separately
under the Creative Commons CC0 1.0 Universal waiver. This does not alter or
supersede the rights, terms, or provenance of the original source datasets.
