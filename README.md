# Sepsis Atlas

Publication-oriented analysis code for construction of the Sepsis Atlas from
publicly available single-cell RNA-sequencing datasets.

## Repository structure

- `code/R/`: reusable preprocessing and quality-control functions.
- `code/loaders/`: study-specific loading procedures.
- `code/canonical_metadata_freeze_v2/`: metadata harmonization.
- `code/gene_feature_mapping_v1/`: common feature-space construction.
- `code/full_atlas_primary_integration_v1/`: primary Atlas integration.
- `code/full_atlas_annotation_v1/`: annotation, validation, reconciliation,
  global QC, and final annotation freezing.
- `workflow/`: portable preprocessing and environment utilities together with
  the canonical downstream reuse workflow.
- `docs/`: documentation for the downstream reuse resources and release
  architecture.
- `environment/`: software-environment information.
- `config/paths.example.env`: example portable path configuration.
- `EXECUTION_ORDER.tsv`: manuscript-associated analysis workflow order.
- `DATA_AVAILABILITY.md`: summary of source and processed-data availability.
- `CODE_SNAPSHOT.txt`: identifier of the frozen source code package from which
  this repository was prepared.

## Input data

Raw sequencing reads are not redistributed with this code repository.
They should be obtained from the original GEO/SRA repositories.

Cell Ranger output locations can be configured with:

    export ATLAS_CELLRANGER_OUTPUT_ROOT=/path/to/cellranger_count/output

## Processed expression data

The public expression matrices contain final-QC, SoupX-corrected integer RNA
counts. They are processed counts, not raw sequencing counts. No
normalized-expression matrix is distributed as the primary public MEX payload.

## Final Atlas

The frozen Atlas contains 665,816 cells from 158 included libraries,
represented in a common 38,606-feature RNA space across nine source studies.

Libraries with fewer than 200 cells remaining after final QC were excluded
from Atlas integration.

## Reproducibility

The manuscript-associated public code snapshot is tagged
`code-v2-20260824`.

The analysis sequence is documented in `EXECUTION_ORDER.tsv`. Paths listed in
that file are relative to `code/`. Software-environment information is
provided in `environment/`, and portable path configuration is illustrated in
`config/paths.example.env`.

The frozen source package used to prepare this repository is identified in
`CODE_SNAPSHOT.txt`.

## Downstream reuse resources

The repository includes compartment-wise batch-corrected reference workflows
for annotation and mapping, rich metadata construction, and preparation of a
study-balanced deconvolution reference. See
`docs/P05E_REUSE_DATA_USE_GUIDE.txt` and `docs/RELEASE_ARCHITECTURE.md`.

Corrected expression is not intended for differential-expression inference;
native final-QC SoupX-corrected integer counts remain authoritative for
differential-expression and pseudobulk analyses.

## License

The analysis code in this repository is released under the MIT License.

The harmonized and derived Sepsis Atlas data release is distributed separately
under the Creative Commons CC0 1.0 Universal waiver. This does not alter or
supersede the rights, terms, or provenance of the original source datasets.
