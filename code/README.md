# Analysis code

Analysis scripts used to construct the Sepsis Atlas.

## Components

- `loaders/`: source-specific data import
- `R/`: shared quality-control and metadata functions
- `canonical_metadata_freeze_v2/`: metadata harmonization
- `gene_feature_mapping_v1/`: common gene-feature mapping
- `full_atlas_primary_integration_v1/`: primary integration
- `full_atlas_annotation_v1/`: cell annotation and validation

The manuscript-associated execution sequence is listed in
`../EXECUTION_ORDER.tsv`.

File names are retained from the executed analysis to preserve provenance.

Code for construction of downstream reference resources is provided under
`../workflow/reuse/p05e_v1/`.
