# P05E downstream reuse workflow v1

This directory contains the canonical scripts used to construct the
downstream reuse resources distributed with the Sepsis Atlas.

## canonical_executed

The following scripts are the executed P05E0-P05E6/finalization workflow:

- `P05E0_reuse_reference_preflight_v1.R`
- `P05E1_define_reuse_reference_plan_v1.R`
- `P05E2_freeze_scmerge2_run_plan_v1.R`
- `P05E2b_refine_reuse_reference_plan_v1.R`
- `P05E3_select_scmerge2_ruvK_v1.R`
- `P05E3b_freeze_scmerge2_selection_v1.R`
- `P05E4_build_final_corrected_reference_v1_fix1.R`
- `P05E4a_build_compact_corrected_reference_v1_fix1.R`
- `P05E4a_finalize_compact_reference_v1_fix2.R`
- `P05E5_build_rich_metadata_v1_fix1.R`
- `P05E6_build_deconvolution_reference_v1.R`
- `P05E6_recover_deconvolution_reference_v1_fix1.R`
- `P05E_FINALIZE_reuse_package_v1_fix1_no_numpy.py`

Superseded development and maintenance scripts are not included in this
public code snapshot.

## Use of downstream resources

Batch-corrected expression resources are intended for within-compartment
reference mapping and annotation.

They should not be used as the quantitative expression space for
differential-expression analysis. Native final-QC SoupX-corrected integer
counts remain the quantitative expression source for differential-expression
and pseudobulk analyses.

The deconvolution reference is derived from the frozen Atlas and is intended
for downstream reuse with the associated Atlas cell-type definitions.
