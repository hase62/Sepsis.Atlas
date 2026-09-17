# Sepsis Atlas release architecture

The public resources are intentionally separated by analytical purpose.

## Native processed counts

The canonical public processed-count release contains final-QC
SoupX-corrected integer RNA counts with no public normalization. This is the
expression authority for differential expression and pseudobulk analyses.

## Compartment-wise corrected reference

A downstream reuse companion resource contains scMerge2-adjusted expression
for five resolved compartments. Correction was fitted independently within
compartments, so this resource is intended for annotation/reference mapping
rather than differential expression or quantitative cross-compartment
expression comparison.

## Rich metadata

Cell-level, corrected-reference-cell, and library-level metadata are supplied
to support study-aware reuse.

## Deconvolution

The default deconvolution resource is a study-balanced CORE reference derived
from native counts. The EXTENDED resource contains weaker-support resolved
identities and is exploratory.

## Sampling caveat

Source studies used heterogeneous sampling and enrichment strategies. Pooled
Atlas cell fractions must not be interpreted as population prevalence or
naive cross-study abundance estimates.
