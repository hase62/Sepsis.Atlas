#!/usr/bin/env python3

from __future__ import annotations

import csv
import gzip
import hashlib
from pathlib import Path
import shutil
import subprocess
import sys
import time
import ast
import struct

ROOT = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
SD = ROOT / "publication" / "scientific_data"

P05B = SD / "processed_data_release_candidate_v1__20260819_194718"
P06 = SD / "final_public_release_freeze_v1__20260819_231631"

EXPECTED_COMPARTMENTS = [
    "T_NK",
    "Monocyte_DC",
    "B_plasma",
    "Neutrophil",
    "Platelet_megakaryocyte",
]

EXPECTED_NATIVE_CELLS = 665816
EXPECTED_LIBRARIES = 158
EXPECTED_FEATURES = 38606
EXPECTED_CORRECTED_CELLS = 597927
EXPECTED_COMMON_GENES = 2971
EXPECTED_CORE = 20
EXPECTED_EXTENDED = 22


def latest_dir(pattern: str) -> Path:
    hits = sorted(
        p for p in SD.iterdir()
        if p.is_dir() and p.match(pattern)
    )
    if not hits:
        raise RuntimeError(f"No directory matching {pattern}")
    return hits[-1]


def require_pass(summary: Path, label: str) -> dict[str, str]:
    if not summary.is_file():
        raise RuntimeError(f"{label}: missing {summary}")
    lines = summary.read_text(encoding="utf-8").splitlines()
    if "status=PASS" not in lines:
        raise RuntimeError(f"{label}: status != PASS")
    kv = {}
    for line in lines:
        if "=" in line:
            k, v = line.split("=", 1)
            kv[k.strip()] = v.strip()
    return kv


def gzip_test(path: Path) -> None:
    cp = subprocess.run(
        ["gzip", "-t", str(path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    if cp.returncode != 0:
        raise RuntimeError(
            f"gzip -t failed: {path}\n{cp.stdout}\n{cp.stderr}"
        )


def sha256_file(path: Path, block: int = 64 * 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while True:
            x = fh.read(block)
            if not x:
                break
            h.update(x)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        w.writeheader()
        w.writerows(rows)


def copy_verified(src: Path, dst: Path, test_gzip: bool = False) -> None:
    if not src.is_file():
        raise RuntimeError(f"Missing source file: {src}")
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    if src.stat().st_size != dst.stat().st_size:
        raise RuntimeError(f"Copy size mismatch: {src} -> {dst}")
    if test_gzip:
        gzip_test(dst)


def gzip_row_count(path: Path) -> int:
    n = 0
    with gzip.open(path, "rt", encoding="utf-8", newline="") as fh:
        for _ in fh:
            n += 1
    return max(0, n - 1)


def _npy_header_v1(shape: tuple[int, int]) -> bytes:
    """Return a NumPy .npy v1.0 header for little-endian float32 C-order."""
    header_dict = (
        "{'descr': '<f4', 'fortran_order': False, "
        f"'shape': ({shape[0]}, {shape[1]}), }}"
    )
    header = header_dict.encode("latin1")

    # v1.0 layout:
    # magic(6) + version(2) + header_len(2) + header
    # Entire preamble+header must be divisible by 16 and end with newline.
    preamble = 10
    pad = 16 - ((preamble + len(header) + 1) % 16)
    if pad == 16:
        pad = 0

    header = header + (b" " * pad) + b"\n"

    if len(header) >= 65536:
        raise RuntimeError("NPY v1 header unexpectedly too large")

    return (
        b"\x93NUMPY" +
        bytes([1, 0]) +
        struct.pack("<H", len(header)) +
        header
    )


def _read_npy_header(path: Path) -> tuple[dict, int]:
    with path.open("rb") as fh:
        magic = fh.read(6)
        if magic != b"\x93NUMPY":
            raise RuntimeError(f"{path.name}: invalid NPY magic")

        version = fh.read(2)
        if version != bytes([1, 0]):
            raise RuntimeError(f"{path.name}: unexpected NPY version")

        header_len = struct.unpack("<H", fh.read(2))[0]
        header = fh.read(header_len).decode("latin1").strip()
        meta = ast.literal_eval(header)
        data_offset = 10 + header_len

    return meta, data_offset


def convert_raw_f32_to_npy(
    raw_path: Path,
    npy_path: Path,
    n_cells: int,
    n_genes: int,
    block_bytes: int = 64 * 1024 * 1024,
) -> None:
    """
    Stream a raw little-endian float32 C-order matrix into standard .npy
    without requiring NumPy in the build environment.
    """
    raw_bytes = n_cells * n_genes * 4

    if raw_path.stat().st_size != raw_bytes:
        raise RuntimeError(
            f"{raw_path.name}: unexpected raw byte size "
            f"{raw_path.stat().st_size} != {raw_bytes}"
        )

    header = _npy_header_v1((n_cells, n_genes))

    with raw_path.open("rb") as src, npy_path.open("wb") as dst:
        dst.write(header)

        copied = 0
        next_report = 512 * 1024 * 1024

        while True:
            chunk = src.read(block_bytes)
            if not chunk:
                break

            dst.write(chunk)
            copied += len(chunk)

            if copied >= next_report:
                pct = copied / raw_bytes * 100
                print(
                    f"  {raw_path.stem}: copied {pct:.1f}%",
                    flush=True,
                )
                next_report += 512 * 1024 * 1024

    if copied != raw_bytes:
        raise RuntimeError(
            f"{npy_path.name}: copied bytes {copied} != {raw_bytes}"
        )

    meta, data_offset = _read_npy_header(npy_path)

    if meta.get("descr") != "<f4":
        raise RuntimeError(f"{npy_path.name}: dtype header check failed")
    if meta.get("fortran_order") is not False:
        raise RuntimeError(f"{npy_path.name}: order header check failed")
    if tuple(meta.get("shape", ())) != (n_cells, n_genes):
        raise RuntimeError(f"{npy_path.name}: shape header check failed")

    expected_npy_bytes = data_offset + raw_bytes
    if npy_path.stat().st_size != expected_npy_bytes:
        raise RuntimeError(
            f"{npy_path.name}: file-size check failed "
            f"{npy_path.stat().st_size} != {expected_npy_bytes}"
        )

    # Spot-check exact raw float32 bytes at 9 positions.
    test_rows = sorted(set([0, n_cells // 2, n_cells - 1]))
    test_cols = sorted(set([0, n_genes // 2, n_genes - 1]))

    with raw_path.open("rb") as src, npy_path.open("rb") as dst:
        for r in test_rows:
            for c in test_cols:
                raw_offset = (r * n_genes + c) * 4
                src.seek(raw_offset)
                a = src.read(4)

                dst.seek(data_offset + raw_offset)
                b = dst.read(4)

                if len(a) != 4 or len(b) != 4 or a != b:
                    raise RuntimeError(
                        f"{npy_path.name}: byte verification failed "
                        f"at ({r},{c})"
                    )


# ----------------------------------------------------------------------
# Resolve and validate frozen sources before creating output
# ----------------------------------------------------------------------

if not P05B.is_dir():
    raise RuntimeError(f"Missing P05B: {P05B}")
if not P06.is_dir():
    raise RuntimeError(f"Missing P06: {P06}")

p4a = latest_dir("reuse_reference_corrected_compact_v1__????????_??????")
p5 = latest_dir("reuse_reference_rich_metadata_v1__????????_??????")
p6 = latest_dir("reuse_reference_deconvolution_v1__????????_??????")

k4a = require_pass(p4a / "P05E4a_SUMMARY.txt", "P05E4a")
k5 = require_pass(p5 / "P05E5_SUMMARY.txt", "P05E5")
k6 = require_pass(p6 / "P05E6_SUMMARY.txt", "P05E6")

if int(k4a["corrected_reference_cells"]) != EXPECTED_CORRECTED_CELLS:
    raise RuntimeError("P05E4a corrected cell count mismatch")
if int(k4a["common_reference_genes"]) != EXPECTED_COMMON_GENES:
    raise RuntimeError("P05E4a common gene count mismatch")
if int(k5["full_atlas_cells"]) != EXPECTED_NATIVE_CELLS:
    raise RuntimeError("P05E5 native cell count mismatch")
if int(k5["libraries"]) != EXPECTED_LIBRARIES:
    raise RuntimeError("P05E5 library count mismatch")
if int(k6["reference_genes"]) != EXPECTED_FEATURES:
    raise RuntimeError("P05E6 gene count mismatch")
if int(k6["CORE_identities"]) != EXPECTED_CORE:
    raise RuntimeError("P05E6 CORE count mismatch")
if int(k6["EXTENDED_total_identities"]) != EXPECTED_EXTENDED:
    raise RuntimeError("P05E6 EXTENDED count mismatch")

corr_manifest_src = p4a / "compact_corrected_reference_manifest.tsv"
corr_genes_src = (
    p4a / "metadata" /
    "SepsisAtlas_corrected_reference_common_genes_v1.tsv"
)
corr_combined_cells_src = (
    p4a / "metadata" /
    "SepsisAtlas_corrected_reference_combined_cell_order_v1.tsv"
)

corr_rows = read_tsv(corr_manifest_src)
if len(corr_rows) != 5:
    raise RuntimeError("Corrected manifest must have 5 rows")

corr_by_comp = {r["compartment"]: r for r in corr_rows}
if set(corr_by_comp) != set(EXPECTED_COMPARTMENTS):
    raise RuntimeError("Corrected compartment set mismatch")

for comp in EXPECTED_COMPARTMENTS:
    row = corr_by_comp[comp]
    if int(row["n_genes"]) != EXPECTED_COMMON_GENES:
        raise RuntimeError(f"{comp}: n_genes != 2971")
    if not (p4a / row["matrix_file"]).is_file():
        raise RuntimeError(f"{comp}: missing matrix")
    if not (p4a / row["cell_order_file"]).is_file():
        raise RuntimeError(f"{comp}: missing cell order")

p5_required = [
    "SepsisAtlas_rich_cell_metadata_v1.tsv.gz",
    "SepsisAtlas_corrected_reference_rich_cell_metadata_v1.tsv.gz",
    "SepsisAtlas_rich_library_metadata_v1.tsv.gz",
    "SepsisAtlas_rich_metadata_dictionary_v1.tsv",
    "rich_metadata_column_coverage.tsv",
    "corrected_reference_status_counts.tsv",
]

p6_required = [
    "SepsisAtlas_deconvolution_reference_CORE_CPM_v1.tsv.gz",
    "SepsisAtlas_deconvolution_reference_CORE_log2CPM_v1.tsv.gz",
    "SepsisAtlas_deconvolution_reference_EXTENDED_CPM_v1.tsv.gz",
    "SepsisAtlas_deconvolution_reference_EXTENDED_log2CPM_v1.tsv.gz",
    "SepsisAtlas_deconvolution_reference_identities_v1.tsv",
    "SepsisAtlas_deconvolution_reference_support_v1.tsv",
    "SepsisAtlas_deconvolution_excluded_identities_v1.tsv",
    "deconvolution_identity_blood_fraction_support.tsv",
    "SepsisAtlas_deconvolution_sample_pseudobulk_counts_v1.mtx.gz",
    "SepsisAtlas_deconvolution_sample_pseudobulk_metadata_v1.tsv.gz",
    "SepsisAtlas_deconvolution_features_v1.tsv.gz",
    "SepsisAtlas_deconvolution_project_condition_balanced_CPM_v1.tsv.gz",
    "SepsisAtlas_deconvolution_project_profile_metadata_v1.tsv.gz",
    "README_SepsisAtlas_deconvolution_reference_v1.txt",
]

for name in p5_required:
    if not (p5 / name).is_file():
        raise RuntimeError(f"Missing P05E5 file: {name}")

for name in p6_required:
    if not (p6 / name).is_file():
        raise RuntimeError(f"Missing P05E6 file: {name}")

for name in p5_required:
    if name.endswith(".gz"):
        gzip_test(p5 / name)

for name in p6_required:
    if name.endswith(".gz"):
        gzip_test(p6 / name)

# ----------------------------------------------------------------------
# Create final downstream reuse package
# ----------------------------------------------------------------------

tag = time.strftime("%Y%m%d_%H%M%S")
out = SD / f"reuse_reference_release_v1__{tag}"

corr_out = out / "corrected_reference"
rich_out = out / "rich_metadata"
deconv_out = out / "deconvolution"
prov_out = out / "provenance"
docs_out = out / "docs"

for d in [corr_out, rich_out, deconv_out, prov_out, docs_out]:
    d.mkdir(parents=True, exist_ok=True)

public_corr_rows = []

for comp in EXPECTED_COMPARTMENTS:
    row = corr_by_comp[comp]
    n_cells = int(row["n_cells"])
    n_genes = int(row["n_genes"])

    raw_src = p4a / row["matrix_file"]
    cells_src = p4a / row["cell_order_file"]

    npy_name = f"{comp}__scMerge2_corrected_expression_common2971_v1.npy"
    npy_path = corr_out / npy_name

    print(f"CONVERT {comp}: {n_cells} x {n_genes}", flush=True)

    convert_raw_f32_to_npy(
        raw_src,
        npy_path,
        n_cells,
        n_genes,
    )

    cells_name = f"{comp}__corrected_reference_cells_v1.tsv"
    copy_verified(
        cells_src,
        corr_out / cells_name,
    )

    public_corr_rows.append({
        "compartment": comp,
        "n_cells": n_cells,
        "n_genes": n_genes,
        "dtype": "float32",
        "array_order": "C",
        "matrix_file": npy_name,
        "cell_order_file": cells_name,
        "correction_scope": "within_compartment",
        "intended_use": "reference_mapping_annotation",
        "DE_use": "NO",
    })

copy_verified(
    corr_genes_src,
    corr_out / "SepsisAtlas_corrected_reference_common_genes_v1.tsv",
)

copy_verified(
    corr_combined_cells_src,
    corr_out / "SepsisAtlas_corrected_reference_combined_cell_order_v1.tsv",
)

write_tsv(
    corr_out / "SepsisAtlas_corrected_reference_manifest_v1.tsv",
    public_corr_rows,
    [
        "compartment",
        "n_cells",
        "n_genes",
        "dtype",
        "array_order",
        "matrix_file",
        "cell_order_file",
        "correction_scope",
        "intended_use",
        "DE_use",
    ],
)

loader_text = """#!/usr/bin/env python3
import csv
import numpy as np
from pathlib import Path

root = Path(__file__).resolve().parent

with open(
    root / "SepsisAtlas_corrected_reference_manifest_v1.tsv",
    encoding="utf-8",
) as fh:
    manifest = list(csv.DictReader(fh, delimiter="\\t"))

for row in manifest:
    x = np.load(root / row["matrix_file"], mmap_mode="r")
    print(row["compartment"], x.shape, x.dtype)

# IMPORTANT:
# Each compartment was batch-corrected independently.
# Use for within-compartment annotation/reference mapping.
# Do not use corrected values for DE or quantitative cross-compartment
# expression comparisons.
"""
(corr_out / "load_corrected_reference_python.py").write_text(
    loader_text,
    encoding="utf-8",
)

corr_readme = """Sepsis Atlas compartment-wise corrected reference v1

The five matrices were fitted independently within:
T_NK, Monocyte_DC, B_plasma, Neutrophil, Platelet_megakaryocyte.

Format
------
Each .npy file is a cells x 2971 genes float32 NumPy array.
Gene order:
  SepsisAtlas_corrected_reference_common_genes_v1.tsv
Cell order:
  matching compartment-specific *_cells_v1.tsv

Intended use
------------
Reference mapping and annotation within the relevant compartment.

Do not use
----------
Do not use corrected matrices for differential-expression inference.
Do not treat corrected values across compartments as one globally fitted
batch-corrected Atlas.

Native final-QC SoupX-corrected integer counts remain authoritative for
DE and pseudobulk analyses.

A single Atlas-wide corrected h5ad is intentionally not created because
that representation would falsely imply one global correction model.
"""
(corr_out / "README_corrected_reference_v1.txt").write_text(
    corr_readme,
    encoding="utf-8",
)

# Rich metadata
for name in p5_required:
    copy_verified(
        p5 / name,
        rich_out / name,
        test_gzip=name.endswith(".gz"),
    )

metadata_expected_rows = {
    "SepsisAtlas_rich_cell_metadata_v1.tsv.gz": 665816,
    "SepsisAtlas_corrected_reference_rich_cell_metadata_v1.tsv.gz": 597927,
    "SepsisAtlas_rich_library_metadata_v1.tsv.gz": 158,
}

for name, expected in metadata_expected_rows.items():
    observed = gzip_row_count(rich_out / name)
    if observed != expected:
        raise RuntimeError(
            f"{name}: row count {observed} != {expected}"
        )

# Deconvolution
for name in p6_required:
    copy_verified(
        p6 / name,
        deconv_out / name,
        test_gzip=name.endswith(".gz"),
    )

# Provenance summaries
for src, name in [
    (p4a / "P05E4a_SUMMARY.txt", "P05E4a_SUMMARY.txt"),
    (p5 / "P05E5_SUMMARY.txt", "P05E5_SUMMARY.txt"),
    (p6 / "P05E6_SUMMARY.txt", "P05E6_SUMMARY.txt"),
]:
    copy_verified(src, prov_out / name)

native_authority = """Sepsis Atlas native processed-count authority

Canonical native-count release:
  processed_data_release_candidate_v1__20260819_194718

Final public release freeze:
  final_public_release_freeze_v1__20260819_231631

Semantics:
  final-QC SoupX-corrected integer RNA counts
  public normalization: none
  665816 cells
  158 libraries
  38606 features

The frozen P05B/P06 native release is not duplicated inside this downstream
reuse companion package. It remains the expression source of truth for
differential expression and pseudobulk analysis.
"""
(prov_out / "NATIVE_COUNTS_AUTHORITY.txt").write_text(
    native_authority,
    encoding="utf-8",
)

data_use = """Sepsis Atlas reuse data-use guide

1. Differential expression / pseudobulk
   Use frozen P05B/P06 native processed counts plus rich metadata.
   Expression authority is final-QC SoupX-corrected integer RNA counts.
   Prefer sample/patient-aware pseudobulk models with study/covariates.
   Never use scMerge2-adjusted expression for DE.

2. Query-cell annotation / reference mapping
   Use corrected_reference/<compartment>__*.npy.
   First assign a coarse compartment, then map within that compartment.
   Correction scope is within_compartment.

3. Rich metadata
   SepsisAtlas_rich_cell_metadata_v1.tsv.gz:
     all 665816 Atlas cells.
   SepsisAtlas_corrected_reference_rich_cell_metadata_v1.tsv.gz:
     597927 corrected-reference cells.
   SepsisAtlas_rich_library_metadata_v1.tsv.gz:
     158-library clinical/technical/source metadata.

4. Bulk-RNA deconvolution
   Default:
     SepsisAtlas_deconvolution_reference_CORE_CPM_v1.tsv.gz
   CORE contains 20 identities and is study-balanced.
   EXTENDED contains 22 identities total and is exploratory only.
   No separate whole-blood-specific reference is frozen because
   cross-study whole-blood support is insufficient across all identities.

5. Sampling / abundance
   Source studies have heterogeneous enrichment and sampling strategies.
   Do not interpret pooled Atlas cell fractions as population prevalence
   or direct cross-study abundance estimates without study-aware modeling.
"""
(docs_out / "DATA_USE_GUIDE.txt").write_text(
    data_use,
    encoding="utf-8",
)

readme = """Sepsis Atlas downstream reuse companion package v1

This is a downstream reuse derivative and does not replace or alter the
frozen P05B/P06 public processed-count release.

Four-pillar design
------------------
1. Native processed counts
   Frozen P05B/P06; not duplicated here.
   Use for DE/pseudobulk.

2. Compartment-wise corrected reference
   597927 cells, 2971 common genes, five scMerge2 models.
   Distributed as per-compartment float32 NumPy arrays.
   Use for annotation/reference mapping, not DE.

3. Rich metadata
   665816-cell rich metadata, 597927-cell corrected-reference metadata,
   and 158-library metadata.

4. Deconvolution reference
   Native-count-derived, study-balanced CORE 20 reference plus
   exploratory EXTENDED 22 reference.

Why there is no single corrected h5ad
-------------------------------------
The five corrected matrices were fitted independently within compartments.
A single X matrix would make the data look like one globally corrected
Atlas, which is not the scientific model that was fitted.
"""
(out / "README.txt").write_text(readme, encoding="utf-8")

# ----------------------------------------------------------------------
# File manifest and SHA256
# ----------------------------------------------------------------------

manifest_rows = []

for path in sorted(p for p in out.rglob("*") if p.is_file()):
    rel = path.relative_to(out).as_posix()

    if rel.startswith("corrected_reference/"):
        pillar = "corrected_reference"
    elif rel.startswith("rich_metadata/"):
        pillar = "rich_metadata"
    elif rel.startswith("deconvolution/"):
        pillar = "deconvolution"
    elif rel.startswith("provenance/"):
        pillar = "provenance"
    elif rel.startswith("docs/") or rel == "README.txt":
        pillar = "documentation"
    else:
        pillar = "other"

    manifest_rows.append({
        "release_relative_path": rel,
        "pillar": pillar,
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path),
        "public_scope": "YES",
    })

manifest_path = out / "P05E_REUSE_RELEASE_MANIFEST.tsv"
write_tsv(
    manifest_path,
    manifest_rows,
    [
        "release_relative_path",
        "pillar",
        "bytes",
        "sha256",
        "public_scope",
    ],
)

npy_files = list(corr_out.glob("*.npy"))
npy_bytes = sum(p.stat().st_size for p in npy_files)

summary_lines = [
    "===== P05E REUSE PACKAGE FINAL FREEZE =====",
    "status=PASS",
    "native_count_authority=P05B_P06",
    "native_cells=665816",
    "native_libraries=158",
    "native_features=38606",
    "native_counts_duplicated_in_reuse_package=NO",
    "native_expression_semantics=final_QC_SoupX_corrected_integer_RNA_counts",
    "native_DE_use=YES",
    "corrected_reference_cells=597927",
    "corrected_reference_common_genes=2971",
    "corrected_reference_compartments=5",
    "corrected_reference_format=NumPy_float32_per_compartment",
    "numpy_required_for_package_build=NO",
    "corrected_reference_scope=within_compartment",
    "corrected_reference_DE_use=NO",
    "single_global_corrected_h5ad=NOT_CREATED_BY_DESIGN",
    f"corrected_reference_npy_bytes={npy_bytes}",
    "rich_cell_metadata_rows=665816",
    "rich_corrected_metadata_rows=597927",
    "rich_library_metadata_rows=158",
    "rich_cell_metadata_columns=109",
    "rich_library_metadata_columns=486",
    "deconvolution_expression_source=native_counts",
    "deconvolution_CORE_identities=20",
    "deconvolution_EXTENDED_total_identities=22",
    "deconvolution_CORE_is_default=YES",
    "whole_blood_specific_reference=NOT_FROZEN_INSUFFICIENT_CROSS_STUDY_SUPPORT",
    "npy_header_shape_dtype_verification=PASS",
    "npy_payload_byte_spotcheck=PASS",
    "metadata_row_count_verification=PASS",
    "gzip_integrity=PASS",
    "expression_modified=NO",
    "annotation_modified=NO",
    "P05B_modified=NO",
    "P06_modified=NO",
    "cellranger_count=NOT_ACCESSED",
]

summary_path = out / "P05E_REUSE_RELEASE_SUMMARY.txt"
summary_path.write_text(
    "\n".join(summary_lines) + "\n",
    encoding="utf-8",
)

sha_path = out / "SHA256SUMS_P05E_REUSE_RELEASE.txt"

def write_and_verify_sha():
    lines = []
    for path in sorted(
        p for p in out.rglob("*")
        if p.is_file() and p.name != sha_path.name
    ):
        digest = sha256_file(path)
        rel = path.relative_to(out).as_posix()
        lines.append(f"{digest}  {rel}")

    sha_path.write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )

    for line in lines:
        expected, rel = line.split("  ", 1)
        if sha256_file(out / rel) != expected:
            raise RuntimeError(f"SHA verification failed: {rel}")

write_and_verify_sha()

summary_lines.append("SHA256=PASS")
summary_path.write_text(
    "\n".join(summary_lines) + "\n",
    encoding="utf-8",
)

write_and_verify_sha()

# Root handoff
handoffs = {
    summary_path: ROOT / "P05E_FINAL_SUMMARY.txt",
    manifest_path: ROOT / "P05E_FINAL_RELEASE_MANIFEST.tsv",
    docs_out / "DATA_USE_GUIDE.txt": ROOT / "P05E_FINAL_DATA_USE_GUIDE.txt",
}

for src, dst in handoffs.items():
    shutil.copy2(src, dst)

print("\n".join(summary_lines))
print(f"\nOUT_DIR={out}")
print("Root handoff files copied: 3")
