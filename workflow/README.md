# Sepsis Atlas workflow utilities

## preprocessing

- atlas_preintegration_template.R
- check_droplet_qc_preqc_outputs_all_projects.R
- check_final_qc_all_projects.R
- run_final_qc_all_libraries.sh
- run_final_qc_one_library.R

## environment

- create_seurat_r43_conda_env.sh
- install_final_qc_packages_r43.R
- seurat_r43_conda_environment.yml

Scripts resolve the Sepsis.Atlas repository root from their installed
workflow location where necessary. Public frozen release files are
independent of this development-side organization.
