#!/bin/bash
#SBATCH -p gpu
#SBATCH -J mean_depth_filt
#SBATCH --mem=5gb
#SBATCH --cpus-per-task=32
#SBATCH --time=01:00:00
#SBATCH -o %x_%j.log

# Print date and time
date
echo

# Load R and bgzip modules
module purge
module load apptainer/1.4.1 rstats/4.5.1
Rscript --version
echo

# R script to run
rscript="mean_depth_filt.R"

# Input
# Scripts directory
scripts_dir="saccharina-speciation/"
rscript="${scripts_dir}${rscript}"
# in_vcf="master_SlaSLCT1FG3_1_AssemblyScaffolds_Repeatmasked.filt_QUAL30_bial.vcf.gz"
in_vcf="$1"
if [[ -f "$in_vcf" ]]
then
  echo "Input VCF: $in_vcf"
else
  echo "Error: Input VCF file ($in_vcf) not found."
  exit 1
fi
echo

# Run Rscript file
cmd=(
    "Rscript"
    "$rscript"
    "$in_vcf"
    "$scripts_dir"
)
echo "${cmd[*]}"
"${cmd[@]}"

# Print time at end
echo
date