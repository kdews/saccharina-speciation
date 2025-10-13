#!/bin/bash
#SBATCH -J depth_filter_vcf
#SBATCH --mem=100gb
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=32
#SBATCH -o %x_%j.log

# Print date and time
date
echo

# Input
in_vcf="$1"
subset_tab="$2"
if [[ -f "$in_vcf" ]] && [[ -f "$subset_tab" ]]
then
  echo "Input VCF: $in_vcf"
  echo "Table of positions to subset: $subset_tab"
else
  echo "Error: Input VCF file ($in_vcf) or table ($subset_tab) not found."
  exit 1
fi
echo

# Output
# Name from input
vcf_basename="$(basename "$in_vcf")"
vcf_base="${vcf_basename%%.*}"
vcf_ext="${vcf_basename#*.}"
out_vcf="${vcf_base}.depth_filt.$vcf_ext"

# Load conda
cond=~/.conda_for_sbatch.sh
if [[ -a "$cond" ]]
then
  source "$cond"
else
  echo "Error on source of $cond"
  exit 1
fi

# Multithread support
if [[ -n "$SLURM_CPUS_PER_TASK" ]]
then
    nthrd="$SLURM_CPUS_PER_TASK"
else
    nthrd="1"
fi

# Activate BCFtools env
module purge
conda activate bcftools
bcftools --version

# Filter VCF for QUAL>30 & biallelic SNPs only
cmd=(
    bcftools view
    --threads "$nthrd"
    -R "$subset_tab"
    -Oz -o "$out_vcf"
    "$in_vcf"
)
echo "${cmd[*]}"
"${cmd[@]}"

# Print date and time
echo "Finished at:"
date