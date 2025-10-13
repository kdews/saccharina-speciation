#!/bin/bash
#SBATCH -J filter_for_PCA
#SBATCH --mem=100gb
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=32
#SBATCH -o %x_%j.log

# Print date and time
date
echo

# Input
in_vcf="$1"
in_vcf="$(realpath -e "$in_vcf")"
if [[ -f "$in_vcf" ]]
then
  echo "Input VCF: $in_vcf"
else
  echo "Error: Input VCF file ($in_vcf) not found."
  exit 1
fi
echo

# Output
# Name from input
vcf_basename="$(basename "$in_vcf")"
vcf_base="${vcf_basename%%.*}"
vcf_ext="${vcf_basename#*.}"
out_vcf="${vcf_base}.filt_QUAL30_bial.$vcf_ext"

# Save current directory path
CURDIR="$(pwd)"
echo "Current directory: $CURDIR"
# Create temporary working directory in /tmp filesystem
# (Avoids mmap errors with BLAST in /scratch1 and /project)
WORKDIR="/tmp/$USER"
mkdir -p "$WORKDIR"
cd "$WORKDIR" || exit 1
echo "Temporary working directory: $WORKDIR"
echo

# # Copy input VCF to /tmp
# echo "Copying input VCF to $WORKDIR"
# cmd="rsync -htu --progress ${in_vcf} ."
# echo "$cmd"
# $cmd

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

# Activate BCFtools module
# module purge
# module load gcc bcftools
# Activate BCFtools env
conda activate bcftools
bcftools --version

# Filter VCF for QUAL>30 & biallelic SNPs only
cmd=(
    bcftools view
    --threads "$nthrd"
    -v snps -m2 -M2
    -i 'QUAL>30'
    -Oz -o "$out_vcf"
    "$in_vcf"
)
echo "${cmd[*]}"
"${cmd[@]}"

# Copy output from /tmp to current directory, preserving timestamps (-t)
echo "Copying output in $WORKDIR/ to $CURDIR/"
rsync -htu --progress "$out_vcf" "$CURDIR/"
# Clean up /tmp
if [[ -d "$WORKDIR" ]]
then
    echo "Cleaning up $WORKDIR"
    rm -rfv "$WORKDIR"
fi

# Print date and time
echo
date