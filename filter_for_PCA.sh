#!/bin/bash
#SBATCH -J filter_for_PCA
#SBATCH --mem=100gb
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=20
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
# R script to run
rscript="mean_depth_filt.R"
# Scripts directory
scripts_dir="saccharina-speciation/"
rscript="${scripts_dir}${rscript}"

# Output
# Name from input
vcf_basename="$(basename "$in_vcf")"
vcf_base="${vcf_basename%%.*}"
vcf_ext="${vcf_basename#*.}"
filt1="QUAL30.BIAL.MISS95"
filt2="MEAN_DP"
filt1_vcf="${vcf_base}.${filt1}.$vcf_ext"
filt2_tab="${vcf_base}.${filt2}.${filt1}.tsv"
filt2_vcf="${vcf_base}.${filt2}.${filt1}.$vcf_ext"

# # Save current directory path
# CURDIR="$(pwd)"
# echo "Current directory: $CURDIR"
# # Create temporary working directory in /tmp filesystem
# # (Avoids mmap errors with BLAST in /scratch1 and /project)
# WORKDIR="/tmp/$USER"
# mkdir -p "$WORKDIR"
# cd "$WORKDIR" || exit 1
# echo "Temporary working directory: $WORKDIR"
# echo

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

echo "Step 1: Biallelic SNP selection and quality + max-missing filtering"

# Activate BCFtools env
conda activate bcftools
bcftools --version

# Filter VCF:
# biallelic SNPs only
# QUAL > 30
# F_MISSING (fraction missing) < 0.05 (95% of samples have call)
cmd=(
    bcftools view
    --threads "$nthrd"
    -v snps -m2 -M2
    -i 'QUAL > 30 && F_MISSING < 0.05'
    -Oz -o "$filt1_vcf"
    "$in_vcf"
)
echo "${cmd[*]}"
"${cmd[@]}"
conda deactivate

# # Copy output from /tmp to current directory, preserving timestamps (-t)
# echo "Copying output in $WORKDIR/ to $CURDIR/"
# rsync -htu --progress "$filt1_vcf" "$CURDIR/"
# # Clean up /tmp
# if [[ -d "$WORKDIR" ]]
# then
#     echo "Cleaning up $WORKDIR"
#     rm -rfv "$WORKDIR"
# fi

echo "Finished step 1."
echo
echo "Step 2: Calculate per-site depth for VCF after first round of filtering"

cmd=(
  bash
  "${scripts_dir}vcf_depth_stats.sh"
  "$filt1_vcf"
)
echo "${cmd[*]}"
"${cmd[@]}"

echo "Finished step 2."
echo
echo "Step 3: Calculating per-site mean depth cutoffs"

# Load R and bgzip modules
module purge
module load apptainer/1.4.1 rstats/4.5.1
Rscript --version

# Run Rscript file
cmd=(
    "Rscript"
    "$rscript"
    "$filt1_vcf"
    "$filt2"
    "$scripts_dir"
)
echo "${cmd[*]}"
"${cmd[@]}"

echo "Finished step 3."
echo
echo "Step 4: Filter VCF with mean depth cutoffs"

# Activate BCFtools env
module purge
conda activate bcftools
bcftools --version

# Filter VCF for only SNPs that pass depth filter set in $filt2_tab
cmd=(
    bcftools view
    --threads "$nthrd"
    -R "$filt2_tab"
    -Oz -o "$filt2_vcf"
    "$in_vcf"
)
echo "${cmd[*]}"
"${cmd[@]}"
conda deactivate

# Print date and time
echo "Finished all steps at: $(date)"