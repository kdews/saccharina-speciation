#!/bin/bash
#SBATCH -J copy_vcf
#SBATCH --mem=100gb
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=32
#SBATCH -o %x_%j.log

# Print date and time
date
echo

# Load conda
cond=~/.conda_for_sbatch.sh
if [[ -a "$cond" ]]
then
  source "$cond"
else
  echo "Error on source of $cond"
  exit 1
fi

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
vcf_last_ext="${vcf_basename##*.}"

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

# Multithread support
if [[ -n "$SLURM_CPUS_PER_TASK" ]]
then
    nthrd="$SLURM_CPUS_PER_TASK"
else
    nthrd="1"
fi

# Activate bcftools conda env
conda activate bcftools
bcftools --version
echo

# Don't compress if already zipped
if [[ $vcf_last_ext == "gz" ]]
then
  out_vcf="$vcf_basename"
  echo "Skipping compression."
  echo "Copying input VCF to $WORKDIR"
  cmd="rsync -htu --progress ${in_vcf} ."
  echo "$cmd"
  $cmd
else
  out_vcf="$vcf_basename.gz"
  # Compress input VCF with bgzip
  cmd=(
        "bgzip"
        "-@"
        "$nthrd"
        "-o"
        "$out_vcf"
        "$in_vcf"
  )
  echo "${cmd[*]}"
  "${cmd[@]}"
  echo
fi

# Index bgzipped VCF with bcftools
# Need tabix format for snpEff/SnpSift
out_idx="$out_vcf.tbi"
# out_idx="$out_vcf.csi"
cmd=(
      "bcftools"
      "index"
      --tbi
      "-f"
      "--threads"
      "$nthrd"
      "-o"
      "$out_idx"
      "$out_vcf"
)
echo "${cmd[*]}"
"${cmd[@]}"
echo

# Copy output from /tmp to current directory, preserving timestamps (-t)
echo "Copying output in $WORKDIR/ to $CURDIR/"
rsync -htu --progress "$out_vcf" "$out_idx" "$CURDIR/"
# Clean up /tmp
if [[ -d "$WORKDIR" ]]; then
    echo "Cleaning up $WORKDIR"
    rm -rfv "$WORKDIR"
fi

# Done
echo
echo "Finished at:"
date