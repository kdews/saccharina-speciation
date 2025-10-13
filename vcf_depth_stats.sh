#!/bin/bash
#SBATCH -J vcf_depth_stats
#SBATCH --mem=100gb
#SBATCH --time=05:00:00
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
vcf_base="${vcf_basename%%\.vcf.*}"
vcf_depth_stats="${vcf_base}.ldepth.mean"
vcf_site_stats="${vcf_base}.site_level.stats"
vcf_sample_stats="${vcf_base}.per_sample.stats"
plot_outdir="stats_plots"

# Clear modules
module purge
# Load conda
cond=~/.conda_for_sbatch.sh
if [[ -a "$cond" ]]
then
  source "$cond"
else
  echo "Error on source of $cond"
  exit 1
fi

# Check for existing result from VCFtools
if [[ -f "$vcf_depth_stats" ]]
then
  echo "Found result file: $vcf_depth_stats"
  echo "Skipping VCFtools step."
else
  # Activate VCFtools env
  module purge
  conda activate vcftools
  vcftools --version
  # Generate table of mean read depth per site (VCFtools)
  cmd=(
      vcftools
      --gzvcf "$in_vcf"
      --site-mean-depth
      --out "$vcf_base"
  )
  echo "${cmd[*]}"
  "${cmd[@]}"
  conda deactivate
fi
echo

# Check for existing result from BCFtools stats
if [[ -f "$vcf_site_stats" ]] && [[ -f "$vcf_sample_stats" ]]
then
  echo "Found result files: $vcf_site_stats & $vcf_sample_stats"
  echo "Skipping BCFtools stats step."
else
  # Activate BCFtools env
  module purge
  conda activate bcftools
  bcftools --version
  # Generate site-level BCFtools stats
  cmd2="bcftools stats $in_vcf"
  echo "$cmd2 > $vcf_site_stats"
  $cmd2 > "$vcf_site_stats"
  # Generate per-sample BCFtools stats
  cmd3="bcftools stats -s - $in_vcf"
  echo "$cmd3 > $vcf_sample_stats"
  $cmd3 > "$vcf_sample_stats"
  conda deactivate
fi
echo

# Check for existing stats plots
if [[ -f "$plot_outdir/summary.pdf" ]]
then
  echo "Found plots: $plot_outdir/summary.pdf"
  echo "Skipping plot-vcfstats step."
else
  # Activate BCFtools env
  module purge
  conda activate bcftools
  bcftools --version
  which python
  # Plot BCFtools stats
  cmd4="plot-vcfstats -p $plot_outdir $vcf_sample_stats"
  echo "$cmd4"
  $cmd4
  conda deactivate
fi
echo

# Done
echo
echo "Finished at:"
date