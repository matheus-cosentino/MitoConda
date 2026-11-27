#!/bin/bash
######################################################################
#                            MitoConda.sh                            #
#                          MSc. Matheus Cosentino                    #
######################################################################
#   \  |  _)   |              ___|                       |           #
#  |\/ |   |   __|    _ \    |        _ \    __ \     _` |    _` |   #
#  |   |   |   |     (   |   |       (   |   |   |   (   |   (   |   #
# _|  _|  _|  \__|  \___/   \____|  \___/   _|  _|  \__,_|  \__,_|   #
#                                                                    #
######################################################################
#                          version: 06.2026                          #
######################################################################

# --- Color Palettes ---
blu="\033[1;34m"  
green="\033[1;32m" 
red="\033[1;31m"   
ylo="\033[1;33m"   
nc="\033[0m"       

# Environment setup
ENV_NAME="mitoconda"
ENV_FILE="workflow/envs/mitoconda.yaml"

#####################
# --- Functions --- #
#####################

# --- Spinner (Safe Version) --- #
run_with_spinner() {
    local pid
    ("$@" > /dev/null 2>&1) &
    pid=$!
    disown $pid 2>/dev/null
    local sequence="GGACTATCCTAAGTGTGACCGTGCTTTGCCGAGCATGATTAGGATGATTTCTGCCATGATACTTGGCTCTAAGCACACAACTTGCTGCACAAATAGTGATAGGTATTACAGATTGTGCAATGAGTTGGCACAAGTGCTCACTGAAGTTGTTTATTCCAATGGTGGTTTTTATTTTAAACCAGGAGGTACAACTTCAGGTGATGCAACTACAGCATATGCCAATTCTGTTTTCAACATATTCCAGGCTGTCAGTGCTAACATTAACCGTTTGCTCACTGTTGACAGTTATGCTATTCATAATGATTCTGTCAAGAGTTTGCAGAGGCAGTTGTATGACAATTGCTACCGTGCCACTTCTGTA"
    local seq_len=${#sequence}
    local width=15
    local pos=0
    
    local cA="\033[1;32m"
    local cT="\033[1;31m"
    local cC="\033[1;34m"
    local cG="\033[1;33m"
    local nc="\033[0m"

    while kill -0 $pid 2>/dev/null; do
        local chunk=""
        for (( j=0; j<width; j++ )); do
            local idx=$(( (pos + j) % seq_len ))
            local base="${sequence:$idx:1}"
            case "$base" in
                A) chunk="${chunk}${cA}A${nc}" ;;
                T) chunk="${chunk}${cT}T${nc}" ;;
                C) chunk="${chunk}${cC}C${nc}" ;;
                G) chunk="${chunk}${cG}G${nc}" ;;
                *) chunk="${chunk}${base}" ;;
            esac
        done
        printf "\r\033[K [ ${chunk} ] Processing..."
        pos=$(( (pos + 1) % seq_len ))
        sleep 0.1
    done
    wait $pid
    local exit_code=$?
    
    printf "\r\033[K"
    if [ $exit_code -eq 0 ]; then
        echo -e "${green}✔ Job done!${nc}"
    else
        echo -e "${red}✖ Job failed with exit code $exit_code.${nc}"
        exit $exit_code
    fi
}

# --- Help --- #
help(){
 echo -e "
 ${green}
 MitoConda${nc}: Amplicon Pipeline & Phylogeny

 ${green}Author${nc}: MSc. Matheus Cosentino 
 ${green}Version${nc}: 06.2026

 ${ylo}Usage:${nc}
  bash MitoConda.sh --output <DIR> [OPTIONS]

 ${ylo}Required Arguments:${nc}
   --output <DIR>       Directory where results will be saved
 
 ${ylo}Input Arguments:${nc}
   --input <DIR>        Directory containing raw reads (.fastq.gz) (Default: data/)

 ${ylo}Module Toggles (Enable Analysis):${nc}
   --diversity          Enable Taxonomic Diversity Analysis (BLAST + LCA Abundance)
   --phylo              Enable Phylogenetic Analysis (Per-sample + Integrated Tree)
   (Note: Quality Control, Reference DB build, and MultiQC always run automatically)
 
 ${ylo}Parameter Overrides:${nc}
   --taxid <INT>        Phylo target taxID (e.g., 9415 for phyllostomidae)
   --gene <STR>         Phylo gene target (e.g., 12s or MTCYB)
   --lca-rank <STR>     Taxonomic rank for LCA analysis (e.g., class, order, family)

 ${ylo}Optional Arguments:${nc}
   -h, --help           Show this help message
   -v, --version        Show version
   --jobs <INT>         Number of total submitted jobs at same time (default: 15)
   --temp-dir <DIR>     Temporary directory (default: system tmp)
 "
}

version(){
    echo "MitoConda v.06.2026"
}

###################################
# --- Environment Management --- #
###################################

manage_environment(){
    echo -e "${blu}[INFO]${nc} Checking Conda environment '${ylo}${ENV_NAME}${nc}'..."
    
    if __conda_setup="$('conda' 'shell.bash' 'hook' 2> /dev/null)"; then
        eval "$__conda_setup"
    else
        if [ -f "${HOME}/miniconda3/etc/profile.d/conda.sh" ]; then
            . "${HOME}/miniconda3/etc/profile.d/conda.sh"
        elif [ -f "${HOME}/anaconda3/etc/profile.d/conda.sh" ]; then
            . "${HOME}/anaconda3/etc/profile.d/conda.sh"
        elif [ -f "/usr/local/bioinfo/Miniforge3-24.9.2-0/etc/profile.d/conda.sh" ]; then
            . "/usr/local/bioinfo/Miniforge3-24.9.2-0/etc/profile.d/conda.sh"
        else
            echo -e "${red}[ERROR]${nc} Could not find conda.sh. Ensure Conda is installed."
            exit 1
        fi
    fi

    if [[ ! -f "$ENV_FILE" ]]; then
         echo -e "${red}[ERROR]${nc} Environment file $ENV_FILE not found!"
         exit 1
    fi

    if conda env list | grep -q "^${ENV_NAME} "; then
        echo -e "       > Environment found: ${green}Yes${nc}"
        echo -e "${blu}[INFO]${nc} Activating environment from ${ylo}$ENV_FILE${nc}..."
    else 
        echo -e "       > Environment found: ${red}No${nc}"
        echo -e "${blu}[INFO]${nc} Creating environment from ${ylo}$ENV_FILE${nc}..."
        run_with_spinner conda env create --name "$ENV_NAME" --file "$ENV_FILE" --quiet
    fi

    echo -ne "${blu}[INFO]${nc} Activating environment... "
    conda activate "$ENV_NAME"
    if [[ $? -eq 0 ]]; then 
        echo -e "${green}Active${nc}"
    else 
        echo -e "${red}Failed${nc}"
        exit 1
    fi
}

###########################
# --- Main Execution --- #
###########################
workdir=$(cd "$(dirname "${BASH_SOURCE[0]}" )" && pwd)
cd "$workdir"
set -e -o pipefail

# Defaults
input="data"
jobs=15
temp_dir=""

# Module Defaults (QC and Reference DB build always run)
mod_diversity="false"
mod_phylo="false"
mod_qc="true"
mod_build_db_only="true"

# Param Overrides
override_taxid=""
override_gene=""
override_lca=""

# --- Argument Parsing --- #
while [[ $# -gt 0 ]]; do
    case "$1" in
        --input) input="$2"; shift 2 ;;
        --output) output="$2"; shift 2 ;;
        --jobs) jobs="$2"; shift 2 ;;
        --temp-dir) temp_dir="$2"; shift 2 ;;

        # --- Module Toggles ---
        --diversity) mod_diversity="true"; shift ;;
        --phylo) mod_phylo="true"; shift ;;

        # --- Parameter Overrides ---
        --taxid) override_taxid="$2"; shift 2 ;;
        --gene) override_gene="$2"; shift 2 ;;
        --lca-rank) override_lca="$2"; shift 2 ;;

        # --- Others --- #
        -h|--help) help; exit 0 ;;
        -v|--version) version; exit 0 ;;
        *) echo -e "${red}[ERROR]${nc} Unknown argument: $1"; help; exit 1 ;;
    esac
done

if [[ -z "$output" ]]; then 
    echo -e "${red}[ERROR]${nc} Output directory is required. Use --output <DIR>"
    exit 1
fi

mkdir -p "$output"

manage_environment

echo -e "${blu}[INFO]${nc} Initializing MitoConda Workflow..."

# 1. Locate the STATIC Main Config
main_config="$workdir/config/config.yaml"
if [[ ! -f "$main_config" ]]; then
    if [[ -f "config.yaml" ]]; then
        main_config="config.yaml"
    else
        echo -e "${red}[ERROR]${nc} Could not find config.yaml in config/ or current dir."
        exit 1
    fi
fi

# 2. Setup Temp Dir
if [[ -z "$temp_dir" ]]; then
    temp_dir="/tmp/${USER}_cavicor"
fi
mkdir -p "$temp_dir"

# 3. Generate the DYNAMIC Override Config
run_overrides="$workdir/run_overrides.yaml"
echo -e "${blu}[INFO]${nc} Generating run overrides: ${ylo}$run_overrides${nc}"

cat <<EOF > "$run_overrides"
# Dynamic Overrides generated by MitoConda.sh
output_dir: "$output"
data_dir: "$input"

modules:
  12s_diversity: $mod_diversity
  quality_control: $mod_qc
  phylogeny: $mod_phylo
  build_db_phylo_only: $mod_build_db_only
EOF

if [[ -n "$override_taxid" ]]; then
    echo "phylo_target:" >> "$run_overrides"
    echo "  - $override_taxid" >> "$run_overrides"
fi

if [[ -n "$override_gene" ]]; then
    echo "phylo_genes:" >> "$run_overrides"
    echo "  - \"$override_gene\"" >> "$run_overrides"
fi

if [[ -n "$override_lca" ]]; then
    echo "taxonkit:" >> "$run_overrides"
    echo "  lca_rank:" >> "$run_overrides"
    echo "    - \"$override_lca\"" >> "$run_overrides"
fi

export TMPDIR="$temp_dir"
export TEMP="$temp_dir"
export TMP="$temp_dir"

# --- Workflow Execution Steps ---
echo -e "\n${green}> Snakemake: Unlocking working directory...${nc}"
snakemake --configfile "$main_config" "$run_overrides" --unlock --quiet || true

echo -e "\n${green}> Snakemake: Creating conda environments (if needed)...${nc}"
snakemake --configfile "$main_config" "$run_overrides" --use-conda --conda-create-envs-only --quiet

echo -e "\n${green}> Snakemake: Performing a dry-run...${nc}"
snakemake --jobs $jobs --use-conda --configfile "$main_config" "$run_overrides" --rerun-incomplete --dry-run
echo "---------------------------------------------------"

echo -e "\n${green}> Snakemake: Starting main execution...${nc}"

JOB_ID=${SLURM_JOB_ID:-local}
SHADOW_DIR="${temp_dir}/cavicor_shadow/${JOB_ID}"
mkdir -p "$SHADOW_DIR"
CONDA_DIR="$workdir/.snakemake/conda"

snakemake \
    --jobs $jobs \
    --use-conda \
    --conda-prefix "$CONDA_DIR" \
    --configfile "$main_config" "$run_overrides" \
    --shadow-prefix "$SHADOW_DIR" \
    --keep-going --rerun-incomplete

echo -e "\n${green}> Snakemake: Creating DAG & Report...${nc}"
mkdir -p "$output/all_samples/pipeline_logs"
snakemake --report "$output/all_samples/pipeline_logs/mitoconda_report.html" \
    --jobs $jobs \
    --use-conda \
    --configfile "$main_config" "$run_overrides" \
    --keep-going 

echo -e "\n${blu}[INFO]${nc} Deactivating ${ylo}${ENV_NAME}${nc} conda environment."
conda deactivate
echo -e "${green}✔ Done.${nc}"
echo -e "\n${green} Thank you for using MitoConda${nc}"
