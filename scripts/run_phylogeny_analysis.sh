#!/bin/bash
#SBATCH --time=12:00:00                     # Job run time (hh:mm:ss) - 12 hours should be plenty
#SBATCH --nodes=1                           # Number of nodes
#SBATCH --ntasks-per-node=32                # Increase from 1 to 32 cores
#SBATCH --mem=32G                           # Increased memory to match core count
#SBATCH --job-name=phylogeny-analysis       # Name of job
#SBATCH --account=25sp-cs581a-eng           # Account for CS 581 students
#SBATCH --partition=eng-instruction         # Partition for CS 581
#SBATCH --constraint=AE7713                 # Target the 128-core nodes
#SBATCH --output=phylogeny_analysis_%j.log  # Output file with job ID
#SBATCH --mail-type=BEGIN,END,FAIL          # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=usharma4@illinois.edu   # Replace with your email

# Print job info
echo "Job ID: $SLURM_JOB_ID"
echo "Running on host: $(hostname)"
echo "Running on nodes: $SLURM_NODELIST"
echo "Number of CPUs: $SLURM_NTASKS"
echo "Current working directory: $(pwd)"

# Directory setup
DATASET_DIR="/u/usharma4/scratch/dataset"
RESULTS_DIR="/u/usharma4/scratch/results"

# Create results directory
mkdir -p $RESULTS_DIR

# Load modules and activate conda properly
module load anaconda3
eval "$(conda shell.bash hook)"
conda activate phylo_env

# Make sure we're in the home directory where our Python scripts are located
cd /u/usharma4

# Define the path to Python in the conda environment
PYTHON_PATH=$(which python)
echo "Using Python from: $PYTHON_PATH"

# Check if dendropy is installed
echo "Checking for dendropy:"
$PYTHON_PATH -c "import sys; print(sys.path); import dendropy; print(f'DendroPy version: {dendropy.__version__}')" || echo "Dendropy not found in Python path"

# Run on both model conditions
for MODEL_COND in "1000M1" "1000M4"; do
    echo "Processing model condition: $MODEL_COND"
    mkdir -p $RESULTS_DIR/$MODEL_COND
    
    # Process replicates (R0 through R4)
    for REP in {0..19}; do
        REP_DIR="$DATASET_DIR/$MODEL_COND/R$REP"
        
        # Skip if the replicate directory doesn't exist
        if [ ! -d "$REP_DIR" ]; then
            echo "  Replicate directory $REP_DIR not found, skipping"
            continue
        fi
        
        echo "  Processing replicate R$REP"
        
        # Define input files
        if [ -f "$REP_DIR/rose.aln.true.fasta" ]; then
            ALIGNMENT="$REP_DIR/rose.aln.true.fasta"
        else
            echo "WARNING: Cannot find alignment file in $REP_DIR"
            continue
        fi
        
        if [ -f "$REP_DIR/rose.tt" ]; then
            TRUE_TREE="$REP_DIR/rose.tt"
        elif [ -f "$REP_DIR/random.tree" ]; then
            TRUE_TREE="$REP_DIR/random.tree"
        else
            echo "WARNING: Cannot find true tree file in $REP_DIR"
            continue
        fi
        
        # Define output directories
        OUT_DIR="$RESULTS_DIR/$MODEL_COND/R$REP"
        mkdir -p $OUT_DIR
        
        # Convert FASTA to PHYLIP format for FastME
        PHYLIP_FILE="$OUT_DIR/alignment.phy"
        echo "  Converting FASTA to PHYLIP format for FastME"
        $PYTHON_PATH -c "
from Bio import AlignIO
try:
    AlignIO.convert('$ALIGNMENT', 'fasta', '$PHYLIP_FILE', 'phylip-relaxed')
    print('Conversion successful')
except Exception as e:
    print(f'Error converting file: {e}')
"
        
        # Run FastTree with GTR model (if not already done)
        if [ ! -f "$OUT_DIR/fasttree_gtr.tree" ]; then
            echo "    Running FastTree GTR"
            FastTree -gtr -nt -fastest "$ALIGNMENT" > "$OUT_DIR/fasttree_gtr.tree"
        else
            echo "    FastTree GTR already completed, skipping"
        fi
        
        # Run FastTree with JC model (if not already done)
        if [ ! -f "$OUT_DIR/fasttree_jc.tree" ]; then
            echo "    Running FastTree JC"
            FastTree -nt -fastest "$ALIGNMENT" > "$OUT_DIR/fasttree_jc.tree"
        else
            echo "    FastTree JC already completed, skipping"
        fi
        
        # Run FastME for NJ with JC distances using gap removal and triangular inequality correction
        if [ ! -f "$OUT_DIR/nj_jc.tree" ] || [ ! -s "$OUT_DIR/nj_jc.tree" ]; then
            echo "    Running FastME with JC distances (with gap removal and triangular inequality correction)"
            # Use -r to remove gaps and -q to correct triangular inequality violations
            fastme -i "$PHYLIP_FILE" -o "$OUT_DIR/nj_jc.tree" -m N -d J -T $SLURM_NTASKS
            
            # # If FastME is still failing, try BioNJ which is more robust with problematic distances
            # if [ ! -s "$OUT_DIR/nj_jc.tree" ]; then
            #     echo "    First FastME attempt failed, trying with BioNJ method"
            #     fastme -i "$PHYLIP_FILE" -o "$OUT_DIR/nj_jc.tree" -m I -d J -r -q -T $SLURM_NTASKS
            # fi
            
            # # If still failing, try RapidNJ as last resort
            # if [ ! -s "$OUT_DIR/nj_jc.tree" ] && which rapidnj > /dev/null; then
            #     echo "    FastME attempts failed. Using RapidNJ with JC distances"
            #     rapidnj "$ALIGNMENT" -i fa -d jc -o t > "$OUT_DIR/nj_jc.tree"
            # fi
        else
            echo "    NJ with JC distances already completed, skipping"
        fi
        
        # Run FastME for NJ with LogDet distances using gap removal and triangular inequality correction
        if [ ! -f "$OUT_DIR/nj_logdet.tree" ] || [ ! -s "$OUT_DIR/nj_logdet.tree" ]; then
            echo "    Running FastME with LogDet distances (with gap removal and triangular inequality correction)"
            # Use -r to remove gaps and -q to correct triangular inequality violations
            fastme -i "$PHYLIP_FILE" -o "$OUT_DIR/nj_logdet.tree" -m N -d L -T $SLURM_NTASKS
            
            # # If FastME is still failing, try BioNJ which is more robust with problematic distances
            # if [ ! -s "$OUT_DIR/nj_logdet.tree" ]; then
            #     echo "    First FastME attempt failed, trying with BioNJ method"
            #     fastme -i "$PHYLIP_FILE" -o "$OUT_DIR/nj_logdet.tree" -m I -d L -r -q -T $SLURM_NTASKS
            # fi
            
            # # If still failing, try RapidNJ as last resort
            # if [ ! -s "$OUT_DIR/nj_logdet.tree" ] && which rapidnj > /dev/null; then
            #     echo "    FastME attempts failed. Using RapidNJ with LogDet distances"
            #     rapidnj "$ALIGNMENT" -i fa -d ld -o t > "$OUT_DIR/nj_logdet.tree"
            # fi
        else
            echo "    NJ with LogDet distances already completed, skipping"
        fi
        
        # Run FastME for NJ with p-distances using gap removal and triangular inequality correction
        if [ ! -f "$OUT_DIR/nj_pdist.tree" ] || [ ! -s "$OUT_DIR/nj_pdist.tree" ]; then
            echo "    Running FastME with p-distances (with gap removal and triangular inequality correction)"
            # Use -r to remove gaps and -q to correct triangular inequality violations
            fastme -i "$PHYLIP_FILE" -o "$OUT_DIR/nj_pdist.tree" -m N -d p -T $SLURM_NTASKS
            
            # # If FastME is still failing, try BioNJ which is more robust with problematic distances
            # if [ ! -s "$OUT_DIR/nj_pdist.tree" ]; then
            #     echo "    First FastME attempt failed, trying with BioNJ method"
            #     fastme -i "$PHYLIP_FILE" -o "$OUT_DIR/nj_pdist.tree" -m I -d p -r -q -T $SLURM_NTASKS
            # fi
            
            # # If still failing, try RapidNJ as last resort
            # if [ ! -s "$OUT_DIR/nj_pdist.tree" ] && which rapidnj > /dev/null; then
            #     echo "    FastME attempts failed. Using RapidNJ with p-distances"
            #     rapidnj "$ALIGNMENT" -i fa -d pd -o t > "$OUT_DIR/nj_pdist.tree"
            # fi
        else
            echo "    NJ with p-distances already completed, skipping"
        fi
        
        # Fix tree files before comparison
        for tree_file in "$OUT_DIR/fasttree_gtr.tree" "$OUT_DIR/fasttree_jc.tree" "$OUT_DIR/nj_jc.tree" "$OUT_DIR/nj_logdet.tree" "$OUT_DIR/nj_pdist.tree"; do
            if [ -f "$tree_file" ] && [ -s "$tree_file" ]; then
                # Fix potential format issues with dendropy
                $PYTHON_PATH -c "
import dendropy
try:
    tree = dendropy.Tree.get(path='$tree_file', schema='newick', rooting='force-rooted')
    tree.write(path='$tree_file.fixed', schema='newick')
    print('Successfully fixed tree format for $tree_file')
except Exception as e:
    print(f'Error fixing tree format: {e}')
"
                # Replace original with fixed version if successful
                if [ -f "$tree_file.fixed" ]; then
                    mv "$tree_file.fixed" "$tree_file"
                fi
            fi
        done
        
        # Compare trees to true tree
        if [ ! -f "$OUT_DIR/fasttree_gtr_comparison.txt" ] && [ -f "$OUT_DIR/fasttree_gtr.tree" ] && [ -s "$OUT_DIR/fasttree_gtr.tree" ]; then
            echo "    Comparing FastTree GTR to true tree"
            $PYTHON_PATH compare_trees.py "$OUT_DIR/fasttree_gtr.tree" "$TRUE_TREE" "$OUT_DIR/fasttree_gtr_comparison.txt"
        else
            echo "    FastTree GTR comparison already completed or tree file missing, skipping"
        fi
        
        if [ ! -f "$OUT_DIR/fasttree_jc_comparison.txt" ] && [ -f "$OUT_DIR/fasttree_jc.tree" ] && [ -s "$OUT_DIR/fasttree_jc.tree" ]; then
            echo "    Comparing FastTree JC to true tree"
            $PYTHON_PATH compare_trees.py "$OUT_DIR/fasttree_jc.tree" "$TRUE_TREE" "$OUT_DIR/fasttree_jc_comparison.txt"
        else
            echo "    FastTree JC comparison already completed or tree file missing, skipping"
        fi
        
        if [ ! -f "$OUT_DIR/nj_jc_comparison.txt" ] && [ -f "$OUT_DIR/nj_jc.tree" ] && [ -s "$OUT_DIR/nj_jc.tree" ]; then
            echo "    Comparing NJ JC to true tree"
            $PYTHON_PATH compare_trees.py "$OUT_DIR/nj_jc.tree" "$TRUE_TREE" "$OUT_DIR/nj_jc_comparison.txt"
        else
            echo "    NJ JC comparison already completed or tree file missing, skipping"
        fi
        
        if [ ! -f "$OUT_DIR/nj_logdet_comparison.txt" ] && [ -f "$OUT_DIR/nj_logdet.tree" ] && [ -s "$OUT_DIR/nj_logdet.tree" ]; then
            echo "    Comparing NJ LogDet to true tree"
            $PYTHON_PATH compare_trees.py "$OUT_DIR/nj_logdet.tree" "$TRUE_TREE" "$OUT_DIR/nj_logdet_comparison.txt"
        else
            echo "    NJ LogDet comparison already completed or tree file missing, skipping"
        fi
        
        if [ ! -f "$OUT_DIR/nj_pdist_comparison.txt" ] && [ -f "$OUT_DIR/nj_pdist.tree" ] && [ -s "$OUT_DIR/nj_pdist.tree" ]; then
            echo "    Comparing NJ p-distance to true tree"
            $PYTHON_PATH compare_trees.py "$OUT_DIR/nj_pdist.tree" "$TRUE_TREE" "$OUT_DIR/nj_pdist_comparison.txt"
        else
            echo "    NJ p-distance comparison already completed or tree file missing, skipping"
        fi
    done
done

# Summarize results
echo "Creating summary results"
$PYTHON_PATH summarize_results.py "$RESULTS_DIR"

echo "Analysis completed. Results are in $RESULTS_DIR"