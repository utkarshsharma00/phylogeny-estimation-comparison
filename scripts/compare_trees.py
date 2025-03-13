#!/usr/bin/env python3
"""
Python script for comparing trees
Usage: python compare_trees.py <estimated_tree> <true_tree> <output_file>
"""

import sys
import os
import dendropy
from dendropy.calculate import treecompare

def main():
    if len(sys.argv) != 4:
        print("Usage: python compare_trees.py <estimated_tree> <true_tree> <output_file>")
        sys.exit(1)
    
    estimated_tree_file = sys.argv[1]
    true_tree_file = sys.argv[2]
    output_file = sys.argv[3]
    
    try:
        # Read trees using DendroPy with more flexible format recognition
        taxa = dendropy.TaxonNamespace()
        
        # Try to read true tree (handle both Newick and Nexus formats)
        try:
            true_tree = dendropy.Tree.get(
                path=true_tree_file,
                schema="newick",
                taxon_namespace=taxa,
                preserve_underscores=True,
                rooting="force-rooted"
            )
        except Exception as e:
            print(f"Error reading true tree as Newick, trying Nexus: {e}")
            true_tree = dendropy.Tree.get(
                path=true_tree_file,
                schema="nexus",
                taxon_namespace=taxa,
                preserve_underscores=True,
                rooting="force-rooted"
            )
        
        # Try to read estimated tree
        try:
            estimated_tree = dendropy.Tree.get(
                path=estimated_tree_file,
                schema="newick",
                taxon_namespace=taxa,
                preserve_underscores=True,
                rooting="force-rooted"
            )
        except Exception as e:
            print(f"Error reading estimated tree as Newick, trying Nexus: {e}")
            estimated_tree = dendropy.Tree.get(
                path=estimated_tree_file,
                schema="nexus",
                taxon_namespace=taxa,
                preserve_underscores=True,
                rooting="force-rooted"
            )
        
        # Print some diagnostic information
        print(f"True tree has {len(true_tree.taxon_namespace)} taxa")
        print(f"Estimated tree has {len(estimated_tree.taxon_namespace)} taxa")
        
        # Ensure both trees are using the same taxon namespace
        common_taxa = set(t.label for t in true_tree.taxon_namespace).intersection(
            set(t.label for t in estimated_tree.taxon_namespace))
        
        print(f"Trees share {len(common_taxa)} taxa in common")
        
        if len(common_taxa) < 3:
            print("WARNING: Less than 3 common taxa found. Tree comparison may not be valid.")
            with open(output_file, 'w') as f:
                f.write("Tree comparison results:\n")
                f.write("WARNING: Less than 3 common taxa found. Comparison not valid.\n")
                f.write("Trees have different taxon sets.\n")
                f.write(f"True tree taxa: {len(true_tree.taxon_namespace)}\n")
                f.write(f"Estimated tree taxa: {len(estimated_tree.taxon_namespace)}\n")
                f.write(f"Common taxa: {len(common_taxa)}\n")
                f.write("RF distance: N/A\n")
                f.write("FN rate: N/A\n")
                f.write("FP rate: N/A\n")
            return
        
        # Prune trees to common taxa if needed
        if len(common_taxa) < len(true_tree.taxon_namespace) or len(common_taxa) < len(estimated_tree.taxon_namespace):
            print("Pruning trees to common taxa set")
            common_taxa_labels = list(common_taxa)
            
            # Create a filter for taxa to keep
            def filter_fn(taxon):
                return taxon.label in common_taxa_labels
            
            # Prune trees
            true_tree.retain_taxa_with_labels(common_taxa_labels)
            estimated_tree.retain_taxa_with_labels(common_taxa_labels)
        
        # Encode bipartitions for comparison
        true_tree.encode_bipartitions()
        estimated_tree.encode_bipartitions()
        
        # Calculate Robinson-Foulds distance
        rf_dist = treecompare.symmetric_difference(true_tree, estimated_tree)
        
        # Calculate FN and FP rates
        true_splits = set(b.split_bitmask for b in true_tree.bipartition_encoding)
        est_splits = set(b.split_bitmask for b in estimated_tree.bipartition_encoding)
        
        # Calculate FN rate (proportion of true splits missed)
        missing_splits = true_splits - est_splits
        FN = len(missing_splits) / (len(true_splits) if true_splits else 1)
        
        # Calculate FP rate (proportion of estimated splits not in true tree)
        extra_splits = est_splits - true_splits
        FP = len(extra_splits) / (len(est_splits) if est_splits else 1)
        
        # Write results to file
        with open(output_file, 'w') as f:
            f.write(f"Tree comparison results:\n")
            f.write(f"RF distance: {rf_dist}\n")
            f.write(f"FN rate: {FN}\n")
            f.write(f"FP rate: {FP}\n")
        
        print(f"Tree comparison completed: {output_file}")
        
    except Exception as e:
        print(f"Error in tree comparison: {e}")
        with open(output_file, 'w') as f:
            f.write(f"Tree comparison error: {str(e)}\n")

if __name__ == "__main__":
    main()