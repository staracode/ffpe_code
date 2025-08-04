import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from pathlib import Path
import numpy as np

def create_heatmap(df, title, filename, figsize=(15, 12), cluster_samples=True):
    """
    Create a heatmap from a dataframe and save it
    """
    if cluster_samples and len(df.columns) > 1:
        # Use seaborn's clustermap for automatic dendrogram + heatmap
        g = sns.clustermap(df,
                          cmap="viridis",
                          cbar_kws={'label': 'Mutation Count'},
                          yticklabels=True,  # Show mutation type names
                          xticklabels=True,
                          figsize=figsize,
                          method='ward',
                          metric='euclidean',
                          row_cluster=False,  # Don't cluster rows (mutation types)
                          col_cluster=True,   # Cluster columns (samples)
                          dendrogram_ratio=0.2,  # Dendrogram takes 20% of height
                          cbar_pos=(0.02, 0.32, 0.03, 0.2))  # Position colorbar
         
        plt.suptitle(title)
        # Rotate x-axis labels
        g.ax_heatmap.tick_params(axis='x', rotation=45)
        
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"✓ Saved {filename}")
        
    else:
        # Create regular heatmap without clustering
        plt.figure(figsize=figsize)
        sns.heatmap(df,
                    cmap="viridis",
                    norm=plt.Normalize(vmin=df.min().min(), vmax=df.max().max()),
                    cbar_kws={'label': 'Mutation Count'},
                    yticklabels=True)  # Show mutation type names
        
        plt.title(title)
        plt.xlabel("Sample")
        plt.ylabel("Mutation Context")
        plt.tight_layout()
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"✓ Saved {filename}")

def main():
    # Create results directory if it doesn't exist
    os.makedirs("results", exist_ok=True)
    
    # Define sample selection
    selected = ["BRP_24_G_13", "BRP_6_H_11", "BRP_2_H_17", "BRP_6_G_7",
                "BRP_1_H_20", "BRP_1_G_3", "BRP_2_H_7", "BRP_6_G_18",
                "BRP_6_H_18", "BRP_1_H_31", "BRP_1_H_15", "BRP_1_G_23",
                "BRP_2_G_18", "BRP_6_H_8", "BRP_24_H_20", "BRP_6_G_14",
                "BRP_24_H_16", "BRP_24_H_7", "BRP_2_G_14", "BRP_24_G_17",
                "BRP_24_G_6", "BRP_2_H_13", "BRP_2_G_25", "BRP_1_G_19"]
    

    print("Reading data files...")
    
    # Original mutation matrix
    csvfile = "results/mutation_matrix.csv"
    df_orig = None
    if os.path.exists(csvfile):
        df_orig = pd.read_csv(csvfile, index_col=0)
        print(f"✓ Loaded {csvfile}")
    else:
        print(f"⚠️  {csvfile} not found")
    
    # Unrepaired corrected matrix
    unrepaired_file = "mutation_unrepaired/all_samples_corrected_matrix.csv"
    df_unrepaired = None
    if os.path.exists(unrepaired_file):
        df_unrepaired = pd.read_csv(unrepaired_file, index_col=0)
        print(f"✓ Loaded {unrepaired_file}")
    else:
        print(f"⚠️  {unrepaired_file} not found")
    
    # Repaired corrected matrix
    repaired_file = "mutation_repaired/all_samples_corrected_matrix.csv"
    df_repaired = None
    if os.path.exists(repaired_file):
        df_repaired = pd.read_csv(repaired_file, index_col=0)
        print(f"✓ Loaded {repaired_file}")
    else:
        print(f"⚠️  {repaired_file} not found")
    
    print("Data loading completed.\n")
    
    # Original mutation matrix heatmap
    print("Creating original mutation matrix heatmap...")
    if df_orig is not None:
        df_sel = df_orig[selected]
        create_heatmap(df_sel, 
                      "Original Mutation Context Heatmap", 
                      "results/heatmap_original.png")
    else:
        print(f"⚠️  Skipping original heatmap - data not loaded")
    
    # Unrepaired corrected matrix heatmap
    print("Creating unrepaired corrected matrix heatmap...")
    if df_unrepaired is not None:
        # Select only the samples that exist in the corrected matrix
        available_samples = [s for s in selected if s in df_unrepaired.columns]
        if available_samples:
            df_unrepaired_sel = df_unrepaired[available_samples]
            # Sort columns alphanumerically
            sorted_columns = sorted(df_unrepaired_sel.columns, key=str)
            df_unrepaired_sel = df_unrepaired_sel[sorted_columns]
            
            # Use original row labels if available
            if df_orig is not None and len(df_unrepaired_sel) == len(df_orig):
                df_unrepaired_sel.index = df_orig.index
            
            create_heatmap(df_unrepaired_sel, 
                          "Unrepaired Corrected Mutation Context Heatmap", 
                          "results/heatmap_unrepaired_corrected.png")
        else:
            print(f"⚠️  No selected samples found in unrepaired data")
    else:
        print(f"⚠️  Skipping unrepaired heatmap - data not loaded")
    
    # Repaired corrected matrix heatmap
    print("Creating repaired corrected matrix heatmap...")
    if df_repaired is not None:
        # Select only the samples that exist in the corrected matrix
        available_samples = [s for s in selected if s in df_repaired.columns]
        if available_samples:
            df_repaired_sel = df_repaired[available_samples]
            # Sort columns alphanumerically
            sorted_columns = sorted(df_repaired_sel.columns, key=str)
            df_repaired_sel = df_repaired_sel[sorted_columns]
            
            # Use original row labels if available
            if df_orig is not None and len(df_repaired_sel) == len(df_orig):
                df_repaired_sel.index = df_orig.index
            
            create_heatmap(df_repaired_sel, 
                          "Repaired Corrected Mutation Context Heatmap", 
                          "results/heatmap_repaired_corrected.png")
        else:
            print(f"⚠️  No selected samples found in repaired data")
    else:
        print(f"⚠️  Skipping repaired heatmap - data not loaded")
    
    # Comparison heatmap (Original vs Repaired Corrected)
    print("Creating comparison heatmap (Original vs Repaired Corrected)...")
    if df_orig is not None and df_repaired is not None:
        try:
            # Find common samples
            common_samples = list(set(df_orig.columns) & set(df_repaired.columns) & set(selected))
            if len(common_samples) >= 2:  # Need at least 2 samples for comparison
                # Create comparison dataframe with original and corrected values side by side
                # First, ensure both matrices have the same structure
                df_orig_subset = df_orig[common_samples].astype(float)
                df_repaired_subset = df_repaired[common_samples].astype(float)
                
                # Reset indices to ensure they match
                df_orig_subset = df_orig_subset.reset_index(drop=True)
                df_repaired_subset = df_repaired_subset.reset_index(drop=True)
                
                comparison_data = {}
                for sample in common_samples:
                    comparison_data[f"{sample}_Original"] = df_orig_subset[sample]
                    comparison_data[f"{sample}_Corrected"] = df_repaired_subset[sample]
                
                df_comparison = pd.DataFrame(comparison_data)
                
                # Use original row labels for the comparison matrix
                if len(df_comparison) == len(df_orig):
                    df_comparison.index = df_orig.index
                
                # Normalize the data (Z-score normalization)
                from sklearn.preprocessing import StandardScaler
                scaler = StandardScaler()
                df_comparison_normalized = pd.DataFrame(
                    scaler.fit_transform(df_comparison),
                    index=df_comparison.index,
                    columns=df_comparison.columns
                )
                
                # Create clustermap with normalized data clustering
                g = sns.clustermap(df_comparison_normalized,
                                  cmap="viridis",
                                  cbar_kws={'label': 'Normalized Mutation Count'},
                                  yticklabels=True,  # Show mutation type names
                                  xticklabels=True,
                                  figsize=(20, 14),
                                  method='ward',
                                  metric='euclidean',
                                  row_cluster=False,  # Don't cluster rows (mutation types)
                                  col_cluster=True,   # Cluster columns (samples)
                                  dendrogram_ratio=0.2,  # Dendrogram takes 20% of height
                                  cbar_pos=(0.02, 0.32, 0.03, 0.2))  # Position colorbar left
                
                # Set titles
                g.fig.suptitle('Original vs Repaired Corrected Mutation Contexts (Normalized)\n', 
                              fontsize=16, y=0.95)
                
                # Rotate x-axis labels
                g.ax_heatmap.tick_params(axis='x', rotation=45)
                
                plt.savefig("results/heatmap_comparison.png", dpi=300, bbox_inches='tight')
                plt.close()
                print("✓ Saved results/heatmap_comparison.png")
                
            else:
                print("⚠️  Not enough common samples for comparison heatmap")
        except Exception as e:
            print(f"⚠️  Error creating comparison heatmap: {e}")
    else:
        print("⚠️  One or both data matrices not loaded, skipping comparison heatmap")
    
    # Summary statistics plot
    print("Creating summary statistics plot...")
    try:
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        # Original vs Unrepaired comparison
        if df_orig is not None and df_unrepaired is not None:
            common_samples = list(set(df_orig.columns) & set(df_unrepaired.columns) & set(selected))
            
            if common_samples:
                orig_means = df_orig[common_samples].mean(axis=1)
                unrep_means = df_unrepaired[common_samples].mean(axis=1)
                
                axes[0, 0].scatter(orig_means, unrep_means, alpha=0.6)
                axes[0, 0].plot([orig_means.min(), orig_means.max()], 
                               [orig_means.min(), orig_means.max()], 'r--', alpha=0.8)
                axes[0, 0].set_xlabel('Original Mean Mutation Count')
                axes[0, 0].set_ylabel('Unrepaired Corrected Mean Mutation Count')
                axes[0, 0].set_title('Original vs Unrepaired Corrected')
                axes[0, 0].grid(True, alpha=0.3)
        
        # Original vs Repaired comparison
        if df_orig is not None and df_repaired is not None:
            common_samples = list(set(df_orig.columns) & set(df_repaired.columns) & set(selected))
            
            if common_samples:
                orig_means = df_orig[common_samples].mean(axis=1)
                rep_means = df_repaired[common_samples].mean(axis=1)
                
                axes[0, 1].scatter(orig_means, rep_means, alpha=0.6)
                axes[0, 1].plot([orig_means.min(), orig_means.max()], 
                               [orig_means.min(), orig_means.max()], 'r--', alpha=0.8)
                axes[0, 1].set_xlabel('Original Mean Mutation Count')
                axes[0, 1].set_ylabel('Repaired Corrected Mean Mutation Count')
                axes[0, 1].set_title('Original vs Repaired Corrected')
                axes[0, 1].grid(True, alpha=0.3)
        
        # Unrepaired vs Repaired comparison
        if os.path.exists(unrepaired_file) and os.path.exists(repaired_file):
            df_unrep = pd.read_csv(unrepaired_file, index_col=0)
            df_rep = pd.read_csv(repaired_file, index_col=0)
            common_samples = list(set(df_unrep.columns) & set(df_rep.columns) & set(selected))
            
            if common_samples:
                unrep_means = df_unrep[common_samples].mean(axis=1)
                rep_means = df_rep[common_samples].mean(axis=1)
                
                axes[1, 0].scatter(unrep_means, rep_means, alpha=0.6)
                axes[1, 0].plot([unrep_means.min(), unrep_means.max()], 
                               [unrep_means.min(), unrep_means.max()], 'r--', alpha=0.8)
                axes[1, 0].set_xlabel('Unrepaired Corrected Mean Mutation Count')
                axes[1, 0].set_ylabel('Repaired Corrected Mean Mutation Count')
                axes[1, 0].set_title('Unrepaired vs Repaired Corrected')
                axes[1, 0].grid(True, alpha=0.3)
        
        # Distribution comparison
        if os.path.exists(csvfile):
            df_orig = pd.read_csv(csvfile, index_col=0)
            orig_means = df_orig[selected].mean(axis=1)
            axes[1, 1].hist(orig_means, bins=20, alpha=0.7, label='Original', density=True)
            
            if df_unrepaired is not None:
                common_samples = list(set(df_unrepaired.columns) & set(selected))
                if common_samples:
                    unrep_means = df_unrepaired[common_samples].mean(axis=1)
                    axes[1, 1].hist(unrep_means, bins=20, alpha=0.7, label='Unrepaired Corrected', density=True)
            
            if df_repaired is not None:
                common_samples = list(set(df_repaired.columns) & set(selected))
                if common_samples:
                    rep_means = df_repaired[common_samples].mean(axis=1)
                    axes[1, 1].hist(rep_means, bins=20, alpha=0.7, label='Repaired Corrected', density=True)
            
            axes[1, 1].set_xlabel('Mean Mutation Count')
            axes[1, 1].set_ylabel('Density')
            axes[1, 1].set_title('Distribution Comparison')
            axes[1, 1].legend()
            axes[1, 1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig("results/heatmap_summary_statistics.png", dpi=300, bbox_inches='tight')
        plt.close()
        print("✓ Saved results/heatmap_summary_statistics.png")
        
    except Exception as e:
        print(f"⚠️  Error creating summary statistics plot: {e}")
    
    # Large mutation frequency heatmap (top 20 types)
    print("Creating large mutation frequency heatmap...")
    try:
        # Use the original mutation matrix for frequency calculation
        if df_orig is not None:
            
            # Calculate total mutations per sample
            total_mutations_per_sample = df_orig.sum(axis=0)
            
            # Calculate frequency of each mutation type per sample
            frequencies = df_orig.div(total_mutations_per_sample, axis=1)
            
            # Get top 20 mutation types by total frequency across all samples
            top_mutations = frequencies.sum(axis=1).nlargest(20).index
            heatmap_data = frequencies.loc[top_mutations]
            
            # Use seaborn's clustermap for frequency heatmap with dendrogram
            g = sns.clustermap(heatmap_data, 
                              cmap='YlOrRd', 
                              cbar_kws={'label': 'Mutation Frequency'},
                              yticklabels=True,  # Show mutation type names
                              xticklabels=True,
                              figsize=(24, 18),
                              method='ward',
                              metric='euclidean',
                              row_cluster=False,  # Don't cluster rows (mutation types)
                              col_cluster=True,   # Cluster columns (samples)
                              dendrogram_ratio=0.2,  # Dendrogram takes 20% of height
                              cbar_pos=(0.02, 0.32, 0.03, 0.2))  # Position colorbar
            
            # Set titles
            g.fig.suptitle('Mutation Frequency Heatmap - Top 20 Mutation Types', 
                          fontsize=16, y=0.95)
            
            # Rotate x-axis labels
            g.ax_heatmap.tick_params(axis='x', rotation=45)
            
            # Add frequency statistics as text
            stats_text = f"""
            Total Samples: {len(heatmap_data.columns)}
            Top 20 Mutation Types: {len(heatmap_data.index)}
            Frequency Range: {heatmap_data.values.min():.4f} - {heatmap_data.values.max():.4f}
            Mean Frequency: {heatmap_data.values.mean():.4f}
            """
            
            g.fig.text(0.02, 0.02, stats_text, fontsize=10, 
                      bbox=dict(boxstyle="round,pad=0.5", facecolor="lightgray", alpha=0.8))
            
            plt.savefig("results/mutation_frequency_heatmap_large.png", dpi=300, bbox_inches='tight')
            plt.close()
            print("✓ Saved results/mutation_frequency_heatmap_large.png")
            
            # Also create a version with sample annotations
            plt.figure(figsize=(28, 18))
            
            # Create annotation matrix (show frequencies as percentages)
            annot_matrix = (heatmap_data * 100).round(2).astype(str) + '%'
            
            sns.heatmap(heatmap_data, 
                       cmap='YlOrRd', 
                       cbar_kws={'label': 'Mutation Frequency', 'shrink': 0.8},
                       annot=annot_matrix,
                       fmt='',
                       linewidths=0.5,
                       linecolor='white',
                       square=False,
                       yticklabels=True,  # Show mutation type names
                       annot_kws={'size': 8})
            
            plt.title('Mutation Frequency Heatmap - Top 20 Mutation Types', 
                     fontsize=22, fontweight='bold', pad=20)
            plt.xlabel('Samples', fontsize=18, fontweight='bold')
            plt.ylabel('Mutation Types', fontsize=18, fontweight='bold')
            
            plt.xticks(rotation=45, ha='right', fontsize=12)
            plt.yticks(fontsize=12)
            
            plt.tight_layout()
            plt.savefig("results/mutation_frequency_heatmap_large_annotated.png", dpi=300, bbox_inches='tight')
            plt.close()
            print("✓ Saved results/mutation_frequency_heatmap_large_annotated.png")
            
        else:
            print(f"⚠️  {csvfile} not found, skipping mutation frequency heatmap")
            
    except Exception as e:
        print(f"⚠️  Error creating mutation frequency heatmap: {e}")
    
    # C>T mutation bar plot
    print("\nCreating C>T mutation bar plot...")
    try:
        if df_orig is not None:
            
            # Filter for C>T mutations (rows that start with "C>T")
            # Convert index to string type first
            df_orig.index = df_orig.index.astype(str)
            ct_mutations = df_orig[df_orig.index.str.startswith('C>T')]
            
            if len(ct_mutations) > 0:
                # Calculate sum of C>T mutations per sample
                ct_sums = ct_mutations.sum(axis=0)
                
                # Sort samples by C>T mutation count
                ct_sums_sorted = ct_sums.sort_values(ascending=False)
                
                # Create bar plot
                plt.figure(figsize=(16, 10))
                
                # Create bars with same color
                bars = plt.bar(range(len(ct_sums_sorted)), ct_sums_sorted.values, 
                              color='steelblue', alpha=0.7, edgecolor='black', linewidth=0.5)
                
                # Add value labels on top of bars
                for i, (bar, value) in enumerate(zip(bars, ct_sums_sorted.values)):
                    plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1, 
                            f'{int(value)}', ha='center', va='bottom', fontsize=10, fontweight='bold')
                
                plt.title('C>T Mutation Counts per Sample)', 
                         fontsize=18, fontweight='bold', pad=20)
                plt.xlabel('Samples', fontsize=14, fontweight='bold')
                plt.ylabel('Total C>T Mutations', fontsize=14, fontweight='bold')
                
                # Set x-axis labels
                plt.xticks(range(len(ct_sums_sorted)), ct_sums_sorted.index, 
                          rotation=45, ha='right', fontsize=10)
                plt.yticks(fontsize=12)
                
                # Add grid
                plt.grid(True, alpha=0.3, axis='y')
                
                # Add statistics text
                stats_text = f"""
                Total C>T Mutation Types: {len(ct_mutations)}
                Total C>T Mutations: {ct_sums_sorted.sum():,}
                Mean per Sample: {ct_sums_sorted.mean():.1f}
                Median per Sample: {ct_sums_sorted.median():.1f}
                Min: {ct_sums_sorted.min()}
                Max: {ct_sums_sorted.max()}
                """
                
                plt.figtext(0.98, 0.98, stats_text, fontsize=12, 
                           bbox=dict(boxstyle="round,pad=0.5", facecolor="lightgray", alpha=0.8),
                           ha='right', va='top')
                
                plt.tight_layout()
                plt.savefig("results/ct_mutations_bar_plot.png", dpi=300, bbox_inches='tight')
                plt.close()
                print("✓ Saved results/ct_mutations_bar_plot.png")
                
                # Save C>T mutation data
                ct_data = pd.DataFrame({
                    'Sample': ct_sums_sorted.index,
                    'C>T_Mutations': ct_sums_sorted.values,
                    'Rank': range(1, len(ct_sums_sorted) + 1)
                })
                ct_data.to_csv("results/ct_mutations_per_sample.csv", index=False)
                print("✓ Saved results/ct_mutations_per_sample.csv")
                
                        
            else:
                print("⚠️  No C>T mutations found in the data")
                
        else:
            print(f"⚠️  {csvfile} not found, skipping C>T mutation analysis")
            
    except Exception as e:
        print(f"⚠️  Error creating C>T mutation bar plot: {e}")
    
    print("\n🎉 Heatmap generation completed!")
    print("\n📊 Generated files:")
    print("  • Original heatmap: results/heatmap_original.png")
    print("  • Unrepaired corrected heatmap: results/heatmap_unrepaired_corrected.png")
    print("  • Repaired corrected heatmap: results/heatmap_repaired_corrected.png")
    print("  • Comparison heatmap (Original vs Corrected): results/heatmap_comparison.png")
    print("  • Summary statistics: results/heatmap_summary_statistics.png")
    print("  • Large frequency heatmaps: results/mutation_frequency_heatmap_large*.png")
    print("  • C>T mutation bar plot: results/ct_mutations_bar_plot.png")

if __name__ == "__main__":
    main()
