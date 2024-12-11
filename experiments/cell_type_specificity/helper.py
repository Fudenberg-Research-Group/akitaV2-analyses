import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import seaborn as sns

def plot_hexbin_accumulated(preds, targets, cell_type_idx, cell_type_name, save_path=None):
    # Concatenate all predictions and targets across the 725 test windows
    pred_vec_all = []
    target_vec_all = []
    
    for i in range(preds.shape[0]):
        # Flatten the 512x512 maps and concatenate
        pred_vec_all.append(preds[i, :, :, cell_type_idx].flatten())
        target_vec_all.append(targets[i, :, :, cell_type_idx].flatten())
    
    # Convert the list of arrays into a single array
    pred_vec_all = np.concatenate(pred_vec_all)
    target_vec_all = np.concatenate(target_vec_all)

    # Filter out NaN values
    mask = np.isfinite(pred_vec_all) & np.isfinite(target_vec_all)
    pred_vec_filtered = pred_vec_all[mask]
    target_vec_filtered = target_vec_all[mask]

    # Calculate Pearson correlation only on finite values
    corr, _ = pearsonr(pred_vec_filtered, target_vec_filtered) if len(pred_vec_filtered) > 0 else (np.nan, np.nan)

    # Create hexbin plot
    plt.figure(figsize=(7, 6))
    hb = plt.hexbin(pred_vec_filtered, target_vec_filtered, gridsize=50, cmap='viridis', bins='log')
    plt.colorbar(label='log10(frequency)')
    
    plt.xlabel(f'pred-{cell_type_name}')
    plt.ylabel(f'target-{cell_type_name}')
    plt.title(f'corr: {corr:.2f}')
    plt.grid(True)
    
    if save_path:
        plt.savefig(save_path, bbox_inches='tight', format="pdf")
    
    plt.show()

    
def calculate_cell_type_correlation(matrix):
    num_cell_types = matrix.shape[-1]
    correlation_matrix = np.zeros((num_cell_types, num_cell_types))
    
    # Iterate over all pairs of cell types
    for i in range(num_cell_types):
        for j in range(i, num_cell_types):
            # Flatten the 512x512 maps for both cell types
            cell_type_i = matrix[:, :, :, i].flatten()
            cell_type_j = matrix[:, :, :, j].flatten()
            
            # Create a mask for valid (non-NaN) entries
            mask = ~np.isnan(cell_type_i) & ~np.isnan(cell_type_j)
            
            if np.any(mask):  # Check if there are any valid pairs
                # Calculate correlation between the two cell types using valid entries
                corr, _ = pearsonr(cell_type_i[mask], cell_type_j[mask])
            else:
                corr = np.nan  # Set correlation to NaN if no valid pairs
            
            correlation_matrix[i, j] = corr
            correlation_matrix[j, i] = corr  # Symmetric matrix
    
    return correlation_matrix


# Function to calculate average correlations between cell types
def average_cell_type_correlations(matrix):
    num_cell_types = matrix.shape[-1]
    
    # Initialize matrix to accumulate correlations across all windows
    corr_sum = np.zeros((num_cell_types, num_cell_types))
    
    # Iterate over all windows
    for i in range(matrix.shape[0]):
        # Calculate correlation matrix for current window
        corr = calculate_cell_type_correlation(matrix[i:i+1])
        
        # Accumulate the correlations
        corr_sum += corr
    
    # Average correlations over all windows
    corr_avg = corr_sum / matrix.shape[0]
    
    return corr_avg
    
    
def calculate_cell_type_differences_and_correlations(preds, targets):
    num_cell_types = preds.shape[3]
    correlations = np.zeros((num_cell_types, num_cell_types))  # Store correlations between differences

    for i in range(num_cell_types):
        for j in range(num_cell_types):
            # Calculate differences between cell types i and j for preds and targets
            preds_diff = preds[:, :, :, i] - preds[:, :, :, j]  # Difference between predictions
            targets_diff = targets[:, :, :, i] - targets[:, :, :, j]  # Difference between targets

            # Flatten the differences for correlation calculation
            preds_diff_flat = preds_diff.reshape(preds_diff.shape[0], -1)
            targets_diff_flat = targets_diff.reshape(targets_diff.shape[0], -1)

            # Compute the correlation only for finite values
            mask = np.isfinite(preds_diff_flat) & np.isfinite(targets_diff_flat)
            if np.any(mask):  # Check if there are any valid values
                corr, _ = pearsonr(preds_diff_flat[mask], targets_diff_flat[mask])
                correlations[i, j] = corr
            else:
                correlations[i, j] = np.nan  # No valid values to compute correlation

    return correlations