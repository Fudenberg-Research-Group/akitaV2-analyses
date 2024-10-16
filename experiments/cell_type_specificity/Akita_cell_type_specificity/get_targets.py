import pandas as pd
import numpy as np
import cooler
from cooltools.lib.numutils import observed_over_expected, adaptive_coarsegrain
from cooltools.lib.numutils import interpolate_bad_singletons, set_diag, interp_nan
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve

test_path = "/project/fudenber_735/tensorflow_models/akita/v2/data/mm10/sequences.bed"
sequences_V2 = pd.read_csv(test_path, sep='\t', names=['chr','start','stop','type'])
fold0 = sequences_V2[sequences_V2["type"] == "fold0"]

def get_target(cooler_path, padding, mseq_str, diagonal_offset=2):
    
    genome_hic_cool = cooler.Cooler(cooler_path)
    
    seq_hic_raw = genome_hic_cool.matrix(balance=True).fetch(mseq_str)

    seq_hic_nan = np.isnan(seq_hic_raw)
    # num_filtered_bins = np.sum(np.sum(seq_hic_nan,axis=0) == len(seq_hic_nan))

    # clip first diagonals and high values
    clipval = np.nanmedian(np.diag(seq_hic_raw, diagonal_offset))
    for i in range(-diagonal_offset+1, diagonal_offset):
        set_diag(seq_hic_raw, clipval, i)
    seq_hic_raw = np.clip(seq_hic_raw, 0, clipval)
    seq_hic_raw[seq_hic_nan] = np.nan
    
    # adaptively coarsegrain based on raw counts
    seq_hic_smoothed = adaptive_coarsegrain(
                        seq_hic_raw,
                        genome_hic_cool.matrix(balance=False).fetch(mseq_str),
                        cutoff=2, max_levels=8)
    seq_hic_nan = np.isnan(seq_hic_smoothed)
    
    # local obs/exp
    seq_hic_obsexp = observed_over_expected(seq_hic_smoothed, ~seq_hic_nan)[0]
    log_hic_obsexp = np.log(seq_hic_obsexp)

    # crop
    if padding > 0:
        log_hic_obsexp = log_hic_obsexp[padding:-padding,:]
        log_hic_obsexp = log_hic_obsexp[:,padding:-padding]

    # aplying Gaussian Kernel
    kernel = Gaussian2DKernel(x_stddev=1)
    kernel_log_hic_obsexp = convolve(log_hic_obsexp, kernel)
    
    return kernel_log_hic_obsexp

def average_pooling(mat, pool_size=2):
    """
    Reduces the size of an 1153x1153 matrix to 513x513 using average pooling.
    Handles matrices that are not perfectly divisible by 2.
    """
    # Crop the matrix to 1152x1152 (since 1153 is not divisible by 2)
    mat_cropped = mat[:1152, :1152]

    # Apply average pooling manually by taking the average over 2x2 blocks
    pooled_mat = np.zeros((513, 513))
    for i in range(513):
        for j in range(513):
            pooled_mat[i, j] = np.mean(mat_cropped[i*2:i*2+2, j*2:j*2+2])
    
    return pooled_mat

padding = (640-512) // 2

cool_paths = {"Hsieh2019_mESC_uC_path" : "/project/fudenber_735/GEO/Hsieh2019/4DN/mESC_mm10_4DNFILZ1CPT8.mapq_30.2048.cool",
"Bonev2017_mESC_path" : "/project/fudenber_735/GEO/bonev_2017_GSE96107/distiller-0.3.1_mm10/results/coolers/HiC_ES_all.mm10.mapq_30.2048.cool",
"Bonev2017_CN_path" : "/project/fudenber_735/GEO/bonev_2017_GSE96107/distiller-0.3.1_mm10/results/coolers/HiC_CN_all.mm10.mapq_30.2048.cool",
"Bonev2017_ncx_CN_path" : "/project/fudenber_735/GEO/bonev_2017_GSE96107/distiller-0.3.1_mm10/results/coolers/HiC_ncx_CN_all.mm10.mapq_30.2048.cool",
"Bonev2017_NPC_path" : "/project/fudenber_735/GEO/bonev_2017_GSE96107/distiller-0.3.1_mm10/results/coolers/HiC_NPC_all_mm10.mapq_30.1024.cool",
"Bonev2017_ncx_NPC_path" : "/project/fudenber_735/GEO/bonev_2017_GSE96107/distiller-0.3.1_mm10/results/coolers/HiC_ncx_NPC_all.mm10.mapq_30.2048.cool"}

all_data = []

for index, row in fold0.iterrows():
    print("index: ", index)
    mseq_str = '%s:%d-%d' % (row.chr, row.start, row.stop)
    
    targets = []
    
    for path in cool_paths:
        targets.append(get_target(cool_paths[path], padding, mseq_str))
       
    targets[4] = average_pooling(targets[4])
        
    targets = np.array(targets)
    all_data.append(targets)
    
    if index % 100 == 0:
        np.save(f"target_matrices_{index}.npy", all_data)
        all_data = []
        
all_data = np.array(all_data)

np.save("target_matrices.npy", all_data)

