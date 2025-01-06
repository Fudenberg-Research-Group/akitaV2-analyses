from akita_utils.format_io import h5_to_df
import pandas as pd
import numpy as np
import pandas as pd
import cooler
from cooltools.lib.numutils import observed_over_expected, adaptive_coarsegrain, set_diag
from bioframe.io.fileops import read_bigwig
from akita_utils.dna_utils import dna_1hot
from akita_utils.stats_utils import slide_diagonal_insulation, _extract_centered_window

def read_tad_disruption_data():
    df_1 = h5_to_df("/project/fudenber_735/akitaX1_analyses_data/TAD_boundary_disruption/model0_1.h5", average=False)
    df_2 = h5_to_df("/project/fudenber_735/akitaX1_analyses_data/TAD_boundary_disruption/model0_2.h5", average=False)
    df_3 = h5_to_df("/project/fudenber_735/akitaX1_analyses_data/TAD_boundary_disruption/model0_3.h5", average=False)
    df_4 = h5_to_df("/project/fudenber_735/akitaX1_analyses_data/TAD_boundary_disruption/model0_4.h5", average=False)
    df_strong = h5_to_df("/project/fudenber_735/akitaX1_analyses_data/TAD_boundary_disruption/model0_strong.h5", average=False)
    
    df = pd.concat([df_1, df_2, df_3, df_4, df_strong], axis=0, ignore_index=True)
    return df


# with removed data extrapolation
def get_target_no_extrapolation(cooler_path, padding, mseq_str, diagonal_offset=2):
    
    genome_hic_cool = cooler.Cooler(cooler_path)
    
    seq_hic_raw = genome_hic_cool.matrix(balance=True).fetch(mseq_str)

    seq_hic_nan = np.isnan(seq_hic_raw)

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
    
    return log_hic_obsexp


def simple_seqs_gen(
    seq_coords_df,
    genome_open
):
    for index, row in seq_coords_df.iterrows():

        wt_seq_1hot = dna_1hot(
            genome_open.fetch(row.chrom	, row.window_start, row.window_end).upper()
        )

        yield wt_seq_1hot


def calculate_min_insulation(target_map, window=16, crop_around_center=True, crop_width=2):
    map_size = target_map.shape[0]
    insulation_scores = slide_diagonal_insulation(target_map, window)
    
    if crop_around_center:
        insulation_scores = _extract_centered_window(
            insulation_scores, window=window, width=3
        )
    
    min_score = np.nanmin(insulation_scores)
    
    return min_score