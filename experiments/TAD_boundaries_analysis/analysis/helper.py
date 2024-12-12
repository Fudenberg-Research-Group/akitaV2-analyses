from akita_utils.format_io import h5_to_df
import pandas as pd

def read_tad_disruption_data():
    df_1 = h5_to_df("/project/fudenber_735/akitaX1_analyses_data/TAD_boundary_disruption/model0_1.h5", average=False)
    df_2 = h5_to_df("/project/fudenber_735/akitaX1_analyses_data/TAD_boundary_disruption/model0_2.h5", average=False)
    df_3 = h5_to_df("/project/fudenber_735/akitaX1_analyses_data/TAD_boundary_disruption/model0_3.h5", average=False)
    df_4 = h5_to_df("/project/fudenber_735/akitaX1_analyses_data/TAD_boundary_disruption/model0_4.h5", average=False)
    df_strong = h5_to_df("/project/fudenber_735/akitaX1_analyses_data/TAD_boundary_disruption/model0_strong.h5", average=False)
    
    df = pd.concat([df_1, df_2, df_3, df_4, df_strong], axis=0, ignore_index=True)
    return df