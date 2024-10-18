# Description:
# This script prepares and runs predictions using a trained deep learning model on genomic sequences,
# focusing on motif disruption analysis. It utilizes SLURM for job scheduling and management,
# supports both single and multi-GPU execution, and handles restarts of partially completed jobs.
# It saves various statistics and prediction maps, and allows for ensemble predictions.
#
# Inputs:
# - params_file: JSON file containing model parameters.
# - model_file: Model file to be used for predictions.
# - motifs_file: File containing motif positions (coordinates).
#
# Parameters:
# - genome_fasta: Genome FASTA file for sequences.
# - plot_lim_min: Minimum limit for heatmap plotting (default: 0.1).
# - plot_freq: Frequency of heatmap plotting (default: 100).
# - plot_map: Whether to plot contact maps for each allele (default: False).
# - out_dir: Output directory for tables and plots (default: "./").
# - processes: Number of processes for multi-GPU execution.
# - rc: Average forward and reverse complement predictions (default: False).
# - stats: Comma-separated list of statistics to save (default: "SCD").
# - shifts: Ensemble prediction shifts (default: "0").
# - targets_file: File specifying target indexes and labels in table format.
# - batch_size: Batch size for predictions (default: 4).
#
# Outputs:
# - Output directory containing results, logs, and pickled options.
#
# Example command-line usage:
# python genomic_disruption_sliding_profile.py params.json model_file.h5 motifs_file.csv

from optparse import OptionParser
import json
import os
import pickle
import random
import pandas as pd
import pysam
from basenji import seqnn, stream

from akita_utils.h5_utils import (
    initialize_stat_output_h5,
    write_stat_metrics_to_h5,
)
from akita_utils.df_utils import split_df_equally
from akita_utils.dna_utils import dna_1hot
from akita_utils.dna_utils import permute_seq_k

################################################################################
# main
################################################################################


def main():
    usage = "usage: %prog [options] <params_file> <model_file> <motifs_file>"
    parser = OptionParser(usage)
    parser.add_option(
        "-f",
        dest="genome_fasta",
        default=None,
        help="Genome FASTA for sequences [Default: %default]",
    )
    parser.add_option(
        "-l",
        dest="plot_lim_min",
        default=0.1,
        type="float",
        help="Heatmap plot limit [Default: %default]",
    )
    parser.add_option(
        "--plot-freq",
        dest="plot_freq",
        default=100,
        type="int",
        help="Heatmap plot freq [Default: %default]",
    )
    parser.add_option(
        "-m",
        dest="plot_map",
        default=False,
        action="store_true",
        help="Plot contact map for each allele [Default: %default]",
    )
    parser.add_option(
        "-o",
        dest="out_dir",
        default="./",
        help="Output directory for tables and plots [Default: %default]",
    )
    parser.add_option(
        "-p",
        dest="processes",
        default=None,
        type="int",
        help="Number of processes, passed by multi script",
    )
    parser.add_option(
        "--rc",
        dest="rc",
        default=False,
        action="store_true",
        help="Average forward and reverse complement predictions [Default: %default]",
    )
    parser.add_option(
        "--stats",
        dest="stats",
        default="SCD",
        type="str",
        help="Comma-separated list of stats to save. [Default: %default]",
    )
    parser.add_option(
        "--shifts",
        dest="shifts",
        default="0",
        type="str",
        help="Ensemble prediction shifts [Default: %default]",
    )
    parser.add_option(
        "-t",
        dest="targets_file",
        default=None,
        type="str",
        help="File specifying target indexes and labels in table format",
    )
    parser.add_option(
        "--batch-size",
        dest="batch_size",
        default=4,
        type="int",
        help="Specify batch size",
    )
    (options, args) = parser.parse_args()

    if len(args) == 3:
        # single worker
        params_file = args[0]
        model_file = args[1]
        windows_file = args[2]

    elif len(args) == 5:  # muliti-GPU option
        # multi worker
        options_pkl_file = args[0]
        params_file = args[1]
        model_file = args[2]
        windows_file = args[3]
        worker_index = int(args[4])

        # load options
        options_pkl = open(options_pkl_file, "rb")
        options = pickle.load(options_pkl)
        options_pkl.close()

        # update output directory
        general_out_dir = options.out_dir
        options.out_dir = "%s/job%d" % (options.out_dir, worker_index)

    else:
        parser.error(
            "Must provide parameters and model files and insertion TSV file"
        )

    if not os.path.isdir(options.out_dir):
        os.mkdir(options.out_dir)
    if options.plot_map:
        plot_dir = options.out_dir
    else:
        plot_dir = None

    stats = options.stats.split(",")

    head_index = int(model_file.split("model")[-1][0])
    model_index = int(model_file.split("c0")[0][-1])

    random.seed(44)

    #################################################################
    # read parameters and targets

    # read model parameters
    with open(params_file) as params_open:
        params = json.load(params_open)
    params_train = params["train"]
    params_model = params["model"]

    if options.batch_size is None:
        batch_size = params_train["batch_size"]
    else:
        batch_size = options.batch_size

    if options.targets_file is not None:
        targets_df = pd.read_csv(options.targets_file, sep="\t", index_col=0)
        target_ids = targets_df.identifier
        target_labels = targets_df.description

    #################################################################
    # load model
    seqnn_model = seqnn.SeqNN(params_model)
    seqnn_model.restore(model_file, head_i=head_index)
    seqnn_model.build_ensemble(options.rc, options.shifts)
    seq_length = int(params_model["seq_length"])

    # dummy target info
    if options.targets_file is None:
        num_targets = seqnn_model.num_targets()
        target_ids = [ti for ti in range(num_targets)]
        target_labels = [""] * len(target_ids)

    #################################################################
    # load motifs

    # filter for worker motifs
    if options.processes is not None:  # multi-GPU option
        # determine boundaries from motif file
        seq_coords_df_full = pd.read_csv(windows_file, sep="\t")
        seq_coords_df = split_df_equally(
            seq_coords_df_full, options.processes, worker_index
        )

    else:
        # read motif positions from csv
        seq_coords_df = pd.read_csv(windows_file, sep="\t")

    ##############################################################

    def tad_disruption_seq_gen(seq_coords_df, genome_open):
        for s in seq_coords_df.itertuples():
            list_1hot = []
            window_start = s.window_start
            window_end = s.window_end
            chrom = s.chrom
            rel_disruption_start = s.rel_disruption_start
            rel_disruption_end =s.rel_disruption_end
            
            # wild type
            wt_seq_1hot = dna_1hot(
                        genome_open.fetch(chrom, window_start, window_end).upper()
                    )
            list_1hot.append(wt_seq_1hot)
        
            alt_seq_1hot = wt_seq_1hot.copy()
            permuted_span = permute_seq_k(
                    alt_seq_1hot[rel_disruption_start:rel_disruption_end], k=8
                )
            alt_seq_1hot[rel_disruption_start:rel_disruption_end] = permuted_span
            list_1hot.append(alt_seq_1hot)
                
            # yielding first the reference, then the permuted sequence
            for sequence in list_1hot:
                yield sequence    

    ##############################################################

    num_experiments = len(seq_coords_df)

    print("===================================")
    print(
        "Number of experiements = ", num_experiments
    )  # Warning! It's not number of predictions. Num of predictions is this number x5 or x6

    # open genome FASTA
    genome_open = pysam.Fastafile(
        options.genome_fasta
    )  # needs to be closed at some point

    #################################################################
    # setup output

    # initialize output
    stats_out = initialize_stat_output_h5(
        options.out_dir, model_file, stats, seq_coords_df
    )

    print("stat_h5_outfile initialized")

    preds_stream = stream.PredStreamGen(
        seqnn_model,
        tad_disruption_seq_gen(seq_coords_df, genome_open),
        batch_size,
    )

    num_preds = len(seq_coords_df) * 2
    exp_index = 0

    for pred_index in range(num_preds):
        if pred_index % 2 == 0:
            # wild types
            ref_preds_matrix = preds_stream[pred_index]
        else:
            permuted_preds_matrix = preds_stream[pred_index]

            write_stat_metrics_to_h5(
                permuted_preds_matrix,
                ref_preds_matrix,
                stats_out,
                exp_index,
                head_index,
                model_index,
                diagonal_offset=2,
                stat_metrics=stats,
            )
            exp_index += 1

    stats_out.close()

    genome_open.close()


################################################################################
# __main__
################################################################################

if __name__ == "__main__":
    main()
