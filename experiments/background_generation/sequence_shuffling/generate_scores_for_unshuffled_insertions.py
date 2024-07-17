from optparse import OptionParser
import json
import os
import pickle
import random
import re
import pandas as pd
import pysam
import numpy as np
import tensorflow as tf
from basenji import seqnn, stream

from akita_utils.seq_gens import unshuffled_insertion_gen
from akita_utils.h5_utils import (initialize_stat_output_h5, write_stat_metrics_to_h5)
from akita_utils.tsv_utils import split_df_equally

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
    parser.add_option(
        "--save-maps",
        dest="save_maps",
        default=False,
        action="store_true",
        help="Save all the maps in the h5 file(for all inserts, all backgrounds used, and all targets)",
    )

    (options, args) = parser.parse_args()

    if len(args) == 3:
        # single worker
        params_file = args[0]
        model_file = args[1]
        shuffled_seqs_tsv = args[2]

    elif len(args) == 5:  # muliti-GPU option
        # multi worker
        options_pkl_file = args[0]
        params_file = args[1]
        model_file = args[2]
        shuffled_seqs_tsv = args[3]
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

    options.shifts = [int(shift) for shift in options.shifts.split(",")]
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
    # load shuffled sequences

    # filter for worker motifs
    if options.processes is not None:  # multi-GPU option
        # determine boundaries from motif file
        seq_coords_full = pd.read_csv(shuffled_seqs_tsv, sep="\t")
        seq_coords_df = split_df_equally(
            seq_coords_full, options.processes, worker_index
        )

    else:
        # read motif positions from csv
        seq_coords_df = pd.read_csv(shuffled_seqs_tsv, sep="\t")

    num_experiments = len(seq_coords_df)*2

    print("===================================")
    print(
        "Number of experiements = ", num_experiments
    )  # Warning! It's not number of predictions. Num of predictions is this number x5 or x6

    # open genome FASTA
    genome_open = pysam.Fastafile(
        options.genome_fasta
    )  # needs to be closed at some point

    #################################################################
    # strong CTCF site coordinates
    ctcf_site_coordinates = ("chr4", 87368268, 87368249, "+")

    # initialize output
    stats_out = initialize_stat_output_h5(options.out_dir, model_file, stats, seq_coords_df)

    print("stat_h5_outfile initialized")

    # if options.save_maps:
        # initlize map h5 files

    preds_stream = stream.PredStreamGen(
        seqnn_model,
        unshuffled_insertion_gen(seq_coords_df, genome_open, ctcf_site_coordinates),
        batch_size,
    )

    for exp_index in range(0, num_experiments, 2):
       
        reference_prediction_matrix = preds_stream[exp_index]
        prediction_matrix = preds_stream[exp_index+1]
            
        write_stat_metrics_to_h5(prediction_matrix,
                                    reference_prediction_matrix,
                                    stats_out,
                                    int(exp_index/2),
                                    head_index,
                                    model_index,
                                     diagonal_offset=2,
                                    stat_metrics=stats)  

        # if options.save_maps:
            # write maps

    stats_out.close()

    # if options.save_maps:
    #     maps_h5_outfile.close()

    genome_open.close()


################################################################################
# __main__
################################################################################


if __name__ == "__main__":
    main()
