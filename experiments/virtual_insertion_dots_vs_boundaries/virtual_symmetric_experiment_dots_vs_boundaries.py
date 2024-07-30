# Description:
# This script manages the execution of dot vs. boundary experiment using a deep learning model.
# It processes sequences from specified FASTA and TSV files, performs predictions, and saves the results.
# The script supports both single and multi-GPU configurations and provides options for saving statistical metrics and contact maps.
# It utilizes TensorFlow for model operations and handles background sequence predictions as well as motif-based predictions.
#
# Inputs:
# - params_file: Path to the JSON file containing model and training parameters.
# - model_file: Path to the file containing the trained model.
# - motifs_file: Path to the TSV file with motif positions.
#
# Options:
# - -f, --genome_fasta: Path to the genome FASTA file [Default: %default].
# - -l, --plot_lim_min: Minimum limit for heatmap plots [Default: %default].
# - --plot-freq: Frequency for heatmap plotting [Default: %default].
# - -m, --plot_map: Whether to plot contact maps for each allele [Default: %default].
# - -o, --out_dir: Output directory for tables and plots [Default: %default].
# - -p, --processes: Number of processes for multi-GPU setup.
# - --rc: Whether to average forward and reverse complement predictions [Default: %default].
# - --stats: Comma-separated list of statistics to save [Default: %default].
# - --shifts: Ensemble prediction shifts [Default: %default].
# - -t, --targets_file: File specifying target indexes and labels in table format.
# - --batch-size: Batch size for processing [Default: %default].
# - --background-file: File with insertion sequences in FASTA format [Default: %default].
# - --save-maps: Whether to save all maps in the H5 file [Default: %default].

# Outputs:
# - Creates an output directory with the following files:
#   - Stat H5 file containing statistical metrics for predictions.
#   - Optionally, maps H5 file with saved contact maps for reference and inserted sequences.
# - Generates log files for process execution.
#
# Example command-line usage:
# python virtual_symmetric_experiment_dots_vs_boundaries.py params.json model.h5 motifs.tsv -f /path/to/genome.fa -m -o results --rc --shifts "1,2" -t targets.tsv --batch-size 16 --background-file backgrounds.fa --save-maps -p 4

from __future__ import print_function
from optparse import OptionParser
import json
import os
import pickle
import random
import pandas as pd
import pysam

import tensorflow as tf

if tf.__version__[0] == "1":
    tf.compat.v1.enable_eager_execution()
gpus = tf.config.experimental.list_physical_devices("GPU")

print(gpus)

from basenji import seqnn, stream, dna_io

from akita_utils.seq_gens import symmertic_insertion_seqs_gen
from akita_utils.df_utils import split_df_equally
from akita_utils.h5_utils import (
    initialize_stat_output_h5,
    initialize_maps_output_h5,
    initialize_maps_output_references,
    write_stat_metrics_to_h5,
    write_maps_to_h5,
)


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
    # insertion-specific options
    parser.add_option(
        "--background-file",
        dest="background_file",
        default=None,
        help="file with insertion seqs in fasta format",
    )
    parser.add_option(
        "--save-maps",
        dest="save_maps",
        default=False,
        action="store_true",
        help="Save all the maps in the h5 file(for all inserts, all backgrounds used, and all targets)",
    )

    (options, args) = parser.parse_args()

    print("\n++++++++++++++++++\n")
    print("INPUT")
    print("\n++++++++++++++++++\n")
    print("options")
    print(options)
    print("\n++++++++++++++++++\n")
    print("args")
    print(args)
    print("\n++++++++++++++++++\n")

    if len(args) == 3:
        # single worker
        params_file = args[0]
        model_file = args[1]
        motif_file = args[2]

    elif len(args) == 5:  # muliti-GPU option
        # multi worker
        options_pkl_file = args[0]
        params_file = args[1]
        model_file = args[2]
        motif_file = args[3]
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

    if options.background_file is None:
        background_file = f"/project/fudenber_735/akitaX1_analyses_data/background_generation/background_generation/background_sequences_model_{model_index}.fa"
    else:
        background_file = options.background_file

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
        seq_coords_full = pd.read_csv(motif_file, sep="\t")
        seq_coords_df = split_df_equally(
            seq_coords_full, options.processes, worker_index
        )

    else:
        # read motif positions from csv
        seq_coords_df = pd.read_csv(motif_file, sep="\t")

    num_experiments = len(seq_coords_df)

    print("===================================")
    print(
        "Number of experiements = ", num_experiments
    )  # Warning! It's not number of predictions. Num of predictions is this number x5 or x6

    # open genome FASTA
    genome_open = pysam.Fastafile(
        options.genome_fasta
    )  # needs to be closed at some point

    background_seqs = []
    with open(background_file, "r") as f:
        for line in f.readlines():
            if ">" in line:
                continue
            background_seqs.append(dna_io.dna_1hot(line.strip()))

    num_backgrounds = len(background_seqs)
    max_bg_index = seq_coords_df["background_index"].max()

    if num_backgrounds < max_bg_index:
        raise ValueError(
            "must provide a background file with at least as many"
            + "backgrounds as those specified in the insert seq_coords tsv."
            + "\nThe provided background file has {len(background_seqs)} sequences."
        )

    #################################################################
    # setup output

    stat_h5_outfile = initialize_stat_output_h5(
        options.out_dir, model_file, stats, seq_coords_df
    )

    print("stat_h5_outfile initialized")

    if options.save_maps:
        maps_h5_outfile = initialize_maps_output_h5(
            options.out_dir, model_file, seqnn_model, seq_coords_df
        )
        print("maps_h5_outfile initialized")

        if options.processes is not None:
            if worker_index == 0:
                refmaps_h5_outfile = initialize_maps_output_references(
                    general_out_dir,
                    model_file,
                    seqnn_model,
                    num_backgrounds=num_backgrounds,
                )
                print("refmaps_h5_outfile initialized")

        else:
            refmaps_h5_outfile = initialize_maps_output_references(
                options.out_dir,
                model_file,
                seqnn_model,
                num_backgrounds=num_backgrounds,
            )
            print("refmaps_h5_outfile initialized")

    #################################################################

    # initialize predictions stream for reference (background) sequences
    refs_stream = seqnn_model.predict(np.array(background_seqs), batch_size=batch_size)

    for background_index in range(num_backgrounds):
        bg_prediction = refs_stream[background_index]

        if options.save_maps:
            # saving reference prediction vectors to the maps_h5_outfile
            if (len(args) == 5 and worker_index == 0) or len(args) == 3:
                # save maps for background sequences
                write_maps_to_h5(
                    bg_prediction,
                    refmaps_h5_outfile,
                    background_index,
                    head_index,
                    model_index,
                    reference=True,
                )
                print(f"Maps for reference {background_index} saved")

    # initialize predictions stream for alternate (ctcf-inserted) sequences
    preds_stream = stream.PredStreamGen(
        seqnn_model,
        symmertic_insertion_seqs_gen(
            seq_coords_df, background_seqs, genome_open
        ),
        batch_size,
    )

    # saving stat metrics
    for exp_index in range(num_experiments):
        # get predictions
        preds_matrix = preds_stream[exp_index]
        background_index = seq_coords_df.iloc[exp_index].background_index
        ref_matrix = refs_stream[background_index]

        write_stat_metrics_to_h5(
            preds_matrix,
            ref_matrix,
            stat_h5_outfile,
            exp_index,
            head_index,
            model_index,
            diagonal_offset=2,
            stat_metrics=stats,
        )

        if options.save_maps:
            # save maps for inserted sequences
            write_maps_to_h5(
                preds_matrix,
                maps_h5_outfile,
                exp_index,
                head_index,
                model_index,
                reference=False,
            )

    stat_h5_outfile.close()

    if options.save_maps:
        maps_h5_outfile.close()

    genome_open.close()


################################################################################
# __main__
################################################################################

if __name__ == "__main__":
    main()
