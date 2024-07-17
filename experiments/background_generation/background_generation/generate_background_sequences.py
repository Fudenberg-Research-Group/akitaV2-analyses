#!/usr/bin/python

# Copyright 2017 Calico LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# =========================================================================

###################################################

"""
generate_background_sequences.py
derived from virtual_symmetric_experiment.py

This scripts computes genomic disruption scores for motifs from a tsv file with:
chrom | start | end
where one row represents a single experiment.
The insertion scores are added as next keys in the h5 format file.

The script requires the following input:

Parameters:
-----------
<params_file> - parameters for the akita model
<model_file> - model in the h5 format
<windows_file> - tsv/csv table specifying genomic windows

Options:
-----------
- path to the mouse or human genome in the fasta format
- batch size 
- output directory for fasta file and plots
- flag -s to save background sequences in fasta format
- flag --max_iters with the maximum number of iterations to perform (for each sequence)
- flag --plot_map to plot maps with the lowest SCD for each sequence
- (optional, specific for plotting) heatmap plot limit
- (optional, specific for plotting) heatmap plot frequency
- (optional) add option --rc to average forward and reverse complement predictions
- (optional) adding --shifts k ensembles prediction shifts by k
"""

from optparse import OptionParser
import json
import os
import pickle
import random
import pandas as pd

# import h5py
import numpy as np

import matplotlib.pyplot as plt
from skimage.measure import block_reduce
import seaborn as sns

import tensorflow as tf
from basenji import seqnn, stream

from akita_utils.background_utils import create_flat_seqs_gen
from akita_utils.dna_utils import dna_1hot_to_seq
from akita_utils.utils import ut_dense
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
        default=True,
        action="store_true",
        help="Plot contact map for each background sequence [Default: %default]",
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
        action="store_false",
        help="Average forward and reverse complement predictions [Default: %default]",
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
        "-s",
        dest="save_seqs",
        default=True,
        action="store_true",
        help="Save the final seqs in fasta format",
    )
    parser.add_option(
        "--max_iters",
        dest="max_iters",
        type="int",
        default=20,
        help="maximum iterations",
    )

    (options, args) = parser.parse_args()

    if len(args) == 3:
        # single worker
        params_file = args[0]
        model_file = args[1]
        windows_tsv_file = args[2]

    elif len(args) == 5:  # muliti-GPU option
        # multi worker
        options_pkl_file = args[0]
        params_file = args[1]
        model_file = args[2]
        windows_tsv_file = args[3]
        worker_index = int(args[4])

        # load options
        options_pkl = open(options_pkl_file, "rb")
        options = pickle.load(options_pkl)
        options_pkl.close()

        # update output directory
        general_out_dir = options.out_dir
        options.out_dir = "%s/job%d" % (general_out_dir, worker_index)

    else:
        parser.error(
            "Must provide parameters and model files and TSV file with genomic windows"
        )

    if not os.path.isdir(options.out_dir):
        os.mkdir(options.out_dir)
    else:
        plot_dir = None

    options.shifts = [int(shift) for shift in options.shifts.split(",")]

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
        seq_coords_full = pd.read_csv(windows_tsv_file, sep="\t")
        seq_coords_df = split_df_equally(
            seq_coords_full, options.processes, worker_index
        )

    else:
        # read motif positions from csv
        seq_coords_df = pd.read_csv(windows_tsv_file, sep="\t")

    num_experiments = len(seq_coords_df)

    print("===================================")
    print(
        "Number of experiements = ", num_experiments
    )  # Warning! It's not number of predictions. Num of predictions is this number x5 or x6

    #################################################################
    # setup output

    # create flat sequences
    flat_seqs = create_flat_seqs_gen(
        seqnn_model,
        options.genome_fasta,
        seq_coords_df,
        max_iters=options.max_iters,
        batch_size=batch_size,
    )

    # save flat sequences in fasta format if requested
    if options.save_seqs is not None:
        with open(f"{options.out_dir}/background_seqs.fa", "w") as f:
            for i in range(len(flat_seqs)):
                f.write(
                    ">shuffled_chr"
                    + str(i)
                    + "_SCD"
                    + "{:.4f}".format(float(flat_seqs[i][2]))
                    + "\n"
                )

                f.write(dna_1hot_to_seq(flat_seqs[i][0]) + "\n")

    # plot flat sequences
    if options.plot_map is not None:
        preds = [flat_seqs[i][1] for i in range(len(flat_seqs))]
        hic_diags = params_model["diagonal_offset"]
        for no, pred in enumerate(preds):
            seq_preds = pred
            seq_map = ut_dense(seq_preds, hic_diags)  # convert back to dense
            _, axs = plt.subplots(
                1, seq_preds.shape[-1], figsize=(24, 4), sharey=True
            )

            scd_preds = np.sqrt((seq_preds**2).sum(axis=0))
            max_scd = np.max(scd_preds)

            for ti in range(seq_preds.shape[-1]):
                seq_map_ti = seq_map[..., ti]
                # TEMP: reduce resolution
                seq_map_ti = block_reduce(seq_map_ti, (2, 2), np.mean)

                sns.heatmap(
                    seq_map_ti,
                    ax=axs[ti],
                    center=0,
                    vmin=-0.6,
                    vmax=0.6,
                    cmap="RdBu_r",
                    xticklabels=False,
                    yticklabels=False,
                )

            plt.tight_layout()
            plt.savefig(f"{options.out_dir}/seq{no}_max-SCD{max_scd}.pdf")
            plt.close()


################################################################################
# __main__
################################################################################


if __name__ == "__main__":
    main()
