# Description:
# This script manages the execution of the flank-core compatibility experiment in a SLURM-managed cluster environment.
# It prepares the environment for running a specific Python script in parallel across multiple nodes or GPUs.
# The script handles job submission, manages output directories, and supports job resumption.
#
# Inputs:
# - params_file: Path to the JSON file containing model parameters.
# - model_file: Path to the file containing the trained model.
# - tsv_file: Path to the TSV file with input data.
#
# Options:
# - -f, --genome_fasta: Path to the genome FASTA file [Default: %default].
# - -m, --plot_map: Whether to plot contact maps for each allele [Default: %default].
# - -l, --plot_lim_min: Minimum limit for heatmap plots [Default: %default].
# - --plot-freq: Frequency for heatmap plotting [Default: %default].
# - -o, --out_dir: Output directory for tables and plots [Default: %default].
# - --rc: Whether to average forward and reverse complement predictions [Default: %default].
# - --shifts: Ensemble prediction shifts [Default: %default].
# - --stats: Comma-separated list of statistics to save [Default: %default].
# - -t, --targets_file: File specifying target indexes and labels in table format.
# - --batch-size: Batch size for processing [Default: %default].
# - --save-maps: Whether to save all maps in the H5 file [Default: %default].
# - --background-file: File with insertion sequences in FASTA format [Default: %default].
# - --cpu: Run without a GPU [Default: %default].
# - --num_cpus: Number of CPUs to use [Default: %default].
# - --name: SLURM job name prefix [Default: %default].
# - --max_proc: Maximum number of concurrent processes [Default: %default].
# - -p, --processes: Number of processes for multi-script setup.
# - -q, --queue: SLURM queue to run jobs on [Default: %default].
# - -r, --restart: Whether to restart a partially completed job [Default: %default].
# - --time: Time allocation for each job [Default: %default].
# - --gres: GPU resources required [Default: %default].
# - --constraint: Constraints for CPU to avoid specific GPUs [Default: %default].
#
# Outputs:
# - Creates an output directory.
# - Submits SLURM jobs to execute the specified Python script (`virtual_insertion_flank_core_compatibility.py`).
# - Supports job resumption if partially completed.
#
# Example command-line usage:
# python multiGPU-virtual_insertion_flank_core_compatibility.py params.json model.h5 input_data.tsv -f /path/to/genome.fa -m -o results --rc --shifts "1,2" -t targets.tsv --batch-size 16 --save-maps --cpu --num_cpus 4 --name my_experiment --max_proc 10 -q gpu --time 02:00:00 --gres gpu --constraint "[xeon-6130|xeon-2640v4]"

from optparse import OptionParser
import os
import pickle
import akita_utils.slurm_utils as slurm
from akita_utils.h5_utils import job_started


################################################################################
# main
################################################################################


def main():
    usage = "usage: %prog [options] <params_file> <model_file> <tsv_file>"
    parser = OptionParser(usage)
    parser.add_option(
        "-f",
        dest="genome_fasta",
        default=None,
        help="Genome FASTA for sequences [Default: %default]",
    )
    parser.add_option(
        "-m",
        dest="plot_map",
        default=False,
        action="store_true",
        help="Plot contact map for each allele [Default: %default]",
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
        "-o",
        dest="out_dir",
        default="scd",
        help="Output directory for tables and plots [Default: %default]",
    )
    parser.add_option(
        "--rc",
        dest="rc",
        default=False,
        action="store_true",
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
        "--stats",
        dest="stats",
        default="SCD",
        type="str",
        help="Comma-separated list of stats to save. [Default: %default]",
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
        default=None,
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
    ## insertion-specific options
    parser.add_option(
        "--background-file",
        dest="background_file",
        default=None,
        help="file with insertion seqs in fasta format",
    )

    # multi
    parser.add_option(
        "--cpu",
        dest="cpu",
        default=False,
        action="store_true",
        help="Run without a GPU [Default: %default]",
    )
    parser.add_option(
        "--num_cpus",
        dest="num_cpus",
        default=2,
        type="int",
        help="Number of cpus [Default: %default]",
    )
    parser.add_option(
        "--name",
        dest="name",
        default="exp",
        help="SLURM name prefix [Default: %default]",
    )
    parser.add_option(
        "--max_proc",
        dest="max_proc",
        default=None,
        type="int",
        help="Maximum concurrent processes [Default: %default]",
    )
    parser.add_option(
        "-p",
        dest="processes",
        default=None,
        type="int",
        help="Number of processes, passed by multi script",
    )
    parser.add_option(
        "-q",
        dest="queue",
        default="gpu",
        help="SLURM queue on which to run the jobs [Default: %default]",
    )
    parser.add_option(
        "-r",
        dest="restart",
        default=False,
        action="store_true",
        help="Restart a partially completed job [Default: %default]",
    )
    parser.add_option(
        "--time",
        dest="time",
        default="01:00:00",
        help="time to run job. [Default: %default]",
    )
    parser.add_option(
        "--gres",
        dest="gres",
        default="gpu",
        help="gpu resources. [Default: %default]",
    )
    parser.add_option(
        "--constraint",
        dest="constraint",
        default="[xeon-6130|xeon-2640v4]",
        help="cpu constraints to avoid the a40 gpus. [Default: %default]",
    )

    (options, args) = parser.parse_args()

    if len(args) != 3:
        parser.error("Must provide parameters and model files and TSV file")
    else:
        params_file = args[0]
        model_file = args[1]
        tsv_file = args[2]

    #######################################################
    # prep work

    # output directory
    if not options.restart:
        if os.path.isdir(options.out_dir):
            # print("Please remove %s" % options.out_dir, file=sys.stderr)
            exit(1)
        os.mkdir(options.out_dir)

    # pickle options
    options_pkl_file = "%s/options.pkl" % options.out_dir
    options_pkl = open(options_pkl_file, "wb")
    pickle.dump(options, options_pkl)
    options_pkl.close()

    #######################################################
    # launch worker threads
    jobs = []
    for pi in range(options.processes):
        if not options.restart or not job_started(options, pi):
            cmd = 'eval "$(conda shell.bash hook)";'
            cmd += "conda activate basenji_py3.9_tf2.15;"
            cmd += (
                "python ${SLURM_SUBMIT_DIR}/virtual_insertion_flank_core_compatibility.py %s %s %d;"
                % (
                    options_pkl_file,
                    " ".join(args),
                    pi,
                )
            )
            cmd += "conda deactivate;"

            name = "%s_p%d" % (options.name, pi)
            outf = "%s/job%d.out" % (options.out_dir, pi)
            errf = "%s/job%d.err" % (options.out_dir, pi)

            num_gpu = 1 * (not options.cpu)

            j = slurm.Job(
                cmd,
                name,
                outf,
                errf,
                queue=options.queue,
                gpu=num_gpu,
                gres=options.gres,
                mem=15000,
                time=options.time,
                cpu=options.num_cpus,
                constraint=options.constraint,
            )
            jobs.append(j)

    slurm.multi_run(
        jobs,
        max_proc=options.max_proc,
        verbose=False,
        launch_sleep=10,
        update_sleep=10,
    )


################################################################################
# __main__
################################################################################

if __name__ == "__main__":
    main()
