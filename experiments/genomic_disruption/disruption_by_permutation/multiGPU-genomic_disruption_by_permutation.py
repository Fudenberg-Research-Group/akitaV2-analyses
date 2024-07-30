# This script facilitates launching multiple computational jobs for genomic disruption analysis using permutation strategies.
# It utilizes the SLURM workload manager for job scheduling and management on HPC clusters.
#
# Command-Line Options:
# -f, --genome_fasta: Genome FASTA file for sequences.
# -m, --plot_map: Enable plotting contact maps for each allele.
# -l, --plot_lim_min: Heatmap plot limit (default: 0.1).
# --plot-freq: Frequency of heatmap plotting (default: 100).
# -o, --out_dir: Output directory for tables and plots (default: "scd").
# -c, --chrom_sizes: Table with chromosome sizes.
# --rc: Average forward and reverse complement predictions (default: False).
# --shifts: Ensemble prediction shifts (default: "0").
# --stats: Comma-separated list of stats to save (default: "SCD").
# -t, --targets_file: File specifying target indexes and labels in table format.
# --batch-size: Specify batch size.
#
# Multi-Processing Options:
# --cpu: Run without using a GPU.
# --num_cpus: Number of CPUs to use (default: 2).
# --name: SLURM job name prefix (default: "exp").
# --max_proc: Maximum concurrent processes.
# -p, --processes: Number of processes for multi-script (multi-GPU execution).
# -q, --queue: SLURM queue for job submission (default: "gpu").
# -r, --restart: Restart a partially completed job.
# --time: Time duration for job execution (default: "01:00:00").
# --gres: GPU resources (default: "gpu").
# --constraint: CPU constraints to avoid specific GPU types (default: "[xeon-6130|xeon-2640v4]").
#
# Functionality:
# - Validates input parameters and initializes necessary directories and serialization of options.
# - Constructs SLURM job commands for launching Python scripts in a specified conda environment.
# - Manages job submission and monitors job progress using SLURM utilities (`slurm_utils`).
# - Supports both single-worker and multi-GPU execution with flexible configuration options.
#
# Example Usage:
# python multiGPU-genomic_disruption_by_permutation.py params.json model_file.h5 tsv_file.tsv


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
        "-c",
        dest="chrom_sizes",
        default=None,
        help="Table with chromosome sizes [Default: %default]",
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
                "python ${SLURM_SUBMIT_DIR}/genomic_disruption_by_permutation.py %s %s %d;"
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
