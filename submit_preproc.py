#!/usr/bin/env python

import sys, os
import argparse
import datetime
import csv

import subprocess

amol3416_dir = os.path.dirname(os.path.realpath(__file__))
git_sha_amol3416 = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=amol3416_dir).split("\n")[-2]

def parse_cmdline_args():
    parser = argparse.ArgumentParser(description='Hummingbird pre-processing submission script for Davinci cluster (LCLS experiment amol3416)')

    required_grp = parser.add_argument_group('required arguments')
    
    required_grp.add_argument('-o', '--output-directory', metavar='/output/directory/', required=True,
                        help="output directory", type=str)
    required_grp.add_argument('-r', '--run-number', metavar='run_number', required=True,
                        help="run number, can also be a series of runs, for example: 1,3,5-10,20,22", type=str)
    parser.add_argument('-y', '--run-type', metavar='run_type', required=False,
                        help="run type of the runs that shall be processed, for example: Sample,Buffer,Background", type=str, default="Sample,Buffer,Background")
    parser.add_argument('-n', '--number-of-frames', metavar='number_of_frames',
                        help="number of frames to be processed (optional)", type=int)
    parser.add_argument('-p', '--number-of-processes', metavar='number_of_processes',
                        help="number of MPI processes to be allocated for job", type=int, default=12)
    parser.add_argument('-e', '--env', metavar='env',
                        help="bash environment file that will be sourced before processing (optional)", type=str)
    parser.add_argument('-l', '--output-level', metavar='output_level',
                        help="output level defines how much data per event will be stored (default=3, 0: no output (\"dry-run\"), 1: only scalar output (hitscores, GMD values, etc.), 2: scalar output and TOF data, 3: scalar output, TOF data and images)", type=int, default=3)
    parser.add_argument('-t', '--hitscore-threshold', metavar='hitscore_threshold',
                        help="Hitscore threshold [if not provided read from CSV file]", type=int)
    parser.add_argument('-g', '--photon-energy-ev', metavar='photon_energy_ev',
                         help="Manually set nominal photon energy in unit electron volts (used for example for pnCCD gain calculation and hit finding) [if not provided read from CSV file]",
                         type=float)
    parser.add_argument('-a', '--rescaling-factors-asics', metavar='rescaling_factors_asics',
                        help="Manually set the 4 rescaling factors for upper right quadrant (for example 2.0,1.0,1.0,1.0) [if not provided read from CSV file]",
                        type=str)
    parser.add_argument('-R', '--do-raw', default=0, metavar='do_raw',
                        help="Only do pedestal correction, no other corrections", type=int)
    parser.add_argument('-Q', '--do-quad', default=1, metavar='do_quad',
                        help="Correct artefacts in faulty (\"upper right\") quadrant", type=int)
    parser.add_argument('-C', '--do-cmc', default=1, metavar='do_cmc',
                        help="Subtract common mode in all pixel rows and columns", type=int)
    parser.add_argument('-M', '--do-metrology', default=1, metavar='do_metrology',
                        help="Move pixels to their physical locations", type=int)
    parser.add_argument('--sbatch-exclude', metavar='sbatch_exclude', type=str,
                        help="List of nodes to exclude, for example c001,a001,a003")
    parser.add_argument('--sbatch-partition', metavar='sbatch_exclude', type=str, 
                        help="SLURM partition that shall be used, for example regular", default="fast")
    parser.add_argument('--sha-version-check', default=1, metavar='sha_version_check',
                        help="Perform version check for amol3416 git repo and hummingbird git repo whenever sbatch script is executed", type=int)

    
    if(len(sys.argv) == 1):
        parser.print_help()
        
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_cmdline_args()
    
    if args.output_directory is None:
        print "ERROR: Output directory must to be provided. Abort."
        sys.exit(1)
    if args.run_number is None:
        print "ERROR: Run number must to be provided. Abort."
        sys.exit(1)
    if args.number_of_processes == 2:
        print "ERROR: The current implementation does not allow the number of processes be two. Change your configuration and try again. Abort."
        sys.exit(1)

    if args.env is None:
        env = "%s/source_this_on_davinci" % amol3416_dir
        print "WARNING: Using Max' environment specified in %s" % env
    else:
        env = args.env

    hummingbird_dir = os.path.dirname(os.path.realpath(subprocess.check_output("source %s; which hummingbird.py" % env, executable="/bin/bash", shell=True).split("\n")[-2]))
    git_sha_hummingbird = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=hummingbird_dir).split("\n")[-2]

    if args.sha_version_check:

        if len(subprocess.check_output(['git', 'diff'], cwd=amol3416_dir)) > 0:
            print "There are uncommitted changes in the working tree of the loaded Git repository \'amol3416\' (%s)." % amol3416_dir
            print "For reproducibility reasons please commit first before submitting jobs for running the pre-processing routine."
            sys.exit(1)
        
        if len(subprocess.check_output(['git', 'diff'], cwd=hummingbird_dir)) > 0:
            print "There are uncommitted changes in the working tree of the loaded Git repository \'hummingbird\' (%s)." % hummingbird_dir
            print "For reproducibility reasons please commit first before submitting jobs for running the pre-processing routine."
            sys.exit(1)

    # Read run numbers from CSV table
    data_mode = 'amol3416'
    filename_csv = '%s/%s_run_params.csv' % (amol3416_dir, data_mode)
    with open(filename_csv, "r") as f:
        _reader = csv.DictReader(f, delimiter='\t')
        reader = [r for r in _reader]
        all_runs = [int(row["RunNr"]) for row in reader]
        type_runs = [int(row["RunNr"]) for row in reader if row["RunType"] in args.run_type.split(",")]
                
    # Construct list of runs that shall be processed
    if args.run_number == "all":
        runs = type_runs
    else:
        runs_tmp = []
        tmp = args.run_number
        for s in tmp.split(','):
            if "-" in s:
                rmin, rmax = s.split('-')
                runs_tmp += range(int(rmin), int(rmax)+1)
            else:
                runs_tmp += [int(s)]
        # Check whether the selected run numbers are all listed in the CSV file, if not skip
        runs = []
        for run in runs_tmp:
            if run not in all_runs:
                print "WARNING: Requested run number %i is not listed in %s and is therefore skipped for processing." % (run, filename_csv)
            elif run not in type_runs:
                print "NOTE: Requested run number %i has wrong type and is therefore skipped for processing." % (run)
            else:
                runs.append(run)
            
    for run in runs:
        logfile = "%s/amol3416_r%04i_ol%i.log" % (args.output_directory, run, args.output_level)
        slurm   = "%s/amol3416_r%04i_ol%i.sh" % (args.output_directory, run, args.output_level)

        s = []
        s += "#!/bin/sh\n"
        s += "#SBATCH --job-name=r%04i/%i\n" % (run, args.output_level)
        s += "#SBATCH --ntasks=%i\n" % args.number_of_processes
        s += "#SBATCH --ntasks-per-node=%i\n" % args.number_of_processes
        s += "#SBATCH --cpus-per-task=1\n"
        if args.sbatch_exclude is not None:
            s += "#SBATCH --exclude=%s\n" % args.sbatch_exclude
        s += "#SBATCH --partition=%s\n" % args.sbatch_partition
        s += "#SBATCH --output=%s\n" % logfile
        cmd = "source %s; " % env
        cmd += "echo $HOSTNAME; "
        cmd += "mpirun -n %i -wd %s " % (args.number_of_processes, amol3416_dir)
        cmd += "hummingbird.py -b conf_preproc.py --lcls-run-number %i --batch-mode --out-dir %s" % (run, args.output_directory)
        if args.hitscore_threshold is not None:
            cmd += " --hitscore-threshold %i" % args.hitscore_threshold
        if args.number_of_frames is not None:
            cmd += " --lcls-number-of-frames %i" % args.number_of_frames
        cmd += " --output-level %i" % args.output_level
        cmd += " --do-raw %i --do-quad %i --do-cmc %i --do-metrology %i" % (args.do_raw, args.do_quad, args.do_cmc, args.do_metrology)
        if args.sha_version_check:
            cmd += " --check-sha-amol3416 %s " % git_sha_amol3416
            cmd += " --check-sha-hummingbird %s " % git_sha_hummingbird
        s += cmd + "\n"
        with open(slurm, "w") as f:
            f.writelines(s)           
        cmd = "sbatch %s" % slurm
        print cmd
        os.system(cmd)

    
