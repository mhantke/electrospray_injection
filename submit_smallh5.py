#!/usr/bin/env python

import sys, os
import argparse
import datetime
import csv
import subprocess

sys.path.append("./src")  
from sizing import run_numbers

amol3416_dir = os.path.dirname(os.path.realpath(__file__))
git_sha_amol3416 = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=amol3416_dir).split("\n")[-2]

def parse_cmdline_args():
    parser = argparse.ArgumentParser(description='Sizing submission script for Davinci cluster (LCLS experiment amol3416)')
    required_grp = parser.add_argument_group('required arguments')
    required_grp.add_argument('-o', '--output-directory', metavar='/output/directory/', required=True,
                        help="output directory", type=str)

    decision_grp = parser.add_mutually_exclusive_group(required=True)
    decision_grp.add_argument('-r', '--run-number', metavar='run_number', 
                              help="run number, can also be a series of runs, for example: 1,3,5-10,20,22", type=str)
    decision_grp.add_argument('-s', '--sample', metavar='sample', 
                              help="Sample type, can be Sucrose, carboxysome, TBSV, PBCV, Trinity or Rubisco", type=str)
    parser.add_argument('--sbatch-exclude', metavar='sbatch_exclude', type=str,
                        help="List of nodes to exclude, for example c001,a001,a003")
    parser.add_argument('--sbatch-partition', metavar='sbatch_exclude', type=str, 
                        help="SLURM partition that shall be used, for example regular")

    if(len(sys.argv) == 1):
        parser.print_help()        
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_cmdline_args()
    
    if args.output_directory is None:
        print "ERROR: Output directory must to be provided. Abort."
        sys.exit(1)
    
    if args.sample is not None:
        runs = run_numbers(amol3416_dir + '/sizing.csv', args.sample)
    else:
        tmp = args.run_number
        runs = []
        for s in tmp.split(','):
            if "-" in s:
                rmin, rmax = s.split('-')
                runs += range(int(rmin), int(rmax)+1)
            else:
                runs += [int(s)]

    for run in runs:
        logfile = "%s/amol3416_r%04i_smallh5.log" % (args.output_directory, run)
        slurm   = "%s/amol3416_r%04i_smallh5.sh" % (args.output_directory, run)

        s = []
        s += "#!/bin/sh\n"
        s += "#SBATCH --job-name=r%04i/smallh5\n" % (run)
        s += "#SBATCH --ntasks=1\n"
        s += "#SBATCH --cpus-per-task=1\n"
        if args.sbatch_exclude is not None:
            s += "#SBATCH --exclude=%s\n" % args.sbatch_exclude
        if args.sbatch_partition is not None:
            s += "#SBATCH --partition=%s\n" % args.sbatch_partition
        s += "#SBATCH --output=%s\n" % logfile
        cmd = "cd scripts; "
        cmd += "./smallh5.py %i " %run
        cmd += "-o %s " %(args.output_directory)
        s += cmd + "\n"
        with open(slurm, "w") as f:
            f.writelines(s)           
        cmd = "sbatch %s" % slurm
        print cmd
        os.system(cmd)

    
