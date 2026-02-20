## Script to process a single chromsome from VCF into a whole genome Zarr directory ##

import multiprocessing
from multiprocessing import Pool
import os
import allel
import numpy as np
import sys
import argparse
import numcodecs
import zarr
#from joblib import Parallel, delayed


def parse_commandline(): 
    """Parse the arguments given on the command-line.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf",
                       help="path to vcf",
                       default=None)                  
    parser.add_argument("--out",
                help="path to directory where zarr will be built",
                default=None)
    parser.add_argument("--nproc",
                help="Job number",
                type=int,
                default=1) 

    args = parser.parse_args()

    return args

def pre_zarr_check(vcf_path):
    """
    Reads in VCF to numpy arrays and check number of variants and max alt alleles 
    
    vcf_path: path to tabix indexed vcf.gz
    region: str of region to check ex. 'I' (tabix syntax)
    Returns:a tuple: num_var , max_var
    
    """
    chrom_callset = allel.read_vcf(vcf_path, fields=['numalt', 'CHROM'], log=sys.stdout)
    chrom_num_alt = chrom_callset['variants/numalt']
    
    max_var = np.max(chrom_num_alt)
    #Get the number of variants per chromosme
    chrom_arr = chrom_callset['variants/CHROM']
    chrom_list = ["I", "II", "III", "IV", "V", "X"]

    chrom_var_list = []
    for chrom in chrom_list:
        filtered_chrom_arr = chrom_arr[np.where(chrom_arr == chrom)]
        num_var = len(filtered_chrom_arr)
        chrom_var_list.append(num_var)
    #num_var = len(chrom_num_alt)
    return chrom_var_list , max_var

def check_zarr_num_var(callset, chrom): 
    """
    Checks the number of variants in a zarr callset by chromosome

    callset: an open zarr file - opened with zarr.open_group(zarr_path)
    chrom_list: list of grouping ids in zarr data base (Typically - 'I')
    returns: list containing the number of variants per chromosome, order reflects order in chrom_lists
    """
    #num_varr_post = [] # Commented out code to fun this across multipe chromosomes
    #for chrom in chrom_list: 
    callset_id = f'{chrom}/variants/CHROM'
    chrom_callset = callset[callset_id]
    num_var = chrom_callset.shape[0]
    return(num_var)
    #num_varr_post.append(num_var)
    #return(num_varr_post)

def read_stats(stats_path): 
    """
    pull number of variants per chromosome and max alt from output of get_stats.py

    """
    with open(stats_path, 'r') as f: 
        lines = f.readlines()
        chrom_array = lines[0].strip().split(sep= ",") 
        max_alt = int(lines[1].strip()) 
    return(chrom_array, max_alt)



if __name__ == '__main__':
    

    args = parse_commandline()
    vcf_path = args.vcf
    print("Processing vcf from: " + vcf_path )
    zarr_path = args.out
    print("Processing zarr to: " + zarr_path )

    chroms = ["I", "II", "III", "IV", "V", "X"]
    n_chroms = len(chroms)
    ## Pull Pre-conversion stats ## 
    
    #load stats 
    print("Checking pre-vcf stats")
    pre_num_var , max_alt = pre_zarr_check(vcf_path)

    ## Conversion ##

    def conversion_wrapper(chrom):
        #print("Pulling from", vcf_path)
        res = allel.vcf_to_zarr(
            vcf_path, 
            zarr_path, 
            group= chrom, 
            fields="*", #This needs to be another field this one doesnt exist in our vcf
            alt_number = max_alt, #set this to max alt value
            region=chrom, 
            log=sys.stdout, 
            compressor=numcodecs.Blosc(cname='zstd', clevel=1, shuffle=False))
        return(res)
        
    print("Starting parallele conversion")
    print("Spreading tasks over: " + str(args.nproc))
    with Pool(args.nproc) as p:
        print(p.map(conversion_wrapper, ["I", "II", "III", "IV", "V", "X"]))



    ## QC ##
    print("starting qc")
    callset = zarr.open_group(zarr_path, mode = 'r')
    
    print("Checking number of variants")
    post_var_per_chrom = []
    for chrom in chroms:
        chrom_var = check_zarr_num_var(callset, chrom) #Check variants per chromosome in zarr
        post_var_per_chrom.append(chrom_var)
    print(pre_num_var)
    print(post_var_per_chrom)



        





