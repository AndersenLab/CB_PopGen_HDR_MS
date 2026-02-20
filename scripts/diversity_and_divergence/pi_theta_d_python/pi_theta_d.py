import zarr 
import pandas as pd
import numpy as np
import sys
import importlib
import argparse
import allel
import os
import pyranges as pr
#from scripts.src.vcf_stats.diversity import mask_divergent
#from scripts.src.vcf_stats.diversity import chrom_subregion_diversity, chrom_window_diversity

#From Notebooks level folder
#sys.path.append('../src/')
#From root
sys.path.append('./')

from vcf_stats import vcf_stat
from vcf_stats import diversity
from vcf_stats import helpers

def parse_commandline(): 
    """Parse the arguments given on the command-line.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--zarr",
                       help="path to zarr group",
                       default=None)                  
    parser.add_argument("--chrom_lengths",
                help="path to chrom_lengths pickle",
                default=None)
    parser.add_argument("--samples",
                help="path to one or more samples.txt file",
                default=None)
    parser.add_argument("--mask_divergent",
                help="Wether or not to mask divergent positions from output calculations only some value needs to be passed ex. T",
                default=None)
    parser.add_argument("--divergent_file",
                    help="Path to divergent mask file (pkl format)",
                    default=None)
    parser.add_argument("--mask_missing",
                        help="To mask missing genotypes from output calculations only some value needs to be passed ex. T",
                        default=None)
    parser.add_argument("--mask_repeats",
                    help="Whether to mask repetitive regions; pass a non-null value (e.g., T) to enable",
                    default=None)
    parser.add_argument("--repeats_file",
                    help="Path to repeats mask file (pkl format)",
                    default=None)
    parser.add_argument("--missing_threshold",
                        help = "Proportion of missing genotypes to maske (.80)",
                        default= 0.80)
    parser.add_argument("--window_size",
                        help = "Size of non-ovlerapping window to calculate diversity",
                        default= 10000) 
    parser.add_argument("--out",
                        help = "path to output directory",
                        default= 0.80)
    # New parameters
    parser.add_argument("--filter_monomorphic",
                        help="Filter monomorphic sites (T/F)",
                        default="T")
    parser.add_argument("--maf_threshold",
                        help="MAF threshold (e.g., 0.05)",
                        type=float,
                        default=None)


    args = parser.parse_args()
    

    return args


args = parse_commandline()
zarr_path = args.zarr 
chrom_dict = args.chrom_lengths
# mask_missing = args.mask_missing
window_size = int(args.window_size)
out_folder = args.out

#Check Outdir Status
check = os.path.isdir(out_folder)

if not check:
    os.makedirs(out_folder)
    print("create folder :", out_folder)
else:
    print(out_folder, "folder already exists")

#Load the Zarr & Chromsome data
zarr_file = zarr.open_group(zarr_path, mode= 'r')
chrom_length = pd.read_pickle(chrom_dict)

#Check Sample Filtering
samples_file = args.samples
if samples_file is not None:
    print("Filtering samples")
    sample_list  = np.loadtxt(samples_file, dtype= str)
else:
    sample_list = None
    

## Site Masking ##
mask_divergent = args.mask_divergent
mask_missing = args.mask_missing
mask_repeats = args.mask_repeats

if mask_divergent is not None:
    print("Loading Masking data")
    divergent_bool = pd.read_pickle(args.divergent_file)
else:
    print("Repeats divergent enabled, but no divergent file provided!")
    divergent_bool = None


if mask_missing is not None:
    threshold = args.missing_threshold
    print("masking variants with greater than: " , threshold)


if mask_repeats is not None:
    if args.repeats_file is not None:
        print("Loading repeats masking data")
        repeats_bool = pd.read_pickle(args.repeats_file)
    else:
        print("Repeats masking enabled, but no repeats file provided!")
        repeats_bool = None


chroms = ["I", "II", "III", "IV", "V", "X"]


genome_subregions = pd.DataFrame()
genome_windows = pd.DataFrame()
for chrom in chroms:
    #Pull the data
    gt , ac , pos = helpers.get_chrom_data(zarr_file, chrom, snps = True, numalt = 1, sample_list = sample_list,filter_monomorphic=(args.filter_monomorphic == "T"))
    
    # filter according to MAF (use functions in diversity.py）
    if args.maf_threshold is not None:
        maf_mask = diversity.filter_by_maf(ac, gt, args.maf_threshold)
        gt = gt[maf_mask]
        ac = ac[maf_mask]
        pos = pos[maf_mask]
        print(f"Filtered {np.sum(~maf_mask)} sites with MAF < {args.maf_threshold}")



    ## Perfrom site masking if required ## 
    mask = np.ones(len(pos), dtype=bool)

    #Divergent region mask
    if mask_divergent is not None:
        if divergent_file is not None:
            print("Masking divergent sites on: ", chrom)
            chrom_div_bool = divergent_bool[chrom]
            mask = np.logical_and(mask, chrom_div_bool)
        else:
            print("Repeats mask requested but divergent file is not loaded.")

    #Missingness Mask
    if mask_missing is not None:
        print("Masking missing sites on: ", chrom)
        missing_mask = diversity.create_missing_mask(gt, pos, threshold, chrom, chrom_length)
        mask = missing_mask 
    

    
    # Repeats mask
    if mask_repeats is not None:
        if repeats_bool is not None:
            print("Masking repeats on:", chrom)
            chrom_repeats_mask = repeats_bool[chrom]
            mask = np.logical_and(mask, chrom_repeats_mask)
        else:
            print("Repeats mask requested but repeats_bool is not loaded.")


    ## Calculate Sub-Region Stats ##
    chrom_regions = diversity.chrom_subregion_diversity(chrom, gt, ac, pos, chrom_length, mask=mask)
    genome_subregions = pd.concat([genome_subregions, chrom_regions], ignore_index = True, sort = False)
    
    #Calculate windowed stats
    chrom_windows, window_df = diversity.chrom_window_diversity(chrom, gt, ac, pos, chrom_length, mask = mask, return_div = False, window_size = window_size)
    genome_windows = pd.concat([genome_windows, chrom_windows], ignore_index= True, sort = False)
    

genome_subregions.to_csv(out_folder + "/" + "sub_region_diversity.csv" )
genome_windows.to_csv(out_folder + "/" + "chromosome_windows_diversity.csv")
window_df.to_csv(out_folder + "/" + "windows.csv")