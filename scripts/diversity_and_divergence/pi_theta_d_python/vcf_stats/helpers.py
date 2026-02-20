from os import remove
import allel
import zarr
#from vcf_stats import samples
import numpy as np


def filter_samps(samps, gt, callset):
    """"
    A function to return a filtered GT array containing samples

    samps: an array of samples you are interested in ex) `np.array(["JU1246", "JU1249"])`
    gt: a genotype array
    callset: Loaded Zarr hierarcy group ex) zarr_file["I"] where zarr file is an open zarr group. 
        This is included for portability. If you are using this as a helper function you may have already loaded this data
    """
    vcf_samples = callset['samples'][:]

    #Find indexes of samples 
    indexes = np.in1d(vcf_samples, samps).nonzero()[0]
    print(indexes)
    #Subset the GT array 
    subset_gt = gt.take(indexes, axis = 1)
    return(subset_gt)


def get_chrom_data(zarr, chrom, snps = True, numalt = 1, sample_list=None, filter_monomorphic=True):
    """
    Pulls basic chrom data from full VCF zarr
    ** Add check for numalt you idiot **

    zarr: open zarr file
    chrom: str id for chromosome ex) "I"
    snps: Logical - filter to just snps
    numalt: max number of alt alleles
    samples: a np array of samples to subset
    returns: tuple of:
         GT array,
         allele count array - note this is a raw allel count
         pos array
    New parameter (20250403):
        filter_monomorphic (bool)
        
    """
    callset = zarr[chrom]
    pos = allel.SortedIndex(callset['variants/POS'])
    gt = allel.GenotypeArray(callset['calldata/GT'])
    
    if snps == True: 
        print("Filtering for SNVS")
        snps = callset['variants/is_snp'][:]
        
        pre = gt.shape[0] #num var
        
        gt = gt.compress(snps, axis = 0) #Filter to just snp variants
        
        post = gt.shape[0] #num var
        
        removed_variants = pre - post
        
        print("There were", removed_variants, "non-snp variants removed")
        
        #Remove the non-snps from the Position Index
        pos = pos[snps]

        ac = gt.count_alleles() 
        
        if sample_list is not None:
            #Filter to samples in samples list
            print("Filtering samples")
            gt = filter_samps(samps = sample_list, gt = gt, callset = callset)
            print(gt)
            #Get a new allele count obj
            ac = gt.count_alleles()
            print(ac)

    # filter out monomorphic sites
    if filter_monomorphic:
        is_monomorphic = (ac[:, 0] == 0) | (ac[:, 1] == 0)  # AC = 0 or reference-only / alternate-only alleles
        gt = gt[~is_monomorphic]
        pos = pos[~is_monomorphic]
        ac = ac[~is_monomorphic]
        print(f"Removed {np.sum(is_monomorphic)} monomorphic sites")
    else: 
        print("Returning all variants")
        ac = gt.count_alleles() 
    
    return(gt, ac, pos)

