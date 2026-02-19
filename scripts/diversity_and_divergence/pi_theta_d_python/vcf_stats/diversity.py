import allel
import numpy as np
#import sys
import statistics
import pandas as pd 
import pyranges as pr
### Divergent Region functions ###

def call_divergent_windows(bed_path, windows, chrom): 
    """
    Helper function for windows plots to identify windows that overlap with divergent regions

    bed_path: path to all isotypes divegent bed
    windows: a numpy array returned from scikit-allele windowed statistic *Must be the same windows used to calculate x-axis
    chrom: chromosome to filter the bedfile to
    returns: overlaps bool from a windows object True means window overlaps with with range
    """

    #Load the bed file 
    hd_bed = pr.read_bed(bed_path)
    chrom_bed = hd_bed[chrom]

    #Create a pyranges obj from the window
    starts = windows[:,0]
    ends = windows[:, 1]
    windows_gr = pr.PyRanges(chromosomes= chrom, starts = starts, ends = ends)

    #Check for overlaps and store overlapping windows in a dictionary with start as keys
    overlaps = windows_gr.overlap(chrom_bed)
    overlaps = overlaps.df #Convert to a dataframe
    overlap_windows = overlaps.set_index('Start').T.to_dict() #Convert the overlaps to dict that can be quickly searched

    #Check if window is in divergent region dictionary
    divergent = []
    for i in range(len(windows)):
        start = windows[i, 0]
        if start in overlap_windows:
            div_bool = True
        else: 
            div_bool = False
        divergent.append(div_bool)
    return(divergent)


# New(20250403) Function: filtering according to allele frequency 
def filter_by_maf(ac, gt, maf_threshold=0.05):
    """
    Filter sites based on minor allele frequency (MAF), taking missing data into account.
    Args:
        ac (ndarray)
        gt (GenotypeArray)
        maf_threshold (float)
    Returns:
        mask (ndarray): True (keep the site)
    """
    # Calculate the number of missing samples at each site
    n_missing = gt.count_missing(axis=1)
    n_samples = gt.shape[1]
    n_valid = n_samples * 2 - n_missing  # Effective number of alleles (assuming diploidy)

    # Calculate the minor allele frequency
    maf = np.min([ac[:, 0]/n_valid, ac[:, 1]/n_valid], axis=0)
    return maf >= maf_threshold


# New(20250403) Function: filter monomorphic sites 
def filter_monomorphic_sites(ac):
    """
    Filter out monomorphic sites
    Args:
        ac (ndarray)
    Returns:
        mask (ndarray): True (keep the site)
    """
    is_monomorphic = (ac[:, 0] == 0) | (ac[:, 1] == 0)
    return ~is_monomorphic


# Functions to calculate and plot Diveristy data over a subregion
def chrom_subregion_diversity(chrom, gt, ac, pos, genome_dict, mask=None, maf_threshold=None):
    # New (20250403)
    if maf_threshold is not None:
        maf_mask = filter_by_maf(ac, gt, maf_threshold)
        gt = gt[maf_mask]
        ac = ac[maf_mask]
        pos = pos[maf_mask]

    chrom_results = []
    for key in genome_dict[chrom]: 
        print(key)
        sub_region = genome_dict[chrom][key]
        start = sub_region["Start"]
        stop = sub_region["Stop"]
        
        #Get mid point of the region for plotting
        t = list(range(start, stop))
        mid = statistics.median(t)

        #Subset the chromosome data 
        try:
            loc_region = pos.locate_range(start, stop - 1)
        except KeyError:
            print(f"Warning: chromosome {chrom} region '{key}' ({start}-{stop-1}) contains no variants, skipping.")
            continue
        pos_region = pos[loc_region]
        ac_region = ac[loc_region]
        gt_region = gt[loc_region]
        
        if mask is not None:
            # Subset mask down to just the region of interest
            #mask_region = mask[start+1 : stop + 2 ] # +2 (+1 to get 0 based + 2 since np is not Right inclusive)
           
            #pull out just masked positions
            masked_positions = np.where(mask == False)
            #print(type(masked_positions))
            masked_positions = masked_positions[0] + 1 
            #Check if variants should be masked
            var_stat = np.isin(np.array(pos_region), masked_positions)
            n_masked = np.count_nonzero(var_stat)
            print("Masking: ", n_masked, " variants")
            
            #Create a GT mask that can be applied to GT array 
            samples = int(gt_region.shape[1])
            gt_mask = create_gt_mask(var_stat, samples)
            
            #Mask sites and get tajimas ac 
            masked_gt_region = gt_region
            masked_gt_region.mask = gt_mask
            ac_d_region = masked_gt_region.count_alleles()
        else:
            mask_region = None
            ac_d_region = ac_region
        #Get the number of variants within the subregion
        n_var = pos_region.shape[0]

        #Unmasked stats
        pi = allel.sequence_diversity(pos_region, 
                                     ac_region, 
                                     start = start, 
                                     stop = stop -1, 
                                     is_accessible= mask)
        #print(pi)
        theta = allel.watterson_theta(pos_region,
                                     ac_region,
                                     start = start, 
                                     stop = stop -1, 
                                     is_accessible= mask)
        #rint(theta)
        d = allel.tajima_d(ac_d_region, 
                          pos_region,
                          start = start,
                          stop = stop -1)
        #print(d)
        #Masked stats



        region_dict = {
            'chrom' : chrom,
            'sub_region': key, 
            'pi': pi, 
            'theta': theta, 
            'd': d, 
            'mid': mid,
            'nvar': n_var}

        chrom_results.append(region_dict)
        chrom_sub_region_df = pd.DataFrame(chrom_results)
    return(chrom_sub_region_df)


def chrom_window_diversity(chrom, gt, ac, pos, genome_dict, mask = None, return_div = False, window_size = 1000, maf_threshold=None):
    """
    mask: an array len = len(chrom) True = included in calculation
    
    Calculates windowed statistics across a chromosome
        Current stats:
           - pi 
           - theta 
           - tajimas
    Retuns data frame [chrom, x, stat_type, stat]
    """
    # New (20250403)
    if maf_threshold is not None:
        maf_mask = filter_by_maf(ac, gt, maf_threshold)
        gt = gt[maf_mask]
        ac = ac[maf_mask]
        pos = pos[maf_mask]
    
    # #Check Mask Status
    # if mask is not None:
    #     # pi theta mask
    #     pi_theta = mask[chrom]
    # else:
    #     pi_theta_divergent = None
    
    #Build Mask for Tajimas
    if mask is not None:
        #Find what positions should be excluded from the calculation
        masked_positions = np.where(mask == False) 
        masked_positions = masked_positions[0] + 1 #Add one to convert from one based index to genome postion
        #Check if variants are in masked position 
        var_stat = np.isin(np.array(pos), masked_positions) #Bool True - should be masked
        n_masked = np.count_nonzero(var_stat)
        print("Masking: ", n_masked, " variants")
        
        #Create GT mask that can be applied to GT array
        samples = int(gt.shape[1])
        gt_mask = create_gt_mask(var_stat, n_samples = samples) 
        
        #Mask sites and get tajimas ac
        masked_gt = gt 
        masked_gt.mask = gt_mask
        ac_d = masked_gt.count_alleles()
    else:
        ac_d = ac
    
    #Missing Data
    n_missing = gt.count_missing(axis= 1)
    n_samples = int(gt.shape[1])

    def PropGT_windowed(wv):
        n_var_win = len(wv)
        #print("There are: " + str(n_var_win) +  " in the window")
        total_gt = n_var_win * n_samples
        missing = np.sum(wv) #add the total number of missing samples
        prop = missing / total_gt
        return(prop)

    mis , windows_m, counts = allel.windowed_statistic(pos, n_missing, statistic= PropGT_windowed, size = window_size, fill = 0)


    #Tajimas
    d , windows_d, values = allel.windowed_tajima_d(pos = pos, 
                                                   ac = ac_d, 
                                                   size = window_size, 
                                                   start = genome_dict[chrom]['full']['Start'] , 
                                                   stop = genome_dict[chrom]['full']['Stop'])
    #Theta
    theta, windows_theta, n_bases, counts = allel.windowed_watterson_theta(pos, ac, size = window_size, is_accessible= mask, start = genome_dict[chrom]['full']['Start'], stop = genome_dict[chrom]['full']['Stop'])
    
    #Pi
    pi, windows_pi, n_bases, counts = allel.windowed_diversity(pos, ac, size = window_size, is_accessible= mask, start = genome_dict[chrom]['full']['Start'], stop = genome_dict[chrom]['full']['Stop'])
    
    #Unmasked X-axis - Should be the same so producing them all is redundant
    theta_x = windows_theta.mean(axis = 1)
    pi_x = windows_pi.mean(axis = 1)
    d_x = windows_d.mean(axis = 1)
    mis_x = windows_m.mean(axis = 1)

    if return_div == True:
        #calculate divergent regions 
        print("WARNING: Return div should only be run on C. elegans")
        bed_path = None
        div_pi_theta = call_divergent_windows(bed_path, windows_theta, chrom)
        div_mis = call_divergent_windows(bed_path, windows_m, chrom )
        div_d = call_divergent_windows(bed_path, windows_d, chrom )

        theta_dict = {
        "chrom" : chrom,
        "x" : theta_x,
        "window_start" : windows_theta[...,0],
        "window_stop" : windows_theta[...,1],
        "stat_type" : "theta",
        "stat" : theta, 
        "div_stat": div_pi_theta
        }
        pi_dict = {
        "chrom" : chrom,
        "x" : pi_x,
        "window_start" : windows_pi[...,0],
        "window_stop" : windows_pi[...,1],
        "stat_type" : "pi",
        "stat" : pi, 
        "div_stat": div_pi_theta
        }
        d_dict = {
        "chrom" : chrom,
        "x" : d_x,
        "window_start" : windows_d[...,0],
        "window_stop" : windows_d[...,1],
        "stat_type" : "d",
        "stat" : d, 
        "div_stat": div_d
        }
        mis_dict = {
        "chrom" : chrom,
        "x" : mis_x,
        "window_start" : windows_m[...,0],
        "window_stop" : windows_m[...,1],
        "stat_type" : "mis",
        "stat" : mis, 
        "div_stat" : div_mis
        }
    else:
    #Append the stats to a dictionary
        theta_dict = {
        "chrom" : chrom,
        "x" : theta_x,
        "window_start" : windows_theta[...,0],
        "window_stop" : windows_theta[...,1],
        "stat_type" : "theta",
        "stat" : theta
        }
        pi_dict = {
        "chrom" : chrom,
        "x" : pi_x,
        "window_start" : windows_pi[...,0],
        "window_stop" : windows_pi[...,1],
        "stat_type" : "pi",
        "stat" : pi
        }
        d_dict = {
        "chrom" : chrom,
        "x" : d_x,
        "window_start" : windows_d[...,0],
        "window_stop" : windows_d[...,1],
        "stat_type" : "d",
        "stat" : d
        }
        mis_dict = {
        "chrom" : chrom,
        "x" : mis_x,
        "window_start" : windows_m[...,0],
        "window_stop" : windows_m[...,1],
        "stat_type" : "mis",
        "stat" : mis
        }
    windows_df = pd.DataFrame(windows_d)
    
    theta_df = pd.DataFrame(theta_dict)
    pi_df = pd.DataFrame(pi_dict)
    d_df = pd.DataFrame(d_dict)
    mis_df = pd.DataFrame(mis_dict)
    chrom_df = pd.concat([theta_df, pi_df, d_df, mis_df], ignore_index = True, sort = False)
    return(chrom_df, windows_df)


def calculate_whole_genome_pi(zarr_file, missing_threshold, chrom_lengths):
    """
    Calculate pi across all chromosome
    
    zarr_file: open zarr file (uses helper.get_chrom_data)
    missing_threshold: threshold for missing genotypes per varint 0 - 1 ex) .80 
    chrom_lengths: dictionary containing start and stop values for each chromosome
    """
    
    #add sample list in the future 
    #
    chroms = ["I", "II", "III", "IV", "V", "X"]
    all_chroms = []
    for chrom in chroms:
        print("working on " + chrom)
        #Pull the zar data
        print("pulling data")
        gt , ac, pos = helpers.get_chrom_data(zarr_file, chrom, 
                                        snps = True, 
                                        numalt = 1,
                                        sample_list = None) #Filter for samples in sample list
        #Calculate Mask
        print("calculating mask")
        mask = diversity.mask_missing(gt , pos, missing_threshold, chrom, chrom_length)

        #Pi 
        print("calculating " + chrom + " pi")
        pi = allel.sequence_diversity(pos, 
                                    ac, 
                                    start = 1, 
                                    stop = chrom_lengths[chrom]['stop'])

        pi_masked = allel.sequence_diversity(pos, 
                                           ac, 
                                           start = 1, 
                                           stop = chrom_lengths[chrom]['stop'], 
                                           is_accessible = mask)
        chrom_data = [chrom, pi, pi_masked]
        all_chroms.append(chrom_data) #Add to all data list
    return(all_chroms)


def calculate_whole_genome_theta(zarr_file, missing_threshold, chrom_lengths):
    """
    Calculate pi across all chromosome
    
    zarr_file: open zarr file (uses helper.get_chrom_data)
    missing_threshold: threshold for missing genotypes per varint 0 - 1 ex) .80 
    chrom_lengths: dictionary containing start and stop values for each chromosome
    """
    
    #add sample list in the future 
    #
    chroms = ["I", "II", "III", "IV", "V", "X"]
    all_chroms = []
    for chrom in chroms:
        print("working on " + chrom)
        #Pull the zar data
        print("pulling data")
        gt , ac, pos = helpers.get_chrom_data(zarr_file, chrom, 
                                        snps = True, 
                                        numalt = 1,
                                        sample_list = None) #Filter for samples in sample list
        #Calculate Mask
        print("calculating mask")
        mask = diversity.mask_missing(gt , pos, missing_threshold, chrom, chrom_length)

        #Pi 
        print("calculating " + chrom + " pi")
        theta = allel.watterson_theta(pos,
                                   ac, 
                                   start = 1, 
                                   stop = chrom_length[chrom]['stop'])

        theta_masked = allel.watterson_theta(pos, 
                                           ac, 
                                           start = 1, 
                                           stop = chrom_length[chrom]['stop'], 
                                           is_accessible = mask)
        chrom_data = [chrom, theta, theta_masked]
        all_chroms.append(chrom_data) #Add to all data list
    return(all_chroms)

#Masking Functions


def create_missing_mask(gt_array, pos, threshold, chrom, genome_dict, start=None, stop=None): 
    """
    Create boolean array if postion passes missing GT threshold. First checks all
    variants if they pass a missing gt threshold. Array contains a value for each 
    position.

    threshold: proportion of samples missing GT information
    returns: bool vector if position passing missing thereshold
    """
    #Get array stats
    pre_n_snps = gt_array.shape[0]
    n_samps = gt_array.shape[1]
    max_missing = round(n_samps * threshold)
    
    #Count missing Gts for each variant
    n_missing_gts = gt_array.count_missing(axis=1) 

    #Make vector of variants that pass filtering threshold
    below_missing_threshold = n_missing_gts < max_missing 
    post_n_snps = np.count_nonzero(below_missing_threshold)

    removed_snps = pre_n_snps - post_n_snps
    print("There were", removed_snps, "variants masked")
    
    #Create an array of all positons from the chromosome length
    if start is not None or stop is not None:
        print("Using usr defined start and stops")
        pos_all = np.arange(start = start , stop = stop +1)
    else:
        start = 1 
        stop = genome_dict[chrom]['full']['Stop']
        pos_all = np.arange(start =1 , stop = stop +1)
   
    
    #Check if the position is a variant
    var_stat = np.in1d(pos_all, pos)
    
    #Check bool array of acessibility 
    # Very inefficent
    
    is_accessible = np.ones_like(pos_all, dtype = bool)
    is_accessible[var_stat] = below_missing_threshold
    
    
    return is_accessible

def create_gt_mask(status_bool, n_samples): 
    """
    Status_bool: vector containing True for variants that pass a filter 

    Returns: a mask for gts based on variant stats bool. Logic is inverted from inital vector
    True - GTs are masked 
    False - GTs are not masked
    """
    mask = []
    for i in range(0, len(status_bool)):
        status = status_bool[i]
        #Inverted logic, mask hides true gts
        if status == True:
            gt_masks = [True] * n_samples
        if status == False:
            gt_masks = [False] * n_samples
        mask.append(gt_masks)
    return(mask)
    
def get_passing_var(pos_array, pos_bool, chrom):
    """
    Writes out a Regions.tsv file (CHROM    POS) containing only 
    variants that are true in the pos_bool
    
    pos_array - sorted index from scikit-allel
    pos_bool - a bool array for all positions on the chromosome
    chrom - "chrom_id"
    returns: np.array 
    """

    def check_var(pos): 
        var_status_index = pos -1
        status = pos_bool[var_status_index]
        if status == True: 
            out =(chrom, int(pos))
            return(out)
        else: 
            out = (None, None)
            return(out)
    all_status = np.array(tuple(map(check_var, pos_array)))
    return(all_status)
