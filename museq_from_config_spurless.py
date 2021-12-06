'''
@author: Dan Levy, Cold Spring Harbor Laboratory
'''

import numpy as np
import os
import time
from support_functions import load_configuration_file, tabprint
from museq_functions import unzip_fastq_pair, get_mutation_rate
from museq_functions import execute_jellyfish
from museq_functions import kmers_to_contigs
from museq_functions import kmers_to_contigs_spurless
from museq_functions import mumdex_map_to_contig, annotate_maps
from museq_functions import extend_contigs_by_align
from museq_functions import unmask_contigs
from museq_functions import align_and_phase_contigs
from museq_functions import resolve_haplotypes
from museq_functions import load_contigs, write_contigs
from museq_functions import write_template_info
import sys


## list of file handles that are common through the code.
CONTIGS_FA             = "contigs.fa"
CONTIG_MAPS            = "contig_maps.txt"
CONTIG_MAPS_PLUS_NPY   = "contig_maps_plus.npy"
CONTIG_MAPS_PLUS_TXT   = "contig_maps_plus.txt"
EXTENDED_CONTIGS_FA    = "extended_contigs.fa"

RAW_TEMPLATES_FA       = "raw_templates.fa"
TEMPLATES_FA           = "templates.fa"
TEMPLATES_UNMASKED_FA  = "templates_unmasked.fa"
TEMPLATE_MASK_FA       = "template_mask.fa"
TEMPLATE_MAPS_TXT      = "template_maps.txt"
TEMPLATE_MAPS_PLUS_NPY = "template_maps_plus.npy"
TEMPLATE_MAPS_PLUS_TXT = "template_maps_plus.txt"
UNMASK_INFO_TXT        = "unmask_info.txt"
TEMPLATE_INFO_TXT      = "template_info.txt"

HAPLOTYPE_SHORTNAME_PATTERN = "templates_unmasked_h%d.fa"

## we initialize the mutation rate to None.
## if it is not computed during unmasking,
## we need to initialize it from the unmasking info
MUT_RATE = None

## check if a file exists
def file_exists(directory, filename):
    return os.path.exists(os.path.join(directory, filename))

TESTING = True
## load command line parameters
if TESTING:
    conf_filename         = '/data/safe/levy/museq/museq.conf'
    BS_0_raw_directory    = "/mnt/wigstore3/user/sili/Muchr18/chr18_unmutated"
    BS_1_raw_directory    = "/mnt/wigstore3/user/sili/Muchr18/chr18pos22"
    parent_data_directory = "/data/safe/levy/museq_pairlibs"
    short_name            = "chr18_pos22"
else:
    conf_filename         = sys.argv[1]
    BS_0_raw_directory    = sys.argv[2]
    BS_1_raw_directory    = sys.argv[3]
    parent_data_directory = sys.argv[4]
    short_name            = sys.argv[5]


## load the configuration file
keywords = load_configuration_file(conf_filename)
## override keywords
# keywords["MIN_DB_COVER"] = 20

keys = np.sort(keywords.keys())
for key in keys:
    print tabprint([key, keywords[key]])


INTEGER_NAMES = True

 
pair_dir       = os.path.join(parent_data_directory, short_name)
bs_0_dir       = os.path.join(pair_dir, "bs_0")
bs_1_dir       = os.path.join(pair_dir, "bs_1")
plot_directory = os.path.join(pair_dir, "plots")

print 
print 
print "********"
print " STEP 1 "
print "********"

print "unzipping the reads and assigning integer read names"
## check if data exist already
bs0r1 = file_exists(bs_0_dir, "r1.fastq")
bs0r2 = file_exists(bs_0_dir, "r2.fastq")
bs1r1 = file_exists(bs_1_dir, "r1.fastq")
bs1r2 = file_exists(bs_1_dir, "r2.fastq")

if bs0r1 and bs0r2 and bs1r1 and bs1r2:
    print "reads are already unpacked" 
else:
    print "unpacking read data"
    print "input data directories:"
    print "unmutated: %s" % BS_0_raw_directory
    print "mutated  : %s" % BS_1_raw_directory
    print
    print "output directories:"
    print "unmutated: %s" % bs_0_dir
    print "mutated  : %s" % bs_1_dir
    print
    sys.stdout.flush()
    t0 = time.time()
    ## for both the BS and no-BS libraries
    ## unzip the initial raw files
    unzip_fastq_pair(BS_0_raw_directory, bs_0_dir, INTEGER_NAMES=INTEGER_NAMES, VERBOSE=True)
    unzip_fastq_pair(BS_1_raw_directory, bs_1_dir, INTEGER_NAMES=INTEGER_NAMES, VERBOSE=True)
    t1 = time.time()
    print "done in %0.2f seconds" % (t1 - t0)
    sys.stdout.flush()


print
print "********"
print " STEP 2 "
print "********"
print "counting kmers using jellyfish"
sys.stdout.flush()
JELLYFISH         = keywords["JELLYFISH"]
JELLYFISH_MEM     = keywords["JELLYFISH_MEM"]
JELLYFISH_THREADS = keywords["JELLYFISH_THREADS"]
KMER_SIZE         = keywords["KMER_SIZE"]
MIN_MER_COUNT     = keywords["MIN_MER_COUNT"]
## check if data already exist
bs0_jf  = file_exists(bs_0_dir, "mer_counts.jf") 
bs0_txt = file_exists(bs_0_dir, "kmer_list.txt")
bs1_jf  = file_exists(bs_1_dir, "mer_counts.jf") 
bs1_txt = file_exists(bs_1_dir, "kmer_list.txt")

if bs0_jf and bs0_txt and bs1_jf and bs1_txt:
    print "already counted the kmers"
    sys.stdout.flush()
else:                        
    t0 = time.time()
    print "kmer counts for unmutated data"
    sys.stdout.flush()
    execute_jellyfish(bs_0_dir, KMER_SIZE, MIN_MER_COUNT, JELLYFISH, JELLYFISH_MEM, JELLYFISH_THREADS, VERBOSE=True)    
    print "kmer counts for mutated data"
    sys.stdout.flush()
    execute_jellyfish(bs_1_dir, KMER_SIZE, MIN_MER_COUNT, JELLYFISH, JELLYFISH_MEM, JELLYFISH_THREADS, VERBOSE=True)    
    t1 = time.time()
    print "done in %0.2f seconds" % (t1 - t0)
    sys.stdout.flush()



print
print "********"
print " STEP 3 "
print "********"
print "For the mutated data, build spur-trimmed contigs"
sys.stdout.flush()
MIN_DB_COVER = keywords['MIN_DB_COVER']
SPUR_TO_ALT_MAX_RATIO= 0.2

if file_exists(bs_1_dir, CONTIGS_FA):
    print "already built contigs"
    sys.stdout.flush()
else:
    t0 = time.time()
    kmers_to_contigs_spurless(bs_1_dir, CONTIGS_FA, MIN_DB_COVER, SPUR_TO_ALT_MAX_RATIO, VERBOSE=True)
    t1 = time.time()
    print "done in %0.2f seconds" % (t1 - t0)
    sys.stdout.flush()



print
print "********"
print " STEP 4 "
print "********"
print "Map mutated reads to the mutated contigs"
sys.stdout.flush()
MUMDEX_FASTQ   = keywords["MUMDEX_FASTQ"]
MIN_MUM_LENGTH = keywords["MIN_MUM_LENGTH"]
if file_exists(bs_1_dir, CONTIG_MAPS):
    print "already mapped contigs"
    sys.stdout.flush()
else:
    t0 = time.time()
    mumdex_map_to_contig(bs_1_dir, CONTIGS_FA, CONTIG_MAPS, MIN_MUM_LENGTH, MUMDEX_FASTQ, VERBOSE=True)
    t1 = time.time()
    print "done in %0.2f seconds" % (t1 - t0)
    sys.stdout.flush()


print
print "********"
print " STEP 5 "
print "********"
print "Annotate the maps with additional information alignment quality beyond the exact match."
print "Saving a pickle of the output as well."
sys.stdout.flush()
ann_txt = file_exists(bs_1_dir, CONTIG_MAPS_PLUS_TXT)
ann_npy = file_exists(bs_1_dir, CONTIG_MAPS_PLUS_NPY)
if ann_txt and ann_npy:
    print "already annotated the maps"
    sys.stdout.flush()
else:
    t0 = time.time()
    annotate_maps(bs_1_dir, CONTIGS_FA, CONTIG_MAPS, PICKLE=True, VERBOSE=True)
    t1 = time.time()
    print "done in %0.2f seconds" % (t1 - t0)
    sys.stdout.flush()
    

print
print "********"
print " STEP 6 "
print "********"
print "Extend contigs with annotated maps."
print "Extensions must be confident and unambiguous." 
print
sys.stdout.flush()
MIN_MAP_RATIO  = keywords["MIN_MAP_RATIO"]
MIN_MAP_LENGTH = keywords["MIN_MAP_LENGTH"]
EXTENT         = keywords["EXTENT"]
MIN_COVER      = keywords["MIN_COVER"]
MIN_RATIO      = keywords["MIN_RATIO"]

if file_exists(bs_1_dir, EXTENDED_CONTIGS_FA):
    print "already extended contigs"
    sys.stdout.flush()
else:
    t0 = time.time()
    extend_contigs_by_align(bs_1_dir, CONTIGS_FA, CONTIG_MAPS_PLUS_NPY, EXTENDED_CONTIGS_FA,
                            KMER_SIZE, MIN_MAP_RATIO, MIN_MAP_LENGTH,
                            EXTENT, MIN_COVER, MIN_RATIO, UNIQUE_ONLY=True, VERBOSE=True)
    t1 = time.time()
    delta_min = (t1 - t0) / 60
    print "done extension in %0.2f minutes" % delta_min
    sys.stdout.flush()

print
print "*********"
print " STEP 6B "
print "*********"
print "Restrict to just the templates we are interested in"
print "Those that are long enough to make something of."
print "Extensions must be confident and unambiguous." 
MIN_ASSEMBLY_LENGTH = keywords['MIN_ASSEMBLY_LENGTH']
## write out just the extended contigs that go into the assembly (better than MIN_ASSEMBLY LENGTH)
contigs, cnames = load_contigs(bs_1_dir, EXTENDED_CONTIGS_FA, return_names=True)
contig_lens = np.array( [len(x) for x in contigs])
max_index   = np.max(np.where(contig_lens >= MIN_ASSEMBLY_LENGTH)[0]) + 1
write_contigs(contigs[:max_index], bs_1_dir, RAW_TEMPLATES_FA, contig_names = cnames[:max_index])

## everything forward will work with the templates, not extended contigs!
## we will change the namespace and then pull out all the common terms as constants
print
print "********"
print " STEP 7 "
print "********"
print "Unmask mutations in extended contigs using unmutated data."
print "Also approximates the mutation rate."
sys.stdout.flush()
MIN_KMER      = keywords['MIN_KMER']
MIN_MAX_RATIO = keywords['MIN_MAX_RATIO']
MUT_RATE_MIN  = keywords['MUT_RATE_MIN']
contigs_unmasked = file_exists(bs_1_dir, TEMPLATES_UNMASKED_FA)
contigs_mask     = file_exists(bs_1_dir, TEMPLATE_MASK_FA)
contigs_info     = file_exists(bs_1_dir, UNMASK_INFO_TXT)
if contigs_unmasked and contigs_mask and contigs_info:
    print "already unmasked extended contigs."
    sys.stdout.flush()
else:
    t0 = time.time()
    MUT_RATE = unmask_contigs(pair_dir, RAW_TEMPLATES_FA, TEMPLATES_FA, TEMPLATES_UNMASKED_FA, TEMPLATE_MASK_FA, UNMASK_INFO_TXT, MIN_KMER, KMER_SIZE, MIN_MAX_RATIO, MUT_RATE_MIN=MUT_RATE_MIN, VERBOSE=True)
    t1 = time.time()
    delta_min = (t1 - t0) / 60
    print "done unmasking in %0.2f minutes" % delta_min
    sys.stdout.flush()


## if we didn't already load the mutation rate, let's do it now
if MUT_RATE == None:
    unmask_info_filename = os.path.join(bs_1_dir, UNMASK_INFO_TXT)
    MUT_RATE = get_mutation_rate(unmask_info_filename)

print
print "********"
print " STEP 8 "
print "********"
print "Align unmasked contigs and phase into two haplotypes."
sys.stdout.flush()
## we take all corrected sequences longer than some minimum length
MIN_ASSEMBLY_LENGTH = keywords["MIN_ASSEMBLY_LENGTH"]
## maximum ratio for considering a possible haplotype polymorphism
BASE_RATIO_CUTOFF = keywords["BASE_RATIO_CUTOFF"]
## the most polymorphisms we are willing to consider at a time (because of exhaustive algorithm)
MAX_POLY_COUNT    = keywords["MAX_POLY_COUNT"]
## for the haplotype splitting
ERR_RATE = keywords["ERR_RATE"] ## error rate
NUMBER_OF_HAPLOTYPES = keywords["NUMBER_OF_HAPLOTYPES"]
## how much h1 better than h2 to get counted as h1 and visa versa
MIN_HAPLO_CONFIDENCE = keywords["MIN_HAPLO_CONFIDENCE"]
MAX_SA_POLY_COUNT    = keywords["MAX_SA_POLY_COUNT"]
ALPHA                = keywords["ALPHA"]
MAX_STUCK            = keywords["MAX_STUCK"]   
MAX_TEMP_STEPS       = keywords["MAX_TEMP_STEPS"]  
## save plots to file?
SAVE_PLOTS          = keywords["SAVE_PLOTS"] == 1
USE_SA              = keywords["USE_SA"] == 1
TRIM_ENDS           = keywords["TRIM_ENDS"]
MOVES_PER_TEMP      = keywords["MOVES_PER_TEMP"]
VERBOSE             = True

if file_exists(bs_1_dir, HAPLOTYPE_SHORTNAME_PATTERN % 1):
    print "already aligned and phased into haplotypes"
    sys.stdout.flush()
else:
    t0 = time.time()
    align_and_phase_contigs(bs_1_dir, TEMPLATES_UNMASKED_FA, 
                                MIN_ASSEMBLY_LENGTH,
                                BASE_RATIO_CUTOFF,
                                MAX_POLY_COUNT,
                                MAX_SA_POLY_COUNT,
                                ERR_RATE,
                                MUT_RATE,
                                MIN_HAPLO_CONFIDENCE,
                                USE_SA=USE_SA,
                                TRIM_ENDS = TRIM_ENDS,
                                MAX_TEMP_STEPS = MAX_TEMP_STEPS,
                                MOVES_PER_TEMP = MOVES_PER_TEMP,
                                ALPHA = ALPHA,
                                MAX_STUCK = MAX_STUCK,
                                NUMBER_OF_HAPLOTYPES = NUMBER_OF_HAPLOTYPES,
                                SAVE_PLOTS=SAVE_PLOTS,
                                VERBOSE=VERBOSE,
                                plot_directory=plot_directory,
                                seq1_ind=0)
    t1 = time.time()
    delta_min = (t1 - t0) / 60
    print "done align and haplo-split in %0.2f minutes" % delta_min
    sys.stdout.flush()

print
print "********"
print " STEP 9 "
print "********"
print "Take consensus over haplotypes" 
sys.stdout.flush()
if file_exists(pair_dir, "haplotype_assemblies.fastq"):
    print "already done."
else:
    resolve_haplotypes(short_name, pair_dir, bs_1_dir, 
                       TEMPLATE_MASK_FA, HAPLOTYPE_SHORTNAME_PATTERN, 
                       plot_directory, ERR_RATE, MUT_RATE, 
                       NUMBER_OF_HAPLOTYPES=NUMBER_OF_HAPLOTYPES,
                       SAVE_PLOTS=True, VERBOSE=True, MAX_PHRED = 91)

print
print "*********"
print " STEP 10 "
print "*********"
print "Using only the extended contigs greater than %d bases" % MIN_ASSEMBLY_LENGTH
print "Remap all reads" 
sys.stdout.flush()

MUMDEX_FASTQ   = keywords["MUMDEX_FASTQ"]
MIN_MUM_LENGTH = keywords["MIN_MUM_LENGTH"]
if file_exists(bs_1_dir, TEMPLATE_MAPS_TXT):
    print "already mapped contigs"
    sys.stdout.flush()
else:
    t0 = time.time()
    mumdex_map_to_contig(bs_1_dir, TEMPLATES_FA, TEMPLATE_MAPS_TXT, MIN_MUM_LENGTH, MUMDEX_FASTQ, VERBOSE=True)
    t1 = time.time()
    print "done in %0.2f seconds" % (t1 - t0)
    sys.stdout.flush()

print
print "*********"
print " STEP 11 "
print "*********"
print "Annotate the template maps with additional information alignment quality beyond the exact match."
print "Saving a pickle of the output as well."
sys.stdout.flush()
ex_ann_txt = file_exists(bs_1_dir, TEMPLATE_MAPS_PLUS_TXT)
ex_ann_npy = file_exists(bs_1_dir, TEMPLATE_MAPS_PLUS_NPY)
if ex_ann_txt and ex_ann_npy:
    print "already annotated the extended maps"
    sys.stdout.flush()
else:
    t0 = time.time()
    annotate_maps(bs_1_dir, TEMPLATES_FA, TEMPLATE_MAPS_TXT, PICKLE=True, VERBOSE=True)
    t1 = time.time()
    print "done in %0.2f seconds" % (t1 - t0)
    sys.stdout.flush()

print "*********"
print " STEP 12 "
print "*********"
print "Generate output file and plots for the templates."
tinfo_txt = file_exists(bs_1_dir, TEMPLATE_INFO_TXT)
if tinfo_txt:
    print "already summarized the output..."
else:
    t0 = time.time()
    write_template_info(pair_dir, bs_1_dir, TEMPLATES_FA, TEMPLATE_MAPS_PLUS_NPY, HAPLOTYPE_SHORTNAME_PATTERN, UNMASK_INFO_TXT, TEMPLATE_INFO_TXT,
                        NUMBER_OF_HAPLOTYPES, MIN_MAP_RATIO, MIN_MAP_LENGTH,
                            VERBOSE=VERBOSE, SAVE_PLOTS=SAVE_PLOTS,plot_directory=plot_directory, plot_sub_directory="coverage")
    t1 = time.time()
    print "done in %0.2f seconds" % (t1 - t0) 
