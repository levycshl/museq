'''
Created on Jun 4, 2014

@author: levy
'''

import numpy as np
from itertools import imap, product
import itertools
import operator
import os
import gzip
import shutil
from collections import defaultdict, Counter
from support_functions import tabprint, path_summary, file2RA
import matplotlib.pyplot as plt
import AlignmentExt
import sys
import matplotlib as mpl
import time

## useful function for unpacking a list of lists into a single list
flatten = lambda l: [item for sublist in l for item in sublist]

BASES = ['A', 'C', 'G', 'T', 'N']
## dictionary used to invert read sequences
INDICT = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N',
          'a':'t', 'c':'g', 'g':'c', 't':'a', 'n':'n'}

BASE2INT = dict([x[::-1] for x in enumerate(BASES)])


# BASES_COLOR = ["blue", "yellow", "green", "red"]
## headache colorscheme from colourlover
DEFAULT_COLOR_LIST = ["#333333",  
                      '#80BCA3', 
                      '#655643',
                      '#E6AC27', 
                      '#BF4D28',
                      '#F6F7BD']

## alternative colors, closer match to IGV
ALT_COLOR_LIST = ["#333333",   
                  'forestgreen',
                  'dodgerblue', 
                  'gold', 
                  'crimson',
                  'gray']


def flip(i):
    return (i+1)%2


def contig_forward_path(cgraph, node):
    path_forward = [node]
    while True:
        last_element = path_forward[-1]
        candidates = cgraph[last_element[0]][flip(last_element[1])]
        ## if there's no choice or more than once choice, get out!
        if len(candidates) != 1:
            break
        ## otherwise, get your choice
        candidate      = candidates[0]
        candidate_flip = (candidate[0], flip(candidate[1]))
        ## for now, they make me feel safe and warm and fuzzy 
        if candidate in path_forward:
            break  ## get out of loops
        if candidate_flip in path_forward:
            break  ## watch out for mobius contigs and hairpins
        ## check there is only one path into candidate...
        into_candidate = cgraph[candidate[0]][candidate[1]]
        if len(into_candidate) != 1:
            break
        path_forward.append(candidate)
    return path_forward

def get_path(cgraph, db_contigs, ind, K):
    forward_path  = contig_forward_path(cgraph, (ind, 0))
    backward_path = contig_forward_path(cgraph, (ind, 1))
    forward_end = forward_path[-1]
    ## if it is a cycle starting at the index:
    if (ind, 0) in cgraph[forward_end[0]][forward_end[1]]:
        total_path = forward_path
    else:
        total_path = [(x[0], flip(x[1])) for x in backward_path[-1:0:-1]] + forward_path
    ind0, rc0 = total_path[0]
    first_contig = db_contigs[ind0] if rc0 == 0 else reverse_complement(db_contigs[ind0])
    contig_list = [first_contig]
    for next_ind, next_rc in total_path[1:]:
        ## get the next contig and lop off the first K-1 bases
        next_contig = db_contigs[next_ind][(K-1):] if next_rc == 0 else reverse_complement(db_contigs[next_ind])[(K-1):]
        contig_list.append(next_contig)
    total_contig = "".join(contig_list)
    return total_contig, total_path

def simplify_cgraph(cgraph, db_contigs, K):
    done = set()
    r = []
    for node in np.arange(len(db_contigs)):
        if node not in done:
            new_string, new_path = get_path(cgraph, db_contigs, node, K)
            for cnode in new_path:
                done.add(cnode[0])
            r.append(new_string)
    contigs  = np.array(r) 
    seq_lens = np.array([len(x) for x in contigs])
    order    = np.argsort(seq_lens)[::-1]
    ans = contigs[order]
    return ans


'''
Given a working directory, convert the kmers into spurless contigs
and write the resulting contigs to file
'''
def kmers_to_contigs_spurless(outdir, contig_shortname, MIN_COUNT, SPUR_TO_ALT_MAX_RATIO=0.2, VERBOSE=False ):
    kmer_filename = os.path.join(outdir, "kmer_list.txt")
    contig_filename = os.path.join(outdir, contig_shortname)
    if VERBOSE:
        print "loading kmers"
    '''PART I: get the kmers'''
    ## load only kmers with count >= MIN_COUNT
    db_graph = kmer_file_to_db_graph(kmer_filename, MIN_COUNT=MIN_COUNT) 
    K = len( db_graph.keys()[0] )
    ## iterate through the db_graph until we get db_contigs
    '''PART II: get the debruijn contigs'''
    done = set()
    r = []
    for node in db_graph:
        if node not in done:
            contig_string, contig = get_contig(db_graph, node)
            for cnode in contig:
                done.add(cnode)
                done.add(reverse_complement(cnode))
            r.append(contig_string)        
    db_contigs  = np.array(r)
    contig_lens = np.array( [len(x) for x in db_contigs] )
    ## now, we've got all the things we need to go spurless, short some parameters
    ## first the contig graph built from the db_contigs..
    '''PART III: build the contig graph, with edge counts'''
    cgraph = contig_graph(db_contigs, K)
    ## then the edge coverage
    edge_cover = []
    for contig in db_contigs:
        edge_cover.append( [db_graph[contig[:K]], db_graph[contig[-K:]]] )
    edge_cover = np.array( edge_cover )
    '''PART IV: pick out possible spurs '''
    ## shorter than 2K -1
    spur_length_filter = contig_lens < (2*K -1)
    ## with exactly one link
    neighbor_count       = np.array( [[len(side) for side in neighbor_pair] for neighbor_pair in cgraph ] )
    neighbor_count_total = np.sum(neighbor_count, axis = 1) 
    spur_one_link_filter = neighbor_count_total == 1     
    spur_candidates = spur_length_filter * spur_one_link_filter
    possible_spur_list = np.where(spur_candidates)[0]
    if VERBOSE:
        print "%d possible spurs" % len(possible_spur_list)
    '''PART V: gather statistics on spurs''' 
    SOURCE = 0
    ALT    = 1
    SPUR   = 2
    all_triples = []
    for ind in possible_spur_list:
        ## the source is the only neighbor to this contig
        source = (cgraph [ind][0] + cgraph [ind][1] )[0]
        source_neighbors = cgraph [source[0]][source[1]]
        source_cover = edge_cover[source]
        ## the triple scores coverage as [SOURCE, ALT, SPUR]
        triple = [source_cover, 0, 0]
        for target in source_neighbors:
            if target[0] == ind:
                triple[2] = edge_cover[target]            
            else:
                if edge_cover[target] > triple[1]:
                    triple[1] = edge_cover[target]
        all_triples.append(triple)
    all_triples  = np.array(all_triples )
    '''PART VI: remove spurs with low ratios compared to the alternative'''
    ratio =  all_triples[:, SPUR] / (all_triples[:, ALT] ).astype(float)
    ratio_filter = ratio < SPUR_TO_ALT_MAX_RATIO
    ## WE USE THE RATIO OF SPUR TO ALT TO FILTER SPURS
    spur_list = possible_spur_list[ratio_filter]
    if VERBOSE:
        x = len(spur_list)
        y = len(db_contigs)
        print "%d spurs identified of %d contigs: %0.3f%% of total" %(x, y , x/float(y))
    spur_filter = np.ones(len(db_contigs), dtype=bool)
    spur_filter[spur_list] = False
    ## MAKE A NEW CONTIG GRAPH THAT EXCLUDES SPURS
    contigs_spurless = [contig for (ind, contig) in enumerate( db_contigs  ) if spur_filter[ind]]
    db_graph_spurless = contig_graph(contigs_spurless, K)
    '''PART VII: eliminate spurs from the spurless graph'''
    simple_contigs = simplify_cgraph(db_graph_spurless, contigs_spurless, K)
    '''PART VIII: write to file '''
    if VERBOSE:    
        print "writing %d contigs to file" % len(simple_contigs)
    outfile     = file(contig_filename, 'w')
    for ind, contig in enumerate( simple_contigs ):
        outfile.write(">%d\n" % ind) 
        outfile.write("%s\n" % contig)
    outfile.close()
    if VERBOSE:
        print "done"


## count how many reads there are
def read_pair_counter(read_dir):
    read1_filename = os.path.join(read_dir, 'r1.fastq')
    infile = file(read1_filename, 'r')
    counter = 0
    for line in infile:
        counter += 1
    infile.close()
    return (counter / 4)


'''
A handy function for iterative downsampling.
Samples uniformly from the first "length_of_input" elements to emit "number to get."
'''
def downsample(input_iterator, length_of_input, number_to_get):
    ## how many are left to get
    left_to_get = number_to_get
    for current_ind, entry in enumerate(input_iterator):
        draw = np.random.randint(0, length_of_input - current_ind)
        if draw < left_to_get:
            left_to_get -= 1
            yield entry
            if left_to_get == 0:
                break

'''
Simple function to iterate through a pair of r1, r2.fastq files.
Returns four lines at a time from each of the pair of files.
'''
def read_pair_iterator(read_dir):
    read1_filename = os.path.join(read_dir, "r1.fastq")
    read2_filename = os.path.join(read_dir, "r2.fastq")
    f1_in = file(read1_filename, 'r')
    f2_in = file(read2_filename, 'r')
    read_buffer1 = []
    read_buffer2 = []
    line_counter = 0
    for line1, line2 in zip(f1_in, f2_in):
        read_buffer1.append(line1)
        read_buffer2.append(line2)
        line_counter += 1
        if line_counter == 4:
            line_counter = 0
            yield (read_buffer1, read_buffer2)
            read_buffer1 = []
            read_buffer2 = []
    f1_in.close()
    f2_in.close()

'''
Takes all the read1 and read2 files from old dirs and merges them into a single read1 and read2 in new_dir
Makes new_dir if needed.
Will also keep an index file showing which reads loaded to which index block
'''
def merge_read_data(old_dirs, new_dir, VERBOSE=False):
    if VERBOSE:
        print "making directory", new_dir
    ## make output directory
    try:
        os.makedirs(new_dir)
    except:
        if VERBOSE:
            print new_dir, "already exists"
            
    ## make new data files
    read1_filename = os.path.join(new_dir, "r1.fastq")
    read2_filename = os.path.join(new_dir, "r2.fastq")
    index_filename = os.path.join(new_dir, "index_starts.txt")
    f1_out = file(read1_filename, 'w')
    f2_out = file(read2_filename, 'w')
    findex = file(index_filename, 'w')
    ## for each of the old directories, make a read pair iterator
    read_index = 0
    for old_dir in old_dirs:
        read_index_start = read_index
        if VERBOSE:
            print old_dir, read_index
        rpi      = read_pair_iterator(old_dir)
        ## iterate through rpi and update to reflect the new counting
        for entry in rpi:
            if VERBOSE and read_index % 100000 == 0:
                print read_index
            ## modify the names so that they are sequential
            new_name_line = "@%d\n" % read_index
            entry[0][0]   = new_name_line
            entry[1][0]   = new_name_line
            ## write the entries to the outfiles
            f1_out.write("".join(entry[0]))
            f2_out.write("".join(entry[1]))
            read_index += 1
        findex.write(tabprint([old_dir, read_index_start, read_index]) + "\n")
    f1_out.close()
    f2_out.close()
    findex.close()


'''
Given an input fastq directory with r1/r2.fastq, 
saves a downsampled version of r1/r2.fastq from old_dir into new_dir
'''
def downsample_fasta(old_dir, new_dir, downsample_count, VERBOSE=False):
    increment = downsample_count / 20
    if VERBOSE:
            print "making directory", new_dir
        ## make output directory
    try:
        os.makedirs(new_dir)
    except:
        if VERBOSE:
            print new_dir, "already exists"                
    ## get total count from old data by iteration
    total_read_count = 0
    rpi = read_pair_iterator(old_dir)    
    for _ in rpi:
        total_read_count += 1
    if VERBOSE:
        print "total read count =", total_read_count
    ## make new data files
    read1_filename = os.path.join(new_dir, "r1.fastq")
    read2_filename = os.path.join(new_dir, "r2.fastq")
    f1_out = file(read1_filename, 'w')
    f2_out = file(read2_filename, 'w')
    ## make a new iterator coupled with the downsampler
    rpi      = read_pair_iterator(old_dir)
    rpi_down = downsample(rpi, total_read_count, downsample_count)
    ## iterate through the downsample
    for ind, entry in enumerate( rpi_down ):
        if VERBOSE and ind % increment == 0:
            print ind
        ## modify the names so that they are sequential
        new_name_line = "@%d\n" % ind
        entry[0][0]   = new_name_line
        entry[1][0]   = new_name_line
        ## write the entries to the outfiles
        f1_out.write("".join(entry[0]))
        f2_out.write("".join(entry[1]))
    f1_out.close()
    f2_out.close()

'''
Given a set of paired-read alignments to a contig of the form:
0 = rindex
1 = cindex
2 = read1_itercept
3 = read1_slope
4 = read2_intercept
5 = read2_slope
Looks up the sequences and puts it into a list along with the map.
Reverses the reads as necessary if FLIP_FORWARD is set to True
Returns the sequences and maps
This is the input to taking the consensus sequence for extension
'''
def read_alignments_seq_and_map(paired_read_alignments, reads, FLIP_FORWARD=False): 
    seq_list     = []
    map_list     = []
    for entry in paired_read_alignments:
        rindex = entry[0]
        for r12, eind in ((0, 2), (1,4)):
            ## load the map and if it is None (no map for that read in the pair),do not add to list.
            rmap = entry[eind:(eind+2)]
            if rmap[0] ==None:
                continue
            else:
                rseq = reads[r12][rindex]
                if FLIP_FORWARD and rmap[1] == -1:
                    rmap = reverse_A(rmap[0], rmap[1], len(rseq))
                    rseq = reverse_complement(rseq)
                seq_list.append(rseq)
                map_list.append(rmap)
    return seq_list, map_list

'''
Given a set of contigs, a kmer size K, and an extent for db graph extensions,
makes a dictionary of mapping extensions between contigs.
Will travel EXTENT bases along DBgraph before stopping.
If a value is specified for contig_index_bound,
will only analyze that many of the contigs.
For this code, contigs are typically in length sorted order 
and this convenience is used to focus only on those contigs
that will require the use of an extension,
specifically, those that are long enough to admit a well-mapped read. 

DB_ext is keyed by (from contig, to contig) and value is S_XC, F_XC
DB_ext_key has for each from contig a list of to contigs
'''
def get_DB_extension(contigs, K, EXTENT, contig_index_bound =None, VERBOSE=False):
    if contig_index_bound == None:
        contig_index_bound = len(contigs)
    ## build the contig graph
    if VERBOSE:
        print "building dbGraph"
    graph = contig_graph(contigs, K)
    DB_ext     = dict()
    DB_ext_key = defaultdict(list)
    if VERBOSE:
        print "getting DB extensions"
    for contig_index in range(contig_index_bound):
        if VERBOSE:
            if contig_index % 1000 == 0:
                print contig_index, 'of', contig_index_bound
                sys.stdout.flush()
        S_XC, F_XC, in_frame = get_map_extension(contig_index, contigs, graph, K, EXTENT)
        for ind, (sxc, fxc, iframe) in enumerate(zip(S_XC, F_XC, in_frame)):
            if iframe:
                DB_ext[contig_index, ind] = (sxc, fxc)
                DB_ext_key[contig_index].append(ind)
    return DB_ext, DB_ext_key

# '''
# Function iterates through a list of maps that are presumed sorted in read name order.
# It then returns sublists with a common read name, along with the index of that read.
# '''
def read_map_iterator(all_maps):
    return block_iterator(all_maps, READ_NAME)

'''
Function iterates through a list of stuff.
Presumes that the list is ordered along "index to track."
It returns sublists of the list_of_stuff with the common index_to_track.
Also returns that index.
'''
def block_iterator(list_of_stuff, index_to_track):
    ans = None
    for element in list_of_stuff:
        ## if we are just getting started
        if ans == None:
            ## make an empty list
            ans = []
            ans.append(element)
            current_index = element[index_to_track]
        else:
            next_index = element[index_to_track]
            if next_index == current_index:
                ans.append(element)
            else:
                ## return present answer
                yield current_index, ans
                ## re initialize answer
                ans = []
                ans.append(element)
                current_index = next_index
    ## be sure to return the last answer when the list is done.
    yield current_index, ans

'''
Takes a list of good ratio maps and looks for alignments to contig cindex.
If there is a good ratio map to the contig in question, that is taken as the answer, no questions asked.
If there is no such map, we look for appropriate extensions in the DB graph.
If all extensions agree on the answer, that answer is returned.
If there is a disagreement, the answer is set to None.
'''
def map_by_extension(gr_maps, cindex, DB_ext):
    ext_maps   = [None , None ]
    ext_done   = [False, False]
    for rmap in gr_maps:
        ri = rmap[READ_12]
        ## if already done that index, continue
        if ext_done[ri]:
            continue
        else:
            cn = rmap[CONTIG_NAME]
            ## if the map is to the contig in question, record and be done
            if cn == cindex:
                ## only set to done if it is the correct contig
                ext_maps[ri] = (rmap[SRC], rmap[FRC])
                ext_done[ri] = True
                ## and if that finishes things, be completely done
                if ext_done[0] and ext_done[1]:
                    break
            else:
                ## check to see if there we can use the db_extension
                try:
                    ## mapped contig to target contig
                    s_xc, f_xc = DB_ext[cindex, cn]
                    s_rx       = rmap[SRC]               ## read to mapped contig
                    f_rx       = rmap[FRC]
                    s_rc, f_rc = compose_maps(s_rx, f_rx, s_xc, f_xc)   ## compose
                    ## if there is no prior map, use this one
                    if ext_maps[ri] == None:
                        ext_maps[ri] = (s_rc, f_rc)
                    else:
                        ## check if this disagrees with the prior map and if it does,
                        ## mark the read map as None and set to done
                        if ext_maps[ri] != (s_rc, f_rc):
                            ## fail, wipe out this map 
                            ext_maps[ri] = None
                            ## and don't bother looking here again
                            ext_done[ri] = True
                except:
                    ## cn is not in the extention of cindex
                    continue
    return ext_maps

'''
Given two map sets: top quality and good ratio maps, and a DB_extension
returns all well mapped alignments.
It will only map reads that are top quality maps to the contig and then use the good ratio maps and the DB ext
as necessary to complete the alignment.
Will only include a read as aligned if both reads in the pair align properly (knowable and unambiguous)
'''
def well_mapped_alignments(top_quality_maps, good_ratio_maps, DB_ext, VERBOSE=False, BOTH_MUST_ALIGN=True):
    alignments = []
    ## We are going to iterator over both the top quality maps (the only maps that really matter)
    ## and the good ratio maps, that back up the top quality maps when there is no clear alignment in the tq maps
    tq_iterator = read_map_iterator(top_quality_maps)
    gr_iterator = read_map_iterator(good_ratio_maps )
    gr_ind, gr_maps = gr_iterator.next()
    ind = 0
    for tq_ind, tq_maps in tq_iterator:
        if VERBOSE:
            if ind % 100000 == 0:
                print ind
                sys.stdout.flush()
        ind += 1
        while tq_ind != gr_ind:
            gr_ind, gr_maps = gr_iterator.next()
        ## list of contig indices that are possibily well mapped by this read pair
        well_mapped_contigs = np.unique([rmap[CONTIG_NAME] for rmap in tq_maps])
        for cindex in well_mapped_contigs:
            map1, map2 = map_by_extension(gr_maps, cindex, DB_ext)
            both_align = ( map1 != None ) and ( map2 != None )
            one_aligns = ( map1 != None ) or  ( map2 != None )
            if BOTH_MUST_ALIGN and both_align:
                info = [tq_ind, cindex]
                info.extend(list(map1))
                info.extend(list(map2))
                alignments.append(info)
            elif not BOTH_MUST_ALIGN and one_aligns:
                info = [tq_ind, cindex]
                info.extend(list(map1) if map1 != None else [None, None]) 
                info.extend(list(map2) if map2 != None else [None, None])
                alignments.append(info)
    return alignments


'''
For a given contig index (cind) and its sublist of alignments (caligns)
1st break apart read-pair alignments into single read sequence and position information.
2nd convert this to base counts
3rd determine the boundaries of proper extension given MIN_COVER and MIN_RATIO bounds
4th find compatible contigs in the DB extension (limit to extendable contigs)
5th return extended sequence and maps to other extendable contigs
'''
def extend_contig_from_align(cind, caligns, reads, contigs, DB_ext, DB_ext_key, MIN_COVER, MIN_RATIO, K, contig_index_bound=None, VERBOSE=True):
    if contig_index_bound == None:
        contig_index_bound = np.max(DB_ext_key.keys()) + 1
    rseqs, rmaps = read_alignments_seq_and_map(caligns, reads, FLIP_FORWARD=True)
    base_count, contig_start  = get_base_counts(rseqs, rmaps, ALL_FORWARD = True)
    contig_end = contig_start + len(contigs[cind])
    extended_seq, s_ce = get_extension(base_count, 
                                               contig_start, 
                                               contig_end, 
                                               MIN_COVER, 
                                               MIN_RATIO)
    extended_seq_rc   = reverse_complement(extended_seq)
    extended_seq_maps = []
    ## if the extension disagrees with its own contig, DO NOT EXTEND, return with no extension
    match, total = check_match(s_ce, 1, contigs[cind], extended_seq, extended_seq_rc)
    if match != total:
        if VERBOSE:
            print tabprint([cind, "there is a mismatch in the extension! do not extend!"])
        extended_seq = contigs[cind]
    else:
        ## use extension map to align extended seq to other contigs in its dbgraph extension
        for extend_to in DB_ext_key[cind]:
            if extend_to < contig_index_bound:
                s_xc, f_xc = DB_ext[cind, extend_to]
                ## compose with map between this contig and its extension
                s_xe, f_xe   = compose_maps(s_xc, f_xc, s_ce, 1)
                ## check the match
                match, total = check_match(s_xe, f_xe, contigs[extend_to], extended_seq, extended_seq_rc)
                ## if the match is perfect, AND NON-TRIVIAL, add to our list of extended maps
                if match == total >= K:
                    map_data = [extend_to, s_xe, f_xe, len(contigs[extend_to]), match]
                    extended_seq_maps.append(map_data)
    return extended_seq, extended_seq_maps


def extensions_and_maps_from_align(alignments, reads, contigs, DB_ext, DB_ext_key, MIN_COVER, MIN_RATIO, K, contig_index_bound=None, VERBOSE=False):
    if contig_index_bound == None:
        contig_index_bound = np.max(DB_ext_key.keys()) + 1
    ## initialize output over extendable sequences
    ## initialized in case cind is not visited by the block iterator (i.e. no good alignments)
    extension     = [contig for contig in contigs[:contig_index_bound]]
    extension_map = [[]     for _      in range  ( contig_index_bound)]    
    ## go through the different contigs and their alignments
    for cind, caligns in block_iterator(alignments, 1):
        extended_seq, extended_seq_maps = extend_contig_from_align(cind, caligns, reads, contigs, DB_ext, DB_ext_key, MIN_COVER, MIN_RATIO, K, contig_index_bound)
        extension    [cind] = extended_seq
        extension_map[cind] = extended_seq_maps
        if VERBOSE:
            rcount   = len(caligns)
            clen     = len(contigs[cind]) 
            elen     = len(extended_seq)
            extended = int(clen != elen)
            info   = [cind, extended, rcount, clen, elen]
            print tabprint(info)
            sys.stdout.flush()
    return extension, extension_map

'''
build lookup dictionaries from extensible contigs to extensions that are non-trivially compatible (K or more bp)
'''
def build_extension_lookup(extension_map):
    S_yX   = dict()
    F_yX   = dict()
    CINDEX      = 0
    SINDEX      = 1
    FINDEX      = 2
    for to_ind, emaps in enumerate(extension_map):
        ## if there are no maps
        if len(emaps) == 0:
            ## make a map to self
            S_yX  [to_ind, to_ind] = 0
            F_yX  [to_ind, to_ind] = 1
        for emap in emaps:
            from_ind = emap[CINDEX]
            S_yX  [from_ind, to_ind] = emap[SINDEX]
            F_yX  [from_ind, to_ind] = emap[FINDEX]
    return S_yX, F_yX

'''
Given a dictionary of mappings from extensible contigs into extensions,
determines all symmetric relationships (A extends to B and B extends to A) and clusters along those.
Meanwhile, also recording the mappings from each extended contig to the first element in its cluster.
This is important for joining clusters by sequence in the next step.
'''
def cluster_extensions_from_lookup(N, S_yX, F_yX):
    incidence_graph = [set() for _ in range(N)]
    ## we insert all the symmetric relations into the incidence graph
    ## that is we require that (a,b) and (b,a) are keys in S_yx
    for a, b in S_yX.keys():
        try: 
            S_yX[b, a]
            incidence_graph[a].add(b)
        except:
            continue    
    ## what we would like to do is a bit more complicated than just clustering
    ## We want to also store information about the mapping.
    S_XC = np.zeros(N, dtype=int)
    F_XC = np.ones(N, dtype=int)        
    ## initialize empty clusters
    clusters = []    
    ## copy the nodes into a set
    nodes = set(np.arange(N))
    ## iterate while there are still nodes in the set
    while nodes:
        ## get a random node and remove from the set
        node = nodes.pop()
        ## since this node is starting a component, it has a (0, 1) mapping
        ## which are the default values so we dont have to do anything here.
        ## establish a set for the next component
        component = { node }
        ## build a queue with just this node
        queue = [ node ]
        ## go through the queue until it is empty
        ## we assume that if you are in the queue, you already have a valid S_XC and F_XC
        while queue:
            ## pop the next item off the queue
            x = queue.pop(0)
            ## get its neighbors
            neighbors = incidence_graph[x] 
            ## forget about things we've already visited
            neighbors.difference_update(component)
            ## remove the remaining nodes from the global set
            nodes.difference_update(neighbors)
            ## add them to the current component
            for y in neighbors:
                ## get all the relevant maps
                s_yY, f_yY = S_yX[y, y], F_yX[y, y]  ## from y to Y (its extension)
                s_yX, f_yX = S_yX[y, x], F_yX[y, x]  ## from y to X (the extended node already in the cluster)
                s_XC, f_XC = S_XC[x], F_XC[x]        ## from X to C (since X is in the cluster, already labeled.
                ## and now some reversals and compositions.
                s_yC, f_yC = compose_maps(s_yX, f_yX, s_XC, f_XC)
                s_Yy, f_Yy = invert_map(s_yY, f_yY)
                s_YC, f_YC = compose_maps(s_Yy, f_Yy, s_yC, f_yC)
                ## and store the answer into the S_XC and F_XC arrays
                S_XC[y] = s_YC
                F_XC[y] = f_YC
            component.update(neighbors)
            ## and add them to the queue (so we visit them)
            queue.extend(neighbors)
        clusters.append(component)
    return clusters, S_XC, F_XC



'''
Given a set of contigs and reads (already mapped and annotated)
1. extend the contigs according to well-mapped reads and their mates,
    using the dbGraph to map mates in-frame
2. join contigs that extend perfectly into each other
3. drop the redundant subsequences

This is very similar to extend_contigs -- but probably much much faster.
This aligns all the reads to the relevant contigs first using the dbGraph.

FUTURE WORK: DB extension and alignments should be handled outside of this function.
But for now, for ease of plugging into existing pipeline, this is the structure.


Parameters are:
outdir            : working directory
contig_shortname  : shortname for contig file
extended_shortname: shortname for output file of extended contigs
K                 : kmer size
MIN_MAP_RATIO     : percent bases agreement with contig to count a map
MIN_MAP_LENGTH    : minimum length of a read to count as well-mapped
EXTENT            : how far to go in the dbGraph for mapping mates (~max value for insert size)
MIN_COVER         : minimum coverage to call an extension
MIN_RATIO         : minimum ratio to call an extension
UNIQUE_ONLY       : if True, remove proper subsequences from final extension (default=True)
'''
def extend_contigs_by_align(outdir, contig_shortname, maps_pickle, extended_shortname,
                            K, MIN_MAP_RATIO, MIN_MAP_LENGTH,
                            EXTENT, MIN_COVER, MIN_RATIO, UNIQUE_ONLY=True, VERBOSE=False):
    if VERBOSE:
        print "loading contigs"
        sys.stdout.flush()
    contigs = load_contigs  (outdir, contig_shortname)
    if VERBOSE:
        print "loading read pairs"
        sys.stdout.flush()
    reads   = load_readpairs(outdir)
    if VERBOSE:
        print "loading maps"
        sys.stdout.flush()
    maps    = load_maps     (outdir, maps_pickle)
    if VERBOSE:
        print "building filters"
        sys.stdout.flush()
    ## establish filters
    map_ratio      = maps[:, MATCH] / maps[:, OVERLAP].astype(float)
    ratio_filter   = map_ratio   >= MIN_MAP_RATIO
    length_filter  = maps[:, OVERLAP] >= MIN_MAP_LENGTH   
    quality_filter = ratio_filter * length_filter
    if VERBOSE:
        print "selecting maps"
        sys.stdout.flush()
    ## get best maps
    top_quality_maps = maps[quality_filter]
    good_ratio_maps  = maps[ratio_filter]
    ## find the biggest index we will have to think about for extensions
    ## this is determined by the MIN_MAP_LENGTH since anything shorter cannot have top quality maps
    contig_index_bound = np.sum([len(x) >= MIN_MAP_LENGTH for x in contigs])
    ## get DB extensions
    if VERBOSE:
        print "getting %d extensions" % contig_index_bound
        sys.stdout.flush()
    DB_ext, DB_ext_key = get_DB_extension(contigs, K, EXTENT, contig_index_bound=contig_index_bound, VERBOSE=VERBOSE)
    if VERBOSE:
        print "getting all alignments"
        sys.stdout.flush()
    alignments = well_mapped_alignments(top_quality_maps, good_ratio_maps, DB_ext, VERBOSE=VERBOSE)
    ## compute contig index and read index for each alignment
    cind = [align[1] for align in alignments]
    rind = [align[0] for align in alignments]
    ## sorting first by contig index and then by read index
    order = np.lexsort([rind, cind])
    ## then re-order alignments
    alignments = [alignments[ind] for ind in order]
    if VERBOSE:
        print "done with alignments, starting to extend"
        sys.stdout.flush()
    # alignments = np.array(alignments)
    ## the following is sure to have a single entry for each index up to contig_index_bound
    extension, extension_map = extensions_and_maps_from_align(alignments, reads, contigs, 
                                                              DB_ext, DB_ext_key, 
                                                              MIN_COVER, MIN_RATIO, K, 
                                                              contig_index_bound=contig_index_bound, 
                                                              VERBOSE=True)
    S_yX, F_yX           = build_extension_lookup(extension_map)
    clusters, S_XC, F_XC = cluster_extensions_from_lookup(contig_index_bound, S_yX, F_yX)    
    cluster_seqs = [cluster_to_sequence(list(cluster_set), S_XC, F_XC, extension) for cluster_set in clusters]
    ## remove matching subsequences
    if UNIQUE_ONLY:
        print 'removing proper sub-sequences'
        sys.stdout.flush()
        cluster_seqs = remove_proper_subsequences(cluster_seqs)
    ## order by length of contig
    order = np.argsort([len(x) for x in cluster_seqs])[::-1]
    cluster_seqs = [cluster_seqs[x] for x in order]
    if VERBOSE:
        print "from %d initial contigs" % len(contigs)
        print "writing %d extended contigs to file" % len(cluster_seqs)
        sys.stdout.flush()
    write_contigs(cluster_seqs, outdir, extended_shortname)


'''
converts a list of cluster indices into a consensus sequence from their extensions
'''
def cluster_to_sequence(cluster_list, S_XC, F_XC, extension):
    clust_seqs                     = [ extension[x]       for x in cluster_list ]
    clust_maps    = [ (S_XC[x], F_XC[x]) for x in cluster_list ] 
    base_count, _ = get_base_counts(clust_seqs, clust_maps, ALL_FORWARD = False)
    max_int       = np.argmax(base_count, axis=1)
    max_seq       = np.array(int_to_seq(max_int))
    return "".join(max_seq)


'''
Given an incidence array and mappings from contigs into extensions,
returns connected components from the incidence array.
Also, returns mappings from each extended contig to the first element in its cluster (S_XC, F_XC)
'''
def get_extension_clusters(incidence_array, S_yX, F_yX):
    N = len(incidence_array)
    ## we convert to a set-graph notation to facilitate the clustering algorithm
    incidence_graph = [set() for _ in range(N)]
    for x, y in zip(*np.where(incidence_array)):
        incidence_graph[x].add(y)
    
    ## what we would like to do is a bit more complicated than just clustering
    ## We want to also store information about the mapping.
    S_XC = np.zeros(N, dtype=int)
    F_XC = np.ones(N, dtype=int)
    
    clusters = []
    
    nodes = np.arange(N)    
    ## copy the nodes into a set
    nodes = set(nodes)
    ## iterate while there are still nodes in the set
    while nodes:
        ## get a random node and remove from the set
        node = nodes.pop()
        ## since this node is starting a component, it has a (0, 1) mapping
        ## which are the default values so we dont have to do anything here.
        ## establish a set for the next component
        component = { node }
        ## build a queue with just this node
        queue = [ node ]
        ## go through the queue until it is empty
        ## we assume that if you are in the queue, you already have a valid S_XC and F_XC
        while queue:
            ## pop the next item off the queue
            x = queue.pop(0)
            ## get its neighbors
            neighbors = incidence_graph[x] 
            ## forget about things we've already visited
            neighbors.difference_update(component)
            ## remove the remaining nodes from the global set
            nodes.difference_update(neighbors)
            ## add them to the current component
            for y in neighbors:
                ## get all the relevant maps
                s_yY, f_yY = S_yX[y, y], F_yX[y, y]  ## from y to Y (its extension)
                s_yX, f_yX = S_yX[y, x], F_yX[y, x]  ## from y to X (the extended node already in the cluster)
                s_XC, f_XC = S_XC[x], F_XC[x]        ## from X to C (since X is in the cluster, already labeled.
                ## and now some reversals and compositions.
                s_yC, f_yC = compose_maps(s_yX, f_yX, s_XC, f_XC)
                s_Yy, f_Yy = invert_map(s_yY, f_yY)
                s_YC, f_YC = compose_maps(s_Yy, f_Yy, s_yC, f_yC)
                ## and store the answer into the S_XC and F_XC arrays
                S_XC[y] = s_YC
                F_XC[y] = f_YC
            component.update(neighbors)
            ## and add them to the queue (so we visit them)
            queue.extend(neighbors)
        clusters.append(component)
    return clusters, S_XC, F_XC


'''
converts an extension map array into a matrix
also returning the mappings from contig y to extension X.
'''
def get_incidence_and_maps(extension_map):
    N = len(extension_map)
    CINDEX      = 0
    SINDEX      = 1
    FINDEX      = 2
    ## build extension array, i, j = 1 if i extends to j
    ## also S_yX, F_yX is the mapping from contig y into extension X.
    extend_array = np.zeros(shape=(N, N), dtype=int)
    S_yX =  np.zeros(shape=(N, N), dtype=int)
    F_yX =  np.ones(shape=(N, N), dtype=int)
    for i in range(N):
        for emap in extension_map[i]:
            j = emap[CINDEX]
            extend_array[i, j] = 1
            S_yX        [j, i] = emap[SINDEX]
            F_yX        [j, i] = emap[FINDEX]
    return extend_array, S_yX, F_yX


'''
Given:
    a list of contigs
    top quality maps of reads to contigs (good length and ratio)
    good ratio maps of reads to contigs (not worried about length)
    a deBruijn graph on the contigs
    the set of reads
    parameters for the minimum cover and minimum ratio to continue an extension
Will return:
    for each contig its extension (according to extend to contgi)
    for each contig extension, a set of maps from other contigs (y) to this one (x)
               provided that this extension is a perfect match overlap with y and at least K bases
'''
def get_extensions_and_maps(contigs, top_quality_maps, good_ratio_maps, graph, reads, MIN_COVER, MIN_RATIO, EXTENT, K, VERBOSE=False):
    N = len(contigs)
    extension     = []
    extension_map = []
    for contig_index in range(N):
        ## find all top quality maps in this contig
        tqm_contig_filter = top_quality_maps[:, CONTIG_NAME] == contig_index
        ## collect their read names (integers)
        well_mapped_reads = np.unique(top_quality_maps[tqm_contig_filter][:, READ_NAME])
        ## generate extended sequence and maps
        extended_seq, extended_seq_maps = extend_contig(contig_index,
                                                        contigs,
                                                        well_mapped_reads, 
                                                        good_ratio_maps, 
                                                        graph, 
                                                        reads, 
                                                        MIN_COVER, MIN_RATIO, 
                                                        EXTENT, K)
        if VERBOSE:
            rcount = len(well_mapped_reads) 
            clen   = len(contigs[contig_index])
            elen   = len(extended_seq)
            extended = int(clen != elen)
            info = [contig_index, extended, rcount, clen, elen]
            ## could modify code to not evaluate when the clen is too short, but they also don't have well mapped reads...
            ## also, this function does not know about MIN_MAP_LEN
            if rcount > 0:
                print tabprint(info)
                sys.stdout.flush()
        extension.append(extended_seq)
        extension_map.append(extended_seq_maps)
    return extension, extension_map


'''
Given a set of well mapped reads for a contig index,
will return  
'''
def extend_contig(contig_index, contigs, well_mapped_reads, good_ratio_maps, graph, reads, MIN_COVER, MIN_RATIO, EXTENT, K):
    if len(well_mapped_reads) == 0:
        ## this is an empty return: your own sequence and an identity map to self 
        return contigs[contig_index], [[contig_index, 0, 1, len(contigs[contig_index]), len(contigs[contig_index])]]
    ## then select those from reads that are well mapped to this contig
    ## crucially, this includes maps to other contigs
    grm_well_mapped_filter = sorted_filter(good_ratio_maps[:, READ_NAME], well_mapped_reads)
    filtered_maps          = good_ratio_maps[grm_well_mapped_filter]
    ## we get the map extension --
    ## i.e. how other contigs align to this contig according to the DB graph
    ## centered at cindex and extended by EXTENT bases in either direction
    S_XC, F_XC, in_frame = get_map_extension(contig_index, contigs, graph, K, EXTENT)
    ## we use the extension to align filtered_maps
    ## converting to a common frame centered on the target contig
    read_maps = align_reads_by_extension(filtered_maps, S_XC, F_XC, in_frame)
    ## collect all reads with unique maps and forward-align them to the target contig 
    rseqs, rmaps = get_singleton_maps(well_mapped_reads, read_maps, reads, FLIP_FORWARD=True)
    ## all sequences are properly oriented 
    ## and the map intercept is the start position for the read
    ## we aggregate information counting bases at each position
    base_count, contig_start  = get_base_counts(rseqs, rmaps, ALL_FORWARD = True)
    contig_end = contig_start + len(contigs[contig_index])
    extended_seq, s_ce = get_extension(base_count, 
                                               contig_start, 
                                               contig_end, 
                                               MIN_COVER, 
                                               MIN_RATIO)    
    extended_seq_rc   = reverse_complement(extended_seq)
    extended_seq_maps = []
    ## use extension map to align extended seq to other contigs in its dbgraph extension
    for cind in np.where(in_frame)[0]:
        ## load map between this contig and cind
        s_xc, f_xc   = S_XC[cind], F_XC[cind]
        ## compose with map between this contig and its extension
        s_xe, f_xe   = compose_maps(s_xc, f_xc, s_ce, 1)
        ## check the match
        match, total = check_match(s_xe, f_xe, contigs[cind], extended_seq, extended_seq_rc)
        ## if the match is perfect, AND NON-TRIVIAL, add to our list of extended maps
        if match == total >= K:
            map_data = [cind, s_xe, f_xe, len(contigs[cind]), match]
            extended_seq_maps.append(map_data)
    return extended_seq, extended_seq_maps

'''
Given a base count array and the contig start and end
and parameters for the minimum coverage and ratio required to call a sticky end,
returns the whole extended contig
'''
def get_extension(base_count, contig_start, contig_end, MIN_COVER, MIN_RATIO):
    ## define a denominator that is the base_sum per position
    ## if coverage goes to zero (shouldn't but could) clipping sets to 1
    coverage = np.clip(np.sum(base_count, axis=1), 1, np.inf)
    base_ratio = base_count / coverage[:, np.newaxis]
    top_ratio = np.max(base_ratio, axis=1)
    ## we are only going to extend the map if there is good sequence agreement
    ## and adequate sequence coverage to justify extending the contig
    ## any disruption in the extension ends the extension. 
    sticky_filter =  (coverage >= MIN_COVER) * (top_ratio > MIN_RATIO)
    touchy_spots  = np.where(~sticky_filter)[0]
    ## (+1) because start at first base after touchy spot
    touchy_start  = touchy_spots[touchy_spots < contig_start]
    sticky_start  = 0 if len(touchy_start) == 0 else np.max(touchy_start) + 1
    touchy_end    = touchy_spots[touchy_spots >= contig_end]
    sticky_end    = len(coverage) if len(touchy_end) == 0 else np.min(touchy_end)
    ## find the sequence by maximum ratio
    max_int     = np.argmax(base_count, axis=1)
    max_seq     = np.array(int_to_seq(max_int))
    ## get full extended sequence
    return "".join(max_seq[sticky_start:sticky_end]), contig_start - sticky_start

'''
Given a list of sequences and forward maps,
returns an array summarizing the read data over the extent covered
(each read contributing a count of one for the base observed)
and the "zero-map" position in the array.
'''
def get_base_counts(rseqs, rmaps, ALL_FORWARD = False, fix_zero = False, max_pos = None, multi=None):
    if not ALL_FORWARD:
        rseqs, rmaps = set_all_forward(rseqs, rmaps)
        
    map_intercept = np.array([rmap[0] for rmap in rmaps])
    read_len      = np.array([len(rseq) for rseq in rseqs])
    ## find min and max positions
    min_pos      = np.min(map_intercept)
    if max_pos == None:
        max_pos      = np.max(map_intercept + read_len)
    if fix_zero:
        min_pos = 0
    try:
        multi[0]
    except:
        multi = np.ones(len(rseqs), dtype=int)
    ## determine the total length of the extended contig
    extended_len = max_pos - min_pos 
    ## subtract min_pos to convert from contig mapping to array mapping
    array_intercept = map_intercept - min_pos
    ## keep track of bases observed in the data
    base_count      = np.zeros(shape = (extended_len, 5), dtype=int)
    ## add them in
    for x, r, m in zip(array_intercept, rseqs, multi):
        pos_list = np.arange(x, x+len(r))
        seq_list = seq_to_int(r)
        np.add.at(base_count, (pos_list, seq_list), m)
    return base_count, -min_pos

def set_all_forward(rseqs, rmaps):
    fseqs, fmaps = [], []
    for rseq, rmap in zip(rseqs, rmaps):
        if rmap[1] == -1:
            rmap = reverse_A(rmap[0], rmap[1], len(rseq))
            rseq = reverse_complement(rseq)
        fseqs.append(rseq)
        fmaps.append(rmap)
    return fseqs, fmaps

'''
given a sorted list of keys (key_list) 
and sorted list of good keys (good_key_list)
returns a boolean numpy vector True if the key_list value is in good_key_list, False otherwise
'''
def sorted_filter(key_list, good_key_list):
    N = len(key_list)
    dfilter = np.zeros(shape=N, dtype=bool)
    if len(good_key_list) == 0:
        return dfilter    
    good_key_iter = iter(good_key_list)
    current_key = good_key_iter.next()
    ind = 0
    for ind, key in enumerate(key_list):
        if key < current_key:
            continue
        elif key == current_key:
            dfilter[ind] = True
        else:
            while current_key < key:
                try:
                    current_key = good_key_iter.next()
                except:
                    break ## we are done (rest of the vector is set to false)
            if key == current_key:
                dfilter[ind] = True        
    return dfilter


'''
Given a list of maps to contigs -AND- a map extension
returns a dictionary keyed by the read_index
such for that for read1 [0] and read2 [1],
contains a Counter of S_RC and F_RC maps to the target contig
'''
def align_reads_by_extension(map_list, S_XC, F_XC, in_frame):
    read_maps = defaultdict(lambda: (set(), set()))
    
    ## for each filtered map
    for rmap in map_list:
        rind = rmap[READ_NAME]      ## which read
        r12  = rmap[READ_12]
        cind = rmap[CONTIG_NAME]        ## to which mapped contig
        if in_frame[cind]:              ## if that contig is in the frame
            s_rx = rmap[SRC]               ## read to mapped contig
            f_rx = rmap[FRC]
            s_xc = S_XC[cind]              ## mapped contig to target contig
            f_xc = F_XC[cind]
            s_rc, f_rc = compose_maps(s_rx, f_rx, s_xc, f_xc)   ## compose 
            read_maps[rind][r12].add( (s_rc, f_rc) )           ## and save
    return read_maps 

'''
Given a read_maps object (defaultdicts on the reads, etc) returns just those with a single map.
If FLIP_FORWARD is true, the sequences returned are forward matched to the reference contig --
if the map is a reverse map, the read is reverse complemented and the map is reversed too.
'''
def get_singleton_maps(well_mapped_reads, read_maps, reads, FLIP_FORWARD=False):
    ## collect all singleton maps
    singleton_seq, singleton_map = [], []
    for rindex in well_mapped_reads:
        for read12 in [0, 1]:
            read_map = read_maps[rindex][read12]
            if len(read_map) == 1:
                rseq = reads[read12][rindex]
                rmap = next(iter(read_map))
                if FLIP_FORWARD and rmap[1] == -1:
                    rmap = reverse_A(rmap[0], rmap[1], len(rseq))
                    rseq = reverse_complement(rseq)
                singleton_seq.append(rseq)
                singleton_map.append(rmap)
    return singleton_seq, singleton_map


'''
Load annotated maps from directory pickle
Annotation is:
0  read_name
1  read1 (0) or read2 (1)
2  contig_name
3  contig position
4  read position
5  match length
6  reverse complement to match (1) or not(0)
7  S_RC (zero intercept of read in contig)
8  F_RC +1 if forward map -1 otherwise
9  bases of agreement (match)
10 bases of overlap   (total)
'''
READ_NAME   = 0
READ_12     = 1
CONTIG_NAME = 2
CONTIG_POS  = 3
READ_POS    = 4
MATCH_LEN   = 5
RC          = 6
SRC         = 7
FRC         = 8
MATCH       = 9
OVERLAP     = 10
def load_maps(outdir, maps_pickle):
    maps_picklefile = os.path.join(outdir, maps_pickle)
    return np.load(maps_picklefile)


'''
Given a list of sequences, returns only those that do not occur as proper substrings in the list
'''
def remove_proper_subsequences(seq_list):
    seq_list_rc = [reverse_complement(seq) for seq in seq_list]
    ans = []
    for target in seq_list:
        forward_count = np.sum([seq.find(target) != -1 for seq in seq_list])
        reverse_count = np.sum([seq.find(target) != -1 for seq in seq_list_rc])
        ## if unique, should only find itself in forward and nothing in reverse
        if (forward_count == 1) and (reverse_count == 0) :
            ans.append(target)
    return ans

'''
Given a set of contigs and reads (already mapped and annotated)
1. extend the contigs according to well-mapped reads and their mates,
    using the dbGraph to map mates in-frame
2. join contigs that extend perfectly into each other
3. drop the redundant subsequences

Parameters are:
outdir            : working directory
contig_shortname  : shortname for contig file
extended_shortname: shortname for output file of extended contigs
K                 : kmer size
MIN_MAP_RATIO     : percent bases agreement with contig to count a map
MIN_MAP_LENGTH    : minimum length of a read to count as well-mapped
EXTENT            : how far to go in the dbGraph for mapping mates (~max value for insert size)
MIN_COVER         : minimum coverage to call an extension
MIN_RATIO         : minimum ratio to call an extension
UNIQUE_ONLY       : if True, remove proper subsequences from final extension (default=True)
'''
def extend_contigs(outdir, contig_shortname, maps_pickle, extended_shortname,
                   K, MIN_MAP_RATIO, MIN_MAP_LENGTH,
                   EXTENT, MIN_COVER, MIN_RATIO, UNIQUE_ONLY=True, VERBOSE=False):
    contigs = load_contigs  (outdir, contig_shortname)
    reads   = load_readpairs(outdir)
    maps    = load_maps     (outdir, maps_pickle)
    
    ## establish filters
    map_ratio      = maps[:, MATCH] / maps[:, OVERLAP].astype(float)
    ratio_filter   = map_ratio   >= MIN_MAP_RATIO
    length_filter  = maps[:, OVERLAP] >= MIN_MAP_LENGTH   
    quality_filter = ratio_filter * length_filter
    
    ## get best maps
    top_quality_maps = maps[quality_filter]
    good_ratio_maps  = maps[ratio_filter]
    
    ## build the contig graph
    print "building dbGraph"
    graph = contig_graph(contigs, K)
    print "getting extensions"
    headings = ['index', 'ext', 'rcount', 'clen', 'elen']
    print tabprint(headings)
    ## get extensions and maps from contigs to each extension (if the contig is a proper match)
    extension, extension_map = get_extensions_and_maps(contigs, 
                                                       top_quality_maps, good_ratio_maps, 
                                                       graph, reads, 
                                                       MIN_COVER, MIN_RATIO, 
                                                       EXTENT, K, VERBOSE)
    
    ## convert the extension maps into a set of arrays
    extend_array, S_yX, F_yX = get_incidence_and_maps(extension_map)
    ## we add an edge if A extends to B and B extends to A.
    ## this is the element-wise AND of the extend_array and its transpose
    incidence_array = extend_array * extend_array.T
    ## use incidence and mapping to get clusters with common cluster mappings
    clusters, S_XC, F_XC = get_extension_clusters(incidence_array, S_yX, F_yX)
    ## convert cluster to their consensus sequence
    cluster_seqs = [cluster_to_sequence(list(cluster_set), S_XC, F_XC, extension) for cluster_set in clusters]
    ## remove matching subsequences
    if UNIQUE_ONLY:
        print 'removing proper sub-sequences'
        cluster_seqs = remove_proper_subsequences(cluster_seqs)
    ## order by length of contig
    order = np.argsort([len(x) for x in cluster_seqs])[::-1]
    cluster_seqs = [cluster_seqs[x] for x in order]
    if VERBOSE:
        print "from %d initial contigs" % len(contigs)
        print "writing %d extended contigs to file" % len(cluster_seqs)
    write_contigs(cluster_seqs, outdir, extended_shortname)

'''
Converts a sequence of letters into the integer equivalents
'''
def seq_to_int(seq):
    return [BASE2INT[x] for x in seq]

'''
Convert a sequence of integers into letter equivalents
'''
def int_to_seq(int_array):
    return [BASES[x] for x in int_array]

        
'''
Returns the reverse complement of the given sequence.
'''
def reverse_complement(seq):
    return "".join([INDICT[base] for base in seq[::-1]])

# '''
# Load kmer count data from text file
# '''
# def get_kmer_count(kmer_filename, MIN_MER = 0):
#     counter = 0
#     with open(kmer_filename, 'r') as kmer_file:
#         for line in kmer_file:
#             count = int(line.strip().split()[1])
#             if count >= MIN_MER:
#                 counter += 1
#     return counter

'''
Load kmer data from a file.
Puts kmers into sorted count order.
'''
def load_kmer_file(kmer_filename):
    with open(kmer_filename, 'r') as kmer_file:
        kmers  = []
        counts = []
        for line in kmer_file:
            info = line.strip().split()
            kmers.append (     info[0] )
            counts.append(int( info[1] ))
    
    kmers  = np.array(kmers )
    counts = np.array(counts)
    order  = np.argsort(counts)[::-1]
    
    return kmers [order], counts[order]



'''
Load kmer data from a file.
Puts kmers into sorted count order.
'''
def load_kmer_file_count_filtered(kmer_filename, MIN_COUNT=0):
    with open(kmer_filename, 'r') as kmer_file:
        kmers  = []
        counts = []
        for line in kmer_file:
            info = line.strip().split()
            count = int( info[1] )
            if count >= MIN_COUNT:
                kmers.append ( info[0] )
                counts.append( count   )
    kmers  = np.array(kmers )
    counts = np.array(counts)
    order  = np.argsort(counts)[::-1]
    return kmers [order], counts[order]

'''
Load kmer data from a file.
Puts kmers into sorted count order.
'''
def kmer_file_to_db_graph(kmer_filename, MIN_COUNT=0):
    db_graph = defaultdict(int)
    with open(kmer_filename, 'r') as kmer_file:
        for line in kmer_file:
            info = line.strip().split()
            count = int( info[1] )
            if count >= MIN_COUNT:
                kmer = info[0]
                db_graph[ kmer                     ] = count
                db_graph[ reverse_complement(kmer) ] = count
    return db_graph


'''
Returns the hamming difference between two strings.
'''
def hamming(str1, str2):
    #ne = str.__ne__  ## this is surprisingly slow
    ne = operator.ne
    return sum(imap(ne, str1, str2))

# 
# def get_key(kmer):
#     kmer_rc = reverse_complement(kmer)
#     return kmer if kmer < kmer_rc else kmer_rc

'''
Get kmers forward as an iterator for the given kmer km.
'''
def fw(km):
    for x in "ACGT":
        yield km[1:] + x

'''
Get kmers backwards as an iterator for the given kmer km.
'''
def bw(km):
    for x in "ACGT":
        yield x + km[:-1]

'''
Convert a sequence of kmers into the string represenetation
'''
def contig_to_string(contig):
    return contig[0] + "".join([x[-1] for x in contig[1:]])

'''
Gets forward contigs.
Used with reverse complement to get reverse (see get_contig)
'''
def get_contig_forward(db_graph, km):
    path_forward = [km]
    while True:
        last_element = path_forward[-1]
        candidates = [x for x in fw(last_element) if x in db_graph]
        ## if there's no choice or more than once choice, get out!
        if len(candidates) != 1:
            break
        ## otherwise, get your choice
        candidate    = candidates[0]
        candidate_rc = reverse_complement(candidate)
        ## these are probably overkill and we can trim them back
        ## for now, they make me feel safe and warm and fuzzy 
        if candidate in path_forward:
            break  ## get out of loops
        if candidate_rc in path_forward:
            break  ## watch out for mobius contigs and hairpins
        ## check that there's only one path into that place
        into_candidate = [x for x in bw(candidate) if x in db_graph]
        if len(into_candidate) != 1:
            break
        path_forward.append(candidate)
    return path_forward

'''
Identify the contig that contains the given node
'''
def get_contig(db_graph, node):
    forward_path  = get_contig_forward(db_graph, node)
    backward_path = get_contig_forward(db_graph, reverse_complement(node))
    ## if it is a cycle starting at node
    if node in fw(forward_path[-1]):
        contig = forward_path
    else:
        contig = [reverse_complement(x) for x in backward_path[-1:0:-1]] + forward_path
    return contig_to_string(contig), contig

'''
Gets all dbGraph contigs by exhaustion of the dbGraph nodes using get_contig
Results are ordered by contig length (longest first)
'''
def all_contigs(db_graph):
    done = set()
    r = []
    for node in db_graph:
        if node not in done:
            contig_string, contig = get_contig(db_graph, node)
            for cnode in contig:
                done.add(cnode)
                done.add(reverse_complement(cnode))
            r.append(contig_string)
    
    contigs  = np.array(r) 
    seq_lens = np.array([len(x) for x in contigs])
    order    = np.argsort(seq_lens)[::-1]
    
    contigs  = contigs[order]
    return contigs

'''
Returns all the kmers in the input string sequentially (iteratively).
'''
def get_kmers(seq, K):
    for i in range(len(seq) - K + 1):
        yield seq[i:i+K]

'''
Given an input path containing files r1.fastq.gz and r2.fastq.gz
generates r1.fastq and r2.fastq in the outpath
If INTEGER_NAMES = True, then it will replace the full read name with an integer
corresponding to its position in the file.
Returns the number of reads unzipped and True/False if it finished (names matched up)
'''
def unzip_fastq_pair(inpath, outpath, INTEGER_NAMES=False, VERBOSE=False):
    if VERBOSE:
        print "making directory", outpath
    ## make output directory
    try:
        os.makedirs(outpath)
    except:
        if VERBOSE:
            print outpath, "already exists"
    
    finished = True
    ## convert read1 into text in outdir
    if VERBOSE:
        print "writing read1"
    read_file = os.path.join(inpath , "r1.fastq.gz")
    out_file  = os.path.join(outpath, "r1.fastq")
    read_name_list = []
    with gzip.open(read_file, 'rb') as f_in:
        with open(out_file, 'w') as f_out:
            if INTEGER_NAMES:                
                for line_counter, line in enumerate(f_in):
                    if line_counter % 4 == 0:
                        rindex = line_counter / 4
                        read_name_list.append(line.split("/")[0])
                        f_out.write( "@%d\n" % rindex)
                    else:
                        f_out.write(line)
                    
            else:
                shutil.copyfileobj(f_in, f_out)
    ## convert read2 into text in outdir
    if VERBOSE:
        print "writing read2"
    read_file = os.path.join(inpath , "r2.fastq.gz")
    out_file  = os.path.join(outpath, "r2.fastq")
    with gzip.open(read_file, 'rb') as f_in:
        with open(out_file, 'w') as f_out:
            if INTEGER_NAMES:
                for line_counter, line in enumerate(f_in):
                    if line_counter % 4 == 0:
                        rindex = line_counter / 4
                        name_check = line.split("/")[0]
                        if name_check != read_name_list[rindex]:
                            print "ERROR! read names do not agree!"
                            finished = False
                            break
                        f_out.write( "@%d\n" % rindex)
                    else:
                        f_out.write(line)
            else:
                shutil.copyfileobj(f_in, f_out)
    return len(read_name_list), finished

'''
This function is useful when combining sequencing runs or datasets.

Given a list of input path containing files r1.fastq.gz and r2.fastq.gz
generates r1.fastq and r2.fastq in the outpath
Will rename everything with integer labels 
and also generate index_starts.txt which gives start and end index for each input path.
'''
def unzip_fastq_pairs(inpaths, outpath, VERBOSE=False):
    if VERBOSE:
        print "making directory", outpath
    ## make output directory
    try:
        os.makedirs(outpath)
    except:
        if VERBOSE:
            print outpath, "already exists"
    
    start_ind = 0
    out_file1  = os.path.join(outpath, "r1.fastq")
    out_file2  = os.path.join(outpath, "r2.fastq")
    out_file3  = os.path.join(outpath, "index_starts.txt")
    
    f_out1     = open(out_file1, 'w') 
    f_out2     = open(out_file2, 'w')
    f_out3     = file(out_file3, 'w')
    for inpath in inpaths:
        if VERBOSE:
            print inpath, start_ind
        ## convert read1 into text in outdir
        if VERBOSE:
            print "writing read1"
        read_file = os.path.join(inpath , "r1.fastq.gz")
        read_name_list = []
        with gzip.open(read_file, 'rb') as f_in:
            for line_counter, line in enumerate(f_in):
                if line_counter % 4 == 0:
                    rindex = line_counter / 4
                    read_name_list.append(line.split("/")[0])
                    f_out1.write( "@%d\n" % (rindex + start_ind))
                else:
                    f_out1.write(line)
        ## convert read2 into text in outdir
        if VERBOSE:
            print "writing read2"
        read_file = os.path.join(inpath , "r2.fastq.gz")

        with gzip.open(read_file, 'rb') as f_in:
            for line_counter, line in enumerate(f_in):
                if line_counter % 4 == 0:
                    rindex = line_counter / 4
                    name_check = line.split("/")[0]
                    if name_check != read_name_list[rindex]:
                        print "ERROR! read names do not agree!"
                        break
                    f_out2.write( "@%d\n" % (rindex + start_ind))
                else:
                    f_out2.write(line)
        f_out3.write(tabprint([inpath, start_ind, rindex+start_ind]) + "\n")
        start_ind = rindex + start_ind + 1
    f_out1.close()
    f_out2.close()
    f_out3.close()

'''
Runs jellyfish on the selected directory (r1.fastq and r2.fastq)
for a given kmer size and minimum mer count.
Also writes the output to a text readable file.
'''
def execute_jellyfish(outdir, KMER_SIZE, MIN_MER_COUNT, JELLYFISH, JELLYFISH_MEM, JELLYFISH_THREADS, VERBOSE=False):
    ## set outdir to current working directory 
    os.chdir(outdir)    
    command_count = "%s count -m %d -s %s -t %d -C -F 2 r1.fastq r2.fastq" % (JELLYFISH, KMER_SIZE, JELLYFISH_MEM, JELLYFISH_THREADS)
    command_dump  = "%s dump -L %d -c mer_counts.jf > kmer_list.txt" % (JELLYFISH, MIN_MER_COUNT)
    if VERBOSE:
        print "generating jellyfish counts"
        print command_count
    os.system(command_count)
    if VERBOSE:
        print "dumping counts to file"
        print command_dump
    os.system(command_dump)


'''
For a set of kmers with counts, we return all contigs in the db graph, provided count > 
'''
def get_contigs(kmers, counts, MIN_COUNT):
    ## initialize db_graph with all mers that meet minimum count
    db_graph = defaultdict(int)
    ## include also reverse complement
    for mer, count in zip(kmers, counts):
        if count >= MIN_COUNT:
            db_graph[ mer                       ] = count
            db_graph[ reverse_complement( mer )  ] = count
    
    return all_contigs(db_graph)

'''
Given a working directory, convert the kmers into contigs
and write the resulting contigs to file
'''
def kmers_to_contigs(outdir, contig_shortname, MIN_COUNT, VERBOSE=False):
    kmer_filename = os.path.join(outdir, "kmer_list.txt")
    contig_filename = os.path.join(outdir, contig_shortname)
    if VERBOSE:
        print "loading kmers"
    kmers, counts   = load_kmer_file(kmer_filename)
    if VERBOSE:
        print np.sum(counts >= MIN_COUNT), "kmers >= ", MIN_COUNT, "depth"
        print "getting contigs"
    contigs         = get_contigs(kmers, counts, MIN_COUNT=MIN_COUNT)
    if VERBOSE:    
        print "writing %d contigs to file" % len(contigs)
    contig_file = file(contig_filename, 'w')
    for ind, contig in enumerate(contigs):
        contig_file.write( ">%d\n" % (ind))
        contig_file.write(contig)
        contig_file.write("\n")
    contig_file.close()
    if VERBOSE:
        print "done"



''' 
Use mumdex to map reads to contigs.

mumdex maps have an output in the form:
    read name
    read1 (0) or read2 (1)
    contig name
    contig position
    read position
    match length
    reverse complement to match (-) or not(+)
'''
def mumdex_map_to_contig(outdir, contig_shortname, maps_shortname, MIN_MUM_LENGTH, MUMDEX_FASTQ, VERBOSE=False):
    contig_filename = os.path.join(outdir, contig_shortname)
    maps_filename   = os.path.join(outdir, maps_shortname)
    reference_directory = contig_filename + ".bin"
    ## remove the target reference directory
    if os.path.isdir(reference_directory):
        os.system('rm -R %s' % reference_directory)
        ## MEMDEX mapping command
    MUMDEX_MAPPER = "%s -l %d %s %s/r1.fastq %s/r2.fastq > %s"
    mumdex_command = MUMDEX_MAPPER % (MUMDEX_FASTQ, MIN_MUM_LENGTH, contig_filename, outdir, outdir, maps_filename)
    if VERBOSE:
        print "mumdexing"
        print mumdex_command
    os.system(mumdex_command)
    if VERBOSE:
        print "and done"

'''
Load just the sequences from the read pairs, and puts into a paired array
to get read names, set INDEX to 0
to get a list of "+" signs, set INDEX to 2
to get quality scores, set INDEX to 3
'''
def load_readpairs(read_dir, INDEX=1):
    read1_filename = os.path.join(read_dir, "r1.fastq")
    read2_filename = os.path.join(read_dir, "r2.fastq")
    
    reads = [[],
             []]
    
    with open(read1_filename, 'r') as f_in:
        for line_counter, line in enumerate(f_in):
            if line_counter % 4 == INDEX:
                reads[0].append(line.strip())
    
    with open(read2_filename, 'r') as f_in:
        for line_counter, line in enumerate(f_in):
            if line_counter % 4 == INDEX:
                reads[1].append(line.strip())
    return reads

''' 
load the contigs into a vector
'''
def load_contigs(outdir, contig_shortname, return_names=False):
    if return_names:
        names = []
    contigs = []
    contig_filename  = os.path.join(outdir, contig_shortname) 
    with open(contig_filename, 'r') as f_in:
        for line_counter, line in enumerate(f_in):
            if return_names:
                if line_counter % 2 == 0:
                    names.append(line.strip()[1:])
            if line_counter % 2 == 1:
                contigs.append(line.strip())
    if return_names:
        return contigs, names
    else:
        return contigs

'''
writes a vector of contigs to file
'''
def write_contigs(contigs, outdir, contig_shortname, contig_names = None):
    contig_filename  = os.path.join(outdir, contig_shortname)
    contig_file = file(contig_filename, 'w')
    for ind, contig in enumerate(contigs):
        if contig_names == None:
            contig_file.write( ">%d\n" % (ind))
        else:
            contig_file.write( ">%s\n" % (contig_names[ind]))
        contig_file.write(contig)
        contig_file.write("\n")
    contig_file.close()


#### SECTION CONTAINS MATH FOR MAPPING READS TO CONTIGS
'''
Given an intercept mapping (S_AB, F_AB), and B' = reverseComplement(B),
Then S_AB' = len_B - 1 - S_AB, and F_AB' = -F_AB
'''
def reverse_B(S_AB, F_AB, len_B):
    return len_B - 1 - S_AB, -F_AB

'''
Given an intercept mapping (S_AB, F_AB), and A' = reverseComplement(A),
Then S_A'B = S_AB + (F_AB*len_A), and F_A'B = -F_AB
'''
def reverse_A(S_AB, F_AB, len_A):
    return S_AB + (F_AB*(len_A - 1)), -F_AB

'''
Given an intercept mapping (S_AB, F_AB) returns (S_BA, F_BA)
'''
def invert_map(S_AB, F_AB):
    S_BA = -F_AB*S_AB
    return S_BA, F_AB

'''
Given a mapping from A->B and B->C returns map from A->C
'''
def compose_maps(S_AB, F_AB, S_BC, F_BC):
    S_AC = S_BC + F_BC*S_AB
    F_AC = F_AB*F_BC
    return S_AC, F_AC


def evaluate_map(S_AB, F_AB, a): 
    return S_AB + F_AB*a

'''
Given a mapping from A->X and B->X returns map from A->B
'''
def compose_from_common_map(S_AX, F_AX, S_BX, F_BX):
    S_XB, F_XB = invert_map(S_BX, F_BX)
    return compose_maps(S_AX, F_AX, S_XB, F_XB)

''' 
return number of matching bases and total number of bases between seq_A and seq_B
given the intercept mapping (presumed forward match, i.e. F_AB = 1)
see check_match for generic version.
'''
def check_match_forward(S_AB, seq_A, seq_B):
    if S_AB < 0:
        A = seq_A[-S_AB:]
        B = seq_B
    else:
        A = seq_A
        
        B = seq_B[S_AB:]
    tot = min(len(A), len(B))
    return tot - hamming(A, B),tot 

''' 
return number of matching bases and total number of bases between seq_A and seq_B
given the intercept mapping
'''
def check_match(S_AB, F_AB, seq_A, seq_B, seq_B_rc):
    if F_AB == 1:
        return check_match_forward(S_AB, seq_A, seq_B)
    else:
        S_ABc, F_ABc = reverse_B(S_AB, F_AB, len(seq_B))
        return check_match_forward(S_ABc, seq_A, seq_B_rc)



'''
given a working directory and a list of contigs and mumdex maps
returns an annotated version of the mumdex maps that includes the intercept mappings,
the total number of matched bases and the total number of overlapping bases.
If PICKLE = True, will also pickle the numpy integer array.

Annotation is:
0  read_name
1  read1 (0) or read2 (1)
2  contig_name
3  contig position
4  read position
5  match length
6  reverse complement to match (1) or not(0)
7  S_RC (zero intercept of read in contig)
8  F_RC +1 if forward map -1 otherwise
9  bases of agreement (match)
10 bases of overlap   (total)
'''
def annotate_maps(outdir, contig_shortname, maps_shortname, PICKLE=True, VERBOSE=False):
    name_head          = maps_shortname.split(".")[0]
    extended_head      = name_head + "_plus"
    extended_shortname = extended_head + ".txt"  
    extended_pickle    = extended_head + ".npy"
    
    contigs = load_contigs  (outdir, contig_shortname)
    reads   = load_readpairs(outdir)
    
    ## so we only have to reverse_complement one time when annotating maps:
    contigs_rc  = [reverse_complement(contig) for contig in contigs]
    contigs_len = [len(contig) for contig in contigs] 
    
    ## load mumdex maps
    maps_filename          = os.path.join(outdir, maps_shortname)
    extended_maps_filename = os.path.join(outdir, extended_shortname)
    extended_maps_pickle   = os.path.join(outdir, extended_pickle)
    
    ## vector to save all maps
    all_maps = []
    print "annotating maps"
    with open(maps_filename, 'r') as f_in:
        for line_counter, line in enumerate(f_in):
            entry = line.strip().split(" ")
            ## convert from string to int
            entry[-1] = entry[-1] == "-"
            entry     = map(int, entry)
            ## extract info
            read_i, read_pair_i, contig_i, contig_x, read_x, match_len, is_rc = entry          
            ## extract read sequence
            read_seq    = reads[read_pair_i][read_i]
            ## the mumdex program identifies the stretch of agreement between the read and the sequence.
            ## S_RC is the zero position of the read in the contig co-ordinates.
            ## F_RC is 1 if forward mapping, -1 otherwise.
            ## if the orientations are the same, S_RC is straightforward (contig_x- read_x)
            ## if the orientations are reversed, S_RC is (contig_x + read_x + match_len - 1)
            ## this can be calculated explicitly from composed maps over the common MUM
            S_RC = contig_x + read_x + match_len - 1 if is_rc else contig_x - read_x
            F_RC = (-1)**is_rc
            ## to extend the annotation, we will check the quality of the match
            ## if the orientations are different, we match to the reverse complement of the contig 
            ## if the orientations are matched, we proceed with the contig sequence
            if is_rc:
                contig_seq   = contigs_rc[contig_i]
                S_RCr, F_RCr = reverse_B(S_RC, F_RC, contigs_len[contig_i])   # F_RCr = 1 by design
                match, total = check_match_forward(S_RCr, read_seq, contig_seq)
            else:
                contig_seq   = contigs[contig_i]
                match, total = check_match_forward(S_RC, read_seq, contig_seq)
            ## add annotations to entry (mapping in S,F notation and match/total)
            entry.extend([S_RC, F_RC, match, total])
            all_maps.append(entry)
            if line_counter % 100000 == 0:
                if VERBOSE:
                    print ".",
    if VERBOSE:
        print
        print "done"
        print "writing to text" 
    
    ## write to text
    with open(extended_maps_filename, 'w') as f_out:
        for entry in all_maps:
            f_out.write(tabprint(entry) + "\n")
    
    if PICKLE:
        if VERBOSE:
            print "writing to pickle"
        ## write to pickle
        contig_maps_extended = np.array(all_maps)
        np.save(extended_maps_pickle, contig_maps_extended)
    if VERBOSE:
        print "done"


'''
Function builds the rest of a dbGraph from the contigs.
For each contig 0...N-1, graph[ind] contains two lists.
The first list are nodes incident on the head of the contig.
The second list are nodes that are incident on the tail of the contig.
A node is encoded as (ind, 0/1) where ind is the index of the contig,
0 is the head of the contig (start) 1 is the tail. 
'''
def contig_graph(contigs, K):
    graph = [([], []) for _ in contigs]
    head_in = dict()
    tail_in = dict()
    ## a unique lookup for the head and tail sequences for each contig
    for ind, seq in enumerate(contigs):
        head_in[seq[:K]] = (ind, 0)
        tail_in[reverse_complement(seq[-K:])] = (ind, 1)
    
    for ind, seq in enumerate(contigs):
        head_out = reverse_complement(seq[:K])
        tail_out = seq[-K:]
        for next_mer in fw(head_out):
            if next_mer in head_in:
                graph[ind][0].append(head_in[next_mer])
            if next_mer in tail_in:
                graph[ind][0].append(tail_in[next_mer])
        for next_mer in fw(tail_out):
            if next_mer in head_in:
                graph[ind][1].append(head_in[next_mer])
            if next_mer in tail_in:
                graph[ind][1].append(tail_in[next_mer])
    return graph


## if given a head node (x, 0) returns the tail node (x, 1)
## and the reverse.
def flip_node(node):
        return node[0], (node[1] + 1) % 2

'''
Given a DB graph, works out the mapping from other contigs
to a target contig by consider all extensions along the DBgraph.
Proceeds breadth first from each end of the contig,
assigning an alignment to each contig within a given extent from the target.

Return contains S_XC, F_XC and in_frame
recording the slope and intercept for contigs within the extent 
(those for which in_frame is True)
'''
def get_map_extension(cindex, contigs, graph, K, extent):
    N = len(contigs)
    contig_len = [len(x) for x in contigs]
    S_XC       = np.zeros(N, dtype=int) 
    F_XC       = np.ones (N, dtype=int)
    in_frame   = np.zeros(N, dtype=bool) 
    
    ## returns the position of the head or tail as indicated
    def posX(node):
        return 0 if node[1] == 0 else contig_len[node[0]] - K + 1
    
    def get_DBG_intercept(A, B):
        posA  = posX(A)
        posB  = posX(B)
        is_rc = A[1] == B[1]
        S_AB  = (posA + posB + K - 2) if is_rc else posB - posA
        F_AB  = (-1)**is_rc
        return S_AB, F_AB
    
    ## given a pair of edges update the data
    def update_data(old_node, new_node):
        nindex = new_node[0]
        ## if new contig already there, done
        if in_frame[nindex]:
            return False
        else:
            ## otherwise
            oindex = old_node[0]
            ## get the map from old_node to cindex
            S_OC, F_OC = S_XC[oindex], F_XC[oindex]
            ## compute map from new_node to old_node for DBGraph
            S_NO, F_NO = get_DBG_intercept(new_node, old_node)
            ## and compose the maps to get the new map
            S_NC, F_NC = compose_maps(S_NO, F_NO, S_OC, F_OC)
            S_XC[nindex] = S_NC
            F_XC[nindex] = F_NC
            in_frame[nindex] = True
            return True
    
    ## moves forward from one node updating data from alignment data
    ## as necessary and building a list of future nodes to visit
    def build_one_out(node, dist, max_dist):
        ## keep a list of new nodes and distances
        ans = []
        ## if we have exceed the boundary, we are done
        if dist > max_dist:
            return ans
        else:
            ## otherwise, find your neighbors
            neighbors = graph[node[0]][node[1]]
            ## and for each neighbor
            for new_node in neighbors:
                ## update data (if false, already inframe)
                if update_data(node, new_node):
                    ## if not already inframe, compute new distance and track new nodes
                    new_dist = dist + contig_len[new_node[0]] - K + 1
                    ans.append((flip_node(new_node), new_dist))
        return ans
    
    ## as above but works on a list
    def build_list_out(node_dist_list, max_dist):
        ans = []
        for node, dist in node_dist_list:
            ans.extend(build_one_out(node, dist, max_dist))
        return ans
    
    ## build back from head and forward from tail
    node_dist_list =[ ((cindex, 0), 0), ((cindex, 1), 0) ]
    while len(node_dist_list) > 0:
        node_dist_list = build_list_out(node_dist_list, extent)

    ## include contig
    in_frame[cindex] = True
    return S_XC, F_XC, in_frame

## check if Istring is a compatible preimage for Mstring under C->T mutation
## fails at the first sign of trouble
def is_35_convert(Istring, Mstring):
    for i, m in zip(Istring, Mstring):
        if i != m and not ((i == 'C') and (m == 'T')):
            return False
    return True

'''
For each contig, will try to reverse C>T conversions that are unambiguous in the unmutated data.
Also, places the sequences in the orientation (forward or reverse complement)
that is best explained by the unmutated data --
such that the conversions are C>T.

inputs:
pair_dir      : parent_directory for the pair bs_0 and bs_1
MIN_KMER      : smallest count for kmers in the unmutated library
KMER_SIZE     : as always the size of the kmer for dbGraphing
MIN_MAX_RATIO : if corrected base exceeds MIN_MAX_RATIO, use corrected base (unambiguous in pre-image)

Returns the observed mutation rate aggregated over all the data
'''
def unmask_contigs(pair_dir, input_contigs_shortname, output_contigs_shortname, unmasked_shortname, mask_shortname, logout_shortname, MIN_KMER, KMER_SIZE, MIN_MAX_RATIO, MUT_RATE_MIN=0.0, VERBOSE=False):
    bs_0_dir   = os.path.join(pair_dir, "bs_0")
    bs_1_dir   = os.path.join(pair_dir, "bs_1")
    ## output filenames
    unmasked_filename = os.path.join(bs_1_dir, unmasked_shortname)
    mask_filename     = os.path.join(bs_1_dir, mask_shortname)
    logout_filename   = os.path.join(bs_1_dir, logout_shortname)
    
    ## load bs1 extended contigs"
    contigs_bs1 = load_contigs(bs_1_dir, input_contigs_shortname)
    
    ## load bs0 kmers
    kmer_filename = os.path.join(bs_0_dir, "kmer_list.txt")
    kmers, counts   = load_kmer_file(kmer_filename)
    
    ## build conversion dictionary
    ## for each kmer (and its rc) of sufficient count (>=MIN_KMER) 
    ## first , keeps a count of unmutated kmers (bs0_count)
    ## second, convert C>T and store unmuated kmer in a list
    if VERBOSE:
        print "building conversion dictionary"
    bs0_count = Counter()
    ct_dict   = defaultdict(list)
    for kmer, count in zip(kmers, counts):
        if count >= MIN_KMER:
            kmer_rc    = reverse_complement(kmer)
            bs0_count[kmer   ] = count
            bs0_count[kmer_rc] = count
            kmer_ct    = kmer   .replace("C", "T") 
            kmer_rc_ct = kmer_rc.replace("C", "T")
            ct_dict[kmer_ct]   .append(kmer)
            ct_dict[kmer_rc_ct].append(kmer_rc)
    
    ## Performs C>T correction on a given sequence
    ## returns the total number of kmers in the preimage and the correction
    ## CHECKS THAT THE KMER IS A 3.5 MATCH
    ## and the correction mask
    def unmask_sequence(seq):
        ## first, convert C>T
        extended_ct  = seq.replace("C", "T")
        cseqs = []      ## list of unmasked kmers
        cmaps = []      ## list of their maps into this sequence
        ## go through the kmers in the sequence and look them up
        ## also keeps along the original kmer to check 3.5 compatibility (no "C" unmasked to "T")
        for base_ind, (kmer, kmer_CT) in enumerate(zip(get_kmers(seq, KMER_SIZE), get_kmers(extended_ct, KMER_SIZE))):
            new_mers = ct_dict[kmer_CT]
            new_mers35 = [new_mer for new_mer in new_mers if is_35_convert(new_mer, kmer)]
            cseqs.extend(new_mers35)
            cmaps.extend([(base_ind, 1)] * len(new_mers35))
        ## look up the count for each map
        multi = np.array([bs0_count[kmer] for kmer in cseqs])
        ## take the total count of compatible prior kmers
        multi_sum = np.sum(multi)
        ## if there are none, return at this point
        if multi_sum == 0:
            return 0, seq, np.zeros(shape=len(seq), dtype=bool)
        ## otherwise, compute consensus sequence from all those maps
        base_counts, intercept = get_base_counts(cseqs, cmaps, ALL_FORWARD = True, fix_zero=True, max_pos=len(seq), multi=multi)
        base_coverage = np.sum(base_counts, axis=1)
        base_ratio    = base_counts.astype(float) / np.maximum(base_coverage[:, np.newaxis], 1)
        max_ratio     = np.max(base_ratio, axis=1)    
        bs0_seq = np.argmax(base_ratio, axis=1)
        bs1_seq = np.array(seq_to_int(seq))
        ## if the pre-image exists and is (pretty much) unambiguous, we use that base
        use_bs0_base = (max_ratio >=  MIN_MAX_RATIO) * (base_coverage > 0)
        ## make a copy of the sequence
        corrected = bs1_seq.copy()
        corrected [use_bs0_base] = bs0_seq[use_bs0_base]
        corrected_seq = "".join(int_to_seq(corrected))
        return multi_sum, corrected_seq, use_bs0_base
    
    total_bits  = 0
    total_flips = 0
    N = len(contigs_bs1)
    unmask_info = []
    seqs_initial   = []
    seqs_corrected = []
    seqs_reversed  = []
    seqs_mask      = []
    if VERBOSE:
        headings = ['index', 'rc', 'length', 'diff', 'fcount', 'rcount', 'prop_off', 'bits', 'flips', "flip_rate", 'is_mutated']
        print tabprint(headings)
    ## for each sequences
    for cindex in range(N):
        ## get seq and reverse complement
        extended_seq    = contigs_bs1[cindex]
        extended_seq_rc = reverse_complement(extended_seq)
        f, seqf, maskf = unmask_sequence(extended_seq   )
        r, seqr, maskr = unmask_sequence(extended_seq_rc)
        if f >= r:
            si = extended_seq
            sc = seqf
            sr = 0
            sm = maskf
        else:
            si = extended_seq_rc
            sc = seqr
            sr = 1
            sm = maskr
        diff  = hamming(si, sc)
        bits  = Counter(sc)["C"]
        flips = Counter(zip(sc, si))[("C", "T")]
        bit_rate = float(flips) / float(max(bits, 1))
        is_mutated = bit_rate >= MUT_RATE_MIN
        is_rc = int(f < r)
        ratio_in_favor = min(f, r)/max(1, float(f+r))
        info = [cindex, is_rc, len(extended_seq), diff, f, r, ratio_in_favor, bits, flips, bit_rate, int(is_mutated)]
        if VERBOSE:
            print tabprint(info)
            sys.stdout.flush()
        if is_mutated:
            total_bits += bits
            total_flips += flips
            seqs_initial.append(extended_seq)
            seqs_corrected.append(sc)
            seqs_reversed.append(sr)
            seqs_mask.append(sm)
            unmask_info.append(info)
        
    unmasked_contig_file = file(unmasked_filename, 'w')
    mask_file            = file(mask_filename, 'w')
    info_file            = file(logout_filename, 'w')
    
    if VERBOSE:
        print "writing %d unmasked contigs (and masks) to file" % len(seqs_corrected)
    for ind, seq in enumerate(seqs_corrected):
        ## also record if the sequence was reversed for unmasking
        unmasked_contig_file.write(">%d_%d\n" % (ind, seqs_reversed[ind]))
        unmasked_contig_file.write(seq)
        unmasked_contig_file.write("\n")
        ##and record the mask used (which positions were corrected)
        mask_file.write(">%d_%d\n" % (ind, seqs_reversed[ind]))
        mask_string = "".join( map(str, seqs_mask[ind].astype(int)) )
        mask_file.write(mask_string )
        mask_file.write("\n")
    unmasked_contig_file.close()
    mask_file.close()
    ## write out the info file
    info_file.write(tabprint(headings) + "\n")
    ## remove the previous index and put in the new one
    for ind, entry in enumerate( unmask_info ):
        info_file.write(tabprint([ind] + entry[1:]) + "\n")
    info_file.close()
    ## finally, write out the templates
    contig_names = map(str, range(len(seqs_initial)))
    write_contigs(seqs_initial, bs_1_dir, output_contigs_shortname, contig_names)
    return float(total_flips) / float(total_bits)

def plot_base_array(base_array, ax, show_letters=False, fontsize=12, cmap=None):
    if cmap == None:
        cmap = mpl.colors.ListedColormap(DEFAULT_COLOR_LIST)
    cax = ax.imshow(base_array, interpolation="nearest", aspect="auto", vmin=-1, vmax=4, cmap=cmap)
    ## and if we show letters
    if show_letters:
        ## go through the array
        for x in range(base_array.shape[0]):
            for y in range(base_array.shape[1]):
                base_new  = base_array[x,y]
                ## if base is not a gap, add it to the axis
                if base_new != -1:
                    text = BASES[base_new]
                    ax.text(y - 0.25, x + 0.25, text, color="white", fontweight="bold", fontsize=fontsize)
    return cax

def add_base_legend(cax, fig):
    cb1 = fig.colorbar(cax, fraction=0.02, pad=0.01)
    ## to align the ticks properly
    half_step = 5/12.
    tick_locs = np.linspace(-1, 4, 7)[:-1] + half_step
    cb1.set_ticks(tick_locs)
    cb1.set_ticklabels(["no base", "A", "C", "G", "T", "N"])
    ticklabs = cb1.ax.get_yticklabels()
    cb1.ax.set_yticklabels(ticklabs, fontsize=20)

'''
Use Needleman-Wunsch algorithm (hereafter NW) to align seq1 and seq2.

NW is a dynamic programming method that identifies a sequence alignment with the best score.
Score depends on base agreement (score matrix) and gap costs (start and extension costs).

Establishes two (N1+1) x (N2+1) matrices, one for scores and one for moves (trace).
Use NW algorithm (written in C as a numpy friendly plugin) to fill in the matrices.
Modifications to the trace matrix allow a free gap at the beginning and end of the sequence.

Returns the score of the alignment and the two aligned sequences (gaps marked with "-").
Has an option (return_aligned_bases) to also return two numpy integer vectors:
the first  vector is the same length as seq1 and has the (integer) matched bases in seq2 (where applicable)
the second vector is the same length as seq2 and has the (integer) matched bases in seq1 (where applicable)

'''
def align_pair(seq1, seq2, score_matrix=None, gap_start=-50, gap_extend=-2, return_aligned_bases=False):
    ## if there is no score matrix, then use the default
    ## -1 for mismatch +1 for match
    try: 
        score_matrix[0]
    except:        
        score_matrix = -np.ones(shape=(4,4), dtype=int)
        score_matrix[np.arange(4), np.arange(4)] = 1
    ## put gap costs into array
    gap_costs = np.array([gap_start, gap_extend])
    ## get sequence lengths
    N1 = len(seq1)
    N2 = len(seq2)
    ## convert to integer array
    S1 = np.array( seq_to_int(seq1) )
    S2 = np.array( seq_to_int(seq2) )
    ## constant assignments for NW algorithm
    LEFT = 0
    DIAG = 1
    UP   = 2
    ## initialize NW matrix for trace and score
    nw_trace = np.zeros(shape=(N1 + 1, N2 + 1), dtype=int)
    nw_score = np.zeros(shape=(N1 + 1, N2 + 1), dtype=int)
    ## allows for a free gap in the beginning of the sequence
    nw_trace[0, :] = LEFT
    nw_trace[:, 0] = UP
    nw_trace[0, 0] = -1
    ## fill in values for nw_score and nw_trace using the C needleman code
    ## is WAY FASTER than using python alone
    _ = AlignmentExt.needleman(S1, S2, 
                               score_matrix, 
                               gap_costs,
                               nw_score, 
                               nw_trace)    
    ## allow for a free gap at the end of the sequence
    ## creates a corridor from the bottom right corridor (start position)
    ## to the best answer on the bottom or right
    x1_bottom, x2_bottom, bottom_score = N1, np.argmax(nw_score[N1, :]), np.max(nw_score[N1, :])
    x1_right , x2_right, right_score   = np.argmax(nw_score[:, N2]), N2, np.max(nw_score[:, N2])    
    if bottom_score >= right_score:
        nw_trace[x1_bottom, (x2_bottom+1):] = LEFT
        total_score = bottom_score
    else:
        nw_trace[(x1_right+1):, x2_right] = UP
        total_score = right_score    
    ## now for the trace-back
    ## also storing aligned bases as integers (if requested)
    if return_aligned_bases:
        aligned_base_1 = -np.ones(N1, dtype=int)  ## s2 bases that align against s1 
        aligned_base_2 = -np.ones(N2, dtype=int)  ## s1 bases that align against s2    
    ## align stores the sequence of bases and gaps visited for seq1 and seq2
    align1, align2 = [], []
    ## x1 and x2 are the present X and Y coordinate in the traceback
    ## initializes to the lower right corner
    x1, x2         = N1, N2
    ## trace_value tells you where to go next (LEFT, DIAG, UP)
    ## where -1 means DONE (and happens at the 0, 0 coordinate)
    trace_value    = nw_trace[x1, x2]    
    ## while we are not done (upper left corner)
    while trace_value != -1:
        ## if the best move is diag, that's a match
        if trace_value == DIAG:
            x1 -= 1 ## move both
            x2 -= 1
            if return_aligned_bases:
                ## store aligned bases
                aligned_base_1[x1] = BASE2INT[seq2[x2]]
                aligned_base_2[x2] = BASE2INT[seq1[x1]]
            align1.append(seq1[x1])
            align2.append(seq2[x2])
        elif trace_value == LEFT:
            x2 -= 1 
            align1.append("-")
            align2.append(seq2[x2])
        elif trace_value == UP:
            x1 -= 1  
            align1.append(seq1[x1])
            align2.append("-")
        ## update trace value
        trace_value = nw_trace[x1, x2]
    seq_align1 = "".join(align1[::-1])
    seq_align2 = "".join(align2[::-1])
    if return_aligned_bases:
        return total_score, seq_align1, seq_align2, aligned_base_1, aligned_base_2
    else:
        return total_score, seq_align1, seq_align2

'''
String representation showing the matching/mismatching of two aligned sequences
'''
def match_strings(seq1, seq2):
    match1 = "".join(["_" if a == b else ("=" if b == "-" else a) for (a,b) in zip(seq1, seq2)])
    match2 = "".join(["_" if a == b else ("=" if b == "-" else a) for (a,b) in zip(seq2, seq1)])
    return (match1, match2)


'''
For each seq2 in seq_list:
1) align seq2                       to seq1
2) align reverse_complement of seq2 to seq1
3) determine which one matches more bases (reverse =? True)
4) record the aligned bases to seq1 using the best of seq2 or seq2_rc

Returns a matrix of aligned bases to seq1
and a vector recording choices to reverse
'''
def get_aligned_bases(seq1, seq_list, VERBOSE=True):
    L    = len(seq1)
    ## how many sequences will we be considering
    M = len(seq_list)
    ## and we make a matrix to store the alignment bases
    ## and whether we match to the reverse complement
    align_seq = np.zeros(shape=(M, L), dtype=int)
    reverse   = np.zeros(shape=(M), dtype=bool)
    ## go through the sequences in seq_list
    for m, seq2 in enumerate(seq_list):
        ## take reverse complement of seq2
        seq2_rc = reverse_complement(seq2)
        ## align and score the forward and reverse matches and store the aligned bases
        total_score   , seq_align1   , seq_align2   , aligned_base_1   , aligned_base_2    = align_pair(seq1, seq2, return_aligned_bases=True)
        total_score_rc, seq_align1_rc, seq_align2_rc, aligned_base_1_rc, aligned_base_2_rc = align_pair(seq1, seq2_rc, return_aligned_bases=True)
        ## compute the number of unmatched base counts in forward and reverse
        forward_blanks = np.sum(aligned_base_1    == -1)
        reverse_blanks = np.sum(aligned_base_1_rc == -1)
        ## who is better? and store that value
        reverse  [m]    = reverse_blanks < forward_blanks
        align_seq[m]    = aligned_base_1_rc if reverse[m] else aligned_base_1
        if VERBOSE:
            if reverse[m]:
                indels = seq1_aligned_indels(seq1, seq_align1_rc, seq_align2_rc, aligned_base_1_rc)
#                 print indels
#                 print "\n".join(match_strings(seq_align1_rc, seq_align2_rc))
            else:
                indels = seq1_aligned_indels(seq1, seq_align1, seq_align2, aligned_base_1)
            info = ["%d of %d" % (m, M), 
                    "seq2 len = %d" % len(seq2), 
                    "rc = %s" % (str(reverse[m])), 
                    indels]
            print tabprint(info)
            sys.stdout.flush()
#                 print indels
#                 print "\n".join(match_strings(seq_align1, seq_align2))
    return align_seq, reverse


'''
Use mutation model (C->T for forward, G->A for reverse) to compute emission probabilities
'''
def compute_emission_probabilities(X, reverse, err, mut):
    M, N = X.shape
    ## count up ACGT base in X
    ## (NOTE, because -1 indicates nothing there and bincount does not like negatives,
    ##        we add one and then ignore the first row.)
    ## has dimensions bases (B) by positions (N)
    base_counts_BN = np.apply_along_axis(np.bincount, 0, X+1, minlength=5)[1:]
    ## take the top two entries (those are the bases with the most coverage)
    top2_counts_BN = np.argsort(base_counts_BN, axis=0)[-2:]
    ## sort those alphabetically and set the A allele to lower and the B allele to the higher.
    allelesAB = np.sort(top2_counts_BN, axis=0)
    # use_emitCT = (A_allele == 1) * (B_allele == 3)
    emit_CT = (allelesAB[0] == 1) * (allelesAB[1] == 3)
    emit_AG = (allelesAB[0] == 0) * (allelesAB[1] == 2)
    ## find emits for each ground truth
    ## value for m, n, i is the emission probability
    ## for observing X[m, n] for sample M, position N 
    ## if the true allele is A (i=0), B (i=1).
    emit_MN2 = np.ones(shape=(M, N, 2))
    for n in range(N):
        A, B = allelesAB[:, n]
        for m in range(M):
            obs = X[m, n]
            ## special cases are:
            # 1) reverse is true and  emit_AG (then G-> ) [B allele] is a special case
            # 2) reverse is false and emit_CT (then C-> ) [A allele] is a special case
            if reverse[m] and emit_AG[n]:
                ## if truth is A (base = A), standard error model
                emit_MN2[m, n, 0] = (obs == A)*(1 - err) + (obs == B)*err
                ## if truth is B (base = G), mutation error model 
                ## NOTE: since truth is B, obs==B is paired with not mutated in contrast to next NOTE.
                emit_MN2[m, n, 1] = (obs == B)*(1 - mut) + (obs == A)*mut
            elif not reverse[m] and emit_CT[n]:
                ## if truth is A (base = C), mutation error model
                ## NOTE: since truth is A, obs==A is paired with not mutated (see above note).
                emit_MN2[m, n, 0] = (obs == A)*(1 - mut) + (obs == B)*mut
                ## if truth is B (base = T), standard error model
                emit_MN2[m, n, 1] = (obs == B)*(1 - err) + (obs == A)*err
            else:
                ## only error models
                emit_MN2[m, n, 0] = (obs == A)*(1 - err) + (obs == B)*err
                emit_MN2[m, n, 1] = (obs == B)*(1 - err) + (obs == A)*err
    ## set unset entries to 1.
    emit_MN2[emit_MN2 == 0] = 1
    return emit_MN2, allelesAB

'''
Exhaustive methnod for finding the best pair of haplotypes.
Given an emission matrix M (samples) by N (positions) by 2 (A and B allele).
First , determines the probability of observing sample m from each possible haplotype.
Second, for each pair of haplotypes, compute the probability of observing all the samples assuming 
        each one derived from the most likely haplotype of the pair.
Third, find the pair that maximizes the total probability.
Report those two haplotypes and then the sample emission (log) probabilities for each.
'''
def find_best_haplopair(emit_MN2):
    log_emit = np.log(emit_MN2)
    M, N, _ = emit_MN2.shape
    ## we define a haplotype as a series of binary choices of alleles
    ## equivalently, a zero-one vector of length N.
    ## we use itertools to generate all of the haplotypes
    haplotypes  = np.array([x for x in product([0, 1], repeat=N)])
    H           = len(haplotypes)
    ## and we keep a matrix, for every haplotype (H), probability of each observation (M)
    hapscore_MH = np.zeros(shape=(M, H), dtype=float)    
    ## the log probability that each observations derives from this haplotype
    for h, haplotype in enumerate(haplotypes):
        hapscore_MH[:, h] = np.sum(log_emit[:, np.arange(N), haplotype], axis=1)
    ## for every pair of haplotypes, compute their best common score
    haplo_pairs = -np.inf * np.ones(shape=(H, H), dtype=float)
    for i in range(H):
        for j in range(i, H):
            ## take the element-wise maximum of the i and j columns and sum
            haplo_pairs[i, j] = np.sum(np.maximum(hapscore_MH[:, i], hapscore_MH[:, j]))
    ## find a value that maximizes the likelihood.
    bestx, besty = np.where(haplo_pairs == np.max(haplo_pairs))
    ## call those H1 and H2
    h1_ind = bestx[0]
    h2_ind = besty[0]
    H1 = haplotypes[h1_ind]
    H2 = haplotypes[h2_ind]
    print "\n".join(map(str, [H1, H2]))
    h1_score = hapscore_MH[:, h1_ind]
    h2_score = hapscore_MH[:, h2_ind]
    return H1, H2, h1_score, h2_score


'''
Use simulated annealing to phase haplotypes.
Needed for determining emission probabilities:
1) het data, observations which are N contigs by P positions.
2) reverse, was the observation N on the reverse strand (G->A mutations)
3) mutation and error rates

Then for the simulated annealing:
1) how many haplotypes
2) most temperatures we are willing to consider (like max iterations)
3) how many moves to try at each temperature
4) temperature reduction parameter (T_next = alhpa * T_current)
5) if we haven't moved in MAX_STUCK temperatures, we are done.

Returns:
The solution with the best probability ever observed in the algorithm.

Haplotype structure (K haplotypes by P positions)
Haplotype probabilities per contig (N contigs by K haplotypes)
Total probability
Did the algorithm terminate before exhausting MAX_TEMP_STEPS?
'''
def apply_SA(het_data, reverse, ERR_RATE, MUT_RATE, NUMBER_OF_HAPLOTYPES, 
             MAX_TEMP_STEPS=10000, STEPS_PER_TEMP=1000, ALPHA=0.95, MAX_STUCK=100, VERBOSE=False):
    ## let het_data be the observation matrix
    N, P = het_data.shape
    ## compute emission probabilities based on observations, error rate and if the sequence was reversed.
    emit_NP2, allelesAB = compute_emission_probabilities(het_data, reverse, ERR_RATE, MUT_RATE)
    ## take the log of the emission
    logemit_NP2 = np.log(emit_NP2)
    ## select number of haplotypes
    K    = NUMBER_OF_HAPLOTYPES
    ## build random haplotypes
    H_KP = np.random.randint(0, 2, size=(K, P))
    ## and take a sum to compute P_NK = P(Contig n |Haplotype k)
    P_NK      = np.sum(logemit_NP2[:, np.arange(P), H_KP], axis=2)
    ## make a copy for swapping
    P_NK_swap = P_NK.copy()
    ## compute switch P from (0 to 1) or (1 to 0), cost per N
    move_NPfrom = logemit_NP2[:, :, ::-1] - logemit_NP2
    ## biggest positive change to a move: 
    ## if there is a negative change -x, opposite column in "from" is +x
    biggest_move = np.max(np.sum(np.maximum(move_NPfrom, 0), axis=0))
    current_temp = biggest_move
    prob         = np.sum(np.max(P_NK, axis=1))
    ## keep track of the best answer
    max_prob          = prob
    max_solution = H_KP.copy()
    ## keep track if you reach the termination conditions
    terminated = False
    ## track how many times we accept no moves at the present temperature
    zero_accept_count = 0
    for ind in range(MAX_TEMP_STEPS):
        ## how many moves are accepted at this temp
        accept_count = 0 
        for iterations in range(STEPS_PER_TEMP ):
            ## pick a random position and haplotype to perturb
            p, k = np.random.randint(P), np.random.randint(K)
            ## alter P_NK_swap
            P_NK_swap[:, k] += move_NPfrom[:, p, H_KP[k, p]]
            ## compute new probability
            new_prob = np.sum(np.max(P_NK_swap, axis=1))
            ## if the move is better or beats a random draw given the temperature...
            if (new_prob > prob) or ( np.exp((new_prob - prob) / current_temp) > np.random.random() ):
                accept_count += 1
                ## change the value in H_KP
                H_KP[k, p] = (H_KP[k, p] + 1) % 2
                ## update P_NK, values in column K
                P_NK[:, k] = P_NK_swap[:, k]
                ## update total probability
                prob = new_prob
                ## check if this is the best ever
                if prob > max_prob:
                    ## and if so, save the score
                    max_prob     = prob
                    max_solution = H_KP.copy()
            else:
                ## reset the values in swap
                P_NK_swap[:, k] = P_NK[:, k]
        ## now, check the accept count. if zero, increment zero_accept_count
        ## otherwise, reset zero_accept_count to zero.
        if accept_count == 0:
            zero_accept_count += 1
            if zero_accept_count >= MAX_STUCK:
                terminated = True
                break
        else:
            zero_accept_count = 0
        if VERBOSE:
            info = [ind, "%0.2f" % current_temp, "%0.2f" % prob, "%0.2f" % max_prob, accept_count, zero_accept_count]
            print tabprint(info)
            sys.stdout.flush()
        ## update the temperature
        current_temp = current_temp * ALPHA
    ## set final answer
    H_KP = max_solution
    P_NK = np.sum(logemit_NP2[:, np.arange(P), H_KP], axis=2)
    return H_KP, P_NK, max_prob, terminated


def get_base_plot(align_seq, show_letters=False):
    fig = plt.figure(figsize=(30, 20))
    ax = fig.add_subplot(1,1,1)
    cax = plot_base_array(align_seq, ax, show_letters=show_letters)
    add_base_legend(cax, fig)
    fig.tight_layout()
    return fig

## double figure showing base call summaries and ratios
## also which points are selected by the BR_CUTOFF
def double_ratio_plot(base_ratios, base_counts, max_br_filter=None, ADJ_BR_CUTOFF=None):
        L = len(base_ratios)
        try:
            max_br_filter[0]
        except:
            max_br_filter = np.zeros(L, dtype=bool)
        fig = plt.figure(figsize=(30, 20))
        ax1 = fig.add_subplot(2, 1, 1)
        ax2 = fig.add_subplot(2, 1, 2, sharex=ax1)
        ## plot sorted lines to get 1st, 2nd and so on...
        ax1.plot(np.sort(base_counts, axis=1), '-', color='gray')
        ax2.plot(np.sort(base_ratios, axis=1), '-', color='gray')
        ## for the first four bases plot their values
        base_colors = DEFAULT_COLOR_LIST[1:5]    
        for b, (BASE, color) in enumerate(zip(BASES[:4], base_colors)):
            ax1.plot(base_counts[:, b], '.', label=BASE, ms=20,color=color)
            ax2.plot(base_ratios[:, b], '.', label=BASE, ms=20,color=color)
        ax1.plot(np.arange(L)[max_br_filter], base_counts[max_br_filter], '.', mec='black', mfc="None", ms=21, mew=2)
        ax2.plot(np.arange(L)[max_br_filter], base_ratios[max_br_filter], '.', mec='black', mfc="None", ms=21, mew=2)
        if ADJ_BR_CUTOFF != None:
            ax2.axhline(ADJ_BR_CUTOFF, ls='--', color ='k')
        ax1.legend(fontsize=20, markerscale=2, loc="upper left")
        return fig
'''
Function for aligning and phasing unmasked contigs.
First trims all raw alignments by TRIM_ENDS bases.
This is important if there are wobble bases or tags in the primer sequence
Also if there is off target priming.
This trimming is not carried through to the final reported sequence, only for calling het loci
as needed to split haplotypes.
Uses NW algorithm to align all contigs exceeding MIN_ASSEMBLY_LENGTH bases against the longest contig.
That might not always be the best choice.
Can specify an alternative index if one desires (seq1_ind=0 by default).

Alignment is computed both the forward and reverse with the best alignment selected from the pair.
After alignment, positions that are observed as hets (subject to a BASE_RATIO_CUTOFF)
and making sure that there aren't too many (MAX_POLY_COUNT) because the algorithms are exhaustive.
Then based on het positions, find haplotypes to explain all the sequences.
That part of the algorithm is sensitive to C->T and G->A conversions.
Can use Simulated annealing (USE_SA) or compute exhaustively.
Exhaustive solution only available for 2 haplotypes and MAX_POLY_COUNT (in practice, 10ish).
Simulated annealing can handle higher haplotype counts and many het loci.

Finally, assigns contigs to haplotypes if the confidence exceeds 1 - MIN_HAPLO_CONFIDENCE.
Unassigned contigs are saved to a separate file
Those assignments are written to .fa files that extend the contig_shortname with h1, h2, other.

If SAVE_PLOTS, then a series of plots will be generated and saved to plot_directory (which should not be None)
'''
def align_and_phase_contigs(contig_directory, contig_shortname, 
                            MIN_ASSEMBLY_LENGTH,
                            BASE_RATIO_CUTOFF,
                            MAX_POLY_COUNT,
                            MAX_SA_POLY_COUNT,
                            ERR_RATE,
                            MUT_RATE,
                            MIN_HAPLO_CONFIDENCE,
                            USE_SA=True,
                            TRIM_ENDS = 0,
                            MAX_TEMP_STEPS = 2000,
                            MOVES_PER_TEMP = 1000,
                            ALPHA = 0.99,
                            MAX_STUCK = 100,
                            NUMBER_OF_HAPLOTYPES = 2,
                            SAVE_PLOTS=False,
                            VERBOSE=False,
                            plot_directory=None,
                            seq1_ind=0):
    corrected_seqs = load_contigs  (contig_directory , contig_shortname)
    ## restrict sequence list based on MIN_ASSEMBLY_LENGTH    
    seq_list_init = [x for x in corrected_seqs if len(x) >= MIN_ASSEMBLY_LENGTH]
    ## and trim TRIM_ENDS if needed
    if TRIM_ENDS > 0:
        seq_list = [x[TRIM_ENDS:-TRIM_ENDS] for x in seq_list_init]
    else:
        seq_list = seq_list_init
    ## define target contig and record its length
    seq1 = seq_list[seq1_ind]
    if VERBOSE:
        print "getting aligned bases"
    align_seq, reverse = get_aligned_bases(seq1, seq_list, VERBOSE=VERBOSE)
    if VERBOSE:
        print "getting base counts from alignments"
    base_counts        = np.array([ np.sum(align_seq  == x, axis=0) for x in range(4)]).T
    ## add and compute ratio
    base_counts_sum    = np.sum(base_counts, axis=1)
    base_ratios        = base_counts.astype(float) / base_counts_sum[:, np.newaxis]
    ## select those locations where the max base ratio is below a pre-defined threshold (0.85?)
    ## these are the likely (but not necessarily) polymorphisms
    max_br        = np.max(base_ratios, axis=1)
    ## depending on if we are using Simulated annealing or not determines the MAX count and so the cutoff
    if USE_SA:
        ADJ_BR_CUTOFF = min(np.sort(max_br)[MAX_SA_POLY_COUNT], BASE_RATIO_CUTOFF)
    else:
        ADJ_BR_CUTOFF = min(np.sort(max_br)[MAX_POLY_COUNT], BASE_RATIO_CUTOFF)
    ## set filter
    max_br_filter = max_br < ADJ_BR_CUTOFF
    ## select just the het loci to make a new matrix
    align_hets_only = align_seq[:, max_br_filter]
    ## saves out the plot of the unsorted haplotypes
    if VERBOSE:
        print "testing polymorphism at %d loci" % np.sum(max_br_filter)
        print 'phasing to sort out % d haplotypes' % NUMBER_OF_HAPLOTYPES
        print "using simulated annealing" if USE_SA == 1 else "exhaustive search"
        
    if USE_SA:
        H_KP, P_NK, max_prob, terminated = apply_SA(align_hets_only, reverse, ERR_RATE, MUT_RATE, NUMBER_OF_HAPLOTYPES, 
                                                    MAX_TEMP_STEPS, MOVES_PER_TEMP, ALPHA, MAX_STUCK, VERBOSE=True)
    else:
        ## we are using an exhaustive algorithm so we can only do HAPLOTYPES = 2.
        NUMBER_OF_HAPLOTYPES = 2
        emit_MN2, allelesAB        = compute_emission_probabilities(align_hets_only, reverse, ERR_RATE, MUT_RATE)
        H1, H2, h1_score, h2_score = find_best_haplopair(emit_MN2)
        H_KP = np.array([H1, H2])                   ## haplotype by position
        P_NK = np.array([h1_score, h2_score]).T     ## contig by haplotype
        max_prob = np.sum(np.max(P_NK, axis=1))     ## total probability
        terminated = True
    if VERBOSE:
        for ind, haplo in enumerate(H_KP):
            print tabprint([ind, haplo])
    ## decide haplotype assignment
    assignment   = np.argmax(P_NK, axis=1)
    ## sort values, sum all and sum all but biggest
    ## that difference is the log_strength 
    P_NK_sort = np.sort(P_NK, axis=1)
    P_NK_sum  = np.logaddexp.reduce(P_NK, axis=1)
    P_NK_num  = np.sum(P_NK_sort[:, :-1], axis=1)
    ## (small is good...)
    log_strength = P_NK_num - P_NK_sum
    ## determine a threshold from MIN_HAPLO_CONFIDENCE parameter 
    log_thresh = np.log(MIN_HAPLO_CONFIDENCE)
    ## set assignment to -1 if less than log_threshold
    assignment[ log_strength > log_thresh ] = -1         
    ## write results to file
    first_name    = contig_shortname.split(".")[0]    
    outname_other = os.path.join(contig_directory, "%s_unassigned.fa" % first_name)
    outfile_other = file(outname_other, 'w')
    outfile = [] 
    for k in range(NUMBER_OF_HAPLOTYPES):
        outname = os.path.join(contig_directory, "%s_h%d.fa" % (first_name, k+1))
        outfile.append( file(outname, 'w'))
    ## write to output
    ## if the alignment requires a reverse complement, that is noted in the name (#, rc)
    ## and the sequence reported as reverse complemented as needed.
    ## NOTE: we are returning the full contigs before TRIM_ENDS
    for n, seq in enumerate(seq_list_init):
        seq_name = "%d_%d" % (n, reverse[n])
        if reverse[n]:
            seq = reverse_complement(seq)
        outline  = ">%s\n%s\n" % (seq_name, seq)
        if assignment[n] == -1:
            outfile_other.write(outline)
        else:
            outfile[assignment[n]].write(outline)
    ## close everything
    outfile_other.close()
    for ofile in outfile:
        ofile.close()
    ## if there are plots to save do it all here:
    if SAVE_PLOTS and plot_directory == None:
        print "no plot directory specified!"
        print "not saving plots!"
        SAVE_PLOTS = False
    if SAVE_PLOTS:
        if not os.path.exists(plot_directory):
            os.makedirs(plot_directory)
        unsorted_allpos_plotname    = os.path.join(plot_directory, "unsorted_all_pos.png"  )
        sorted_allpos_plotname    = os.path.join(plot_directory  , "sorted_all_pos.png"    )
        base_summary_plotname       = os.path.join(plot_directory, "base_summary.png"      )
        hets_only_unsorted_plotname = os.path.join(plot_directory, "unsorted_hets_only.png")
        hets_only_sorted_plotname   = os.path.join(plot_directory, "sorted_hets_only.png"  )
        ## first plot of all unsorted positions
        fig = get_base_plot(align_seq)
        plt.savefig(unsorted_allpos_plotname, dpi=200, facecolor='w', edgecolor='w',
            papertype=None, format=None, transparent=False)
        plt.close()
        ## second double plot showing selection of het positions
        fig = double_ratio_plot(base_ratios, base_counts, max_br_filter, ADJ_BR_CUTOFF)
        plt.savefig(base_summary_plotname, dpi=200, facecolor='w', edgecolor='w',
            papertype=None, format=None, transparent=False)
        plt.close()
        ## third plot showing unsorted calls at het positions
        letter_count = np.sum(align_hets_only >= 0)
        show_letters = letter_count <= 2000
        fig = get_base_plot(align_hets_only, show_letters=show_letters)
        ## save unaligned_plot
        plt.savefig(hets_only_unsorted_plotname, dpi=200, facecolor='w', edgecolor='w',
            papertype=None, format=None, transparent=False)
        plt.close()
        ## gets splits and labels for the plot
        order = np.lexsort([log_strength, reverse, assignment])
        assignment_ordered = assignment[order]
        reverse_ordered    = reverse[order]
        ylocs   = []
        ylabels = []
        main_split = []
        sub_split  = []
        for k in range(-1, NUMBER_OF_HAPLOTYPES):
            kmatch     = np.where(assignment_ordered == k)[0]
            kmatch_for = np.where((assignment_ordered == k) * ~reverse_ordered)[0]
            kmatch_rev = np.where((assignment_ordered == k) * reverse_ordered)[0]
            if len(kmatch_rev) > 0:
                if k == -1:
                    ylabels.append("unassigned\nreverse")
                else:
                    ylabels.append("haplo %d\nreverse" % (k+1))        
                ylocs.append( np.mean(kmatch_rev) )
            if len(kmatch_for) > 0:
                if k == -1:
                    ylabels.append("unassigned\nforward")
                else:
                    ylabels.append("haplo %d\nforward" % (k+1))        
                ylocs.append( np.mean(kmatch_for) )
            if len(kmatch) > 0:
                main_split.append(np.max(kmatch) + 0.5)
            if len(kmatch_rev) > 0 and len(kmatch_for) > 0:
                sub_split.append(min(np.max(kmatch_rev), np.max(kmatch_for)) + 0.5)
        ## make sorted het only plot
        fig = plt.figure(figsize=(30, 20))
        ax = fig.add_subplot(1,1,1)
        cax = plot_base_array(align_hets_only[order] , ax, show_letters=show_letters)
        add_base_legend(cax, fig)
        for yval in main_split:
            plt.axhline(yval, color="lime", lw=4)
        for yval in sub_split:
            plt.axhline(yval, color="lime", lw=2, ls='--')    
        ax.set_yticks(ylocs)
        ax.set_yticklabels(ylabels, fontsize=15)
        plt.tight_layout()
        plt.savefig(hets_only_sorted_plotname, dpi=200, facecolor='w', edgecolor='w',
                papertype=None, format=None, transparent=False)
        plt.close()
        ## make sorted all of everything plot
        fig = plt.figure(figsize=(30, 20))
        ax = fig.add_subplot(1,1,1)
        cax = plot_base_array(align_seq[order] , ax, show_letters=False)
        add_base_legend(cax, fig)
        for yval in main_split:
            plt.axhline(yval, color="lime", lw=4)
        for yval in sub_split:
            plt.axhline(yval, color="lime", lw=2, ls='--')    
        ax.set_yticks(ylocs)
        ax.set_yticklabels(ylabels, fontsize=15)
        plt.tight_layout()
        plt.savefig(sorted_allpos_plotname, dpi=200, facecolor='w', edgecolor='w',
                papertype=None, format=None, transparent=False)
        plt.close()


'''
Uses the results of the pairwise alignment method "align_pair" to decode the alignment gaps.
This is done with respect to seq1, that is as though seq1 were the reference.
Then gaps are interpreted as insertions or deletions into seq1.
Base position (relative to seq1), length of indel and content.

Returns a list of all indels denoted by a tuple with the following entries:
0: +1/-1 insertion or deletion
1: position in seq1
2: length of indel
3: bases deleted or inserted

'''
def seq1_aligned_indels(seq1, seq_align1, seq_align2, aligned_base_1):
    indels = []
    ## we are treating seq1 as though the reference.
    ## therefore:
    ##
    ## gap locations in seq1_align correspond to insertions
    ## we record the location of the insertions against seq1, their length and content 
    gap_list = path_summary(np.array([x == "-" for x in seq_align1]))
    ## only those positions which are gaps (==1)
    gap_list = gap_list[gap_list[:, 2] == 1]
    for entry in gap_list:
        start, end, _ = entry
        ## want to ignore starting and ending gaps
        if start == 0 or end == len(seq_align1) - 1:
            continue
        else:
            indels.append(( 1, start, end - start + 1, seq_align2[start:(end+1)]) )
    ## gaps locations in aligned_base_1 correspond to deletions relative to seq1
    ## NOTE: THIS IS NOT SYMMETRIC WITH THE ABOVE INSERTIONS
    ## to keep frame with seq1, we use aligned_base_1
    gap_list = path_summary(aligned_base_1 == -1)
    ## only those positions which are gaps (==1)
    gap_list = gap_list[gap_list[:, 2] == 1]
    for entry in gap_list:
        start, end, _ = entry
        ## want to ignore starting and ending gaps
        if start == 0 or end == len(aligned_base_1) - 1:
            continue
        else:
            indels.append((-1, start, end - start + 1, seq1[start:(end+1)]))
    return indels

'''
Returns the consensus sequence from a set of haplotypes
Uses the mask information to model error and mutation when computing the consensus
Along with the "best" answer, emission probabilities for each position are returned.
Watch the output in VERBOSE to note recurrent indels which may be a sign of a poor choice in the seq1_ind.
The method iterates MAX_INDEL_CLEARING_ITERATIONS until there are no common recurrent indels in the reference.
Each iteration will clear at most one erroneous indel at a time.
Also return a sequence of "phred" quality scores which summarizes a normalized version of the emission probabilities
if TRIM_ENDS > 0, will lop off the first and last TRIM_ENDS bases from each contig before assembly.
'''
def get_consensus(contig_directory, haplotype_shortname, contigs_mask_shortname, err, mut, seq1_ind=0, TRIM_ENDS = 0, VERBOSE=True, MAX_INDEL_CLEARING_ITERATIONS=20):
    ## load contig masks
    contigs_mask     = load_contigs(contig_directory, contigs_mask_shortname)
    ## load haplotype contigs and names
    haplo_contigs, names = load_contigs  (contig_directory , haplotype_shortname, return_names=True)
    if len(haplo_contigs) == 0:
        if VERBOSE:
            print "we got no contigs!"
        return None
    ## track if the contig was reversed (G->A conversions, not C->T)
    reverse        = np.array([int(x.split("_")[1]) for x in names], dtype=bool)
    ## get its initial_index in the mask list
    initial_ind    = np.array([int(x.split("_")[0]) for x in names], dtype=int)
    
    
    seq_list_init = haplo_contigs
    if VERBOSE:
        print "TRIMMING ENDS by %d bp on either side" % TRIM_ENDS
    if TRIM_ENDS > 0:
        seq_list = [x[TRIM_ENDS:-TRIM_ENDS] for x in seq_list_init]
    else:
        seq_list = seq_list_init
    
    ## load relevant sequence masks, reversing as needed
    seq_masks = []
    for ind, rev in zip(initial_ind, reverse):
        corrected = np.array(map(int, contigs_mask[ind]))
        if rev:
            corrected = corrected[::-1]
        ## if we trim ends, we have to take them off the masks too.
        if TRIM_ENDS > 0:
            corrected = corrected[TRIM_ENDS:-TRIM_ENDS]
        seq_masks.append(corrected)
    
    ## establish reference sequence
    seq1 = seq_list[seq1_ind]
    for trial_number in range(MAX_INDEL_CLEARING_ITERATIONS):
        aligned_bases      = []
        aligned_bases_mask = []
        all_indels         = []
        all_intervals      = []
        ## give it a try to resolve this sequence
        ## some work to get the masks aligned too.
        if VERBOSE:
            print "If you see many repeated indels here that is a problem."
            print "If they exactly match and are more than half the coverage"
            print "a correction will be made and this code will rerun"
            sys.stdout.flush()
        for ind, (seq2, seq2_mask) in enumerate(zip(seq_list, seq_masks)):
            ## align the pair
            total_score, seq_align1, seq_align2, aligned_base_1, aligned_base_2 = align_pair(seq1, seq2, return_aligned_bases=True)
            ## apply the alignment to the mask (to follow corrected bases)
            ## initialize the seq_array2m alignment to -1 everywhere (will denote gaps)
            seq_array2m  = -np.ones(len(seq_align2), dtype=int)
            seq2_covered = np.array([x != '-' for x in seq_align2], dtype=bool)
            ## at covered positions, substitute the seq2_mask for seq2_init (1 corrected, 0 not corrected)
            seq_array2m[seq2_covered] = seq2_mask
            ## along the alignment, where are the ungapped positions in seq1?
            seq1_covered = np.array([x != "-" for x in seq_align1], dtype=bool) 
            ## restricting to those gets the aligned_bases for the mask
            aligned_base_1m = seq_array2m[seq1_covered]
            ## get indels
            indels = seq1_aligned_indels(seq1, seq_align1, seq_align2, aligned_base_1)
            ## where is there a match
            has_match = np.where ( [x != '-' for x in seq_align2] )[0]
            ## identify start and end
            start, end = np.min(has_match), np.max(has_match)
            ## and save results
            if VERBOSE:
                print names[ind], indels
                sys.stdout.flush()
            aligned_bases.append(aligned_base_1)
            aligned_bases_mask.append(aligned_base_1m)
            all_indels.append(indels)
            all_intervals.append((start, end))
        interval_array = np.array(all_intervals)
        starts = interval_array[:, 0]
        ends   = interval_array[:, 1]
        ## test if we've got lots of the exact same indel
        ## might need to work on improving this for near-matches
        indel_counts = Counter(flatten(all_indels))
        useful_indels = []
        rates         = []
        for indel in indel_counts.keys():
            count = indel_counts[indel]
            updown, pos, length, seq = indel
            cover = np.sum( (pos > starts) * (pos < ends) )
            rate = float(count) / cover
            ## if more than half the sequences appear to have this event...
            if rate > 0.5:
                useful_indels.append(indel)
                rates.append(rate)
        ## if there's nothing in useful indels, we are done with this loop
        if len(useful_indels) == 0:
            break
        else:
            ## take the one with the best rate...
            top_index = np.argmax(rates)
            top_indel = useful_indels[top_index]
            if VERBOSE:
                print "Well, seems we had a very common indel:"
                print top_indel, rates[top_index]
                print "correcting and rerunning...."
                sys.stdout.flush()
            updown, pos, length, seq = top_indel
            ## and correct seq1
            if updown == -1:
                seq1 = seq1[:pos] + seq1[(pos+length):]
            else:
                seq1 = seq1[:pos] + seq + seq1[pos:]
    
    ## convert alignments into arrays
    align_seq      = np.array(aligned_bases)
    align_seq_mask = np.array(aligned_bases_mask)
    base_counts    = np.array([ np.sum(align_seq  == x, axis=0) for x in range(4)]).T
    
    ## EMISSION ERROR MODEL governed by error and emission.
    log_err    = np.log(err)
    log_no_err = np.log(1 - err)
    log_mut    = np.log(mut)
    log_no_mut = np.log(1 - mut)
    ## NOTE, in what follows, we shift the base assignments by 1 so that -1 (gap) = 0, 0 (A) = 1, etc.
    ## this is because -1 is not useful as an index
    ## -1 will not contribute to observation probabilities.
    ## Assumption is that indels are sparse.
    ## for those cases that are driven only by error
    p0_scores = np.array([0., log_no_err, log_err   , log_err   , log_err   ])
    p1_scores = np.array([0., log_err   , log_no_err, log_err   , log_err   ])
    p2_scores = np.array([0., log_err   , log_err   , log_no_err, log_err   ])
    p3_scores = np.array([0., log_err   , log_err   , log_err   , log_no_err])
    ## for those cases that are driven by mutation too
    p1_scores_mut = np.array([0., log_err, log_no_mut, log_err   , log_mut])
    p2_scores_mut = np.array([0., log_mut, log_err   , log_no_mut, log_err])
    ## there may be a more elegant solution but this works well enough
    N, L = align_seq.shape
    base_probs = []
    ## for each position, we note the observations and masks
    for pos in range(L):
        pdata       = align_seq     [:, pos]
        pdata_m     = align_seq_mask[:, pos]
        ## noting if the particular value was uncorrected
        uncorrected = pdata_m == 0
        ## look up all the values as needed
        p0  = p0_scores    [pdata + 1]
        p1  = p1_scores    [pdata + 1]
        p2  = p2_scores    [pdata + 1]
        p3  = p3_scores    [pdata + 1]
        p1m = p1_scores_mut[pdata + 1] 
        p2m = p2_scores_mut[pdata + 1]
        ## if it is the right type of mutation and the base is uncorrected
        ## use the mutation models
        use_mut1      = ~reverse * uncorrected
        use_mut2      = reverse* uncorrected
        ## replace by filter
        p1[use_mut1] = p1m[use_mut1]
        p2[use_mut2] = p2m[use_mut2]
        ## finally, sum to get the total (log) probability
        probs = np.sum(np.array([p0, p1, p2, p3]), axis=1)
        base_probs.append(probs)
    
    prob_array = np.array(base_probs)
    best_call = np.argmax(prob_array, axis=1)
    ## sort the probability array for calculating phred score
    prob_array_sort = np.sort(prob_array, axis=1)
    ## want the three smallest values, log_summed
    numer = np.logaddexp.reduce(prob_array_sort[:, :3], axis=1)
    ## compared to everything log_summed
    denom = np.logaddexp.reduce(prob_array_sort       , axis=1)
    log_score   = numer - denom
    ## convert to log 10
    log10_score = log_score * np.log10(np.e)
    ## convert to phred
    phred_score = (-10 * log10_score).astype(int)
    ## make final sequence from best_calls
    final_seq = "".join(int_to_seq(best_call))
    return final_seq, prob_array, base_counts, phred_score


## double figure showing base counts (coverage) and probabilities for the haplotype
def coverage_prob_plot(base_counts, prob_array):
    ## put into 10* log10 PHRED space
    prob_norm = 10* np.log10(np.e) * (prob_array - np.logaddexp.reduce(prob_array, axis=1)[:, np.newaxis])
    fig = plt.figure(figsize=(30, 20))
    ax1 = fig.add_subplot(2, 1, 1)
    ax2 = fig.add_subplot(2, 1, 2, sharex=ax1)
    ## plot sorted lines to get 1st, 2nd and so on...
    ax1.plot(np.sort(base_counts, axis=1), '-', color='gray')
    ax2.plot(np.sort(prob_norm, axis=1), '-', color='gray')
    ## for the first four bases plot their values
    base_colors = DEFAULT_COLOR_LIST[1:5]    
    for b, (BASE, color) in enumerate(zip(BASES[:4], base_colors)):
        ax1.plot(base_counts[:, b], '.', label=BASE, ms=20,color=color)
        ax2.plot(prob_norm[:, b], '.', label=BASE, ms=20,color=color)
    ax1.legend(fontsize=20, markerscale=2, loc="upper left")
    return fig


'''
Does the final big of aligning and taking a consensus of each subset of haplotype strings.
Is strand-aware and mutation aware, leveraging information from the BS0 (unmutated) data to model error
Also records information about the probability summarized as a phred score and recorded as an integer string (PHRED+33)
Most of the heavy lifting is in get_consensus.
Really just repeats that for all haplotypes and saves the answer to file.
If TRIM_ENDS > 0, then it will lop off the first and last TRIM_ENDS base pairs from each template before final assembly.
'''
def resolve_haplotypes(short_name, pair_dir, contig_directory, contigs_mask_shortname, haplotype_shortname_pattern, plot_directory, ERR_RATE, MUT_RATE, TRIM_ENDS=0, MAX_PHRED=41, NUMBER_OF_HAPLOTYPES=2, SAVE_PLOTS=False, VERBOSE=False):
    final_haplotype_fastq = os.path.join(pair_dir, "haplotype_assemblies.fastq")
    if SAVE_PLOTS and plot_directory == None:
            print "no plot directory specified!"
            print "not saving plots!"
            SAVE_PLOTS = False
    elif SAVE_PLOTS:
        haplotype_probs_plotname_pattern = os.path.join(plot_directory, "haplotype_%d_base_prob.png")
        if not os.path.exists(plot_directory):
            os.makedirs(plot_directory)
    
    haplo_file = file(final_haplotype_fastq, 'w')
    for hindex in range(NUMBER_OF_HAPLOTYPES):
        if VERBOSE:
            print "resolving haplotype %d" % (hindex + 1)
        haplotype_shortname = haplotype_shortname_pattern % (hindex + 1)
        ## we write the answer to ans first in case get_consensus comes up empty
        ans = get_consensus(contig_directory, haplotype_shortname, contigs_mask_shortname, ERR_RATE, MUT_RATE, TRIM_ENDS=TRIM_ENDS, seq1_ind=0, VERBOSE=VERBOSE)
        if ans == None:
            continue
        else:
            final_seq, prob_array, base_counts, phred_int = ans
        phred_int_bound = np.minimum(phred_int, MAX_PHRED)
        phred_seq       = "".join([chr(x + 33) for x in phred_int_bound])
        seq_name        = "@%s_h%d" % (short_name, hindex+1)
        haplo_string    = "\n".join([seq_name, final_seq, "+", phred_seq])
        haplo_file.write(haplo_string + "\n")
        ## save a plot if needed
        if SAVE_PLOTS:
            fig = coverage_prob_plot(base_counts, prob_array)
            fig.suptitle("Haplotype %d consensus probability" % (hindex + 1), fontsize=20)
            plt.savefig(haplotype_probs_plotname_pattern % (hindex + 1), dpi=200, facecolor='w', edgecolor='w',
                papertype=None, format=None, transparent=False)
            plt.close()
    haplo_file.close()


'''
Function loads: 
    a set of contigs
    a set of reads
    a set of annotated maps from those reads to those contigs
Then, picking the best maps (good ratio and overlap length), find well mapped read pairs,
assign them to the contig and examine the consensus:
First, as a plot (if SAVE_PLOT)
Second, as statistics for length of contig, number of read pairs with HQ maps and median coverage over contig.
Aggregate statistics are reported.
'''
def establish_coverage_profile(outdir, contig_shortname, maps_pickle, 
                               MIN_MAP_RATIO, MIN_MAP_LENGTH, 
                               VERBOSE=False, SAVE_PLOTS=False, 
                               plot_directory=None, plot_sub_directory=None):
    contigs = load_contigs  (outdir, contig_shortname)
    reads   = load_readpairs(outdir)
    maps    = load_maps     (outdir, maps_pickle)
    ## establish filters
    map_ratio      = maps[:, MATCH] / maps[:, OVERLAP].astype(float)
    ratio_filter   = map_ratio   >= MIN_MAP_RATIO
    length_filter  = maps[:, OVERLAP] >= MIN_MAP_LENGTH   
    quality_filter = ratio_filter * length_filter
    ## get best maps
    top_quality_maps = maps[quality_filter]
    good_ratio_maps  = maps[ratio_filter]
    ## we are not using the extentions for mapping here so we make an empty DB extension.
    DB_ext = dict()
    alignments = well_mapped_alignments(top_quality_maps, good_ratio_maps, DB_ext, VERBOSE=VERBOSE, BOTH_MUST_ALIGN=True)
    ## re-order the alignemnts to reflect contig assignment
    cind = [align[1] for align in alignments]
    rind = [align[0] for align in alignments]
    ## sorting first by contig index and then by read index
    order = np.lexsort([rind, cind])
    ## then re-order alignments
    alignments = [alignments[ind] for ind in order]    
    ## make plots and generate coverage data
    if plot_directory == None:
        SAVE_PLOTS = False
    if plot_sub_directory == None:
        plot_sub_directory = "extended_plots"
    if SAVE_PLOTS:
        extended_plot_directory = os.path.join(plot_directory, plot_sub_directory)
        extended_contig_plot_pattern = os.path.join(extended_plot_directory, "extended_contig_%d.png")
        if not os.path.exists(extended_plot_directory):
                    os.makedirs(extended_plot_directory)    
    contig_lens = np.array( [len(x) for x in contigs] )
    contig_info = []
    for cindex, caligns in block_iterator(alignments, 1):
        ## how many read pairs?
        num_reads = len(caligns)
        ## otherwise convert paired reads to reads and maps
        rseqs, rmaps = read_alignments_seq_and_map(caligns, reads, FLIP_FORWARD=True)
        ## make that into a base count array
        base_count, contig_start  = get_base_counts(rseqs, rmaps, ALL_FORWARD = True)
        ## restrict base counts to get rid of N
        base_counts   = base_count[:, :4]
        ## figure out where the contig starts and ends
        contig_end = contig_start + len(contigs[cindex])
        ## calculate hte total base coverage and get ratio
        base_coverage = np.sum(base_count, axis=1) 
        base_ratios   = base_counts.astype(float)/ np.maximum(base_coverage[:, np.newaxis], 1)
        ## want to also compute coverage restricted to the contig
        contig_coverage = np.zeros(contig_lens[cindex], dtype=int)
        ## figure out where to start and end the "base_coverage" info that gets added to contig_coverage
        low_index  = max(0, contig_start)
        high_index = min(contig_end, len(base_coverage))
        cov_start  = max(0, -contig_start) 
        cov_end    = cov_start + high_index - low_index
        ## update contig coverage
        contig_coverage[cov_start:cov_end] = base_coverage[low_index:high_index]
        ## compute median
        mcov  = np.median(contig_coverage)
        ## record the relevant information
        info = [cindex, contig_lens[cindex], num_reads, mcov]
        contig_info .append(info)
        if VERBOSE:
            print tabprint(info)
        if SAVE_PLOTS:
            ## make the plot to output
            fig = plt.figure(figsize=(30, 20))
            fig.suptitle("Extended contig %d" % cindex, fontsize=25)
            ax1 = fig.add_subplot(2, 1, 1)
            ax2 = fig.add_subplot(2, 1, 2, sharex=ax1)
            ax1.set_title("coverage", fontsize=20)
            ax2.set_title("ratio", fontsize = 20)
            ax1.axvline(contig_start - 0.5, color='red', lw=2)
            ax2.axvline(contig_start - 0.5, color='red', lw=2)
            ax1.axvline(contig_end - 0.5, color='red', lw=2)
            ax2.axvline(contig_end - 0.5, color='red', lw=2)
            ## plot sorted lines to get 1st, 2nd and so on...
            ax1.plot(np.sort(base_counts, axis=1), '-', color='gray')
            ax2.plot(np.sort(base_ratios, axis=1), '-', color='gray')
            ## for the first four bases plot their values
            base_colors = DEFAULT_COLOR_LIST[1:5]    
            for b, (BASE, color) in enumerate(zip(BASES[:4], base_colors)):
                ax1.plot(base_counts[:, b], '.', label=BASE, ms=20,color=color)
                ax2.plot(base_ratios[:, b], '.', label=BASE, ms=20,color=color)
            ax1.plot(np.arange(contig_start, contig_end), contig_coverage, '-k', lw=2)
            ax1.axhline(mcov , linestyle='-' , color='black', lw=2)
            ax1.legend(fontsize=20, markerscale=2, loc="upper left")
            plt.savefig(extended_contig_plot_pattern % (cindex), dpi=200, facecolor='w', edgecolor='w',
                            papertype=None, format=None, transparent=False)
            plt.close()
    return contig_info

'''
Load the template file and establishes some important information such as:
1. number of conservative read pairs mapped to a contig
2. median base coverage over the contig
3. convertable positions (approximately, up to unmasking)
4. corrected positions (apprxoimately # of mutations)
5. assignment to haplotype 1 or 2 (or more) 0 if unassigned
6. was the sequence reversed to coalign with the rest of the haplotype?
Can also save the template coverage plots if interested.
'''
def write_template_info(pair_dir, bs_1_dir, template_shortname, maps_pickle, haplotype_shortname_pattern, UNMASK_INFO_TXT, outtable_shortname,
                        NUMBER_OF_HAPLOTYPES, MIN_MAP_RATIO, MIN_MAP_LENGTH,
                        VERBOSE=False, SAVE_PLOTS=False,plot_directory=None, plot_sub_directory=None):
    ## load the contig information based on well mapped reads to the contigs
    contig_info = establish_coverage_profile(bs_1_dir, template_shortname, maps_pickle,
                               MIN_MAP_RATIO, MIN_MAP_LENGTH, 
                               VERBOSE=True, SAVE_PLOTS=SAVE_PLOTS,plot_directory=plot_directory, plot_sub_directory=plot_sub_directory)
    ## now get the haplotype assignments
    contigs = load_contigs(bs_1_dir, template_shortname)
    N = len(contigs)
    ## pull the haplotypes
    final_haplotype_fastq = os.path.join(pair_dir, "haplotype_assemblies.fastq")
    infile = file(final_haplotype_fastq, 'r')
    haplotypes = []
    for ind, line in enumerate(infile):
        if ind % 4 == 1:
            haplotypes.append(line.strip())
    infile.close()
    ## establish haplotype length (haplo 0 is unassigned and so length 0)
    hlen = [0] + [len(x) for x in haplotypes]
    contig_assignment = np.zeros(N, dtype=int)
    contig_reverse    = np.zeros(N, dtype=int)
    contig_start      = np.zeros(N, dtype=int)
    contig_end        = np.zeros(N, dtype=int) 
    ## for each haplotype
    for hindex in range(NUMBER_OF_HAPLOTYPES):
        haplotype_shortname = haplotype_shortname_pattern % (hindex + 1)
        ## load haplotype contigs and names
        haplo_contigs, names = load_contigs  (bs_1_dir , haplotype_shortname, return_names=True)
        ## track if the contig was reversed (G->A conversions, not C->T)
        reverse_name     = np.array([int(x.split("_")[1]) for x in names], dtype=bool)
        reverse = np.zeros(len(names))
        ## get its initial_index in the mask list
        initial_ind    = np.array([int(x.split("_")[0]) for x in names], dtype=int)
        seq1 = haplotypes[hindex]
        for m, seq2 in enumerate(haplo_contigs):
            ## take reverse complement of seq2
            seq2_rc = reverse_complement(seq2)
            ## align and score the forward and reverse matches and store the aligned bases
            total_score   , seq_align1   , seq_align2   , aligned_base_1   , aligned_base_2    = align_pair(seq1, seq2   , return_aligned_bases=True)
            total_score_rc, seq_align1_rc, seq_align2_rc, aligned_base_1_rc, aligned_base_2_rc = align_pair(seq1, seq2_rc, return_aligned_bases=True)
            ## compute the number of unmatched base counts in forward and reverse
            forward_blanks = np.sum(aligned_base_1    == -1)
            reverse_blanks = np.sum(aligned_base_1_rc == -1)
            ## who is better? and store that value
            reverse  [m]    = reverse_blanks < forward_blanks
            if reverse[m]:
                indels = seq1_aligned_indels(seq1, seq_align1_rc, seq_align2_rc, aligned_base_1_rc)
                has_match = np.where ( [x != '-' for x in seq_align2_rc] )[0]
            else:
                indels = seq1_aligned_indels(seq1, seq_align1, seq_align2, aligned_base_1)
                has_match = np.where ( [x != '-' for x in seq_align2] )[0]
            start, end = np.min(has_match), np.max(has_match)
            info = [initial_ind[m], int(reverse[m]), hindex+1, start, end, indels]
            ind = initial_ind[m]
            contig_assignment[ind] = hindex + 1
            contig_reverse   [ind] = int( reverse_name[m] )
            contig_start     [ind] = start
            contig_end       [ind] = end
            print tabprint(info)
            sys.stdout.flush()
    ## now get the parts to do with conversion (mostly)
    unmask_info_filename = os.path.join(bs_1_dir, UNMASK_INFO_TXT)
    unmask_info = file2RA(unmask_info_filename)
    ## and put it all together
    headings = ['cindex', 'length', 'pair_count', 'median_coverage', 'bits', 'flips', 'haplotype', 'reverse', 'start', 'end', 'haplotype_len']
    outtable_filename  = os.path.join(bs_1_dir, outtable_shortname)
    outtable = file(outtable_filename, 'w')
    outtable.write(tabprint(headings) + "\n")
    for cindex in range(N):
        info  = contig_info[cindex] + [unmask_info['bits' ][cindex], unmask_info['flips'][cindex]]
        info += [contig_assignment[cindex],contig_reverse[cindex],contig_start[cindex], contig_end[cindex], hlen[contig_assignment[cindex]]]
        outtable.write(tabprint(info) + "\n")
    outtable.close()

'''
Given the directory wherein live the bs_1 data
and the unmask-info-txt filename, load the data and compute average mutation rate.
'''
def get_mutation_rate(unmask_info_filename):
    unmask_info = file2RA(unmask_info_filename)
    bits = unmask_info['bits']
    flips = unmask_info['flips']
    return np.sum(flips).astype(float) / np.sum(bits)


## checks if lookup[node][tail] has a value, and if it does, does it return to here?
## if either is no, return None, None.
## otherwise, return the next_node, and whether the match is to the head=0 or tail=1.
def get_strong_match(node, tail, lookup):
    next_node = lookup[node][tail]
    if next_node == None:
        return None, None
    try:
        next_tail = lookup[next_node].index(node)
    except:
        return None, None
    return next_node, next_tail

##
## Some code to build on the clique finding code taken from Bron-Kerbosch implementation below.
## builds the dictionary graph based on a list of nodes and legal pairs from the compatible pair data.
##
def clique_set(nodes, pairs):
    graph = dict()
    for ind in nodes:
        graph[ind] = list()
    ordered = np.sort(nodes)
    for opair in itertools.combinations(ordered, 2):
        if opair in pairs:
            a, b = opair
            graph[a].append(b)
            graph[b].append(a)
    cliques = find_cliques(graph)
    return cliques

##
## Some code to build on the clique finding code taken from Bron-Kerbosch implementation below.
## builds the dictionary graph based on a list of nodes and legal pairs from the compatible pair data.
## node_filter is a dictionary that returns True/False for each node.
##
def clique_set_filtered(nodes, pairs, node_filter):
    graph = dict()
    for ind in nodes:
        if node_filter[ind]:
            graph[ind] = list()
    ordered = np.sort(nodes)
    for opair in itertools.combinations(ordered, 2):
        if opair in pairs:
            a, b = opair
            if node_filter[a] and node_filter[b]:
                graph[a].append(b)
                graph[b].append(a)
    cliques = find_cliques(graph)
    return cliques


##
## for each clique, returns the total observations over all elements in the clique.
## it also returns the lowest ACU score and a node that has that score.
## the final results are sorted by count (largest first) [jitter to break ties].
def score_cliques(cliques, counter, ACU_score):
    ans = []
    counts = []
    for clique in cliques:
        min_acu = 0
        total_count = 0
        min_node = None
        for node in clique:
            total_count += counter[node]
            if ACU_score[node] < min_acu:
                min_acu = ACU_score[node]
                min_node = node
        ans.append( [min_node, total_count, min_acu, clique] )
        counts.append(total_count)
    ## add a little jitter to counts to randomly break ties...
    counts = np.array(counts) + 0.001*np.random.random(len(counts))
    order = np.argsort(counts)[::-1]
    return [ans[x] for x in order]

# This code taken in toto from github abhin4v/maximal_cliques.py
#
# Finds all maximal cliques in a graph using the Bron-Kerbosch algorithm. The input graph here is 
# in the adjacency list format, a dict with vertexes as keys and lists of their neighbors as values.
# https://en.wikipedia.org/wiki/Bron-Kerbosch_algorithm
def find_cliques(graph):
    p = set(graph.keys())
    r = set()
    x = set()
    cliques = []
    for v in degeneracy_ordering(graph):
        neighs = graph[v]
        find_cliques_pivot(graph, r.union([v]), p.intersection(neighs), x.intersection(neighs), cliques)
        p.remove(v)
        x.add(v)
    return cliques

def find_cliques_pivot(graph, r, p, x, cliques):
    if len(p) == 0 and len(x) == 0:
        cliques.append(r)
    else:
        u = iter(p.union(x)).next()
        for v in p.difference(graph[u]):
            neighs = graph[v]
            find_cliques_pivot(graph, r.union([v]), p.intersection(neighs), x.intersection(neighs), cliques)
            p.remove(v)
            x.add(v)

def degeneracy_ordering(graph):
    ordering = []
    ordering_set = set()
    degrees = defaultdict(lambda : 0)
    degen = defaultdict(list)
    max_deg = -1
    for v in graph:
        deg = len(graph[v])
        degen[deg].append(v)
        degrees[v] = deg
        if deg > max_deg:
            max_deg = deg
    while True:
        i = 0
        while i <= max_deg:
            if len(degen[i]) != 0:
                break
            i += 1
        else:
            break
        v = degen[i].pop()
        ordering.append(v)
        ordering_set.add(v)
        for w in graph[v]:
            if w not in ordering_set:
                deg = degrees[w]
                degen[deg].remove(w)
                if deg > 0:
                    degrees[w] -= 1
                    degen[deg - 1].append(w)
    ordering.reverse()
    return ordering

## returns a,b as a sorted pair
def pair(a, b):
    if a < b: 
        return a, b 
    else: 
        return b, a

class forwardBranching_pairinfo():
    def __init__(self, cindex, tail, db_graph, contig_lens, K, MAX_DISTANCE, MAX_DEPTH):
        self.db_graph      = db_graph
        self.contig_lens   = contig_lens  
        self.choice_list   = [[(cindex, tail)]]
        self.choice_index  = [0]
        self.K             = K
        self.distance_list = [-contig_lens[cindex] + self.K - 1]
        self.choice_counts = [1]
        self.MAX_DISTANCE  = MAX_DISTANCE
        self.MAX_DEPTH     = MAX_DEPTH
        self.forward_node  = set()
        self.forward_pair  = set()
        self.extend_forward()    
    ## add the next set of neighbors and distances, if possible.
    def extend_one_forward(self):
        cindex, tail = self.choice_list[-1][self.choice_index[-1]]
        ## get path to here:
        path_nodes = self.current_path_nodeonly()[:-1]
        ## marking for node visited
        self.forward_node.add(cindex)
        ## marking all pairs to this node:
        self.forward_pair.update( [pair(x, cindex) for x in path_nodes] )
        distance = self.distance_list[-1]    
        tail_flip = (tail + 1) % 2
        neighbors = self.db_graph[cindex][tail_flip]
        if len(neighbors) == 0:
            return False
        else:
            self.choice_list.append(neighbors)
            self.choice_counts.append(len(neighbors))
            self.choice_index.append(0)
            new_distance = distance - self.K + 1 + self.contig_lens[cindex]
            self.distance_list.append(new_distance)
            return True
    ## extend all steps forward up to MAX_DISTANCE, MAX_DEPTH, or no more extending.
    def extend_forward(self):
        while True:
            if self.distance_list[-1] > self.MAX_DISTANCE:
                break
            if len(self.choice_list)  > self.MAX_DEPTH:
                break
            if not self.extend_one_forward():
                break
    ## increment the choice index
    def next_index(self):
        while(len(self.choice_index) > 0):
            ## increment the last bit
            self.choice_index[-1] += 1
            if self.choice_index[-1] == self.choice_counts[-1]:
                ## ran out of index... pop off the last place...
                self.choice_list.pop()
                self.choice_index.pop()
                self.distance_list.pop()
                self.choice_counts.pop()
            else:
                break
    ## get path
    def current_path_nodeonly(self):
        return [choices[index][0] for (index, choices) in zip(self.choice_index, self.choice_list)]
    ## default iter function for acting as an iterator
    def __iter__(self):
        return self
    ## return the current path, increment index, and re-extend
    def next(self):
        if len(self.choice_index) == 0:
            raise StopIteration
        else:
            self.next_index()
            if len(self.choice_index) > 0:
                self.extend_forward()
            return True

class forwardBranching_to():
    def __init__(self, cindex, tail, target_ind, db_graph, contig_lens, K, MAX_DISTANCE, MAX_DEPTH):
        self.db_graph      = db_graph
        self.contig_lens   = contig_lens  
        self.choice_list   = [[(cindex, tail)]]
        self.target_ind    = target_ind
        self.choice_index  = [0]
        self.K             = K 
        self.distance_list = [-contig_lens[cindex] + K - 1]
        self.choice_counts = [1]
        self.MAX_DISTANCE  = MAX_DISTANCE
        self.MAX_DEPTH     = MAX_DEPTH
        self.extend_forward()    
    ## add the next set of neighbors and distances, if possible.
    def extend_one_forward(self):
        cindex, tail = self.choice_list[-1][self.choice_index[-1]]
        distance = self.distance_list[-1]    
        tail_flip = (tail + 1) % 2
        neighbors = self.db_graph[cindex][tail_flip]
        ## ALSO STOP IF YOU HIT THE TARGET INDEX,WE'RE DONE HERE!
        if len(neighbors) == 0 or cindex == self.target_ind:
            return False
        else:
            self.choice_list.append(neighbors)
            self.choice_counts.append(len(neighbors))
            self.choice_index.append(0)
            new_distance = distance - self.K + 1 + self.contig_lens[cindex]
            self.distance_list.append(new_distance)
            return True
    ## extend all steps forward up to MAX_DISTANCE, MAX_DEPTH, or no more extending.
    def extend_forward(self):
        while True:
            if self.distance_list[-1] > self.MAX_DISTANCE:
                break
            if len(self.choice_list)  > self.MAX_DEPTH:
                break
            if not self.extend_one_forward():
                break
    ## increment the choice index
    def next_index(self):
        while(len(self.choice_index) > 0):
            ## increment the last bit
            self.choice_index[-1] += 1
            if self.choice_index[-1] == self.choice_counts[-1]:
                ## ran out of index... pop off the last place...
                self.choice_list.pop()
                self.choice_index.pop()
                self.distance_list.pop()
                self.choice_counts.pop()
            else:
                break
    ## get path
    def current_path(self):
        return [choices[index] for (index, choices) in zip(self.choice_index, self.choice_list)]
    ## default iter function for acting as an iterator
    def __iter__(self):
        return self
    ## return the current path, increment index, and re-extend
    def next(self):
        if len(self.choice_index) == 0:
            raise StopIteration
        else:
            next_path = self.current_path()
            distance  = self.distance_list[-1] 
            self.next_index()
            if len(self.choice_index) > 0:
                self.extend_forward()
            return distance, next_path


def offset_seq(offset, seq, max_len):
    ans = [" "]*max_len
    ans[offset:(offset+len(seq))] = list(seq)
    return "".join(ans)

'''
STEP 1: identify bridge reads

1. load all the read maps 
2. find all reads with "good ratio"
3. find "bridge reads": those with good ratio maps to more than one contig.
4. store bridge information at the read and the contig.
'''
def load_bridge_read_data(bs_1_dir, contig_maps_plus_npy, MIN_MAP_RATIO, VERBOSE=False):
    ## A (shortened) guide to annotated maps
    CONTIG_NAME = 2
    MATCH       = 9
    OVERLAP     = 10
    ## load contig maps
    if VERBOSE:
        print "loading contig maps"
    contig_maps = load_maps(bs_1_dir, contig_maps_plus_npy)
    if VERBOSE:
        print len(contig_maps), "total maps"
    ## determine map ratio
    map_ratio  = contig_maps[:, MATCH] / contig_maps[:, OVERLAP].astype(float)
    ## filter for maps that have a "good ratio," better than MIN_MAP_RATIO
    good_ratio_filter = map_ratio >= MIN_MAP_RATIO
    ## those are the good ratio maps
    gr_maps = contig_maps[good_ratio_filter]
    if VERBOSE:
        print len(gr_maps), "good ratio maps"
    ## a bridge read is a read with good ratio maps to more than one contig.
    ## for each bridge read, we store all contigs with good ratio maps
    bridge_read_contigs    = dict()
    ## for each contig, we keep a list of bridge reads.
    contig_to_bridge_reads = defaultdict(list)
    ## some counters to keep the pace.
    counter  = 0
    mapcount = 0
    t0 = time.time()
    ## for each read pair in good maps
    for rind, rmaps in read_map_iterator(gr_maps):
        counter += 1
        ## turn the read maps into a numpy array
        rmaps = np.array(rmaps)
        mapcount += len( rmaps )
        ## list the unique contig targets for this read pair
        contig_list = np.unique(rmaps[:, CONTIG_NAME])
        ## if more than one, this is a "bridge-read"
        if len(contig_list) > 1:
            ## we store the list of contigs connected in bridge_read_contigs
            bridge_read_contigs[rind] = contig_list
            ## for each target contig, remember the name of this read.
            for contig in contig_list:
                contig_to_bridge_reads[contig].append(rind)
        if VERBOSE and counter % 100000 == 0:
            t1 = time.time()
            delta = t1 - t0
            part_done = mapcount / float(len(gr_maps))
            info = ["%0.2f done"%part_done, "%0.2f elapsed"%delta, "%0.2f est" % (delta / part_done) , "%0.2f left" % ((delta / part_done) - delta)]
            print tabprint(info)
    return bridge_read_contigs, contig_to_bridge_reads         

'''
STEP 2: load and score the contigs
1. load the contigs from file
2. get lengths and reverse complements as needed
3. build deBruijn graph (contig graph)
4. load unmasking data to score contigs by how "almost certainly unique" they are. (ACU score).
'''
def load_contigs_and_scores(bs_1_dir, contigs_fa, unmask_info_txt, K):
    ## load the contigs
    contigs     = load_contigs(bs_1_dir, contigs_fa)
    ## reverse complements
    contigs_rc  = [reverse_complement(contig) for contig in contigs]
    ## and lengths
    contig_lens = np.array([len(x) for x in contigs])
    ## for the given K, build the contig edge graph.
    dbgraph = contig_graph(contigs, K)
    ## load mask data
    unmask_info_filename = os.path.join(bs_1_dir, unmask_info_txt)
    unmask_info = file2RA(unmask_info_filename)
    bits  = unmask_info['bits']
    flips = unmask_info['flips']
    ## compute flip rate per contig
    rate  = flips.astype(float) / np.maximum(1, bits)
    logmr = np.log10(np.maximum(1e-6, rate))
    log1mr =  np.log10(np.maximum(1e-6, 1 - rate))
    not_flips = bits - flips
    ## determine the almost certainly unique (ACU) score
    ## i.e. how unlikely to observe this sequence of flips given the observed flip rate.
    log10_prob = flips*logmr + not_flips*log1mr
    ACU_score  = log10_prob
    return contigs, contigs_rc, contig_lens, dbgraph, ACU_score

'''
STEP 3: gather forward/reverse path/read information.
For each contig, we get
1. list of contigs that appear on forward (reverse) paths from here
2. which contigs co-occur in forward (reverse) paths (i.e. are "compatible")
3. assign best bridge-to contig for each read based on ACU score
4. count the best assigned reads for the forward (reverse) contigs.
'''
def forward_reverse_contig_info(dbgraph, contig_lens, bridge_read_contigs, contig_to_bridge_reads, ACU_score, 
                                K, MAX_DISTANCE, MAX_DEPTH, VERBOSE=False):
    N = len(dbgraph)
    t0 = time.time()
    forward_counts  = dict()    # which contigs are forward  and their counts in the bridge reads
    backward_counts = dict()    # which contigs are backward and ...
    forward_pairs  = dict()     # which contig pairs co-occur on forward  paths
    backward_pairs = dict()     # which contig pairs co-occur on backward paths
    ## for each index
    if VERBOSE:
        print "finding and counting branches"
    for cindex in range(N):
        ## find forward and backward branching paths
        forward  = forwardBranching_pairinfo(cindex, 0, dbgraph, contig_lens, K, MAX_DISTANCE, MAX_DEPTH)
        backward = forwardBranching_pairinfo(cindex, 1, dbgraph, contig_lens, K, MAX_DISTANCE, MAX_DEPTH)
        ## these are a bit hacked from iterators, so to fill the fields we have to iterate to completion.
        for _ in forward:
            continue
        for _ in backward:
            continue
        ## pull the forward and backward nodes (removing this index)
        fnode = forward .forward_node
        bnode = backward.forward_node
        fnode.remove(cindex)
        bnode.remove(cindex)
        ## keep track of all compatible pairs
        forward_pairs[cindex]  = forward .forward_pair
        backward_pairs[cindex] = backward.forward_pair    
        ## pull all of the bridge reads for this contig
        bridge_reads = contig_to_bridge_reads[cindex]
        ## keep counts for the forward and reverse
        forward_counter  = Counter()
        backward_counter = Counter()
        ## for each read
        for bindex in bridge_reads:
            ## get the best scoring nodes, if they exist
            best_forward   = None
            best_backward  = None
            forward_score  = 0
            backward_score = 0
            for ind in bridge_read_contigs[bindex]:
                if ind in fnode:
                    if ACU_score[ind] < forward_score:
                        best_forward  = ind
                        forward_score = ACU_score[ind]
                if ind in bnode:
                    if ACU_score[ind] < backward_score:
                        best_backward  = ind
                        backward_score = ACU_score[ind]
            ## if they exist, increment the forward and backward counters
            if best_forward != None:
                forward_counter[best_forward] += 1
            if best_backward != None:
                backward_counter[best_backward] += 1
        ## when done with all reads, store the answers to forward/backward counts
        forward_counts [cindex] = forward_counter
        backward_counts[cindex] = backward_counter
        if VERBOSE and (cindex% 1000) == 0:
            t1 = time.time()
            delta = t1- t0
            print cindex, "%0.2f seconds" % delta
    return forward_counts, backward_counts, forward_pairs, backward_pairs


'''
STEP 4: identify and score compatible node sets.
For each contig, we get
1. find the forward (reverse) maximal cliques of compatible sets.
2. score each clique: i.e. total read count, best ACU, and node with best ACU.
'''
def find_and_score_compatible_sets(forward_counts, forward_pairs, backward_counts, backward_pairs, ACU_score, acu_filter, VERBOSE=False):
    N = len(forward_counts)
    all_scores = []
    for ind in range(N):
        if VERBOSE and (ind % 1000 == 0):
            print ind
        if acu_filter[ind]:
            ## forward
            counter = forward_counts[ind]
            pairs  = forward_pairs[ind]
            nodes  = counter.keys()
            cliques = clique_set_filtered(nodes, pairs, acu_filter)
            forward_scores = score_cliques(cliques, counter, ACU_score)
            ## backward
            counter = backward_counts[ind]
            pairs  = backward_pairs[ind]
            nodes  = counter.keys()
            cliques = clique_set_filtered(nodes, pairs, acu_filter)
            backward_scores = score_cliques(cliques, counter, ACU_score)
            scores = (forward_scores, backward_scores)
        else:
            scores = [[], []]
        all_scores.append(scores)
    return all_scores

'''
STEP 5: find the best forward and reverse (by count)
'''
def find_favorites(all_scores):
    favorites = []
    for scores in all_scores:
        fscores, bscores = scores
        try:
            best_forward = fscores[0][0]
        except:
            best_forward = None
        try:
            best_backward = bscores[0][0]
        except:
            best_backward= None        
        favorites.append( [best_backward, best_forward] )
    return favorites

'''
STEP 6: find connected components with strong matches.
'''
def strong_match_connected_components(favorites, acu_filter):
    N = len(acu_filter)
    components = []
    node_seen = np.zeros(N, dtype=bool)
    for start_node in range(N):
        ## if you're not in the filter, forget about it...
        if not acu_filter[start_node]:
            continue
        ## if we haven't seen this node yet
        if not node_seen[start_node]:
            ## store paths for forward and reverse
            paths = [[],[]]
            for start_tail in range(0, 2):
                path = paths[start_tail]
                node = start_node
                tail = start_tail
                ## build path as long as there is a strong match
                while True:
                    path.append((node, tail))
                    flip = (tail + 1) % 2
                    node, tail = get_strong_match(node, flip, favorites)
                    if node == None:
                        break
            ## reverse the second path
            rev = [(ind, (tail + 1)%2) for (ind, tail) in paths[1]][::-1]
            ## append it to the first path (without the first element which is already in rev)
            total_path = rev + paths[0][1:]
            ## update what's been seen
            for ind, tail in total_path:
                node_seen[ind] = True
            ## append this to the components
            components.append(total_path)
    return components

'''
STEP 7: convert connected components (paths) into full paths.
Component paths move from good ACU score nodes.
Need to fill in.
To do so, we find all compatible paths that span and count how often they are realized in the reads.
The one most often realized (with jitter to randomize equality) is selected.
'''
def complete_component_paths(components, 
                             bridge_read_contigs, contig_to_bridge_reads, 
                             dbgraph, contig_lens, K, MAX_DISTANCE, MAX_DEPTH):
    full_paths = []
    for path in components:
        L = len(path)
        if L == 1:
            full_path = path
        else:
            bridges = []
            ## iterating over the pairs of junctions
            for ind in range(L - 1):
                ## determine the current node and next
                current_node, next_node = path[ind], path[ind+1]
                cindex = current_node[0]
                tail   = current_node[1]
                target_ind = next_node[0]
                ## find all forward branching paths that terminate in target ind
                fpaths = forwardBranching_to(cindex, tail, target_ind, dbgraph, contig_lens, K, MAX_DISTANCE, MAX_DEPTH)
                good_paths = []
                ## the forward branching to just stops searching when the target is found, 
                ## it doesn't exclude other paths and so we do that here.
                ## SHOULD CLEAN UP FORWARD_BRANCHING classes to handle this.
                for dist, fpath in fpaths:
                    if fpath[-1][0] == target_ind:
                        good_paths.append(fpath)
                ## from the good paths, we want to check bridge reads for "perfect coverage"
                P = len(good_paths)
                perfect_score = np.array( [len(x) for x in good_paths] )
                ## we make a lookup table saying which good_path each index appears in
                gp_lookup = defaultdict(list)
                for ind, gpath in enumerate(good_paths):
                    for index, tail in gpath:
                        gp_lookup[index].append(ind)
                ## now we find the bridge reads:
                ## reads from the cindex and to the cindex
                br_from = contig_to_bridge_reads[cindex]
                br_to   = contig_to_bridge_reads[target_ind]
                ## and then the intersect
                br_both = np.intersect1d(br_from, br_to, assume_unique=True, return_indices=False)
                ## this is going to keep track of how often is each path perfectly covered.
                path_read_count = np.zeros(P, dtype=int)
                ## for reads that bridge A and B,
                for rind in br_both:
                    ## how many of the good_path contigs did you see?
                    gp_score = np.zeros(P, dtype=int)
                    ## these are the contigs this read bridges
                    bcontigs = bridge_read_contigs[rind]
                    ## increment the counter according to the lookup
                    for cind in bcontigs:
                        np.add.at(gp_score, gp_lookup[cind], 1)
                    ## wherever you hit every single contig, increment that counter.
                    path_read_count[gp_score == perfect_score] += 1
                ## return the best path
                best_path = good_paths[np.argmax(path_read_count + 0.001*np.random.random(P))]
                bridges.append(best_path)
            ## combine bridges to get full path
            full_path = [] + bridges[0]
            for bridge in bridges[1:]:
                full_path += bridge[1:]
    #     print full_path
        full_paths.append(full_path)
    return full_paths

'''
STEP 8: convert paths back into strings
pretty straightfoward.
flip if it says flip.
clip off the K-1 overlap when joining.
'''
def paths_to_seqs(full_paths, contigs, contigs_rc, K):
    final_seqs = []
    for path in full_paths:
        ans = []
        ## store sequence or reverse according to ans.
        for node, reverse in  path:
            node_str = contigs[node] if reverse == 0 else contigs_rc[node]
            ans.append(node_str)
        ## start with the zero entry
        path_string = ans[0]
        for contig_string in ans[1:]:
            ## and append everything but the first K-1 bases for each subseuqence string
            path_string += contig_string[(K-1):]
        final_seqs.append(path_string)
        
        final_lens = [len(seq) for seq in final_seqs]
        ans = [final_seqs[x] for x in np.argsort(final_lens)[::-1]]
    return ans


'''
This function uses the dbGraph, the unmasking data, and reads that map to multiple contigs
to assemble paths through the dbgraph.

We use the unmasking data to identify contigs that are almost certainly unique.
This is determined by their ACU_score, computed from # of bits and observed mutation rate.
The score is log10 and lower is "more likely unique"

We examine each read-pair to see if it connects contigs.
Connects between ACU contigs are counted and confident connections are identified.
ACU contigs that are unambiguously joined connect into chains.

These "ACU paths" are then filled in with intermediate nodes, preferring the most observed path.
These filled in paths are converted into sequences and saved to file.
'''
def acu_contig_extension(bs_1_dir, contig_maps_plus_npy, contigs_fa, unmask_info_txt, acu_extension_fa,
                         MIN_MAP_RATIO, MAX_DISTANCE, MAX_DEPTH, KMER_SIZE, MAX_ACU, VERBOSE=False):
    K = KMER_SIZE
    bridge_read_contigs, contig_to_bridge_reads = load_bridge_read_data(bs_1_dir, 
                                                                        contig_maps_plus_npy, 
                                                                        MIN_MAP_RATIO, 
                                                                        VERBOSE)
    contigs, contigs_rc, contig_lens, dbgraph, ACU_score = load_contigs_and_scores(bs_1_dir,
                                                                                   contigs_fa, 
                                                                                   unmask_info_txt,
                                                                                   KMER_SIZE)
    forward_counts, backward_counts, forward_pairs, backward_pairs = forward_reverse_contig_info(dbgraph, 
                                                                                                 contig_lens,
                                                                                                 bridge_read_contigs, 
                                                                                                 contig_to_bridge_reads, 
                                                                                                 ACU_score,
                                                                                                 K,
                                                                                                 MAX_DISTANCE,
                                                                                                 MAX_DEPTH,
                                                                                                 VERBOSE)
    acu_filter = ACU_score <= MAX_ACU    
    all_scores = find_and_score_compatible_sets(forward_counts, 
                                                forward_pairs, 
                                                backward_counts, 
                                                backward_pairs,
                                                ACU_score, 
                                                acu_filter)
    favorites  = find_favorites(all_scores)
    components = strong_match_connected_components(favorites, acu_filter)
    full_paths = complete_component_paths(components, 
                                          bridge_read_contigs, contig_to_bridge_reads, 
                                          dbgraph, contig_lens, K, MAX_DISTANCE, MAX_DEPTH)
    acu_seqs   = paths_to_seqs(full_paths, contigs, contigs_rc, K)
    write_contigs(acu_seqs, bs_1_dir, acu_extension_fa)
    
    

