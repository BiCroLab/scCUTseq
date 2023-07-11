#!/usr/bin/env python

import argparse
import gzip
import re
import os
from multiprocessing import Pool
from functools import partial
from collections import Counter
from itertools import groupby
from operator import itemgetter
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from distance import hamming

# Function to append a suffix to the base name 
# Return a string
def fastq_suffix(filename, suffix, directory=None, extension='.fq.gz'):
    dirname = directory if directory else os.path.dirname(filename)
    dirname = '.' if dirname == '' else dirname
    basename = os.path.basename(filename)
    basename = basename.replace('.gz', '')
    basename = basename.replace('.fastq', '.fq')
    basename = basename.replace('.fq', ''.join(['_', suffix, extension]))
    name = os.path.join(dirname, basename)
    return name

# Function to read the FastQ file in batches
# Return list of tuples
def batch_iterator(iterator, batch_size):
    '''Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    '''
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = iterator.__next__()
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch

# Function to read the batch list and pass the entries to the map function 
# Return a tuple
def read_batch(in_handle):
    for title, seq, qual in in_handle:
        yield title, seq, qual

# Function to extract the barcode from the read
# Return a tuple
def process_barcode(read, barPos, umiPos, keepPos, patLen):
    id, seq, qual = read
    id = id.split()
    id_extra = []

    bc_seq = ''.join([seq[beg:end] for beg, end in barPos])
    bc_tag = ':'.join(['BC', 'Z'] + [bc_seq])
    qt_tag = ':'.join(['QT', 'Z'] + [''.join([qual[beg:end] for beg, end in barPos])])
    id_extra = id_extra + [bc_tag, qt_tag]
    if umiPos:
        rx_tag = ':'.join(['RX', 'Z'] + [''.join([seq[beg:end] for beg, end in umiPos])])
        qx_tag = ':'.join(['QX', 'Z'] + [''.join([qual[beg:end] for beg, end in umiPos])])
        id_extra = id_extra + [rx_tag, qx_tag]
    if keepPos:
        seq = ''.join([seq[beg:end] for beg, end in keepPos] + [seq[patLen:]])
        qual = ''.join([qual[beg:end] for beg, end in keepPos] + [qual[patLen:]])
    else:
        seq = seq[patLen:]
        qual = qual[patLen:]

    record = [ ' '.join(['_'.join(['@'+id[0]] + id_extra)] + id[1:]), seq, '+', qual ]
    return bc_seq, record

def process_barcode_paired(read1, read2, barPos, umiPos, keepPos, patLen):
    # Process read 1
    id, seq, qual = read1
    id = id.split()
    id_extra = []

    bc_seq = ''.join([seq[beg:end] for beg, end in barPos])
    bc_tag = ':'.join(['BC', 'Z'] + [bc_seq])
    qt_tag = ':'.join(['QT', 'Z'] + [''.join([qual[beg:end] for beg, end in barPos])])
    id_extra = id_extra + [bc_tag, qt_tag]
    if umiPos:
        rx_tag = ':'.join(['RX', 'Z'] + [''.join([seq[beg:end] for beg, end in umiPos])])
        qx_tag = ':'.join(['QX', 'Z'] + [''.join([qual[beg:end] for beg, end in umiPos])])
        id_extra = id_extra + [rx_tag, qx_tag]
    if keepPos:
        seq = ''.join([seq[beg:end] for beg, end in keepPos] + [seq[patLen:]])
        qual = ''.join([qual[beg:end] for beg, end in keepPos] + [qual[patLen:]])
    else:
        seq = seq[patLen:]
        qual = qual[patLen:]

    record1 = [ ' '.join(['_'.join(['@'+id[0]] + id_extra)] + id[1:]), seq, '+', qual ]

    # Process read 2
    id, seq, qual = read2
    id = id.split()
    record2 = [ ' '.join(['_'.join(['@'+id[0]] + id_extra)] + id[1:]), seq, '+', qual ]

    return bc_seq, (record1, record2)

if __name__ == '__main__':

    # Parse arguments
    parser = argparse.ArgumentParser(description='Split reads by barcode into multiple files', prefix_chars='-')
    parser.add_argument('--version', action='version', version='Version: %(prog)s 1.0')

    parser.add_argument('-i', '--input', type=str, dest='infile', nargs='+', help='the input fastq file', required=True)
    parser.add_argument('-i2', '--input2', type=str, dest='infile2', nargs='+', help='the input fastq file of R2', required=False)
    parser.add_argument('-o', '--outdir', type=str, default=[None], dest='outdir', nargs=1, help='the output directory', required=True)
    parser.add_argument('-b', '--barcode', type=str, default=[None], dest='barcode', nargs=1, help='the file containing the barcode sequences', required=True)
    parser.add_argument('-d', '--distance', type=int, dest='distance', nargs=1, choices=range(3), default=[2], help='Hamming distance (default: 2)')
    parser.add_argument('-p', '--pattern', type=str, default=None, dest='pattern', nargs=1, help='the sequence pattern to process (can contain B, U, K and D)', required=True)
    parser.add_argument('-t', '--threads', type=int, default=[1], dest='threads', nargs=1, help='the number of parallel threads to use during the splitting')
    parser.add_argument('-s', '--size', type=int, default=[5000000], dest='batchsize', nargs=1, help='larger numbers will result in higher memory consumption')
    parser.add_argument('-k', '--keep', help='write discarded reads to a separate file', action='store_true', default=False)
    parser.add_argument('-r', '--run', help='run type', dest='run', nargs=1, choices=['single', 'paired'], required=True)
    parser.add_argument('-v', '--verbose', help='increase output verbosity', action='store_true', default=False)
    args = parser.parse_args()

    run = args.run[0]
    inFile = args.infile
    if run == "paired":
        inFile2 = args.infile2
    distance = args.distance[0]
    verbose = args.verbose

    output_dir = args.outdir[0] if args.outdir[0][-1] == '/' else args.outdir[0] + '/'
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    barcode = args.barcode[0]
    if barcode == None:
        print('Barcode(s) not specified!')
        exit(1)
    pattern = args.pattern[0]
    threads = args.threads[0]
    batch_size = args.batchsize[0]
    keep = args.keep

    if re.fullmatch(r'[BUKD]+', pattern) == None:
        print('Error: pattern contains invalid character(s):', set(re.sub(r'[BUKD]+', '', pattern)))
        exit(1)
    else:
        len_pattern = len(pattern)
        pos_keep = [ i.span() for i in re.finditer('[K]+', pattern) ]
        pos_keep = pos_keep if len(pos_keep) > 0 else False
        pos_discard = [ i.span() for i in re.finditer('[D]+', pattern) ]
        pos_discard = pos_discard if len(pos_discard) > 0 else False
        pos_barcode = [ i.span() for i in re.finditer('[B]+', pattern) ]
        if sum([b-a for a, b in pos_barcode]) == 0:
            print('No barcode (B) found in pattern', pattern)
            exit(1)
        pos_umi = [ i.span() for i in re.finditer('[U]+', pattern) ]
        pos_umi = pos_umi if len(pos_umi) > 0 else False

    # Output command information if verbosity is set
    if verbose:
        if run == "single":
            print('Input file(s):', inFile)
        else:
            print('Input file(s) Read 1:', inFile)
            print('Input file(s) Read 2:', inFile2)
        print('Output directory:', output_dir)
        print('Refenence barcodes file:', barcode)
        print('Defined UMIs coordinates offset:', str(pos_umi))
        print('Number of parallel processes:', str(threads))
        print('Input/output batch size:', str(batch_size))
        print('Write reads with ambiguous or incorrect pattern to separate file:', keep)
        print('Pattern:', pattern, ''.join(['(', str(len(pattern)), ')']))
        if pos_barcode:
            print('Barcode position(s):', [ (a+1, b) for a,b in pos_barcode ])
            print('Barcode distance (Hamming):', distance)
        if pos_umi:
            print('UMI position(s):', [ (a+1, b) for a,b in pos_umi ])
        if pos_keep:
            print('Pattern position(s) to keep:', [ (a+1, b) for a,b in pos_keep ])
        if pos_discard:
            print('Pattern position(s) to discard:', [ (a+1, b) for a,b in pos_discard ])

    # Start loop through input files
    # single-end processing
    if run == "single":
        # Read in barcodes and perform barcode checks
        with open(barcode, 'rt') as bc_handle:
            barcode_files = { bc.strip() : output_dir + bc.strip() + '.fq.gz' for bc in bc_handle }
            
            # Barcode distance check
            for key1 in barcode_files.keys():
                for key2 in barcode_files.keys():
                    if hamming(key1, key2) <= distance and key1 != key2:
                        print(f'Warning: two barcodes ({key1} and {key2}) have Hamming distance <= {distance}')

            for key in barcode_files.keys():
                if os.path.isfile(barcode_files[key]):
                    os.remove(barcode_files[key])
            if os.path.isfile(output_dir + 'discarded.fq.gz') and keep:
                os.remove(output_dir + 'discarded.fq.gz')

            len_barcode = len(list(barcode_files.keys())[0])
            if len_barcode != sum([b-a for a, b in pos_barcode]):
                print("Error: Pattern (" + str(sum([b-a for a, b in pos_barcode])) + ") and file (" + str(len_barcode) + ") barcode lengths differ!")
                exit(1)

            global_counter = Counter([ key for key in barcode_files.keys()] + ['discarded'])
            for input_file in inFile:
                if os.path.splitext(input_file)[1] == '.gz':
                    file_iterator = FastqGeneralIterator(gzip.open(input_file, 'rt'))
                else:
                    file_iterator = FastqGeneralIterator(open(input_file, 'rt'))
                
                # Process reads
                processed = 0
                local_counter = Counter([ key for key in barcode_files.keys()] + ['discarded'])
                with Pool(processes = threads) as pool:
                    for i, batch in enumerate(batch_iterator(file_iterator, batch_size)):
                        g = pool.map(partial(process_barcode, barPos=pos_barcode, umiPos=pos_umi, keepPos=pos_keep, patLen=len_pattern), read_batch(batch))
                        processed += len(g)
                        # Group reads by barcode
                        sortkeyfn = itemgetter(0)
                        g.sort(key=sortkeyfn)
                        g = {key : list(v[1] for v in valuesiter) for key, valuesiter in groupby(g, key=sortkeyfn)}
                        # Write (append) to files
                        for bc in g.keys():
                            if bc in barcode_files.keys():
                                with gzip.open(barcode_files[bc], 'at') as out_handle:
                                    [ out_handle.write(field+'\n') for read in g[bc] for field in read ]
                                global_counter[bc] += len(g[bc])
                                local_counter[bc] += len(g[bc])
                            else:
                                found = False
                                for key in barcode_files.keys():
                                    if hamming(bc, key) <= distance:
                                        with gzip.open(barcode_files[key], 'at') as out_handle:
                                            [ out_handle.write(field+'\n') for read in g[bc] for field in read ]
                                        global_counter[key] += len(g[bc])
                                        local_counter[key] += len(g[bc])
                                        found = True
                                        break
                                if not found:
                                    with gzip.open(output_dir + 'discarded.fq.gz', 'at') as out_handle:
                                        [ out_handle.write(field+'\n') for read in g[bc] for field in read ]
                                    global_counter['discarded'] += len(g[bc])
                                    local_counter['discarded'] += len(g[bc])
                        if verbose:
                            print(os.path.basename(input_file) + ':', 'processed', processed, 'reads')
                    file_iterator.close()
                out_stats = re.sub('(fq|fastq).*$', 'stats.txt', os.path.basename(input_file))
                with open(output_dir + out_stats, 'wt') as out_handle:
                    for key in barcode_files.keys():
                        out_handle.write('\t'.join(['Barcode', key, str(local_counter[key]-1), str((local_counter[key]-1)/processed) ]) + '\n')
                    out_handle.write('\t'.join(['Barcode', 'Discarded', str(local_counter['discarded']-1), str((local_counter['discarded']-1)/processed)]) + '\n')
            if len(inFile) > 1:
                with open(output_dir + 'Cumulative.stats.txt', 'wt') as out_handle:
                    for key in barcode_files.keys():
                        out_handle.write('\t'.join(['Barcode', key, str(local_counter[key]-1), str((local_counter[key]-1)/processed)]) + '\n')
                    out_handle.write('\t'.join(['Barcode', 'Discarded', str(local_counter['discarded']-1), str((local_counter['discarded']-1)/processed)]) + '\n')
        
    # Start loop through input files
    # paired-end processing
    if run == "paired":
        # Read in barcodes and perform barcode checks
        with open(barcode, 'rt') as bc_handle:
            barcode_files_r1 = { bc.strip() : output_dir + bc.strip() + '_R1.fq.gz' for bc in bc_handle }
        barcode_files_r2 = { bc : output_dir + bc + '_R2.fq.gz' for bc in barcode_files_r1.keys() }
            
        # Check if barcode file keys are identical
        if not barcode_files_r1.keys() == barcode_files_r2.keys():
            print("Error: Barcode R1 and R2 files are not identical")
            exit(1)
        
        # Barcode distance check
        for key1 in barcode_files_r1.keys():
            for key2 in barcode_files_r1.keys():
                if hamming(key1, key2) <= distance and key1 != key2:
                    print(f'Warning: two barcodes ({key1} and {key2}) have Hamming distance <= {distance}')

        # Remove files if present
        for key in barcode_files_r1.keys():
            if os.path.isfile(barcode_files_r1[key]):
                os.remove(barcode_files_r1[key])
        for key in barcode_files_r2.keys():
            if os.path.isfile(barcode_files_r2[key]):
                os.remove(barcode_files_r2[key])

        # Remove discarded file if present
        if os.path.isfile(output_dir + 'discarded_R1.fq.gz') and keep:
            os.remove(output_dir + 'discarded_R1.fq.gz')
        if os.path.isfile(output_dir + 'discarded_R2.fq.gz') and keep:
            os.remove(output_dir + 'discarded_R2.fq.gz')

        len_barcode = len(list(barcode_files_r1.keys())[0])
        if len_barcode != sum([b-a for a, b in pos_barcode]):
            print("Error: Pattern (" + str(sum([b-a for a, b in pos_barcode])) + ") and file (" + str(len_barcode) + ") barcode lengths differ!")
            exit(1)
    
        global_counter = Counter([ key for key in barcode_files_r1.keys()] + ['discarded'])
        for input_file, input_file2 in zip(inFile, inFile2):
            if os.path.splitext(input_file)[1] == '.gz':
                file_iterator_r1 = FastqGeneralIterator(gzip.open(input_file, 'rt'))
                file_iterator_r2 = FastqGeneralIterator(gzip.open(input_file2, 'rt'))
            else:
                file_iterator_r1 = FastqGeneralIterator(open(input_file, 'rt'))
                file_iterator_r2 = FastqGeneralIterator(open(input_file2, 'rt'))
            
            # Process reads
            processed = 0
            local_counter = Counter([ key for key in barcode_files_r1.keys()] + ['discarded'])
            with Pool(processes = threads) as pool:
                for i, (batch1, batch2) in enumerate(zip(batch_iterator(file_iterator_r1, batch_size), 
                batch_iterator(file_iterator_r2, batch_size))):
                    g = pool.starmap(partial(process_barcode_paired, barPos=pos_barcode, umiPos=pos_umi, keepPos=pos_keep, patLen=len_pattern), 
                    zip(read_batch(batch1), read_batch(batch2)))
                    processed += len(g)
                    # Group reads by barcode
                    sortkeyfn = itemgetter(0)
                    g.sort(key=sortkeyfn)
                    g = {key : list(v[1] for v in valuesiter) for key, valuesiter in groupby(g, key=sortkeyfn)}
                    # Write (append) to files
                    for bc in g.keys():
                        if bc in barcode_files_r1.keys(): # Doesn't matter if _r1 or _r2 is used here
                            with gzip.open(barcode_files_r1[bc], 'at') as out_handle:
                                [ out_handle.write(field+'\n') for read in g[bc] for field in read[0] ]
                            with gzip.open(barcode_files_r2[bc], 'at') as out_handle:
                                [ out_handle.write(field+'\n') for read in g[bc] for field in read[1] ]
                            global_counter[bc] += len(g[bc])
                            local_counter[bc] += len(g[bc])
                        else:
                            found = False
                            for key in barcode_files_r1.keys():
                                if hamming(bc, key) <= distance:
                                    with gzip.open(barcode_files_r1[key], 'at') as out_handle:
                                        [ out_handle.write(field+'\n') for read in g[bc] for field in read[0] ]
                                    with gzip.open(barcode_files_r2[key], 'at') as out_handle:
                                        [ out_handle.write(field+'\n') for read in g[bc] for field in read[1] ]
                                    global_counter[key] += len(g[bc])
                                    local_counter[key] += len(g[bc])
                                    found = True
                                    break
                            if not found:
                                with gzip.open(output_dir + 'discarded_R1.fq.gz', 'at') as out_handle:
                                    [ out_handle.write(field+'\n') for read in g[bc] for field in read[0] ]
                                with gzip.open(output_dir + 'discarded_R2.fq.gz', 'at') as out_handle:
                                    [ out_handle.write(field+'\n') for read in g[bc] for field in read[1] ]
                                global_counter['discarded'] += len(g[bc])
                                local_counter['discarded'] += len(g[bc])
                    if verbose:
                        print(os.path.basename(input_file) + ':', 'processed', processed, 'read pairs')
                file_iterator_r1.close()
                file_iterator_r2.close()
            out_stats = re.sub('(_R1.fq|_R1.fastq).*$', '.stats.txt', input_file)
            with open(out_stats, 'wt') as out_handle:
                for key in barcode_files_r1.keys():
                    out_handle.write('\t'.join(['Barcode', key, str(local_counter[key]-1), str((local_counter[key]-1)/processed) ]) + '\n')
                out_handle.write('\t'.join(['Barcode', 'Discarded', str(local_counter['discarded']-1), str((local_counter['discarded']-1)/processed)]) + '\n')
        if len(inFile) > 1:
            with open('Comulative.stats.txt', 'wt') as out_handle:
                for key in barcode_files_r1.keys():
                    out_handle.write('\t'.join(['Barcode', key, str(local_counter[key]-1), str((local_counter[key]-1)/processed)]) + '\n')
                out_handle.write('\t'.join(['Barcode', 'Discarded', str(local_counter['discarded']-1), str((local_counter['discarded']-1)/processed)]) + '\n')

            