#!/usr/bin/env python

import simplesam
import argparse
import re
import sys

if __name__ == '__main__':
    
    # Parse arguments
    parser = argparse.ArgumentParser(description='Move tags from header to sam tags in bam file', prefix_chars='-')
    parser.add_argument('-i', '--infile', type=str, dest='infile', nargs=1, help='the input bam file', required=True)
    args = parser.parse_args()    

    infile = args.infile[0]

    with simplesam.Reader(open(infile)) as in_bam:
        with simplesam.Writer(sys.stdout, in_bam.header) as out_sam:
            for read in in_bam:
                # Get header name and split by "_"
                header = read.qname.split("_")

                # Set new header and new tags
                read.qname = header[0]
                read.tags['BC'] = re.sub(".*:", "", header[1])
                read.tags['QT'] = re.sub(".*:", "", header[2])
                read.tags['RX'] = re.sub(".*:", "", header[3])
                read.tags['QX'] = re.sub(".*:", "", header[4])

                # Write new reads
                out_sam.write(read)
