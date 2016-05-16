#! /usr/bin/env python
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2014-2015, Michigan State University.
# Copyright (C) 2015-2016, The Regents of the University of California.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#
#     * Neither the name of the Michigan State University nor the names
#       of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written
#       permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Contact: khmer-project@idyll.org
from __future__ import print_function
import argparse
import screed
import sys
import khmer
from itertools import izip, imap
from multiprocessing import Pool
from multiprocessing.managers import BaseManager


def output_single(read):
    if hasattr(read, 'quality'):
        return "@%s\n%s\n+\n%s\n" % (read.name, read.sequence, read.quality)
    else:
        return ">%s\n%s\n" % (read.name, read.sequence)

def output_txt(read):
    name, seq, qual = read
    return "@%s\n%s\n+\n%s\n" % (name, seq, qual)

def txt_generator(r1,r2):
    return (r1.name, r1.sequence, r1.quality), (r2.name, r2.sequence, r2.quality)

class HT_Manager(BaseManager):
    pass    

def thread_init(ht, mn, mx, sf):
    global htable
    global min_coverage
    global max_coverage
    global single_file
    htable = ht
    min_coverage=mn
    max_coverage=mx
    single_file = sf

def find_cov(x):
    read1, read2 = x
    name1, seq1, qual1 = read1
    name2, seq2, qual2 = read2
    seq = seq1.upper()
    seq = seq.replace('N', 'A')
    pseq = seq2.upper()
    pseq = pseq.replace('N','A')

    try:
        med, _, _ = htable.get_median_count(seq)
        pmed, _, _ = htable.get_median_count(pseq)
    except ValueError:
        return (None, None)

    keep = True
    pkeep = True

    if min_coverage:
        if med < min_coverage:
            keep = False
        if pmed < min_coverage:
            pkeep = False

    if max_coverage:
        if med > max_coverage:
            keep = False
        if pmed > max_coverage:
            pkeep = False

    if keep and pkeep:
        return (read1,read2)

    elif single_file:
        if keep:
            return (read1,None)
        elif pkeep:
            return (None,read2)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--min-coverage', type=int, default=None)
    parser.add_argument('-M', '--max-coverage', type=int, default=None)
    parser.add_argument('-T', '--threads', type=int, default=None)
    parser.add_argument('input_count_graph')
    parser.add_argument('input_readfile')
    parser.add_argument('input_pairfile')
    parser.add_argument('output_readfile')
    parser.add_argument('output_pairfile')
    parser.add_argument('output_singlefile', default=None, nargs='?')
    args = parser.parse_args()

    print('min_coverage: %s' % args.min_coverage, file=sys.stderr)
    print('max_coverage: %s' % args.max_coverage, file=sys.stderr)

    if not (args.min_coverage or args.max_coverage):
        print("neither min nor max coverage specified!? exiting!", file=sys.stderr)
        sys.exit(1)

    if args.min_coverage and args.max_coverage and \
       args.max_coverage < args.min_coverage:
        print("min_coverage > max_coverage!? exiting!", file=sys.stderr)
        sys.exit(1)

    output_file = args.output_readfile
    output_fp = open(output_file, 'w')
    output_pairfile = args.output_pairfile
    output_pfp = open(output_pairfile,'w')

    if args.output_singlefile:
        output_singlefile = args.output_singlefile
        output_sfp = open(output_singlefile,'w')

    n_kept = 0
    n = 0

    if args.threads:
        HT_Manager.register('get_ht',khmer.load_countgraph)
        manager = HT_Manager()
        manager.start()
        htable = manager.get_ht(args.input_count_graph)
        pool = Pool(processes=args.threads, initializer=thread_init,
                    initargs=[htable,args.min_coverage,args.max_coverage,args.output_singlefile])
        it1 = iter(screed.open(args.input_readfile))
        it2 = iter(screed.open(args.input_pairfile))
        passable_generator = imap(txt_generator,it1,it2)
        it3 = pool.imap(find_cov,passable_generator,chunksize=1000)
        for read1, read2 in it3:
            if n % 100000 == 0:
                print('...', n, n_kept, file=sys.stderr)
            n += 2
            if read1 and read2:
                n_kept += 2
                output_fp.write(output_txt(read1))
                output_pfp.write(output_txt(read2))
            elif args.output_singlefile:
                if read1:
                    n_kept += 1
                    output_sfp.write(output_txt(read1))
                elif read2:
                    n_kept += 1
                    output_sfp.write(output_txt(read2))

    else:
        htable = khmer.load_countgraph(args.input_count_graph)
        pair_iter = iter(screed.open(args.input_pairfile))

        for n, record in enumerate(screed.open(args.input_readfile)):
            if n % 100000 == 0:
                print('...', n, n_kept, file=sys.stderr)

            seq = record.sequence.upper()
            seq = seq.replace('N', 'A')
            pair_record = pair_iter.next()
            pseq = pair_record.sequence.upper()
            pseq = pseq.replace('N','A')

            try:
                med, _, _ = htable.get_median_count(seq)
                pmed, _, _ = htable.get_median_count(pseq)
            except ValueError:
                continue

            keep = True
            pkeep = True
            if args.min_coverage:
                if med < args.min_coverage:
                    keep = False
                if pmed < args.min_coverage:
                    pkeep = False

            if args.max_coverage:
                if med > args.max_coverage:
                    keep = False
                if pmed > args.max_coverage:
                    pkeep = False

            if keep and pkeep:
                n_kept += 2
                output_fp.write(output_single(record))
                output_pfp.write(output_single(pair_record))

            elif args.output_singlefile:
                if keep:
                    n_kept += 1
                    output_sfp.write(output_single(record))
                elif pkeep:
                    n_kept += 1
                    output_sfp.write(output_single(pair_record))

    print('consumed %d reads; kept %d' % (2 * (n + 1), n_kept), file=sys.stderr)

if __name__ == '__main__':
    main()
