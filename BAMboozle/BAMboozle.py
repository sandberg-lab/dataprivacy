# BAMboozle.py - de-identify sequencing data
# Author: Christoph Ziegenhain / christoph.ziegenhain@ki.se
# Last update: 11-01-2021

import os
import pysam
import argparse
import multiprocessing as mp
import itertools

def makeBAMheader(args, v):
    bam = pysam.AlignmentFile(args.bam, 'rb')
    hdr = bam.header.to_dict()
    bam.close()

    cmdlinecall = 'BAMboozle --bam '+args.bam+' --out '+args.out+' --fa '+args.fa+' --p '+str(args.p)
    if args.strict:
        cmdlinecall = cmdlinecall+' --strict'
    if args.keepunmapped:
        cmdlinecall = cmdlinecall+' --keepunmapped'
    if args.keepsecondary:
        cmdlinecall = cmdlinecall+' --keepsecondary'

    pg = {'ID': 'BAMboozle', 'PN': 'BAMboozle',
          'CL': cmdlinecall, 'VN': v}

    if 'PG' in hdr:
        pglines = hdr['PG']
        pglines.append(pg)
    else:
        pglines = [pg]
    hdr['PG'] = pglines

    return hdr

def idx_bam(bam, threads):
    threads = str(threads)
    try:
        pysam.index("-@"+threads,bam)
    except:
        outcome = 'idxerror'
    else:
        outcome = 'idxsuccess'
    if outcome == 'idxerror':
        print("indexing failed, trying to sort bam file...")
        inbam = bam
        bam = bam+".sorted.bam"
        pysam.sort("-@"+threads,"-o", bam, inbam)
        print("indexing bam file...")
        pysam.index(bam)
    return bam

def collect_bam_chunks(inpath, chrs, outpath, unmapped):
    allpaths = [inpath+".tmp."+c+".bam" for c in chrs]
    if unmapped:
        allpaths = [inpath+".tmp."+c+".bam" for c in chrs[:-1]]
        allpaths.append(inpath+".tmp."+"unmapped"+".bam")
    cat_args = ['-o', outpath]+allpaths
    pysam.cat(*cat_args)
    x = [os.remove(f) for f in allpaths]

def remove_tag(read, rtag):
    all_tags = read.get_tags()
    to_keep = [t[0] != rtag for t in all_tags]
    kept_tags = [tag for tag, keep in zip(all_tags, to_keep) if keep]
    read.set_tags(kept_tags)
    return read

def count_ref_consuming_bases(cigartuples):
    bases = 0
    for cig in cigartuples:
        if cig[0] in [0,2,7,8]:
            bases = bases+cig[1]
    return bases

def clean_bam(inpath, threads, fastapath, chr, strict, keepunmapped, keepsecondary, anonheader):
    fa = pysam.FastaFile(fastapath)

    if chr == '*':
        chrlabel = 'unmapped'
    else:
        chrlabel = chr

    #open in/out files
    outpath = inpath+".tmp."+chrlabel+".bam"
    tmppath = inpath+".tmp."+chrlabel+".tmp.bam"
    inp = pysam.AlignmentFile(inpath, 'rb', threads = threads)
    out = pysam.AlignmentFile(tmppath, 'wb', header = anonheader, threads = threads)
    for read in inp.fetch(chr):
        # deal with unmapped reads
        if chrlabel == 'unmapped':
            trim_tags = ['uT', 'nM', 'NM', 'XN', 'XM', 'XO', 'XG']
            if strict:
                trim_tags = trim_tags+['NH','HI','IH','AS','MQ','H1','H2','OA','OC','OP','OQ','SA','SM','XA','XS']
            for t in trim_tags:
                if read.has_tag(t):
                    read = remove_tag(read, t)
            out.write(read)
            continue

        #only use primary alignments
        if not keepsecondary and read.is_secondary:
            continue

        #determine some basics
        readlen = read.query_length
        qual = read.query_qualities
        if read.is_paired:
            readtype = 'PE'
            readtype_int = 2
        else:
            readtype = 'SE'
            readtype_int = 1

        #modify tags
        trim_tags = ['MC','XN','XM','XO','XG']
        for t in ['NM', 'nM']:
            if read.has_tag(t):
                read.set_tag(tag = t, value_type = 'I', value = 0)
        #if read.has_tag('MC'):
            #read = remove_tag(read,'MC')
            #read.set_tag(tag = 'MC', value_type = 'Z', value = str(readlen)+'M')
        if read.has_tag('MD'): #do we fix it or remove it?
            #read = read.remove_tag('MD')
            read.set_tag(tag = 'MD', value_type = 'Z', value = str(readlen))

        if strict:
            if read.has_tag('AS'):
                read.set_tag(tag = 'AS', value_type = 'I', value = readtype_int*readlen)
            if read.has_tag('MQ'):
                read.set_tag(tag = 'MQ', value_type = 'I', value = readtype_int*readlen)
            if read.has_tag('NH'):
                read.set_tag(tag = 'NH', value_type = 'I', value = 1)
            trim_tags = trim_tags+['HI','IH','H1','H2','OA','OC','OP','OQ','SA','SM','XA','XS']
            read.mapping_quality = 255

        for t in trim_tags:
            if read.has_tag(t):
                read = remove_tag(read, t)

        #some aligners like to keep the unmapped reads at the same posiiton as their mate, deal with them
        if read.is_unmapped:
            if keepunmapped:
                out.write(read)
            continue

        #look at cigar value
        incigar = read.cigartuples
        present_cigar_types = [x[0] for x in incigar]

        #if reads start with clipping or deletion, change the start position too
        #however, in the case of PE, changing the mapping position will mess up SAM fields RNEXT, PNEXT, TLEN
        #so this should only happen for SE
        if readtype == 'SE' and present_cigar_types[0] in [2, 4, 5]:
            fill_len = incigar[0][1]
            new_start = read.reference_start-(fill_len+1)
            if new_start < 1: #take care that things dont go out of range
                new_start = 1
            read.reference_start = new_start
            if present_cigar_types[1] == 0: #next segment is mapped, so add the mapped length to cigar
                incigar[1] = (incigar[1][0], incigar[1][1] + fill_len)
                present_cigar_types, incigar = present_cigar_types[1:], incigar[1:]
            else: #otherwise set first segment to M
                present_cigar_types[0], incigar[0] = 0, (0, incigar[0][1])

        if 3 not in present_cigar_types: #unspliced alignment, simply fix sequence
            final_outseq = fa.fetch(chr, read.reference_start, read.reference_start+readlen)
            final_cigar = [(0, readlen)]
        else: #spliced alignment
            splice_fields = [x == 3 for x in present_cigar_types]
            splice_field_idx = [idx for idx, splice in enumerate(splice_fields) if splice]
            num_splices = sum(splice_fields)
            splicecigar = [incigar[idx] for idx in splice_field_idx]

            #save info for each aligned exon piece in a list
            outsegments = {}
            for segment in range(num_splices+1):
                if segment == 0: #first segment
                    remaininglen = readlen
                    seglength = count_ref_consuming_bases(incigar[:splice_field_idx[segment]])
                    segstartpos = read.reference_start
                elif segment == num_splices: #last segment, just put the
                    seglength = remaininglen
                    segstartpos = read.reference_start+sum([l[1][2] for l in outsegments.items()])+sum([l[1] for l in splicecigar]) #sum of all so far mapped lenghts and skipped intron lengths
                else: #segments between splices
                    seglength = count_ref_consuming_bases(incigar[splice_field_idx[segment-1]+1:splice_field_idx[segment]])
                    segstartpos = read.reference_start+sum([l[1][2] for l in outsegments.items()])+sum([l[1] for l in splicecigar[:segment]]) #sum of all so far mapped lenghts and so far skipped intron lengths
                segendpos = segstartpos + seglength
                outsegments[segment] = (segstartpos, segendpos, seglength)
                remaininglen = remaininglen - seglength
                if (remaininglen < 0) or (remaininglen == 0 and segment < num_splices): #in case deletions we filled up consumed so much sequence length that we do not make it to the next splice
                    seglength = seglength + remaininglen #substract the overlap
                    segendpos = segendpos + remaininglen #substract the overlap
                    outsegments[segment] = (segstartpos, segendpos, seglength) #overwrite output
                    splicecigar = splicecigar[:segment] #keep only splices dealt with so far
                    break #don't forget to break the loop

            #combine ouput for bam record
            outseq = [fa.fetch(chr, outsegments[idx][0], outsegments[idx][1]) for idx in outsegments]
            outcigar = [(0, outsegments[idx][2]) for idx in outsegments]
            combined_cigar = [field for pair in itertools.zip_longest(outcigar, splicecigar) for field in pair if field is not None]
            #set for output read
            final_outseq = ''.join(outseq)
            final_cigar = combined_cigar

        if len(final_outseq) != len(qual): #the sanitized output sequence cannot be longer than a contig (reason:deletions)
            len_diff = len(qual)-len(final_outseq)
            old_len = final_cigar[-1][1]
            new_len = old_len - len_diff
            final_cigar[-1] = (0,new_len)
            qual = qual[:len(final_outseq)]
        #set alignment record and write
        read.query_sequence = final_outseq
        read.query_qualities = qual
        read.cigartuples = final_cigar
        out.write(read)
    inp.close()
    out.close()

    #resorting SE reads is necessary since we may mess with read start positions
    try:
        test = readtype
    except NameError:
        readtype = 'NA'
    if readtype == 'SE':
        pysam.sort("-o", outpath, tmppath)
        if os.path.exists(tmppath):
            os.remove(tmppath)
    else:
        os.rename(tmppath, outpath)


def main():
    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument('--bam', type=str, metavar='FILENAME',
                        help='Path to input BAM file', required = True)
    parser.add_argument('--out', type=str, metavar='FILENAME',
                        help='Path to output bam file', required = True)
    parser.add_argument('--fa', type=str, metavar='FILENAME',
                        help='Path to genome reference fasta', required = True)
    parser.add_argument('--p', type=int, default = 10,
                        help='Number of processes to use')
    parser.add_argument('--strict', action = 'store_true',
                        help='Strict: also sanitize mapping score & auxiliary tags (eg. AS / NH).')
    parser.add_argument('--keepsecondary', action = 'store_true',
                        help='Keep secondary alignments in output bam file.')
    parser.add_argument('--keepunmapped', action = 'store_true',
                        help='Keep ummapped reads in output bam file.')


    args = parser.parse_args()
    v = '0.5.0'
    print("BAMboozle.py v"+v)
    bampath = args.bam
    try:
        fa = pysam.FastaFile(args.fa)
    except ValueError:
        print("Error: Reference fasta file is not indexed!")
        print("Please run: samtools faidx "+args.fa)
        quit()

    if not os.path.exists(bampath+'.bai'):
        print("input bam index not found, indexing...")
        bampath = idx_bam(bampath,args.p)

    #Construct the new bam header to work with
    bamheader = makeBAMheader(args, v)
    print("Working...")

    chrs = pysam.idxstats(bampath).split('\n')
    chrs = [c.split('\t')[0] for c in chrs[:-1]]
    chrs = [c for c in chrs if c in fa.references] #we can only sanitize contigs that have a reference seq
    if args.keepunmapped:
        chrs.append('*')
    fa.close()

    if args.p > 20:
        pysam_workers = 2
        n_jobs = int(args.p/2)
    else:
        pysam_workers = 1
        n_jobs = args.p

    pool = mp.Pool(n_jobs)
    results = [pool.apply_async(clean_bam, (args.bam,pysam_workers,args.fa,chr,args.strict,args.keepunmapped,args.keepsecondary,bamheader,  )) for chr in chrs]
    x = [r.get() for r in results]
    #single threaded below:
    #[clean_bam(bampath,pysam_workers,args.fa,chr,args.strict) for chr in chrs]

    print("Creating final output .bam file...")
    collect_bam_chunks(inpath = bampath, chrs = chrs, outpath = args.out, unmapped = args.keepunmapped)
    print("Indexing final output .bam file...")
    y = idx_bam(args.out,args.p)

    print("Done!")

if __name__ == "__main__":
    main()
