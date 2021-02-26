#!/usr/bin/env python
# -*- coding: utf-8

__author__ = "Laura G. Macias"
__email__ = "laugmacias@gmail.com"

from argparse import ArgumentParser
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation

# Check whether a gen is correctly annotated:
# this functions checks first and last codon
def gen_ok(nts):
    if "N" not in nts:
        if nts.startswith("ATG"):
            if nts.endswith("TAA") or nts.endswith("TAG") or nts.endswith("TGA"):
                if len(nts) % 3 == 0:
                    return True
    else:
        return False

def main():
    """
    This script extracts a gene sequence from a Genbank file and its promoter region.

    The name of the gene sequence to be extracted must be provided using the -name argument.
    If the gene name is stored in a GBK qualifier different than the "gene" qualifier (e.g.
    locus_tag) must be indicated using the -qual option.

    Output: a fasta file is created (geneName_outName.fna) with two sequences. The first sequence
    correspond to the gene (ORF) and the second to the promoter region (ORF +- 1000 bp)
    
    """
    parser = ArgumentParser(description=main.__doc__)
    parser.add_argument("-gbk", dest="genbank", help="annotated genbank file", type=str)
    parser.add_argument("-name", dest="gene_name", help="name of the gene sequence", type=str)
    parser.add_argument("-qual", dest="tag", \
                        help="gene/locustag qualifier of the gene", type=str, default="gene")
    parser.add_argument("-out", dest="out", help="name of the output e.g. spp/strain name", type=str)
    args = parser.parse_args()

    # Parse genbank file
    annot = SeqIO.parse(args.genbank,"genbank")
    found = False
    for rec in annot:
        for feat in rec.features:
            if feat.type == "CDS":
                bases = feat.location.extract(rec.seq)
                if gen_ok(bases):
                    if args.tag == "locustag":
                        if "locus_tag" in feat.qualifiers.keys():
                            gen = feat.qualifiers["locus_tag"][0]
                            if gen == args.gene_name:
                                output = open(args.gene_name+"_"+args.out+".fna", "a")
                                found = True
                                start = int(feat.location.start - 1000)
                                end = int(feat.location.end + 1000)
                                promoter_loc = FeatureLocation(start,end,strand=feat.location.strand)

                                # write gene seq to file
                                sequence_object = Seq(str(bases))
                                record = SeqRecord(sequence_object, id=args.gene_name+"_"+args.out, description="")
                                SeqIO.write(record, output, "fasta")

                                # write gene seq + promoter to file
                                sequence_promoter = Seq(str(promoter_loc.extract(rec.seq)))
                                record_promoter = SeqRecord(sequence_promoter, id=args.gene_name+"+/-1000bp"+"_"+args.out, description="")
                                SeqIO.write(record_promoter, output, "fasta")
                                output.close()
                        

                    else:
                        if "gene" in feat.qualifiers.keys():
                            gen = feat.qualifiers["gene"][0]
                            if gen == args.gene_name:
                                output = open(args.gene_name+"_"+args.out+".fna", "a")
                                found = True
                                start = int(feat.location.start - 1000)
                                end = int(feat.location.end + 1000)
                                promoter_loc = FeatureLocation(start,end,strand=feat.location.strand)

                                # write gene seq to file
                                sequence_object = Seq(str(bases))
                                record = SeqRecord(sequence_object, id=args.gene_name+"_"+args.out, description="")
                                SeqIO.write(record, output, "fasta")

                                # write gene seq + promoter to file
                                sequence_promoter = Seq(str(promoter_loc.extract(rec.seq)))
                                record_promoter = SeqRecord(sequence_promoter, id=args.gene_name+"+/-1000bp"+"_"+args.out, description="")
                                SeqIO.write(record_promoter, output, "fasta")
                                output.close()





    if not found:
        print("gene not found")


if __name__ == "__main__":
    main()
