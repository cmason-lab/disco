#!/usr/bin/env python
import argparse
import pkg_resources
from disco import *
import os


def main():
    parser = argparse.ArgumentParser(prog='disco',
                                     description='DISCO: Distributions of Isoforms in Single Cell Omics',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s '+pkg_resources.require("disco")[0].version)
    parser.add_argument('sampleannfile', metavar='SampleAnnotationFile', type=str,
                        help='filename of tab separated text, no header, with columns: '
                             '<path to miso summary file> <sample name> <group name>')
    parser.add_argument('group1', metavar='Group1', type=str,
                        help='must match a group name in sample annotation file')
    parser.add_argument('group2', metavar='Group2', type=str,
                        help='must match a group name in sample annotation file')
    parser.add_argument('--outdir', metavar="", dest='outdir', type=str,
                        default="./disco_output/",
                        help="Output directory")
    parser.add_argument('--pkldir', metavar="", dest='pkldir', type=str,
                        default="./pkldir",
                        help="Directory to store intermediate data processing files")
    parser.add_argument('--group1color', metavar="", dest='group1color', type=str,
                        default="r",
                        help="Color in plots for group 1; can be {y, m, c, r, g, b, w, k} or html code")
    parser.add_argument('--group2color', metavar="", dest='group2color', type=str,
                        default="b",
                        help="Color in plots for group 2; can be {y, m, c, r, g, b, w, k} or html code")
    parser.add_argument('--group1file', metavar="", dest='group1file', type=str,
                        help="output file for sample group 1. If not specified, "
                             "will save to <outdir>/<group1name>_alldatadf.txt")
    parser.add_argument('--group2file', metavar="", dest='group2file', type=str,
                        help="output file for sample group 2. If not specified, "
                             "will save to <outdir>/<group2name>_alldatadf.txt")
    parser.add_argument('--geneannotationfile', metavar="", dest='geneannotationfile', type=str,
                        default=None,
                        help="Mapping of Ensembl gene IDs to HGNC symbol and gene descriptions")
    parser.add_argument('--transcriptannotationfile', metavar="", dest='transcriptannotationfile', type=str,
                        default=None,
                        help="Mapping of Ensembl transcript IDs to isoform function (ex. protein coding, NMD, etc)")
    parser.add_argument('--maxciwidth', metavar="", dest='maxciw', type=float,
                        default=1.0,
                        help="Maximum width of confidence interval of PSI estimate")
    parser.add_argument('--mininfreads', metavar="", dest='mininfreads', type=int,
                        default=0,
                        help="Minimum number of informative reads to include PSI estimate")
    parser.add_argument('--mindefreads', metavar="", dest='mindefreads', type=int,
                        default=0,
                        help="Minimum number of definitive reads to include PSI estimate")
    parser.add_argument('--minavgpsi', metavar="", dest='minavgpsi', type=float,
                        default=0.0,
                        help="Do not run statistical tests for isoforms with average PSI "
                             "in both groups less than minavgpsi")
    parser.add_argument('--minnumcells', metavar="", dest='minnumcells', type=int,
                        default=0,
                        help="Do not run statistical test for isoform if less than minnumcells have information")
    parser.add_argument('--minavgshift', metavar="", dest='minavgshift', type=float,
                        default=0,
                        help="Do not run statistical test for isoform if shift in mean PSI between the two groups is "
                             "less than minavgshift")
    parser.add_argument('--stattest', metavar="", dest='stattest', choices=["KS", "T"], type=str,
                        default="KS",
                        help="Which test to run? Options: {'KS', 'T'}")  # options kstest or ttest
    parser.add_argument('--alpha', metavar="", dest='alpha', type=float,
                        default=0.05,
                        help="Adjustded p-value threshold for statistical significance")
    parser.add_argument('--multitestmethod', metavar="", dest='multitestmethod',
                        choices=['bonferroni', 'sidak', 'holm-sidak', 'holm', 'simes-hochberg', 'hommel',
                                 'fdr_bh', 'fdr_by', 'fdr_tsbh', 'fdr_tsbky', 'none'], type=str, default='fdr_bh',
                        help="Method for multiple testing correction. Options: "
                             "{'bonferroni', 'sidak', 'holm-sidak', 'holm', 'simes-hochberg', 'hommel', "
                             "'fdr_bh', 'fdr_by', 'fdr_tsbh', 'fdr_tsbky', 'none'}. "
                             "See statsmodels.stats.multitest.multipletests for more info")

    args = parser.parse_args()

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    if not os.path.exists(args.pkldir):
        os.makedirs(args.pkldir)

    group1out = args.outdir+"/"+args.group1+"_alldatadf.txt" if args.group1file is None else args.group1file
    group2out = args.outdir+"/"+args.group2+"_alldatadf.txt" if args.group2file is None else args.group2file
    disco1 = Disco(args.sampleannfile, args.group1, group1out, args.pkldir)
    disco2 = Disco(args.sampleannfile, args.group2, group2out, args.pkldir)

    statres = stat_test(disco1,
                        disco2,
                        args.outdir+"/"+args.group1+"vs"+args.group2+"_"+args.stattest+"_allresults.txt",
                        maxciw=args.maxciw,
                        mininfreads=args.mininfreads,
                        mindefreads=args.mindefreads,
                        minavgpsi=args.minavgpsi,
                        minnumcells=args.minnumcells,
                        minavgshift=args.minavgshift,
                        multitestmethod=args.multitestmethod,
                        alpha=args.alpha,
                        geneannfile=args.geneannotationfile,
                        transcriptannfile=args.transcriptannotationfile,
                        testtype=args.stattest)

    sigks = getsig(statres,
                   outfile=args.outdir+"/"+args.group1+"vs"+args.group2+"_"+args.stattest+"_significantresults.txt")

    plotsig_violin(sigks, statres, disco1, disco2,
                   args.outdir+"/"+args.group1+"vs"+args.group2+"_"+args.stattest+"_violinplots.pdf",
                   args.group1, args.group2, args.group1color, args.group2color)


if __name__ == '__main__':
    main()
