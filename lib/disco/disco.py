def main():
    # from disco import *
    import argparse

    parser = argparse.ArgumentParser(description='Run disco', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('samp1list', metavar='SampleGroup1', type=str,
                        help='comma delimited list of group 1 sample names')
    parser.add_argument('samp2list', metavar='SampleGroup2', type=str,
                        help='comma delimited list of group 2 sample names')
    parser.add_argument('--datapath', metavar="", dest='datapath', type=str,
                        help='path with miso summary files; searching in /datapath/sample*')
    parser.add_argument('-g1', '--group1name', metavar="", dest='group1name', type=str,
                        help="Name for group 1 to use in tables and plots")
    parser.add_argument('-g2', '--group2name', metavar="", dest='group2name', type=str,
                        help="Name for group 2 to use in tables and plots")
    parser.add_argument('-o1', '--outfile1', metavar="", dest='outfile1', type=str,
                        help="output file for sample group 1")
    parser.add_argument('-o2', '--outfile2', metavar="", dest='outfile2', type=str,
                        help="output file for sample group 2")
    parser.add_argument('--comppref', metavar="", dest='comppref', type=str,
                        help="Prefix for naming files for this comparison")
    parser.add_argument('--pkldir', metavar="", dest='pkldir', type=str, default="./",
                        help="Directory to store per sample pickle files")
    parser.add_argument('--annotationfile', metavar="", dest='annotationfile', type=str,
                        help="File with gene annotations")
    parser.add_argument('--maxciwidth', metavar="", dest='maxciw', type=float, default=1.0,
                        help="Filter out PSI estimates with width of confidence greater than maxciwidth")
    parser.add_argument('--mininfreads', metavar="", dest='mininfreads', type=int, default=0,
                        help="Minimum number of informative reads to include PSI estimate")
    parser.add_argument('--mindefreads', metavar="", dest='mindefreads', type=int, default=0,
                        help="Minimum number of definitive reads to include PSI estimate")
    parser.add_argument('--minavgpsi', metavar="", dest='minavgpsi', type=float, default=0.0,
                        help="Do not run statistical tests for isoforms with average PSI "
                             "in both groups less than minavgpsi")
    parser.add_argument('--minnumcells', metavar="", dest='minnumcells', type=int, default=0,
                        help="Do not run statistical test for isoform if less than minnumcells have information")
    parser.add_argument('--minmedianshift', metavar="", dest='minmedianshift', type=float, default=0,
                        help="Do not run statistical test for isoform if shift in median between the two groups is "
                             "less than minmedianshift")
    parser.add_argument('--stattest', metavar="", dest='stattest', choices=["KS", "T"], type=str,
                        default="KS", help="Which statistical test to run? options: {KS, T}")  # options kstest or ttest

    # default, choices, help
    args = parser.parse_args()

    print args
    disco1 = Disco(args.samp1list.split(","), args.datapath, args.outfile1, args.pkldir)
    disco2 = Disco(args.samp2list.split(","), args.datapath, args.outfile2, args.pkldir)
    statres = stat_test(disco1, disco2, args.comppref+"_"+args.stattest+"res.txt",
                        args.maxciw, args.mininfreads, args.mindefreads, args.minavgpsi,
                        args.minnumcells, args.annotationfile, args.stattest)
    sigks = getsig(statres, minmedianshift=args.minmedianshift, outfile=args.comppref+"_sig"+args.stattest+"test.txt")
    plotsig_violin(sigks, statres, disco1, disco2, args.comppref+"_sig"+args.stattest+"test_violinplots.pdf",
                   args.group1name, args.group2name)

if __name__ == '__main__':
    main()
