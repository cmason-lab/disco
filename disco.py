#!/usr/bin/env python
import pandas as pd
import glob
import os
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.backends.backend_pdf as pdf
import argparse

# import multiprocessing as mp

sns.set_context("talk", font_scale=1.5)
sns.set_style("darkgrid")


class Disco:
    def __init__(self, samples, datapath, outfile="alldatadf.txt", pkldir="picklefiles_isocentric"):
        """

        :type datapath: str
        :type samples: list
        """
        self.samples = samples
        self.datapath = datapath
        self.outfile = outfile
        self.filenames = None
        self.cellnames = None
        self.cellsfailed = []
        self.pkldir = pkldir
        if os.path.exists(outfile):
            print "processed data exists, reading from", outfile
            datadf = pd.DataFrame.from_csv(outfile, sep="\t")
            # reading in changes colname gene!isoform to gene!isoform.1 for some reason, change back:
            cols = list(datadf.columns)
            if "gene!isoform.1" in cols:
                cols[cols.index("gene!isoform.1")] = "gene!isoform"
            datadf.columns = cols
            self.alldatadf = datadf
        else:
            self.alldatadf = self._process()
        # todo add transcript type annotation

    def _process(self):
        filenames = []
        for sample in self.samples:
            filenames += glob.glob(self.datapath + "/" + sample + "*")
        # print filenames
        # filenames.remove("../miso/summaryfiles/MDS4-B12_miso-output.miso_summary")  # temp workaround for unknown bug
        self.filenames = filenames
        cellnames = [i.split("/")[len(i.split("/")) - 1].split("_")[0] for i in filenames]
        self.cellnames = cellnames
        pklfiles = [self.readsummfile(i) for i in range(len(filenames))]
        results = []
        for p in pklfiles:
            if p is None:
                continue
            else:
                results.append(pd.read_pickle(p))
        # results = [pd.read_pickle(p) for p in pklfiles]
        resultsdf = pd.concat(results, axis=0)
        resultsdf.index = resultsdf["gene!isoform"]
        # resultsdf["ciwidth_i"] = resultsdf["cihigh_i"] - resultsdf["cilow_i"]
        resultsdf.to_csv(self.outfile, sep="\t")
        return resultsdf

    def readsummfile(self, sampindex):
        filename = self.filenames[sampindex]
        cellname = self.cellnames[sampindex]
        picklefile = self.pkldir + "/" + cellname + ".pkl"
        print cellname
        if os.path.exists(picklefile):
            return picklefile
        # df1 = pd.read_table(filename, sep="\t")
        # print df1["event_name"].value_counts()
        # print df1.shape
        # df1["cellname"] = pd.Series(np.repeat(cellname, df1.shape[0]), index=df1.index)
        # df2 = df1.apply(self._readhelper1, 1)
        # print df2.shape
        # df3 = list(df2.apply(self._readhelper2, 1))
        # df4 = pd.concat(df3)
        # print df4.shape
        # df4.to_pickle(picklefile)
        # todo fix bug that causes some cells to fail and remove try/catch
        try:
            df1 = pd.read_table(filename, sep="\t")
            print df1.shape
            df1["cellname"] = pd.Series(np.repeat(cellname, df1.shape[0]), index=df1.index)
            df2 = df1.apply(self._readhelper1, 1)
            print df2.shape
            df3 = list(df2.apply(self._readhelper2, 1))
            df4 = pd.concat(df3)
            print df4.shape
            df4.to_pickle(picklefile)
        except:
            print cellname, "failed"
            self.cellsfailed.append(cellname)
            return None
        print
        return picklefile

    @staticmethod
    def _readhelper1(x):
        gene = x["event_name"]
        # print gene
        isfs = x["isoforms"].split(",")
        rownames = []
        isfshortnames = []
        for i in range(len(isfs)):
            rownames.append(gene + "!" + isfs[i].strip("\'"))
            isfshortnames.append("isf-"+str(i+1))
        psis = map(float, x["miso_posterior_mean"].split(","))
        cilows = map(float, x["ci_low"].split(","))
        cihighs = map(float, x["ci_high"].split(","))
        # todo add isf chromosome start stops to long format
        if len(psis) == 1:
            psis.append(1 - psis[0])
            # extrapolate confidence interval to psi of next isoform by maintaining width and position within range
            cilows.append(psis[1]-(psis[0]-cilows[0]))
            cihighs.append(psis[1]+(cihighs[0]-psis[0]))
        # calculate num reads informative and definitive
        ciwidths = list(np.array(cihighs) - np.array(cilows))
        readlist = x["counts"].split("(")
        numreadsinf = 0
        numreadsdef = 0
        for seg in readlist[1:]:
            segsplit = seg.strip(",").split(":")
            sumisfcode = sum(map(int, segsplit[0].strip(")").split(",")))
            if sumisfcode > 0:
                numreadsinf += int(segsplit[1])
            if sumisfcode == 1:
                numreadsdef += int(segsplit[1])
        y = []
        for i in range(len(rownames)):
            y.append(x.append(pd.Series([rownames[i], isfshortnames[i], psis[i], cilows[i], cihighs[i], ciwidths[i],
                                         numreadsinf, numreadsdef],
                                        index=["gene!isoform", "isfshortname", "psi_i", "cilow_i", "cihigh_i",
                                               "ciwidth_i", "numreadsinf", "numreadsdef"])))
        return y

    @staticmethod
    def _readhelper2(x):
        # print len(x)
        return pd.concat(x, axis=1).T

    def plotgenehist(self, gene):
        genedf = self.alldatadf[self.alldatadf["event_name"] == gene]
        isomat = genedf.pivot(index="cellname", columns="gene!isoform", values="psi_i")
        newcols = ["isf-"+str(i) for i in range(1, isomat.shape[1]+1)]
        isomat.columns = newcols
        ax = isomat.plot(kind="hist", alpha=0.75, xlim=(0, 1))
        plt.xlabel("PSI (% spliced isoform)")
        plt.ylabel("Number of cells")
        plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
        plt.show()
        # todo only keep top 7
        return ax


# todo put comparison functions in a class

# todo test only 1 isoform when 2 exist

# Functions for filtering and comparison
def stat_test(disco1, disco2, outfile, maxciw=1, mininfreads=0, mindefreads=0, minavgpsi=0, minnumcells=0,
              annotationfile=None, testtype="KS"):
    """
    testtype: {KS, T}
    possible filters:
    1. width of confidence interval < maxciw
    2. number of informative reads >= mininfreads
    3. number of definitve reads >= mindefreads
    4. average psi >= minavgpsi in both groups
    5. number of cells called >= numcells
    :type disco1, disco2: object of Disco class
    """
    samp1datadf = disco1.alldatadf
    print "Sample 1:", samp1datadf.shape
    samp1psis = filterdf(samp1datadf, maxciw, mininfreads, mindefreads, minnumcells)
    print "Filtered Sample 1:", samp1psis.shape
    print
    samp2datadf = disco2.alldatadf
    print "Sample 2:", samp2datadf.shape
    samp2psis = filterdf(samp2datadf, maxciw, mininfreads, mindefreads, minnumcells)
    print "Filtered Sample 2:", samp2psis.shape
    print
    intersectkeys = set(samp1psis.index) & set(samp2psis.index)
    print "Number of shared isoforms:", len(intersectkeys)

    samp1kept = samp1psis.loc[intersectkeys, ]
    samp2kept = samp2psis.loc[intersectkeys, ]

    joineddf = samp1kept.join(samp2kept)
    newkeys = filterminavgpsi(joineddf, minavgpsi)
    print "After filtering for min avg. psi:", len(newkeys)

    samp1kept2 = samp1kept.loc[newkeys, ]
    samp2kept2 = samp2kept.loc[newkeys, ]

    if testtype == "KS":
        print "Running KS test"
        statsres = samp1kept2.apply(lambda x: pd.Series(stats.ks_2samp(x, samp2kept2.loc[x.name, ]),
                                                        index=["K", "pvalue"], name=x.name), axis=1)
    else:
        print "Running t test"
        # Mask NAs
        maskmat1 = np.isfinite(samp1kept2)
        maskmat2 = np.isfinite(samp2kept2)
        # statsres = samp1kept2.apply(_ttest, axis=1, args=(samp2kept2, maskmat1, maskmat2))
        statsres = samp1kept2.apply(lambda x: pd.Series(stats.ttest_ind(x[maskmat1.loc[x.name]],
                                                                        samp2kept2.loc[x.name][maskmat2.loc[x.name]]),
                                                        index=["T", "pvalue"], name=x.name), axis=1)
    #  append median shift to statsres DF
    medianshifts = samp1kept2.median(axis=1) - samp2kept2.median(axis=1)
    statsres2 = statsres.apply(lambda x: x.append(pd.Series([medianshifts.loc[x.name]],
                                                            index=["median_shift"], name=x.name)), axis=1)

    if annotationfile is None:
        print
        print "Final result:", statsres.shape
        return statsres
    annotationdf = pd.DataFrame.from_csv(annotationfile, sep="\t", index_col=9)
    statsres_ann = statsres2.apply(_getannotation, axis=1, args=(annotationdf, "d"))

    # todo add isf shortname to this DF
    # statsres_ann["isfshortname"] = samp1datadf.loc[statsres_ann.index]["isfshortname"]

    statsres_ann.to_csv(outfile, sep="\t")
    print
    print "Final result:", statsres_ann.shape
    return statsres_ann


def _getannotation(x, annotationdf, dummy):
    ensid = x.name.split("!")[0]
    if ensid in annotationdf.index:
        annrow = annotationdf.loc[ensid]
        y = pd.Series([ensid, annrow["Approved Symbol"], annrow["Approved Name"], annrow["Chromosome"]],
                      index=["Ensemble_ID", "Gene_Symbol", "Gene_Name", "Locus"],
                      name=x.name)
    else:
        y = pd.Series([ensid, None, None, None],
                      index=["Ensemble_ID", "Gene_Symbol", "Gene_Name", "Locus"],
                      name=x.name)
    return x.append(y)


def filterdf(alldatadf, maxciw=1, mininfreads=0, mindefreads=0, minnumcells=0):
    """ defaults are no filtering """
    filt1 = alldatadf[alldatadf["ciwidth_i"] <= maxciw]
    print "After filtering for conf. int. width:", filt1.shape
    filt2 = filt1[filt1["numreadsinf"] >= mininfreads]
    print "After filtering for num reads inf.:", filt2.shape
    filt3 = filt2[filt2["numreadsdef"] >= mindefreads]
    print "After filtering for num reads def.:", filt3.shape
    filt3piv = filt3.pivot(index="gene!isoform", columns="cellname", values="psi_i")
    print "# isoforms x # cells:", filt3piv.shape
    numcells = filt3piv.apply(lambda x: pd.Series(x > 0).sum(), axis=1)
    keystokeep2 = numcells[numcells > minnumcells].index
    final = filt3piv.loc[keystokeep2]
    print "After keeping rows with min num cells called:", final.shape
    return final


def filterminavgpsi(joineddf, minavgpsi=0):
    avgpsi = joineddf.mean(axis=1)
    keystokeep = avgpsi[avgpsi > minavgpsi].index
    return keystokeep


def getsig(statsres, pvalue=0.05, bonferroni=True, minmedianshift=0, outfile=None):
    thresh = pvalue/statsres.shape[0] if bonferroni else pvalue
    # significant
    sigres = statsres[statsres["pvalue"] < thresh]
    print "Significant at pvalue", pvalue, "and bonferroni", bonferroni, ":", sigres.shape
    toreturn = sigres
    if sigres.shape[0] > 0:
        sigres2 = sigres[sigres["median_shift"].abs() > minmedianshift]
        print "Number of sites with median shift >", minmedianshift, ":", sigres2.shape
        toreturn = sigres2
    if outfile is not None:
        toreturn.to_csv(outfile, sep="\t", index_label="gene!isoform")
    return toreturn


def _renamer(x, sigkeys, kskeys, enstdf):
    isfshort = x["isfshortname"]
    geneisf = x["gene!isoform"]
    if geneisf in sigkeys:
        # newlabel = isfshort+"*"
        newlabel = isfshort.strip("isf-")+"*"
    elif geneisf in kskeys:
        # newlabel = isfshort
        newlabel = isfshort.strip("isf-")+"-"
    else:
        # newlabel = isfshort+"/"
        newlabel = isfshort.strip("isf-")
    return x.append(pd.Series([newlabel], index=["newlabel"]))


# def _renamer(x, sigkeys, kskeys, enstdf):
#     # isfshort = x["isfshortname"]
#     geneisf = x["gene!isoform"]
#     enst = geneisf.split("!")[1].split(".")[0]
#     if enst in enstdf.index:
#         transcripttype = enstdf.loc[enst]["Transcript type"]
#     else:
#         transcripttype = ""
#     if geneisf in sigkeys:
#         # newlabel = isfshort+"*"
#         # newlabel = isfshort.strip("isf-")+"*"
#         newlabel = enst+"*"
#     elif geneisf in kskeys:
#         # newlabel = isfshort
#         # newlabel = isfshort.strip("isf-")+"-"
#         newlabel = enst+"-"
#     else:
#         # newlabel = isfshort+"/"
#         # newlabel = isfshort.strip("isf-")
#         newlabel = enst
#     print newlabel+" "+transcripttype
#     return x.append(pd.Series([newlabel+" "+transcripttype], index=["newlabel"]))


def plotsig_violin(sigks, statsres, sciso1, sciso2, outfile, sampname1, sampname2, enstfile):
    enstdf = pd.DataFrame.from_csv(enstfile, sep="\t", index_col=1)
    sns.set(style="dark", palette="muted", color_codes=True, font_scale=1.5)
    plotspdf = pdf.PdfPages(outfile)
    alldatadf1 = sciso1.alldatadf
    alldatadf2 = sciso2.alldatadf
    alldatadf1["group"] = pd.Series(np.repeat(sampname1, alldatadf1.shape[0]), index=alldatadf1.index)
    alldatadf2["group"] = pd.Series(np.repeat(sampname2, alldatadf2.shape[0]), index=alldatadf2.index)
    genestoplot = sigks["Ensemble_ID"].unique()
    print "Plotting ", len(genestoplot), " genes"
    for gene in genestoplot:
        # print gene
        genedf1 = alldatadf1[alldatadf1["event_name"] == gene]
        genedf2 = alldatadf2[alldatadf2["event_name"] == gene]
        merged = pd.concat([genedf1, genedf2])
        merged2 = merged.apply(_renamer, axis=1, args=(sigks.index, statsres.index, enstdf))
        print merged2.shape
        sns.violinplot(x="newlabel", y="psi_i", hue="group", data=merged2, split=True,
                       inner="quart", cut=0, palette={sampname1: "r", sampname2: "b"})
        annrow = sigks[sigks["Ensemble_ID"] == gene].iloc[0]
        if annrow["Gene_Symbol"] is None:
            plt.title(gene)
            print gene
        else:
            plt.title(annrow["Gene_Symbol"]+"\n"+annrow["Gene_Name"]+"; "+annrow["Locus"])
            print annrow["Gene_Symbol"]
        plt.ylabel("Distribution of Isoform Expression")
        plt.xlabel("Isoform ( - tested, * significant )")
        plt.legend(loc='upper left', bbox_to_anchor=(1.0, 0.9))
        plt.subplots_adjust(top=0.88, right=0.8, bottom=0.12)
        plotspdf.savefig()
        plt.close()
    plotspdf.close()
    print "Saved plots to", outfile
    sns.set()
    return


def plotsig_hist(sigks, statsres, disco1, disco2, outfile, sampname1, sampname2):
    pass


def _renamer1samp(x, enstdf, dummy):
    isfshort = x["isfshortname"]
    # geneisf = x["gene!isoform"]
    # if geneisf in sigkeys:
    #     # newlabel = isfshort+"*"
    #     newlabel = isfshort.strip("isf-")+"*"
    # elif geneisf in kskeys:
    #     # newlabel = isfshort
    #     newlabel = isfshort.strip("isf-")+"-"
    # else:
        # newlabel = isfshort+"/"
    newlabel = isfshort.strip("isf-")
    return x.append(pd.Series([newlabel], index=["newlabel"]))


def plotviolin_1samp(disco, genestoplot, outfile, enstfile):
    enstdf = pd.DataFrame.from_csv(enstfile, sep="\t", index_col=1)
    sns.set(style="dark", palette="muted", color_codes=True, font_scale=1.5)
    plotspdf = pdf.PdfPages(outfile)
    alldatadf = disco.alldatadf
    dummy=None
    # genestoplot = sigks["Ensemble_ID"].unique()
    print "Plotting ", len(genestoplot), " genes"
    for gene in genestoplot:
        # print gene
        genedf = alldatadf[alldatadf["event_name"] == gene]
        print genedf
        genedf2 = genedf.apply(_renamer1samp, axis=1, args=(enstdf, dummy))
        sns.violinplot(x="newlabel", y="psi_i", data=genedf2,
                       inner="quart", cut=0)
        plt.title(genedf)
        # annrow = sigks[sigks["Ensemble_ID"] == gene].iloc[0]
        # if annrow["Gene_Symbol"] is None:
        #     plt.title(gene)
        #     print gene
        # else:
        #     plt.title(annrow["Gene_Symbol"]+"\n"+annrow["Gene_Name"]+"; "+annrow["Locus"])
        #     print annrow["Gene_Symbol"]
        plt.ylabel("Distribution of Isoform Expression")
        plt.xlabel("Isoform")
        # plt.legend(loc='upper left', bbox_to_anchor=(1.0, 0.9))
        # plt.subplots_adjust(top=0.88, right=0.8, bottom=0.12)
        plotspdf.savefig()
        plt.close()
    plotspdf.close()
    print "Saved plots to", outfile
    sns.set()
    return
    pass


def main():

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
                        default="t-test", help="Which statistical test to run: ks test or t test (default)?")  # options kstest or ttest

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

# def _medianshift(x, samp1psis, samp2psis):
#     shift = samp1psis.loc[x.name].
# def _ttest(x, samp2kept2, maskmat1, maskmat2):
#     totest1 = x[maskmat1.loc[x.name]]
#     totest2 = samp2kept2.loc[x.name][maskmat2.loc[x.name]]
#     return stats.ttest_ind(totest1, totest2, index=["T", "pvalue"], name=x.name)
