#!/usr/bin/env python
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt

# import multiprocessing as mp


class Disco:
    def __init__(self, sampleannfile, groupname, outfile="alldatadf.txt", pkldir="picklefiles_isocentric"):
        """
        """
        sampleann = pd.DataFrame.from_csv(sampleannfile, sep="\t", header=None)
        groupdf = sampleann[sampleann[2] == groupname]
        # print sampleann.head()
        # print groupdf.head()
        self.sampleann = groupdf
        self.outfile = outfile
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

    def _process(self):
        # filenames = list(self.sampleann.index)
        # cellnames = list(self.sampleann[1])
        # print filenames
        # filenames.remove("../miso/summaryfiles/MDS4-B12_miso-output.miso_summary")  # temp workaround for unknown bug

        # self.cellnames = cellnames
        pklfiles = [self.readsummfile(i) for i in range(self.sampleann.shape[0])]
        # print self.sampleann.shape
        results = []
        for p in pklfiles:
            if p is None:
                continue
            else:
                results.append(pd.read_pickle(p))
        # results = [pd.read_pickle(p) for p in pklfiles]
        # if len(results) > 1:
        #     resultsdf = pd.concat(results, axis=0)
        # else:
        #     resultsdf = results[0]
        resultsdf = pd.concat(results, axis=0)
        resultsdf.index = resultsdf["gene!isoform"]
        # resultsdf["ciwidth_i"] = resultsdf["cihigh_i"] - resultsdf["cilow_i"]
        resultsdf.to_csv(self.outfile, sep="\t")
        return resultsdf

    def readsummfile(self, sampindex):
        filename = self.sampleann.index[sampindex]
        cellname = self.sampleann.loc[filename][1]
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
            print df1.columns
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
        starts = map(int, x["mRNA_starts"].split(","))
        ends = map(int, x["mRNA_ends"].split(","))
        # todo add isf chromosome start stops to long format
        if len(psis) == 1:
            psis.append(1 - psis[0])
            # extrapolate confidence interval to psi of next isoform by maintaining width and position within range
            cilows.append(psis[1]-(psis[0]-cilows[0]))
            cihighs.append(psis[1]+(cihighs[0]-psis[0]))

        # calculate num reads informative and definitive
        ciwidths = list(np.array(cihighs) - np.array(cilows))
        lengthdiffs = list(np.array(ends) - np.array(starts))
        isflengths = [abs(i) for i in lengthdiffs]
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
        #
        tmpcounts = x["assigned_counts"].split(",")
        # print tmpcounts
        counts = {}
        for c in tmpcounts:
            csplit = map(int, c.split(":"))
            counts[csplit[0]] = csplit[1]
        # counts = [c.split(":")[1] for c in tmpcounts]
        # print len(counts), len(psis)
        y = []
        for i in range(len(rownames)):
            y.append(x.append(pd.Series([rownames[i], isfshortnames[i], psis[i], cilows[i], cihighs[i], ciwidths[i],
                                         numreadsinf, numreadsdef, counts.get(i, 0), isflengths[i]],
                                        index=["gene!isoform", "isfshortname", "psi_i", "cilow_i", "cihigh_i",
                                               "ciwidth_i", "numreadsinf", "numreadsdef", "assignedcounts_i",
                                               "isflength_i"])))

        # y = []
        # for i in range(len(rownames)):
        #     y.append(x.append(pd.Series([rownames[i], isfshortnames[i], psis[i], cilows[i], cihighs[i], ciwidths[i],
        #                                  numreadsinf, numreadsdef],
        #                                 index=["gene!isoform", "isfshortname", "psi_i", "cilow_i", "cihigh_i",
        #                                        "ciwidth_i", "numreadsinf", "numreadsdef"])))
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

# def _medianshift(x, samp1psis, samp2psis):
#     shift = samp1psis.loc[x.name].
# def _ttest(x, samp2kept2, maskmat1, maskmat2):
#     totest1 = x[maskmat1.loc[x.name]]
#     totest2 = samp2kept2.loc[x.name][maskmat2.loc[x.name]]
#     return stats.ttest_ind(totest1, totest2, index=["T", "pvalue"], name=x.name)
