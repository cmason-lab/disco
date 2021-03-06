import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.stats.multitest as smm


# Functions for filtering and comparison
def stat_test(disco1, disco2, outfile, maxciw=1, mininfreads=0, mindefreads=0, minavgpsi=0, minnumcells=0,
              minavgshift=0, multitestmethod="fdr_bh", alpha=0.05, geneannfile=None, transcriptannfile=None,
              testtype="KS"):
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
    samp1psis.to_csv(disco1.outfile.rstrip(".txt")+"_filtpsimat.txt", sep="\t")
    print
    samp2datadf = disco2.alldatadf
    print "Sample 2:", samp2datadf.shape
    samp2psis = filterdf(samp2datadf, maxciw, mininfreads, mindefreads, minnumcells)
    print "Filtered Sample 2:", samp2psis.shape
    samp2psis.to_csv(disco2.outfile.rstrip(".txt")+"_filtpsimat.txt", sep="\t")
    print
    # intersectkeys = set(samp1psis.index) & set(samp2psis.index)
    # print "Number of shared isoforms:", len(intersectkeys)
    # samp1kept = samp1psis.loc[intersectkeys, ]
    # samp2kept = samp2psis.loc[intersectkeys, ]

    unionkeys = set(samp1psis.index) | set(samp2psis.index)
    print "Number of total isoforms:", len(unionkeys)
    samp1kept = samp1psis.loc[unionkeys, ]
    samp2kept = samp2psis.loc[unionkeys, ]

    joineddf = samp1kept.join(samp2kept)
    newkeys = filterminavgpsi(joineddf, minavgpsi)
    print "After filtering for min avg. psi:", len(newkeys)

    samp1kept2 = samp1kept.loc[newkeys, ]
    samp2kept2 = samp2kept.loc[newkeys, ]
    meanshifts = samp2kept2.mean(axis=1) - samp1kept2.mean(axis=1)
    newkeys2 = meanshifts[meanshifts.abs() > minavgshift].index
    print "After filtering for min avg. shift in psi:", len(newkeys2)
    # medianshifts = samp1kept2.median(axis=1) - samp2kept2.median(axis=1)
    # newkeys2 = medianshifts[medianshifts > minmedianshift].index
    # print "After filtering for min median shift in psi:", len(newkeys2)

    samp1kept2 = samp1kept2.loc[newkeys2, ]
    samp2kept2 = samp2kept2.loc[newkeys2, ]
    meanshifts = meanshifts.loc[newkeys2]
    # todo test only 1 isoform when 2 exist (might already be doing this..)
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
    # medianshifts = samp1kept2.median(axis=1) - samp2kept2.median(axis=1)
    statsres2 = statsres.apply(lambda x: x.append(pd.Series([meanshifts.loc[x.name]],
                                                            index=["mean_shift"], name=x.name)), axis=1)
    statsres = statsres2
    if multitestmethod == "none":
        adjp = pd.Series(statsres['pvalue'], index=statsres.index, name="adjusted_pvalue")
        sig = pd.Series(adjp < alpha, index=statsres.index, name="significant")
    else:
        multitestres = smm.multipletests(statsres['pvalue'], method=multitestmethod, alpha=alpha)
        adjp = pd.Series(multitestres[1], index=statsres.index, name="adjusted_pvalue")
        sig = pd.Series(multitestres[0], index=statsres.index, name="significant")
    statsres2 = pd.concat([statsres, adjp], axis=1)
    statsres = statsres2
    statsres2 = pd.concat([statsres, sig], axis=1)
    statsres = statsres2
    print "Significant at mulitple testing corrected pvalue", alpha, ":", len(statsres[statsres["significant"]])
    if geneannfile is not None:
        # todo clean up hardcoding
        annotationdf = pd.DataFrame.from_csv(geneannfile, sep="\t", index_col=9)
    else:
        annotationdf = None
    statsres_ann = statsres.apply(_getannotation, axis=1, args=(annotationdf, "dummy"))
    statsres = statsres_ann
    if transcriptannfile is not None:
        # todo clean up hardcoding
        enstannotationdf = pd.DataFrame.from_csv(transcriptannfile, sep="\t", index_col=None)
        enstdf2 = enstannotationdf.drop_duplicates(['Ensembl Transcript ID'])
        enstdf2.index = enstdf2['Ensembl Transcript ID']
    else:
        enstdf2 = None
    statsres_ann = statsres.apply(_getenstannotation, axis=1, args=(enstdf2, "dummy"))
    statsres = statsres_ann
    # todo add isf shortname to this DF
    # statsres_ann["isfshortname"] = samp1datadf.loc[statsres_ann.index]["isfshortname"]

    statsres.to_csv(outfile, sep="\t")
    print
    # print "Final result:", statsres.shape
    return statsres


def _getannotation(x, annotationdf, dummy):
    ensid = x.name.split("!")[0]
    if annotationdf is not None and ensid in annotationdf.index:
        annrow = annotationdf.loc[ensid]
        y = pd.Series([ensid, annrow["Approved Symbol"], annrow["Approved Name"], annrow["Chromosome"]],
                      index=["Ensemble_ID", "Gene_Symbol", "Gene_Name", "Locus"],
                      name=x.name)
    else:
        y = pd.Series([ensid, None, None, None],
                      index=["Ensemble_ID", "Gene_Symbol", "Gene_Name", "Locus"],
                      name=x.name)
    return x.append(y)


def _getenstannotation(x, enstdf, dummy):
    # example: ENSG00000162244!ENST00000479017.exon3_ENST00000479017.exon2_ENS
    ensid = x.name.split("!")[1].split(".")[0]
    if enstdf is not None and ensid in enstdf.index:
        annrow = enstdf.loc[ensid]
        y = pd.Series([ensid, annrow["Transcript type"], annrow["Transcript length (including UTRs and CDS)"]],
                      index=["Transcript_ID", "Isoform_Function", "Isoform_Length"],
                      name=x.name)
    else:
        y = pd.Series([ensid, None, None],
                      index=["Transcript_ID", "Isoform_Function", "Isoform_Length"],
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


def getsig(statsres, outfile=None):
    sigres = statsres[statsres["significant"]]
    if outfile is not None:
        sigres.to_csv(outfile, sep="\t", index_label="gene!isoform")
    return sigres

# def getsig(statsres, pvalue=0.05, bonferroni=True, outfile=None): #minmedianshift=0,
#     thresh = pvalue/statsres.shape[0] if bonferroni else pvalue
#     # significant
#     sigres = statsres[statsres["pvalue"] < thresh]
#     print "Significant at pvalue", pvalue, "and bonferroni", bonferroni, ":", sigres.shape
#     toreturn = sigres
#     # if sigres.shape[0] > 0:
#     #     sigres2 = sigres[sigres["median_shift"].abs() > minmedianshift]
#     #     print "Number of sites with median shift >", minmedianshift, ":", sigres2.shape
#     #     toreturn = sigres2
#     if outfile is not None:
#         toreturn.to_csv(outfile, sep="\t", index_label="gene!isoform")
#     return toreturn

