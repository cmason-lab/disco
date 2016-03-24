import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as pdf
import seaborn as sns

sns.set_context("talk", font_scale=1.5)
sns.set_style("darkgrid")


def plotsig_violin(sigks, statsres, sciso1, sciso2, outfile, sampname1, sampname2, color1="r", color2="b"):
    # enstdf = pd.DataFrame.from_csv(enstfile, sep="\t", index_col=1)
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
        merged2 = merged.apply(_renamer, axis=1, args=(sigks.index, statsres.index))
        print merged2.shape
        sns.violinplot(x="newlabel", y="psi_i", hue="group", data=merged2, split=True,
                       inner="quart", cut=0, palette={sampname1: color1, sampname2: color2})
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


def _renamer(x, sigkeys, kskeys):
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


def plotviolin_1samp(disco, genestoplot, outfile):
    # todo this is probably not working, fix it
    # enstdf = pd.DataFrame.from_csv(enstfile, sep="\t", index_col=1)
    sns.set(style="dark", palette="muted", color_codes=True, font_scale=1.5)
    plotspdf = pdf.PdfPages(outfile)
    alldatadf = disco.alldatadf
    # genestoplot = sigks["Ensemble_ID"].unique()
    print "Plotting ", len(genestoplot), " genes"
    for gene in genestoplot:
        # print gene
        genedf = alldatadf[alldatadf["event_name"] == gene]
        print genedf
        genedf2 = genedf.apply(lambda x: x.append(pd.Series([x["isfshortname"].strip("isf-")], index=["newlabel"])),
                               axis=1)
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

# def plotsig_hist(sigks, statsres, disco1, disco2, outfile, sampname1, sampname2):
#     pass


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
