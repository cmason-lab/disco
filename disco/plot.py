import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as pdf
import seaborn as sns

sns.set_context("talk", font_scale=1.5)


def plotsig_violin(sigks, statsres, sciso1, sciso2, outfile, sampname1, sampname2, color1="r", color2="b"):
    # enstdf = pd.DataFrame.from_csv(enstfile, sep="\t", index_col=1)
    sns.set(style="dark", palette="muted", color_codes=True, font_scale=1.5)
    plotspdf = pdf.PdfPages(outfile)
    alldatadf1 = sciso1.alldatadf
    alldatadf2 = sciso2.alldatadf
    alldatadf1["group"] = pd.Series(np.repeat(sampname1, alldatadf1.shape[0]), index=alldatadf1.index)
    alldatadf2["group"] = pd.Series(np.repeat(sampname2, alldatadf2.shape[0]), index=alldatadf2.index)
    genestoplot = sigks.sort('median_shift')["Ensemble_ID"].unique()
    print "Plotting ", len(genestoplot), " genes"
    if len(genestoplot) > 100:
        genestoplot = genestoplot[0:100]
        print "plotting top 100 genes (by median shift) of ", len(genestoplot), "total significant genes"
    for gene in genestoplot:
        # print gene
        genedf1 = alldatadf1[alldatadf1["event_name"] == gene]
        genedf2 = alldatadf2[alldatadf2["event_name"] == gene]
        merged = pd.concat([genedf1, genedf2])
        merged2 = merged.apply(_renamer, axis=1, args=(sigks.index, statsres.index, statsres['Isoform_Function']))
        # print merged2.shape
        sns.violinplot(x="newlabel", y="psi_i", hue="group", data=merged2, split=True,
                       inner="quart", cut=0, palette={sampname1: color1, sampname2: color2})
        annrow = sigks[sigks["Ensemble_ID"] == gene].iloc[0]
        if annrow["Gene_Symbol"] is None:
            plt.title(gene)
            print gene
        else:
            plt.title(annrow["Gene_Symbol"]+"\n"+annrow["Gene_Name"]+"; "+annrow["Locus"])
            print annrow["Gene_Symbol"]
        plt.ylim(0,1)
        plt.ylabel("Distribution of Isoform Expression")
        plt.xlabel("Isoform ( - tested, * significant )")
        plt.legend(loc='upper left', bbox_to_anchor=(1.0, 0.9))
        plt.subplots_adjust(top=0.88, right=0.8, bottom=0.2)
        plotspdf.savefig()
        plt.close()
    plotspdf.close()
    print "Saved plots to", outfile
    sns.set()
    return


def _renamer(x, sigkeys, kskeys, isffunc):
    isfshort = x["isfshortname"]
    geneisf = x["gene!isoform"]
    if geneisf in sigkeys:
        # newlabel = isfshort+"*"
        func = str(isffunc[geneisf])
        func = "" if func is None or func == 'nan' else func
        newlabel = isfshort.strip("isf-")+"*\n"+func
    elif geneisf in kskeys:
        # newlabel = isfshort
        func = str(isffunc[geneisf])
        func = "" if func is None or func == 'nan' else func
        newlabel = isfshort.strip("isf-")+"-"+"\n"+func
    else:
        # newlabel = isfshort+"/"
        newlabel = isfshort.strip("isf-")
    return x.append(pd.Series([newlabel], index=["newlabel"]))


def plotviolin_1samp(disco, sigks, outfile):
    sns.set(style="dark", palette="muted", color_codes=True, font_scale=1.5)
    plotspdf = pdf.PdfPages(outfile)
    alldatadf = disco.alldatadf
    genestoplot = sigks["Ensemble_ID"].unique()
    print "Plotting ", len(genestoplot), " genes"
    for i in range(len(genestoplot)):
        gene = genestoplot[i]
        # print gene
        genedf = alldatadf[alldatadf["event_name"] == gene]
        # print genedf.shape
        genedf2 = genedf.apply(lambda x: x.append(pd.Series([x["isfshortname"].strip("isf-")], index=["newlabel"])),
                               axis=1)
        sns.violinplot(x="newlabel", y="psi_i", data=genedf2, inner="quartile", cut=0)
        plt.ylim(0,1)
        annrow = sigks[sigks["Ensemble_ID"] == gene].iloc[0]
        if annrow["Gene_Symbol"] is None or str(annrow["Gene_Symbol"]) == "nan":
            plt.title(gene)
            print gene
        else:
            plt.title(annrow["Gene_Symbol"]+"\n"+annrow["Gene_Name"]+"; "+annrow["Locus"])
            print annrow["Gene_Symbol"]
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


def plotbulk(sigks, statsres, disco1, disco2, sampname1, sampname2, color1, color2, outfile):
    sns.set(style="dark", palette="muted", color_codes=True, font_scale=1.5)
    plotspdf = pdf.PdfPages(outfile)
    alldatadf1 = disco1.alldatadf
    alldatadf2 = disco2.alldatadf
    alldatadf1["group"] = pd.Series(np.repeat(sampname1, alldatadf1.shape[0]), index=alldatadf1.index)
    alldatadf2["group"] = pd.Series(np.repeat(sampname2, alldatadf2.shape[0]), index=alldatadf2.index)
    genestoplot = sigks["Ensemble_ID"].unique()
    for gene in genestoplot:
        genedf1 = alldatadf1[alldatadf1["event_name"] == gene]
        genedf1.head()
        genedf2 = alldatadf2[alldatadf2["event_name"] == gene]
        genedf2.head()
        merged = pd.concat([genedf1, genedf2])
        merged2 = merged.apply(_renamer, axis=1, args=(sigks.index, statsres.index, statsres['Isoform_Function']))
        sns.barplot(x="newlabel", y="psi_i", hue="group", data=merged2,
                    palette={sampname1: color1, sampname2: color2})
        plt.ylim(0, 1)
        annrow = sigks[sigks["Ensemble_ID"] == gene].iloc[0]
        if annrow["Gene_Symbol"] is None:
            plt.title(gene)
            print gene
        else:
            plt.title(annrow["Gene_Symbol"]+"\n"+annrow["Gene_Name"]+"; "+annrow["Locus"])
            print annrow["Gene_Symbol"]
        plt.ylabel("Distribution of Isoform Expression")
        plt.xlabel("Isoform ( - tested in sc, * significant in sc )")
        plt.legend(loc='upper left', bbox_to_anchor=(1.0, 0.9))
        plt.subplots_adjust(top=0.88, right=0.8, bottom=0.2)
        plotspdf.savefig()
        plt.close()
    plotspdf.close()
    print "Saved plots to", outfile
    sns.set()
    return
#
#
#
# #incomplete:


def plotsinglebulk(sigks, statsres, disco1, sampname1, color1, outfile):
    sns.set(style="dark", palette="muted", color_codes=True, font_scale=1.5)
    plotspdf = pdf.PdfPages(outfile)
    alldatadf1 = disco1.alldatadf
    #alldatadf2 = disco2.alldatadf
    alldatadf1["group"] = pd.Series(np.repeat(sampname1, alldatadf1.shape[0]), index=alldatadf1.index)
    #alldatadf2["group"] = pd.Series(np.repeat(sampname2, alldatadf2.shape[0]), index=alldatadf2.index)
    genestoplot = sigks["Ensemble_ID"].unique()
    for gene in genestoplot:
        genedf1 = alldatadf1[alldatadf1["event_name"] == gene]
        genedf1.head()
        #genedf2 = alldatadf2[alldatadf2["event_name"] == gene]
        #genedf2.head()
        #merged = pd.concat([genedf1, genedf2])
        merged = genedf1
        merged2 = merged.apply(_renamer, axis=1, args=(sigks.index, statsres.index, statsres['Isoform_Function']))
        sns.barplot(x="newlabel", y="psi_i", hue="group", data=merged2,
                    palette={sampname1: color1})
        plt.ylim(0, 1)
        annrow = sigks[sigks["Ensemble_ID"] == gene].iloc[0]
        if annrow["Gene_Symbol"] is None:
            plt.title(gene)
            print gene
        else:
            plt.title(annrow["Gene_Symbol"]+"\n"+annrow["Gene_Name"]+"; "+annrow["Locus"])
            print annrow["Gene_Symbol"]
        plt.ylabel("Distribution of Isoform Expression")
        plt.xlabel("Isoform ( - tested in sc, * significant in sc )")
        plt.legend(loc='upper left', bbox_to_anchor=(1.0, 0.9))
        plt.subplots_adjust(top=0.88, right=0.7, bottom=0.2)
        plotspdf.savefig()
        plt.close()
    plotspdf.close()
    print "Saved plots to", outfile
    sns.set()
    return


def plot3bulk(sigks, statsres, disco1, disco2, disco3, sampname1, sampname2, sampname3,
              color1, color2, color3, outfile):
    sns.set(style="dark", palette="muted", color_codes=True, font_scale=1.5)
    plotspdf = pdf.PdfPages(outfile)
    alldatadf1 = disco1.alldatadf
    alldatadf2 = disco2.alldatadf
    alldatadf3 = disco3.alldatadf
    alldatadf1["group"] = pd.Series(np.repeat(sampname1, alldatadf1.shape[0]), index=alldatadf1.index)
    alldatadf2["group"] = pd.Series(np.repeat(sampname2, alldatadf2.shape[0]), index=alldatadf2.index)
    alldatadf3["group"] = pd.Series(np.repeat(sampname3, alldatadf3.shape[0]), index=alldatadf3.index)
    genestoplot = sigks["Ensemble_ID"].unique()
    for gene in genestoplot:
        genedf1 = alldatadf1[alldatadf1["event_name"] == gene]
        genedf1.head()
        genedf2 = alldatadf2[alldatadf2["event_name"] == gene]
        genedf2.head()
        genedf3 = alldatadf3[alldatadf3["event_name"] == gene]
        genedf3.head()
        merged = pd.concat([genedf1, genedf2, genedf3])
        merged2 = merged.apply(_renamer, axis=1, args=(sigks.index, statsres.index, statsres['Isoform_Function']))
        sns.barplot(x="newlabel", y="psi_i", hue="group", data=merged2,
                    palette={sampname1: color1, sampname2: color2, sampname3: color3})
        plt.ylim(0, 1)
        annrow = sigks[sigks["Ensemble_ID"] == gene].iloc[0]
        if annrow["Gene_Symbol"] is None:
            plt.title(gene)
            print gene
        else:
            plt.title(annrow["Gene_Symbol"]+"\n"+annrow["Gene_Name"]+"; "+annrow["Locus"])
            print annrow["Gene_Symbol"]
        plt.ylabel("Distribution of Isoform Expression")
        plt.xlabel("Isoform ( - tested in sc, * significant in sc )")
        plt.legend(loc='upper left', bbox_to_anchor=(1.0, 0.9))
        plt.subplots_adjust(top=0.88, right=0.7, bottom=0.2)
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
