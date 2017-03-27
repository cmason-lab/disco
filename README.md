# DISCO
Distributions of Isoforms in Single Cell Omics

## Installation

```
python setup.py install
```

## Usage


**To run with default settings:**

```
disco path/to/samplefilelist.txt group1 group2
```

**To see other options, including filtering (default is no filtering):**


```shell
disco -h

usage: disco [-h] [-v] [--outdir] [--pkldir] [--group1color] [--group2color]
         [--group1file] [--group2file] [--geneannotationfile]
         [--transcriptannotationfile] [--maxciwidth] [--mininfreads]
         [--mindefreads] [--minavgpsi] [--minnumcells] [--minavgshift]
         [--stattest] [--alpha] [--multitestmethod]
         SampleAnnotationFile Group1 Group2


DISCO: Distributions of Isoforms in Single Cell Omics


positional arguments:


    SampleAnnotationFile  
                      filename of tab separated text, no header, with
                      columns: <path to miso summary file> <sample name>
                      <group name>
    Group1                
                      must match a group name in sample annotation file

    Group2                
                      must match a group name in sample annotation file


optional arguments:

    -h, --help            show this help message and exit

    -v, --version         show program's version number and exit

    --outdir              Output directory (default: ./disco_output/)

    --pkldir              Directory to store intermediate data processing files
                          (default: ./pkldir)

    --group1color         Color in plots for group 1; can be {y, m, c, r, g, b,
                          w, k} or html code (default: r)

    --group2color         Color in plots for group 2; can be {y, m, c, r, g, b,
                          w, k} or html code (default: b)

    --group1file          output file for sample group 1. If not specified, will
                          save to <outdir>/<group1name>_alldatadf.txt (default:
                          None)

    --group2file          output file for sample group 2. If not specified, will
                          save to <outdir>/<group2name>_alldatadf.txt (default:
                          None)

    --geneannotationfile  
                          Mapping of Ensembl gene IDs to HGNC symbol and gene
                          descriptions (default: None)

    --transcriptannotationfile 
                          Mapping of Ensembl transcript IDs to isoform function
                          (ex. protein coding, NMD, etc) (default: None)

    --maxciwidth          Maximum width of confidence interval of PSI estimate
                          (default: 1.0)

    --mininfreads         Minimum number of informative reads to include PSI
                          estimate (default: 0)

    --mindefreads         Minimum number of definitive reads to include PSI
                          estimate (default: 0)
                    
    --minavgpsi           Do not run statistical tests for isoforms with average
                          PSI in both groups less than minavgpsi (default: 0.0)

    --minnumcells         Do not run statistical test for isoform if less than
                          minnumcells have information (default: 0)

    --minavgshift         Do not run statistical test for isoform if shift in
                          mean PSI between the two groups is less than
                          minavgshift (default: 0)

    --stattest            Which test to run? options: {KS, T} (default: KS)

    --alpha               Adjustded p-value threshold for statistical
                          significance (default: 0.05)

    --multitestmethod     Method for multiple testing correction. Available methods are: 
                          `bonferroni` : one-step correction
                          `sidak` : one-step correction
                          `holm-sidak` : step down method using Sidak adjustments
                          `holm` : step-down method using Bonferroni adjustments
                          `simes-hochberg` : step-up method  (independent)
                          `hommel` : closed method based on Simes tests (non-negative)
                          `fdr_bh` : Benjamini/Hochberg  (non-negative)
                          `fdr_by` : Benjamini/Yekutieli (negative)
                          `fdr_tsbh` : two stage fdr correction (non-negative)
                          `fdr_tsbky` : two stage fdr correction (non-negative)
                          `none` : no multiple testing correction
                          (default: fdr_bh)
```

## Example

```
cd example/
./run.sh >> run.log 2>&1
```

