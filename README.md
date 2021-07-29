# Multiple Rare-variants and Phenotypes Mixed Model

This directory contains code and documentation for the [Multiple Rare-variants and Phenotypes Mixed Model method]() by Venkataraman *et. al.*. The repository is maintained by Guhan Ram Venkataraman (GitHub: guhanrv).

## Script details and options

A full list of options can be obtained by running `python3 mrpmm.py -h`, and the output is replicated here for reference:

```{bash}
(base) user$ python3 mrpmm.py -h
usage: mrpmm.py [-h] --variants VARIANTS --phenotypes PHENOTYPES
                --metadata_path METADATA_PATH
                [--variant_filters {pcv,pav,ptv} [{pcv,pav,ptv} ...]]
                [--out_folder OUT_FOLDER] [--C CLUSTERS [CLUSTERS ...]]
                [--se_thresh SE_THRESH] [--maf_thresh MAF_THRESH]
                [--fout FOUT]

MRPMM takes in several variables that affect how it runs.

optional arguments:
  -h, --help            show this help message and exit
  --variants VARIANTS   path to file containing list of variants to include,
                                 one line per variant. Has a header of "V".
                        
                                 format of file:
                        
                                 V
                                 1:69081:G:C
                                 1:70001:G:A
                                 
  --phenotypes PHENOTYPES
                        path to tab-separated file containing list of: 
                                 summary statistic file paths,
                                 phenotypes, and
                                 whether or not to use the file in R_phen generation.
                               
                                 format:
                                 
                                 path        study     pheno     R_phen
                                 /path/to/file1   study1    pheno1     TRUE
                                 /path/to/file2   study2    pheno2     FALSE
                                 
  --metadata_path METADATA_PATH
                        path to tab-separated file containing:
                                 variants,
                                 gene symbols,
                                 consequences,
                                 and HGVSp annotations.
                               
                                 format:
                                 
                                 V       gene_symbol     most_severe_consequence HGVSp  
                                 1:69081:G:C     OR4F5   5_prime_UTR_variant     ""
                                
  --variant_filters {pcv,pav,ptv} [{pcv,pav,ptv} ...]
                        variant set(s) to consider. 
                                 options: proximal coding [pcv], 
                                          protein-altering [pav], 
                                          protein truncating [ptv],
                                          (default: ptv). can run multiple.
  --out_folder OUT_FOLDER
                        folder to which output(s) will be written (default: current folder).
                                 if folder does not exist, it will be created.
  --C CLUSTERS [CLUSTERS ...]
                        what number of clusters to use. must be valid ints. can input multiple
                                 (default: 1).
  --se_thresh SE_THRESH
                        SE threshold for variant inclusion
  --maf_thresh MAF_THRESH
                        MAF threshold for variant inclusion
  --fout FOUT           file prefix with which output(s) will be written (default: underscore-delimited
                                 phenotypes).
```

## Variant groupings

The following groups are how we assign priors (`sigma_m`) to variants: Protein-truncating variants (PTVs), Protein-altering variants (PAVs), proximal coding variants (PCVs), and intronic variants ("all" setting above).

`ptv = ['frameshift_variant', 'splice_acceptor_variant', 'splice_donor_variant', 'stop_gained', 'start_lost', 'stop_lost']`. By default, PTVs are assigned a spread (sigma) of 0.2.

`pav = ['protein_altering_variant', 'inframe_deletion', 'inframe_insertion', 'splice_region_variant', 'start_retained_variant', 'stop_retained_variant', 'missense_variant']`. By default, PAVs are assigned a spread (sigma) of 0.05.

`proximal_coding = ['synonymous_variant', '5_prime_UTR_variant', '3_prime_UTR_variant', 'coding_sequence_variant', 'incomplete_terminal_codon_variant', 'TF_binding_site_variant']`. By default, PCVs are assigned a spread (sigma) of 0.03.

`intron = ['regulatory_region_variant', 'intron_variant', 'intergenic_variant', 'downstream_gene_variant', 'mature_miRNA_variant', 'non_coding_transcript_exon_variant', 'upstream_gene_variant', 'NA', 'NMD_transcript_variant']`. By default, intronic variants are assigned a spread (sigma) of 0.02.

These groupings and their sigma values can be changed in the `filter_category` method within [`mrpmm.py`](https://github.com/rivas-lab/mrpmm/blob/master/mrpmm.py).

**IMPORTANT NOTE:** If `pav` is selected as the analysis type, then PAVs **and** PTVs are included in the analysis (cascading down). If `pcv` is selected, then PCVs, PAVs **and** PTVs are all included.

## Downloading variant metadata files

We provide variant metadata files for both UK Biobank [array](https://biobankengine.stanford.edu/static/ukb_cal-consequence_wb_maf_gene_ld_indep_mpc_pli.tsv.gz) and [exome](https://biobankengine.stanford.edu/static/ukb_exm_oqfe-consequence_wb_maf_gene_ld_indep_mpc_pli.tsv.gz) for direct download. We generated these files as described in [Venkataraman et. al.](https://www.biorxiv.org/content/10.1101/257162v7).

## Example use cases

We provide for download exome meta-analysis summary statistics for HDL cholesterol, triglycerides, and LDL cholesterol.

|       Trait     |     Summary Statistic    |
| --------------  | ------------ |
| HDL Cholesterol | https://biobankengine.stanford.edu/static/hdl_mrpmm.metal.tsv.gz |
| Triglycerides  | https://biobankengine.stanford.edu/static/tg_mrpmm.metal.tsv.gz |
| LDL Cholesterol | https://biobankengine.stanford.edu/static/ldl_mrpmm.metal.tsv.gz |

### Protein-truncating variant clustering analysis

Say we want to obtain variant clusterings for PTVs in genes signfiicantly associated to these phenotypes (HDL Cholesterol, Triglycerides, LDL Cholesterol) and want an analysis of how well these variants fit into up to 5 clusters of effects. We can find which genes are candidate genes by looking at our exome-wide gene-based results [here](https://biobankengine.stanford.edu/RIVAS_HG38/mrpgene/all) or the same for array [here](https://biobankengine.stanford.edu/RIVAS_HG19/mrpgene/all). Let's say we're interested in investigating the PTVs within the *PCSK9* and *APOB* genes. Using the appropriate metadata file ([array](https://biobankengine.stanford.edu/static/ukb_cal-consequence_wb_maf_gene_ld_indep_mpc_pli.tsv.gz), HG19, or [exome](https://biobankengine.stanford.edu/static/ukb_exm_oqfe-consequence_wb_maf_gene_ld_indep_mpc_pli.tsv.gz), HG38), we can extract those variants that we want to a file. The file (which will be input into the `--variants` parameter) should look like this:

```
V
1:55039741:C:G
1:55039742:A:G
1:55039749:C:G
1:55039750:C:T
1:55039753:C:T
1:55039755:GCAA:G
1:55039769:G:A
1:55039770:C:G
1:55039772:G:C
1:55039774:C:T
...
```

Then, we can create a "map file" (to be input into `--phenotypes`) that looks like the following:

```
path	study	pheno	R_phen
/path/to/hdl_mrpmm.glm.linear.gz	metal	hdl	TRUE
/path/to/tg_mrpmm.glm.linear.gz        metal   tg     TRUE
/path/to/ldl_mrpmm.glm.linear.gz        metal   ldl     TRUE
```

Then, from the command line, we run:

`python3 mrpmm.py --variants /path/to/variant_file --phenotypes /path/to/map_file --metadata_path /path/to/ukb_exm_oqfe-consequence_wb_maf_gene_ld_indep_mpc_pli.tsv.gz --C 1 2 3 4 5 --se_thresh 100 --variant_filters ptv`

Note the use of a non-default higher SE threshold for quantitative trait. The defaults should take care of the rest.

### Other optionalities

As mentioned in the `-h` command, MRPMM is flexible enough to perform clusterings for any valid integer value of clusters, use different standard error and minor allele frequency thresholds, and analyze PCVs, PAVs, or PTVs. We caution against using MRPMM for a large number of genes (containing a large number of variants) as this may result in none of the hypothesized number of clusters fitting well.
