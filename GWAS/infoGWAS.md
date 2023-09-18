# GWAS filter

The Genotyping data I shared with you had done partial of sample QC (exclude heterozygosity rate part), and did not perform QC on SNP. 
Investigators can perform QC with their own criteria. Criteria that they used obtained by the section of “Genome-wide Genotyping” of CANDLE paper at https://pubmed.ncbi.nlm.nih.gov/36055779/.

## Genome-wide genotyping
Among the 450 CANDLE mothers, 351 had available genotype data obtained from the Affymetrix Axiom® Genome-Wide AFR 1 Array Set (Affymetrix) after the quality control (QC)
procedures were performed, including removing individuals with sex discordance, high genotype missing rates (>3.0%), cryptic relatedness (identity by descent >0.1875), or a very high or low
heterozygosity rate (> mean+3SD or < mean–3SD). We also excluded single nucleotide polymorphisms (SNPs) with minor allele frequencies <0.01, genotype call rates <95%, or Hardy–Weinberg equilibrium testing P values <10−6 . A total of 802,649 autosomal SNPs passed QC. Principal component analysis was performed using PLINK, and the first 10 principal components were calculated for subsequent analysis in an attempt to correct for the potential genetic structure of the study subjects.

* Publications about CANDLE: https://candlestudy.uthsc.edu/publications-presentations/
