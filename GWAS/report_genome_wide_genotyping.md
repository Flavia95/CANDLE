# Preparing GWAS vcf file for use in  analysis 

## 1. Prepare candle file

*  Drop off 23 genotype samples for whom we don’t have the transcriptomics data

```
bcftools view -S ^samples_removed.txt -o CANDLE_GWAS_mother_623.vcf CANDLE_GWAS_mother_648_vcf_all.vcf
```

*  Minor Allele Frequency spectrum

```
plink2 --vcf CANDLE_GWAS_mother_623.vcf --freq
```
![MAF](https://github.com/Flavia95/CANDLE/assets/52487106/5225ba3a-4dd2-4bf0-862a-04e3050baa99)
The plot shows what we expected, but let's see what happened when we extracted the rows with ALT_FREQS equal to 0. The result is that we have a lot of missing data.
```
filtered_df <- x %>% filter(ALT_FREQS == 0 )

"KG_1_9876305_f   T   C         0   1238   9953718
13        1           KG_1_11737998_r   A   .         0   1240  11815411
14        1    indel_rs35607549_unw_r   -   .         0   1244  11858936"
```
---

# Analyses on Chr21:

* Working only on chromosome 21

```
bcftools view -r 21 CANDLE_GWAS_mother_623.vcf.gz > CANDLE_GWAS_mother_623.vcf.chr21.vcf
```
## 2. Statistics on vcf file

* How many variants types (SNPs, indels etc)? How many multiallelic sites?

```
bcftools stats CANDLE_GWAS_mother_623.vcf.chr21.vcf.gz
```

Samples           | Records |SNPs | no-ALTs | Multiallelic sites       
--------------| ----|---------|---------|------
623| 12,947 |  12,939 |8| 0 | 


## 4. Check distribution of missing data per site and per individual

### Missing data 

```
vcftools --gzvcf CANDLE_GWAS_mother_623.vcf.chr21.vcf.gz --missing-site --out missing_alleles_chr21
vcftools --gzvcf CANDLE_GWAS_mother_623.vcf.chr21.vcf.gz --missing-indv --out missing_alleles_chr21
```

Column 6 of output file is frequency of missing alleles 
Summarising no. SNPs with more than 5% missing data, more than 10%, and more than 20%

```
cat missing_alleles_chr21.lmiss | awk '$6 > 0.05' | wc -l
cat missing_alleles_chr21.lmiss | awk '$6 > 0.1' | wc -l
cat missing_alleles_chr21.lmiss | awk '$6 > 0.2' | wc -l
```
% missing data           |   result|    
--------------| ----|
more than 5% | 535 |
more than 10% | 267 | 
more than 20% | 55 |

Frequency plot of missing data:  
![frequencyplotmissingdata](https://github.com/Flavia95/CANDLE/assets/52487106/777f73cb-fe06-4c82-94a2-591b2e458dae)

Will remove sites with more than 10% missing data

Plots show the missing data is pretty evenly distributed across the chromosome  
![missingalongchr](https://github.com/Flavia95/CANDLE/assets/52487106/0beb0ac0-7f0c-44e6-841b-81b2ba3bcefe)

Plot show the missing data is pretty evenly distributed across the individuals
![missingindv](https://github.com/Flavia95/CANDLE/assets/52487106/81f62103-d442-4c2c-948d-8ef402c1f6bf)


## 5. HWE

```
plink2 --vcf CANDLE_GWAS_mother_623.vcf.chr21.vcf.gz --hardy
hwe_data %>% filter(V10<0.0001) %>% tally()
    n
1 645
hwe_data %>% filter(V10<0.05) %>% tally()
     n
1 2563
```
ld)  

Only 645 sites have P value minor of 0.0001, it means that they deviated for HW equlibrium. Only 2,563 sites have P value minor of 0.05.

## 6. Linkage Disequilibrium analysis

```
plink2 --vcf CANDLE_GWAS_mother_623.vcf.chr21.vcf.gz --double-id --make-bed --out LD_plink_chr21
plink2 --bfile LD_plink_chr21 --indep-pairwise 1500 150 0.7
plink2 --bfile LD_plink_chr21 --extract plink2.prune.in --make-bed --out ft_ld

```
--ld-window 9999 (all variant pairs with >= 9999 SNPs between them will be ignored)  
--ld-window-kb 500 (all variant pairs with > 500kb between them will be ignored)  
--ld-window-r2 0.4 (print all R2 values above this threshold)  

Plot shows how linkage disequilibrium breaks down over physical distance on a chromosome.Each dot represents the R2 value between a pair of SNPs at a certain physical distance apart.
There is a trend for R2 values to decrease as the physical distance between SNPs increases. LD decays fairly rapidly on this chromosome, with R2 dropping below 0.4 for SNP pairs more than ~100,000 bp apart.
![linkagedacaychr21](https://github.com/Flavia95/CANDLE/assets/52487106/8111f198-8e02-47a1-a8d1-c41db2c38b67)

## 7. Filter data before do the PCA. 

* Keep SNPs with no more than 10% missing alleles.
* Remove low maf variants less than 0.01.
* Remove HWE p-value less than 1e-5.
* Removes correlated pairs of SNPs so that the remaining SNPs are roughly independent.

```
plink2 --vcf CANDLE_GWAS_mother_623.vcf.chr21.vcf.gz --hwe 1e-5 keep-fewhet --make-bed --out filter_hwe
plink --bfile filter_hwe --maf 0.01 --indep-pairwise 50 5 0.2  --make-bed --out chr21-filter
plink2 --bfile chr21-filter --freq --out chr21-filter-freq
plink2 --bfile chr21-filter --read-freq chr21-filter-freq.afreq --pca 
```
The ouput file: leave SNPs with MAF at least 1%, with no pairs remaining with r2>0.2. Snps with no more than 10% of missing alles and keep HWE p-value more than 1e-5.
- 57 variants  (of tot) removed due to Hardy-Weinberg exact test
- 677 variants removed due to minor allele threshold and linkage disequilibrium

PCA on 12213 variants and 623 people that pass filters and QC.



