# Identifying patterns of introgression in complex natural monkeyflower hybrid zones

This code describes the variant calling pipeline (using GATK Best Practices) and downstream analysis of WGS data.

## GATK Variant Calling Pipeline 

This takes raw fastq files, aligns them to the TOLv5 reference genome, and calls variants. Output is hybrid1_all.vcf.
Steps that say complex are for samples sequenced across multiple lanes, simple means across one lane. 

Step1_readQC
Step2_preprocessing_complex_pt1
Step2_preprocessing_complex_pt2
Step2_preprocessing_simple
Step3_complex_pt1
Step3_complex_pt2
Step3.2_temp
  - (this is to run Step 3 chromosome by chromosome, which I had to do here because its a big dataset)
  - sample_map.py is required to make the sample_map for the GATK haplotypeCaller tool in all Steps 3

## Downstream Analysis of WGS Data

The following steps take the output of the GATK pipeline (hybrid1_all.vcf).

### Filtering

Visualize metrics using ```FS_visualization``` input of this is hybrid1_all.vcf and output is hybrid1_all.hf.vcf

### Plink PCA

Input is hybrid1_all.hf.vcf and output is Plink bfiles as well as PCA files (eigenvec and eigenval). 
Run using ```vcf2plink```.
There are also some files for making a pca with smartpca (EIGENSOFT). I didn't end up using this.

### Fast Stucture

This takes the output of ```vcf2plink``` in the prior section (bed/bim/fam). 
Set up the environment using ```fastStructure_setup``` and then run using ```1_faststucture```.

### WinPCA

This is a windowed PCA approach (https://academic.oup.com/bioinformatics/article/41/10/btaf529/8261369) 
for identifying possible structural across the genome. 
Also referenced this tutorial (https://github.com/clairemerot/Tutorial_SV/blob/main/01_pca_haploblocks/README.md)
Dependencies are: numpy pandas numba scikit-allel plotly.
Requires chromosome lengths:
```
#from TOLv5 annotation
chr_lengths <- c("1" = 11879706, "2" = 17967367, "3" = 16906840, "4"=20228821, "5"=19037101, "6"=18017287,
                 "7"=16736285, "8"=24504531, "9"=14486675, "10"=19650840, "11"=16175311, "12"=19059609,
                 "13"= 18954856, "14"=26762607)
```

Set up:
```
module load python
conda create -n winpca
git clone https://github.com/MoritzBlumer/winpca.git  # clone github repository
chmod +x winpca/winpca                                # make excutable
mv winpca /home/dtataru/.conda/envs/winpca/bin/

```
To run: run_winpca.sh
```
#!/bin/bash
#SBATCH --job-name=winpca
#SBATCH --output=/project/dtataru/hybrids/hybrids/logs/winpca.out
#SBATCH --error=/project/dtataru/hybrids/hybrids/logs/winpca.err
#SBATCH --time=0-24:00:00
#SBATCH -N 1
#SBATCH --cpus-per-task=20
#SBATCH -A loni_ferrislac
#SBATCH --partition=single

### LOAD MODULES ###

module load bcftools/1.18
module load python/3.11.5-anaconda  
eval "$(conda shell.bash hook)"
conda activate /home/dtataru/.conda/envs/winpca

cd /project/dtataru/hybrids/winpca

#filter for biallelic SNPs - otherwise we run into an error:
bcftools view -m2 -M2 -v snps /project/dtataru/hybrids/4_GATKvarcall/3_Genotyped_GVCFs/hybrid1_all.hf.vcf > /project/dtataru/hybrids/4_GATKvarcall/3_Genotyped_GVCFs/hybrid1_all.hf.biallelic.vcf

#Then let's run winpca code

#make the PCAs
#-w for window size -i for increment size --np to remove filters creating an error -v GT to precise the type of data.
#then there are three argument "$PREFIX" "$VCF" "$REGION"
winpca pca -w 10000 -i 10000 --np -v GT winpca_out_chr1 /project/dtataru/hybrids/4_GATKvarcall/3_Genotyped_GVCFs/hybrid1_all.hf.biallelic.vcf Chr1:1-11879706

#polarize them
winpca polarize winpca_out_chr1

#plot PCs
winpca chromplot winpca_out_chr1 Chr1:1-11879706

#get out of the env
conda deactivate
```


## Consensus Genome

This was written prior to the release of the M. laciniatus reference genome, it is to create a high quality pseudo-reference from WLF47.
This pseudoreference is used in both downstream ancestry analyses- spIDder and Ancestry-HMM.
