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

## Local PCA to identify structural variants

Taking the  subset_230_biallelic_only_alt_imputed.vcf.gz used as imput for GEMMA and calculating genotype likelihoods from it using workflow from Gompert et al. 20205, Science. 

First, run vcf2gl.pl to convert the vcf to genotype likelihood format, then run estpEM and gl2genest. These scripts are in /project/dtataru/SY2/bams, so I can just call that directory and run a separate script in /project/dtataru/hybrids/lostruct called run_glestimation.sh

```
#!/bin/bash
#SBATCH --job-name=glest
#SBATCH --output=/project/dtataru/hybrids/logs/glest.out
#SBATCH --error=/project/dtataru/hybrids/logs/glest.err
#SBATCH --time=0-24:00:00
#SBATCH -N 1
#SBATCH --cpus-per-task=20
#SBATCH -A loni_ferrislac
#SBATCH --partition=single

SCRIPT_DIR="/project/dtataru/SY2/bams"
VCF="/project/dtataru/hybrids/GEMMA/subset_230_biallelic_only_alt_imputed.vcf"
OUTPUT_DIR="/project/dtataru/hybrids/lostruct"

cd ${OUTPUT_DIR}

#convert vcf to genotype likelihood
perl ${SCRIPT_DIR}/vcf2gl.pl maf ${VCF}

#name of output files and first line of .gl to have #ind #loci instead of 0 0

#run estpEM
#${SCRIPT_DIR}/estpEM -i subset_230_biallelic_only_alt_imputed.gl -o subset_230_biallelic_only_alt_imputed_estpEM.txt -e 0.001 -m 50 -h 1

#run gl2genest
## posterior mode
#${SCRIPT_DIR}/gl2genestMax subset_230_biallelic_only_alt_imputed_estpEM.txt  subset_230_biallelic_only_alt_imputed.gl
## posterior mean
#${SCRIPT_DIR}/gl2genest subset_230_biallelic_only_alt_imputed_estpEM.txt subset_230_biallelic_only_alt_imputed.gl
```
okay there's acutally no GL data for that vcf, just GT, checked using:
```
grep -v "^#" subset_230_biallelic_only_alt_imputed.vcf | head -2 | cut -f9
GT
grep -v "^#" hybrid1_SH.hf.vcf | head -2 | cut -f9
GT:AD:DP:GQ:PGT:PID:PL:PS
GT:AD:DP:GQ:PGT:PID:PL:PS

```
One option I have is to convert my GT to hard calls (0/1/2) which should be fine with my high coverage data like so:

```
#!/bin/bash
#SBATCH --job-name=run_lostruct
#SBATCH --output=/project/dtataru/hybrids/logs/lostruct.out
#SBATCH --error=/project/dtataru/hybrids/logs/lostruct.err
#SBATCH --time=0-24:00:00
#SBATCH -N 1
#SBATCH --cpus-per-task=20
#SBATCH -A loni_ferrislac
#SBATCH --partition=single

INPUT_DIR="/project/dtataru/hybrids/GEMMA"
OUTPUT_DIR="/project/dtataru/hybrids/lostruct"
PREFIX="subset_230_biallelic_only_alt_imputed"
VCF="subset_230_biallelic_only_alt_imputed.vcf"

#load modules
module load python/3.11.5-anaconda
module load r

#activate plink env
#eval "$(conda shell.bash hook)"
#conda activate /home/dtataru/.conda/envs/plink
#export LD_LIBRARY_PATH=${CONDA_PREFIX}/lib:$LD_LIBRARY_PATH
cd ${OUTPUT_DIR}
#with full vcf
#plink --vcf ${INPUT_DIR}/${VCF} --recode A --out geno

#echo "plink pruning"
# Step 1: Calculate LD and identify SNPs to keep
#plink --bfile ${INPUT_DIR}/${PREFIX} \
#      --indep-pairwise 200 20 0.3 \
#      --out ld_pruned

# Step 2: Extract the pruned SNP set
#plink --bfile ${INPUT_DIR}/${PREFIX} \
#      --extract ld_pruned.prune.in \
#      --make-bed \
#      --out "${PREFIX}_pruned"

# Step 3: Recode to .raw for lostruct
#plink --bfile "${PREFIX}_pruned" \
#      --recode A \
#      --out geno_pruned

#7762060 sites remaining after pruning

#echo "extract positions from .bim"
#awk '{print $1, $4}' "${PREFIX}_pruned.bim" > positions.txt

#echo "transpose output to genotype matrix"
#python3 transpose_geno.py geno_pruned.raw genotype_matrix.txt

echo "run lostruct"
#Usage: Rscript ${SCRIPTDIR}/localpca_manyMDSaxes_v2.R <input_file> <output_prefix> <window_size_snps> <n_axes>
Rscript localpca_manyMDSaxes_withfiltering.R genotype_matrix.txt hybrids1_subset_227_biallelic_only_alt_imputed_1k 1000 10

echo "done"

```

## Run Entropy to compare to year 2 & 3

salloc --time=6:00:00 --ntasks=12 --nodes=1 --account=loni_ferrislac --partition=single

want to be able to compare to year 2 & 3, so need to subset the WGS data

```
#!/bin/bash
#SBATCH --job-name=bcftoolsisec
#SBATCH --output=/project/dtataru/hybrids/logs/bcftoolsisec.out
#SBATCH --error=/project/dtataru/hybrids/logs/bcftoolsisec.err
#SBATCH --time=0-72:00:00
#SBATCH -N 1
#SBATCH --cpus-per-task=20
#SBATCH -A loni_ferrislac
#SBATCH --partition=single

module load bcftools/1.18
module load htslib/1.23

cd /project/dtataru/hybrids/4_GATKvarcall/3_Genotyped_GVCFs

#rename ddRad datasets to match WGS
bcftools view -h morefilter_1x_hybrids2.maxdepth6000.vcf.gz | grep "^##contig" | sed -E 's/.*ID=([^,]+).*/\1/' | \
  awk '{old=$0; new=$0; gsub("-","_",new); print old"\t"new}' > rename_chrs.txt
bcftools annotate --rename-chrs rename_chrs.txt morefilter_1x_hybrids2.maxdepth6000.vcf.gz -Oz -o morefilter_1x_hybrids2.maxdepth6000_renamed.vcf.gz
tabix morefilter_1x_hybrids2.maxdepth6000_renamed.vcf.gz

bcftools view -h morefilter_1x_hybrids3.maxdepth6000.allsamples.vcf.gz | grep "^##contig" | sed -E 's/.*ID=([^,]+).*/\1/' | \
  awk '{old=$0; new=$0; gsub("-","_",new); print old"\t"new}' > rename_chrs.txt
bcftools annotate --rename-chrs rename_chrs.txt morefilter_1x_hybrids3.maxdepth6000.allsamples.vcf.gz -Oz -o morefilter_1x_hybrids3.maxdepth6000.allsamples_renamed.vcf.gz
tabix morefilter_1x_hybrids3.maxdepth6000.allsamples_renamed.vcf.gz

#SY2
bcftools isec -n=2 -w1 hybrid1_all.hf_renamed.vcf.gz /project/dtataru/SY2/bams/morefilter_1x_hybrids2.maxdepth6000_renamed.vcf.gz -Oz -o hybrid1_all.hf.subsetSY2.vcf.gz
#bcftools view -H hybrid1_all.hf.subsetSY2.vcf.gz | wc -l
#ended up with 62040/90625 variants

#SY3
bcftools isec -n=2 -w1 hybrid1_all.hf_renamed.vcf.gz /project/dtataru/SY3/bams/morefilter_1x_hybrids3.maxdepth6000.allsamples_renamed.vcf.gz -Oz -o hybrid1_all.hf.subsetSY3.vcf.gz
bcftools view -H hybrid1_all.hf.subsetSY3.vcf.gz | wc -l
#ended up with 1811 variants

#index new vcfs
tabix hybrid1_all.hf.subsetSY2.vcf.gz
tabix hybrid1_all.hf.subsetSY3.vcf.gz

#compute intersection (0002.vcf or 0003.vcf)
bcftools isec -p isec_out hybrid1_all.hf.subsetSY2.vcf.gz hybrid1_all.hf.subsetSY3.vcf.gz

#compute union
bcftools concat isec_out/0000.vcf isec_out/0001.vcf isec_out/0002.vcf -Oz -o union.vcf.gz
bcftools sort union.vcf.gz -Oz -o union_sorted.vcf.gz

#make position lists of each
bcftools query -f '%CHROM\t%POS\n' union_sorted.vcf.gz > union_pos.txt #62802
bcftools query -f '%CHROM\t%POS\n' isec_out/0002.vcf > intersection_pos.txt #1049 sites

mv union_sorted.vcf.gz WGS1_ddRADunion_sorted.vcf.gz

cd /project/dtataru/hybrids/4_GATKvarcall/3_Genotyped_GVCFs/isec_out/

#make genotypelikelihood file
perl /project/dtataru/SY2/bams/vcf2gl.pl 0.0 WGS1_ddRADunion_sorted.vcf
#Number of loci: 62802; number of individuals 297

#make treatComb.txt and KeepInds.txt files to split pops looks like list with sample and pop
yes 1 | head -n 297 > KeepInds.txt
#split pops
perl /project/dtataru/SY2/bams/splitPops.pl WGS1_ddRADunion_sorted.gl

#runestpEM
/project/dtataru/SY2/bams/estpEM -i GB_WGS1_ddRADunion_sorted.gl -o GB_WGS1_ddRADunion_estpEM.txt -e 0.001 -m 50 -h 1 #128
/project/dtataru/SY2/bams/estpEM -i SH_WGS1_ddRADunion_sorted.gl -o SH_WGS1_ddRADunion_estpEM.txt -e 0.001 -m 50 -h 1 #66
/project/dtataru/SY2/bams/estpEM -i HH_WGS1_ddRADunion_sorted.gl -o HH_WGS1_ddRADunion_estpEM.txt -e 0.001 -m 50 -h 1 #103

# split faststructure file by pop and make sure it matches the order
perl splitfaststructure.pl GB_WGS1_ddRADunion_sorted.gl
perl splitfaststructure.pl SH_WGS1_ddRADunion_sorted.gl
perl splitfaststructure.pl HH_WGS1_ddRADunion_sorted.gl
```

I need to check why there were so few shared sites, especially for sample year 1&3 I am guessing that this issue may have to do with multiallelic sites, so I will try to maybe just run position matching:

Aha! I figured out the issue. I had been playing around with what it would have looked like to keep all of the samples (n=221) for SY3 and call variants. this resulted in 2771 variants, too few. So i want to use the filtered dataset with just 159 samples(that is what I used for entropy anyways).
```
salloc --time=6:00:00 --ntasks=12 --nodes=1 --account=loni_ferrislac --partition=single

module load bcftools

#renanme the correct sy3 file
bcftools annotate --rename-chrs rename_chrs.txt morefilter_1x_hybrids3.maxdepth6000.vcf -Oz -o morefilter_1x_hybrids3.maxdepth6000_renamed.vcf.gz
tabix morefilter_1x_hybrids3.maxdepth6000_renamed.vcf.gz

# extract positions from each
bcftools query -f '%CHROM\t%POS\n' hybrid1_all.hf.vcf.gz | sort -u > wgs1_pos.txt
bcftools query -f '%CHROM\t%POS\n' /project/dtataru/SY2/bams/morefilter_1x_hybrids2.maxdepth6000_renamed.vcf.gz | sort -u > sy2_pos.txt
bcftools query -f '%CHROM\t%POS\n' /project/dtataru/SY3/bams/morefilter_1x_hybrids3.maxdepth6000_renamed.vcf.gz | sort -u > sy3_pos.txt

#find shared positions between pairs
comm -12 wgs1_pos.txt sy2_pos.txt > sy1_sy2_shared_pos.txt #73333 sites
comm -12 wgs1_pos.txt sy3_pos.txt > sy1_sy3_shared_pos.txt #155373 sites

#combine all three
cat sy1_sy2_shared_pos.txt sy1_sy3_shared_pos.txt | sort -u -k1,1 -k2,2n > union_shared_pos.txt #190980 sites MUCH BETTER

#turn it into a regions file
bgzip union_shared_pos.txt
tabix -s1 -b2 -e2 union_shared_pos.txt.gz

bcftools view -R union_shared_pos.txt.gz hybrid1_all.hf.vcf.gz -Oz -o hybrid1_all.hf_union_subset.vcf.gz
tabix hybrid1_all.hf_union_subset.vcf.gz

```

now to run entropy:

```
#!/bin/bash
#SBATCH --output=/project/dtataru/hybrids/logs/entropy_%A_%a.out
#SBATCH --error=/project/dtataru/hybrids/logs/entropy_%A_%a.err
#SBATCH --time=3-00:00:00
#SBATCH -p single
#SBATCH -N 1
#SBATCH --cpus-per-task=20
#SBATCH -A loni_ferrislac

eval "$(conda shell.bash hook)"
conda activate entropy

cd /project/dtataru/hybrids/4_GATKvarcall/3_Genotyped_GVCFs/isec_out
population=(SH GB HH)

for POP in "${population[@]}"
do
    # run three chains for each population
    entropy -i /project/dtataru//project/dtataru/hybrids/4_GATKvarcall/3_Genotyped_GVCFs/isec_out/${POP}_WGS1_ddRADunion_sorted.gl -m 1 \
        -n 2 \
        -k 3 -q faststructure_input/${POP}_WGS1_ddRADunion_sorted.meanQ  \
        -l 2000 -b 1500 -t 10 \
        -o ${POP}_output/${POP}_mcmcoutk3chain1.hdf5

    entropy -i /project/dtataru/SY2/bams/${POP}_WGS1_ddRADunion_sorted.gl -m 1 \
        -n 2 \
        -k 3 -q faststructure_input/${POP}_WGS1_ddRADunion_sorted.meanQ   \
        -l 2000 -b 1500 -t 10 \
        -o ${POP}_output/${POP}_mcmcoutk3chain2.hdf5

    entropy -i /project/dtataru/SY2/bams/${POP}_WGS1_ddRADunion_sorted.gl -m 1 \
        -n 2 \
        -k 3 -q faststructure_input/${POP}_WGS1_ddRADunion_sorted.meanQ  \
        -l 2000 -b 1500 -t 10 \
        -o ${POP}_output/${POP}_mcmcoutk3chain3.hdf5

    # Run assessconvergence.R on each chain individually
    for chain in 1 2 3
    do
        Rscript auxfiles/assessconvergence.R ${POP}_output/${POP}_mcmcoutk3chain${chain}.hdf5
    done

	#assess convergence, R ~1 is good
	estpost.entropy -p q -s 4 ${POP}_output/${POP}_mcmcoutk3chain1.hdf5 \
		${POP}_output/${POP}_mcmcoutk3chain2.hdf5 \
		${POP}_output/${POP}_mcmcoutk3chain3.hdf5 q
		-o ${POP}_output/${POP}_qmcmcdiag.txt

done
```

and then run_entropypostprocessing.sh:

```
#!/bin/bash
#SBATCH --output=/project/dtataru/hybrids/logs/entropypost_%A_%a.out
#SBATCH --error=/project/dtataru/hybrids/logs/entropypost_%A_%a.err
#SBATCH --time=1-00:00:00
#SBATCH -p single
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -A loni_ferrislac
#SBATCH --cpus-per-task=12
#SBATCH --mem=200G

eval "$(conda shell.bash hook)"
conda activate entropy
module load r

POP="GB"
cd /project/dtataru/hybrids/4_GATKvarcall/3_Genotyped_GVCFs/isec_out/${POP}_output/

# getting genotype estimates
estpost.entropy -p gprob -s 0 "${POP}_mcmcoutk3chain1.hdf5" "${POP}_mcmcoutk3chain2.hdf5" "${POP}_mcmcoutk3chain3.hdf5" -o genoest.txt

# getting ancestry (q) estimates
estpost.entropy -p q -s 0 "${POP}_mcmcoutk3chain1.hdf5" "${POP}_mcmcoutk3chain2.hdf5" "${POP}_mcmcoutk3chain3.hdf5" -o admixest.txt

# getting WAIC values
estpost.entropy -p deviance -s 3 "${POP}_mcmcoutk3chain1.hdf5"

#plot admixture plots
#Rscript plotadmix.R
```
