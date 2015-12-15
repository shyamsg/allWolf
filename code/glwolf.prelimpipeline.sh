#! /bin/bash

#####################################################################
# File: glwolf.prelimpipeline.sh                                    #
# Author: Shyam Gopalakrishnan                                      #
# Date: 7th December 2015                                           #
# Description: Hastily put together pipeline for the initial results#
# on the gl wolf and dog data - mostly for the poster that Mikkel is#
# presenting on Friday, 11th December, 2015.                        #
#####################################################################

PROJECT_HOME=/home/shyam/projects/allWolf
DATA_HOME=$PROJECT_HOME/data/mappedToWolf

## First, subset the data to set of loci where you want data.
## This set is obtained from the other project - wolfRefGenome
module load bedtools/2.25.0
module load bcftools/1.2

SITESBCF=/home/shyam/projects/wolfRefGenome/data/wolfRef/allAlignedToWolf.miss0.75.maf0.05.minq30.mingq20.forPCA.bcf
cd $DATA_HOME
if [ ! -e $DATA_HOME/`basename $SITESBCF .bcf`.vcf ]; then
    bcftools view $SITESBCF > $DATA_HOME/`basename $SITESBCF .bcf`.vcf
fi
echo "Made chosenSNPs vcf."

for vcf in *scf.vcf.bgz; do
    if [ -e `basename $vcf .vcf.bgz`.snpsonly.vcf -o -e `basename $vcf .vcf.bgz`.snpsonly.vcf.gz ]; then
	echo "$vcf already processed."
    else
	zcat $vcf | head -10000 | grep "^#" > `basename $vcf .vcf.bgz`.snpsonly.vcf 
	bedtools intersect -wa -a $vcf -b allAlignedToWolf.miss0.75.maf0.05.minq30.mingq20.forPCA.vcf >> `basename $vcf .vcf.bgz`.snpsonly.vcf & 
	echo "Done with $vcf."
    fi
done

## bgzip and tabix
for vcf in *.scf.snpsonly.vcf; do
    if [ ! -e $vcf.gz.tbi ]; then
	bgzip $vcf
	tabix -p vcf $vcf.gz
    fi
    echo "Done with $vcf."
done

## merge using bcftools
module load bcftools/1.2
if [ ! -e merged.snpsonly.vcf ]; then
    echo "starting merge."
    bcftools merge -Ou *.scf.snpsonly.vcf.gz | bcftools filter -Ov -i 'MAF[0]>0.01 & TYPE="snp" & (REF!="N" & N_ALT=1)' -o merged.snpsonly.vcf 
fi
echo "Done merging."

## Remove Krummelangso from the set
if [ ! -e merged.noKL.snpsonly.vcf ]; then
    vcftools --vcf merged.snpsonly.vcf --recode-INFO-all --remove-indv Krummelangso --stdout --recode | sed 's/-1\/-1/.\/./g' > merged.noKL.snpsonly.vcf
fi

## filter the files for pca
module load vcftools/0.1.14
if [ ! -e merged.snpsonly.filtered.forPCA.vcf ]; then
    vcftools --vcf merged.snpsonly.vcf --stdout --recode-INFO-all --recode --maf 0.01 --max-missing 0.75 --minQ 50 --minGQ 30 | sed 's/-1\/-1/.\/./g' > merged.snpsonly.filtered.forPCA.vcf
fi 

## filter the files for pca
module load vcftools/0.1.14
if [ ! -e merged.noKL.snpsonly.filtered.forPCA.vcf ]; then
    vcftools --vcf merged.noKL.snpsonly.vcf --stdout --recode-INFO-all --recode --maf 0.01 --max-missing 0.75 --minQ 50 --minGQ 30 | sed 's/-1\/-1/.\/./g' > merged.noKL.snpsonly.filtered.forPCA.vcf
fi 

## PCA
nsites=$(grep -vc "^#" merged.noKL.snpsonly.filtered.forPCA.vcf)
nsamps=$(head -10000 merged.noKL.snpsonly.filtered.forPCA.vcf | grep CHROM | cut -f10- | wc -w)
if [ ! -e merged.noKL.snpsonly.filtered.forPCA.geno ]; then
    angsd -vcf-gl merged.noKL.snpsonly.filtered.forPCA.vcf -fai /home/joseas/data/Wolf/VCFsWolf/data/prefixes/L.Dalen_14_wolf.scf.fasta.fai -nInd $nsamps -doMajorMinor 1 -doMaf 1 -doPost 1 -doGeno 32 -out merged.noKL.snpsonly.filtered.forPCA && gunzip merged.noKL.snpsonly.filtered.forPCA.geno.gz
fi

## PCA with KL
nsites=$(grep -vc "^#" merged.snpsonly.filtered.forPCA.vcf)
nsamps=$(head -10000 merged.snpsonly.filtered.forPCA.vcf | grep CHROM | cut -f10- | wc -w)
if [ ! -e merged.snpsonly.filtered.forPCA.geno ]; then
    angsd -vcf-gl merged.snpsonly.filtered.forPCA.vcf -fai /home/joseas/data/Wolf/VCFsWolf/data/prefixes/L.Dalen_14_wolf.scf.fasta.fai -nInd $nsamps -doMajorMinor 1 -doMaf 1 -doPost 1 -doGeno 32 -out merged.snpsonly.filtered.forPCA && gunzip merged.snpsonly.filtered.forPCA.geno.gz
fi

cd $PROJECT_HOME/analysis/pca
## NGSCovar used to generate pca plots
if [ ! -e merged.noKL.snpsonly.filtered.forPCA.maf0.05.normed.covar ]; then
    /home/fgvieira/data/appz/ngsTools/ngsPopGen/ngsCovar -probfile $DATA_HOME/merged.noKL.snpsonly.filtered.forPCA.geno -outfile merged.noKL.snpsonly.filtered.forPCA.maf0.05.normed.covar -nind $nsamps -nsites $nsites -call 0 -minmaf 0.05 -norm 1
fi
if [ ! -e merged.noKL.snpsonly.filtered.forPCA.maf0.1.normed.covar ]; then
    /home/fgvieira/data/appz/ngsTools/ngsPopGen/ngsCovar -probfile $DATA_HOME/merged.noKL.snpsonly.filtered.forPCA.geno -outfile merged.noKL.snpsonly.filtered.forPCA.maf0.1.normed.covar -nind $nsamps -nsites $nsites -call 0 -minmaf 0.1 -norm 1
fi
if [ ! -e merged.noKL.snpsonly.filtered.forPCA.maf0.2.normed.covar ]; then
    /home/fgvieira/data/appz/ngsTools/ngsPopGen/ngsCovar -probfile $DATA_HOME/merged.noKL.snpsonly.filtered.forPCA.geno -outfile merged.noKL.snpsonly.filtered.forPCA.maf0.2.normed.covar -nind $nsamps -nsites $nsites -call 0 -minmaf 0.2 -norm 1
fi

## NGSCovar used to generate pca plots
if [ ! -e merged.snpsonly.filtered.forPCA.maf0.05.normed.covar ]; then
    /home/fgvieira/data/appz/ngsTools/ngsPopGen/ngsCovar -probfile $DATA_HOME/merged.snpsonly.filtered.forPCA.geno -outfile merged.snpsonly.filtered.forPCA.maf0.05.normed.covar -nind $nsamps -nsites $nsites -call 0 -minmaf 0.05 -norm 1
fi
if [ ! -e merged.snpsonly.filtered.forPCA.maf0.1.normed.covar ]; then
    /home/fgvieira/data/appz/ngsTools/ngsPopGen/ngsCovar -probfile $DATA_HOME/merged.snpsonly.filtered.forPCA.geno -outfile merged.snpsonly.filtered.forPCA.maf0.1.normed.covar -nind $nsamps -nsites $nsites -call 0 -minmaf 0.1 -norm 1
fi
if [ ! -e merged.snpsonly.filtered.forPCA.maf0.2.normed.covar ]; then
    /home/fgvieira/data/appz/ngsTools/ngsPopGen/ngsCovar -probfile $DATA_HOME/merged.snpsonly.filtered.forPCA.geno -outfile merged.snpsonly.filtered.forPCA.maf0.2.normed.covar -nind $nsamps -nsites $nsites -call 0 -minmaf 0.2 -norm 1
fi

PLOTPCA=0
if [ $PLOTPCA == 1 ]; then
    ## Run eigen decomposition on the covar matrix using R
    Rscript /home/shyam/projects/glWolf/code/pca.ngscovar.R merged.noKL.snpsonly.filtered.forPCA.maf0.05.normed.covar glwolf.noKL.samples.txt merged.noKL.snpsonly.filtered.maf0.05.ngscovar.pca.pdf
    Rscript /home/shyam/projects/glWolf/code/pca.ngscovar.R merged.noKL.snpsonly.filtered.forPCA.maf0.1.normed.covar glwolf.noKL.samples.txt merged.noKL.snpsonly.filtered.maf0.1.ngscovar.pca.pdf
    Rscript /home/shyam/projects/glWolf/code/pca.ngscovar.R merged.noKL.snpsonly.filtered.forPCA.maf0.2.normed.covar glwolf.noKL.samples.txt merged.noKL.snpsonly.filtered.maf0.2.ngscovar.pca.pdf
    
    ## Run eigen decomposition on the covar matrix using R
    Rscript /home/shyam/projects/glWolf/code/pca.ngscovar.R merged.snpsonly.filtered.forPCA.maf0.05.normed.covar glwolf.samples.txt merged.snpsonly.filtered.maf0.05.ngscovar.pca.pdf
    Rscript /home/shyam/projects/glWolf/code/pca.ngscovar.R merged.snpsonly.filtered.forPCA.maf0.1.normed.covar glwolf.samples.txt merged.snpsonly.filtered.maf0.1.ngscovar.pca.pdf
    Rscript /home/shyam/projects/glWolf/code/pca.ngscovar.R merged.snpsonly.filtered.forPCA.maf0.2.normed.covar glwolf.samples.txt merged.snpsonly.filtered.maf0.2.ngscovar.pca.pdf
fi

TREEMIX=0
if [ $TREEMIX == 1 ]; then
    ## Treemix on those datasets
    cd $DATA_HOME
    if [ ! -e merged.snpsonly.maf0.05.treemix.gz ]; then
	python /home/shyam/projects/SantasHelpers/misc/vcf2Treemix.py merged.snpsonly.filtered.forPCA.vcf glwolf.samples.txt 0.05 | gzip -c > merged.snpsonly.maf0.05.treemix.gz
    fi
    if [ ! -e merged.snpsonly.maf0.2.treemix.gz ]; then
	python /home/shyam/projects/SantasHelpers/misc/vcf2Treemix.py merged.snpsonly.filtered.forPCA.vcf glwolf.samples.txt 0.2 | gzip -c > merged.snpsonly.maf0.2.treemix.gz
    fi
    cd $PROJECT_HOME/analysis/treemix
    for rep in {1..10}; do 
	treemix -i merged.snpsonly.maf0.05.treemix.gz -o tree0.05_rep$rep &
	sleep 1
	treemix -i merged.snpsonly.maf0.2.treemix.gz -o tree0.2_rep$rep &
	sleep 1
	treemix -i merged.snpsonly.maf0.05.treemix.gz -m 1 -o tree0.05_mig1_rep$rep &
	sleep 1
	treemix -i merged.snpsonly.maf0.2.treemix.gz -m 1 -o tree0.2_mig1_rep$rep &
	sleep 1
	echo "Done with rep $rep."
    done

    Rscript /home/shyam/projects/SantasHelpers/misc/plotTreemix.R tree0.05_rep1 merged.snpsonly.filtered.maf0.05.treemix.mig0.pdf
    Rscript /home/shyam/projects/SantasHelpers/misc/plotTreemix.R tree0.2_rep1 merged.snpsonly.filtered.maf0.2.treemix.mig0.pdf
    Rscript /home/shyam/projects/SantasHelpers/misc/plotTreemix.R tree0.05_mig1_rep1 merged.snpsonly.filtered.maf0.05.treemix.mig1.pdf
    Rscript /home/shyam/projects/SantasHelpers/misc/plotTreemix.R tree0.2_mig1_rep1 merged.snpsonly.filtered.maf0.2.treemix.mig1.pdf
fi

## Run ngsadmix on the dataset.
## Convert the files to beagle format
cd $DATA_HOME
if [ ! -e merged.snpsonly.filtered.forPCA.beagle.pl.gz ]; then
    python /home/shyam/projects/SantasHelpers/pythonScripts/vcf2Beagle.py --vcf merged.snpsonly.filtered.forPCA.vcf --out merged.snpsonly.filtered.forPCA.beagle.pl
    gzip merged.snpsonly.filtered.forPCA.beagle.pl
fi
NGSADMIX=0
if [ $NGSADMIX == 1 ]; then
    cd $PROJECT_HOME/analysis/admix
    for K in {2..5}; do
	for rep in {1..15}; do
	    if [ ! -e merged.K$K.rep$rep.qopt ]; then
		/home/fgvieira/data/appz/ngsTools/angsd/misc/NGSadmix -likes $DATA_HOME/merged.snpsonly.filtered.forPCA.beagle.pl.gz -K $K -minMaf 0.05 -P 2 -minInd 5 -outfiles merged.K$K.rep$rep &
		sleep 1
	    fi
	done
	echo "Launched NGSAdmix for K=$K" 
	wait
	bestKFile=`grep best merged.K$K.rep*.log | sort -k2.6,2f | head -1 | cut -f1 -d:`
	bestKFile=`basename $bestKFile .log`
	python /home/shyam/projects/SantasHelpers/pythonScripts/plotNGSadmix.py -q $bestKFile.qopt -o `echo $bestKFile | cut -f1-2 -d.`.pdf -l $DATA_HOME/glwolf.samples.txt
	echo "Done NGSAdmix for K=$K"
    done

    
fi


