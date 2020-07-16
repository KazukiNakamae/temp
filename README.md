MaChIAto is a third party software of CRISPResso (https://github.com/pinellolab/CRISPResso2) to get stricter classification of the targeted amplicon sequencing (called also deep sequencing) and profile the relationship with genomic context (e.g. Thermodynamics and mutation pattarns) in the target site.

MaChIAto would be helpful for people who want to
- get more accurate classification of knock-out, hmology-based knock-in, and Prime Editing. (Case I)
- quantify the rate of imprecise editing reads including the reads with donor-reversed integration. (Case II)
- obtain the local alignment of target site and knock-in junction. (Case III)
- analyze the microhmology involved with InDels and imprecise editing. (Case IV)
- find the useful features for the efficiency prediction. (Case V)

If you have the same purpose with the above cases, you can junp to [the section of Case Studies](#Case Studies)!

# Quick Start

We demonstrate the simple and basic example. You can download MaChIAto with example files.

![download_zip.png](https://github.com/KazukiNakamae/temp/blob/temp-images/download_zip.png)

# Requirement

# set bash as default shell
chsh -s /bin/tcsh;

# install miniconda3 from cran (https://docs.conda.io/en/latest/miniconda.html)
# install r from cran (https://cran.ism.ac.jp)
# version R-4.0.1
# install xcode
# install dependency
sudo xcodebuild -license;
# install XQuartz from xquartz.macosforge.org


# MaChIAto (MaChIAto Classifier)
conda create --name MaChIAto_env;
source activate MaChIAto_env;
conda install -c anaconda python=3.8
pip install --upgrade pip
pip install regex tqdm argparse biopython numpy matplotlib GPy gpyopt datetime pandas;



# MaChIAto Aligner

cd /Volumes/databank2;
conda update -n base -c defaults conda;
conda create --name MaChIAto_Aligner_env;
source activate MaChIAto_Aligner_env;
conda install -c bioconda python=3.8 bwa=0.7.17 samtools=1.9;


# type agree

５０遺伝子座ノックイン

for fn in `ls /Volumes/databank2/MaChIAto2_v1.7.0_output_200311 | grep "MaChIAto"`;do Rscript /Volumes/databank2/temp/MaChIAto/MaChIAto_Aligner/MaChIAtoAligner.R /Volumes/databank2/MaChIAto2_v1.7.0_output_200311/$fn /Users/nedo01/Desktop/MaChIAtoAlignerv14_output_200501_from_MaChIAto2_v1.7.0_output_200311/ GTTTGG CCAAAC;done;
# when you use external storage for save output, the process can be too slow.


# MaChIAto Analyzer

conda create --name MaChIAto_Analyzer_env;
source activate MaChIAto_Analyzer_env;
conda install -c anaconda python=3.8;
pip install --upgrade pip;
pip install numpy regex tqdm pandas;

conda install -c anaconda wget;
WORKING_DIRECTORY=$PWD;
cd /Volumes/databank2/temp/MaChIAto/MaChIAto_Analyzer;
wget http://unafold.rna.albany.edu/cgi-bin/OligoArrayAux-download.cgi?oligoarrayaux-3.8.tar.gz;
tar zxvf ./OligoArrayAux-download.cgi?oligoarrayaux-3.8.tar.gz;
cd oligoarrayaux-3.8;
./configure --prefix=/Volumes/databank2/temp/MaChIAto/MaChIAto_Analyzer;
make;
make install;
cd $WORKING_DIRECTORY

conda install -c bioconda emboss;

python /Volumes/databank2/temp/MaChIAto/MaChIAto_Analyzer/collect_MaChIAto_data.py \
-i ./MaChIAto2_v1.7.0_output_200311 \
-o ./20200311_MaChIAto2_data_analysis_mode1 \
-sc gttttagagctaggccaacatgaggatcacccatgtctgcagggcctagcaagttaaaataaggctagtccgttatcaacttggccaacatgaggatcacccatgtctgcagggccaagtggcaccgagtcggtgc \
-ul A \
-ol B \
-il C D \
-t bmh \
--ignore_list ./20190916_MaChIAto2_data_analysis/ignore_list.csv;

Rscript /Volumes/databank2/temp/MaChIAto/MaChIAto_Analyzer/MaChIAtoAnalyzer.R \
./20200311_MaChIAto2_data_analysis_mode1 \
./MaChIAtoAnalyzer_beta1.6_output \
InDelphi \
./extra_table/bmh_extra_data_InDelphi.csv \
Accessibility \
./extra_table/bmh_extra_data_accessibility.csv \
Genome_Property \
./extra_table/bmh_extra_data_genomeprop.csv;

# MaChIAto Reviewer

conda create -n MaChIAto_Reviewer_env --clone MaChIAto_Aligner_env;
source activate MaChIAto_Reviewer_env;
Rscript /Volumes/databank2/temp/MaChIAto/MaChIAto_Reviewer/MaChIAtoReviewer.R \
./MaChIAto2_v1.7.0_output_200311 \
./MaChIAtoAlignerv14_output_200501_from_MaChIAto2_v1.7.0_output_200311 \
./20200311_MaChIAto2_data_analysis_mode1 \
./MaChIAtoReviewer_beta2.4_output_200610_from_MaChIAto2_v1.7.0_output_200311;


# Case Studies