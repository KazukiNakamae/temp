MaChIAto is a third party software of CRISPResso (https://github.com/pinellolab/CRISPResso2) to get stricter classification of the targeted amplicon sequencing (called also deep sequencing) and profile the relationship with genomic context (e.g. Thermodynamics and mutation pattarns) in the target site.

MaChIAto would be helpful for people who want to
- get more accurate classification of knock-out, hmology-based knock-in, and Prime Editing. (Case I)
- quantify the rate of imprecise editing reads including the reads with donor-reversed integration. (Case II)
- obtain the local alignment of target site and knock-in junction. (Case III)
- analyze the microhmology involved with InDels and imprecise editing. (Case IV)
- find the useful features for the efficiency prediction. (Case V)

If you have the same purpose with the above cases, you can junp to [the section of Case Studies](#CaseStudies)!

# Quick Start

We demonstrate the simple example.
Before that, you have to build the environment according to [the section of Requirement](#Requirement)

1. Download MaChIAto with example files by

clicking "Download ZIP".
![download_zip.png](https://github.com/KazukiNakamae/temp/blob/temp-images/download_zip.png)

Or, enter the folowing command terminal.

```bash
git clone https://github.com/KazukiNakamae/MaChIAto.git;
```

1. Open terminal and go the directory.

```bash
cd <...>/MaChIAto
```

1. Re-classify the allele frequency table derived from CRISPResso2.

```bash
source activate MaChIAto_env; # *You do not need enter it again once you did it.
# ↓The re-classification command of MaChIAto (Classifier) 
python MaChIAto/MaChIAto.py -ccf <Alleles_frequency_table.zip of the CRISPResso2 output> \
-o <output directory> \
-a <wt amplicon sequence>\
-g <protospaser sequence>\
<other parameters>;

conda deactivate;
```

The other parameters varies per experiment type {knock-out, homology-based knokc-in, Prime Editing}
The example is for the knock-out analysis.

##### The case of knock-out
```bash
python MaChIAto/MaChIAto.py -ccf <Alleles_frequency_table.zip of the CRISPResso2 output> \
-o <output directory> \
-a <wt amplicon sequence>\
-g <protospaser sequence>;
```

*The case of the knock-out and knock-in is shown in the section of [MaChIAto_(MaChIAto Classifier)](#MaChIAto_(MaChIAto Classifier)).

1. Check the classification result.

You can get the classification result of the editing. The result is visualized as the pie chart based on the ALL_dataframe,txt.

We provide the detailed description how to read them in the section of [MaChIAtoClassifier result](#MaChIAtoClassifier_Result).

1. Get local alignment and the mutation profiling.

```bash
source activate MaChIAto_Aligner_env; # *You do not need enter it again once you did it.
# ↓The alignment command of MaChIAto Aligner
Rscript MaChIAto_Aligner/MaChIAtoAligner.R \
<directory of MaChIAto Classifier> \
<output directory>;

conda deactivate;
```

If your donor insert has the extra sequence in the outside of homology arm, you should enter it.
Please refer to the section of [MaChIAtoAligner](#MaChIAtoAligner).

1. Check the alignment result.

You can get the local alignment of the editing. The result is visualized as the map.

We provide the detailed description how to read them in the section of [MaChIAtoAligner result](#MaChIAtoAligner_Result).

1. Aggregate the multiple data.

If you have multiple data (n>3) to profile the characteristics, you should aggregate the data using the following command.

```bash
source activate MaChIAto_Analyzer_env; # *You do not need enter it again once you did it.
python MaChIAto_Analyzer/collect_MaChIAto_data.py \
-i <the prefix of MaChIAto Classifier result> \
-o <output directory> \
-ol <knock-out label> \
-t <calculation target> \

conda deactivate;
```

It is used for the knock-out analysis. The calculation target for MaChIAto Analyzer/Reviewer is protpspacer.
You can target the homology arm and RT template using another command. Please see the section of [MaChIAtoAnalyzer](#MaChIAtoAnalyzer).

1. Investivate the relationship with (epi-)genomic context

You can see the correlation with the >70 (epi-)genomic context using MaChIAto Analyzer.

```bash
Rscript MaChIAto_Analyzer/MaChIAtoAnalyzer.R \
<directory of collect_MaChIAto_data.py> \
<output directory>;

conda deactivate;
```

MaChIAtoAnalyzer basically use > 70 genomic context.
If you want to use epigenomic context or other context, please see the section of [MaChIAtoAnalyzer](#MaChIAtoAnalyzer).

1. Check the correlation.

You can get the correlation between the editing efficacy and the context. The result is visualized as the scatter plot.

We provide the detailed description how to read them in the section of [MaChIAtoAnalyzer result](#MaChIAtoAnalyzer_Result).

1. Profile the mutation tendency.


```bash
source activate MaChIAto_Reviewer_env; # *You do not need enter it again once you did it.
Rscript MaChIAto_Reviewer/MaChIAtoReviewer.R \
<the prefix of MaChIAto Classifier result> \
<the prefix of MaChIAto Aligner result> \
<directory of collect_MaChIAto_data.py> \
<output directory>;

1. Check the mutation profile.

You can get the local alignment of the editing. The result is visualized as the bar plot and pie chart.

We provide the detailed description how to read them in the section of [MaChIAtoReviewer result](#MaChIAtoReviewer_Result).






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


# MaChIAto_(MaChIAto Classifier)
conda create --name MaChIAto_env;
source activate MaChIAto_env;
conda install -c anaconda python=3.8
pip install --upgrade pip
pip install regex tqdm argparse biopython numpy matplotlib GPy gpyopt datetime pandas;

##### The case of knock-out
```bash
python MaChIAto/MaChIAto.py -ccf <Alleles_frequency_table.zip of the CRISPResso2 output> \
-o <output directory> \
-a <wt amplicon sequence>\
-g <protospaser sequence>;
```

##### The case of homology-based knock-in
```bash
python MaChIAto/MaChIAto.py -ccf <Alleles_frequency_table.zip of the CRISPResso2 output> \
-o <output directory> \
-a <wt amplicon sequence>\
-g <protospaser sequence> \
-e <expected editing amplicon seqeunce> \
-d <donor insert sequence> \
-lh <length of 5\' homology arm> \
-rh <length of 3\' homology arm>;
```

##### The case of Prime Editing (Substitution/Deletion editing)
```bash
python MaChIAto/MaChIAto.py \
-ccf <Alleles_frequency_table.zip of the CRISPResso2 output> \
-o <output directory> \
-a <wt amplicon sequence>\
-g <protospaser sequence> \
-e <expected editing amplicon seqeunce> \
-lh <length of prime binding site> \
-rh <length of RT template> \
-cn <distance between two nick sites> \
--primeediting_analysis;
```

##### The case of Prime Editing (Insertion editing)
```bash
python MaChIAto/MaChIAto.py \
-ccf <Alleles_frequency_table.zip of the CRISPResso2 output> \
-o <output directory> \
-a <wt amplicon sequence>\
-g <protospaser sequence> \
-e <expected editing amplicon seqeunce> \
-d <donor insert sequence> \
-lh <length of prime binding site> \
-rh <length of RT template> \
-cn <distance between two nick sites> \
--primeediting_analysis;
```

If you want to use the result of CRISPResso version1, the -ccf should be replased into -cf.
Example for CRISPResso version1
```bash
python MaChIAto/MaChIAto.py -cf <Alleles_frequency_table.txt of the CRISPResso output> \
-o <output directory> \
-a <wt amplicon sequence>\
-g <protospaser sequence>;
```

# MaChIAtoClassifier_Result

# MaChIAtoAligner

cd /Volumes/databank2;
conda update -n base -c defaults conda;
conda create --name MaChIAto_Aligner_env;
source activate MaChIAto_Aligner_env;
conda install -c bioconda python=3.8 bwa=0.7.17 samtools=1.9;

# type agree

```bash
source activate MaChIAto_Aligner_env; # *You do not need enter it again once you did it.
# ↓The alignment command of MaChIAto Aligner
Rscript MaChIAto_Aligner/MaChIAtoAligner.R \
<directory of MaChIAto Classifier> \
<output directory> \
<extra sequence at the 5\' side of donor insert> \
<extra sequence at the 3\' side of donor insert>;

conda deactivate;
```

# when you use external storage for save output, the process can be too slow.

# MaChIAtoAligner_Result

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

```bash
source activate MaChIAto_Analyzer_env; # *You do not need enter it again once you did it.
python MaChIAto_Analyzer/collect_MaChIAto_data.py \
-i <the prefix of MaChIAto Classifier result> \
-o <output directory> \
-sc gttttagagctaggccaacatgaggatcacccatgtctgcagggcctagcaagttaaaataaggctagtccgttatcaacttggccaacatgaggatcacccatgtctgcagggccaagtggcaccgagtcggtgc \ #note
-ul A \
-ol B \
-il C D \
-t bmh \
--ignore_list ./20190916_MaChIAto2_data_analysis/ignore_list.csv;

```


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


# CaseStudies