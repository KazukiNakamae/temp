MaChIAto is a third party software of CRISPResso (https://github.com/pinellolab/CRISPResso2) to get stricter classification of the targeted amplicon sequencing (also called deep sequencing) and profile the relationship with genomic context (e.g. Thermodynamics and mutation pattarns) in the target site.

MaChIAto would be helpful for people who want to
- get more accurate classification of knock-out, hmology-based knock-in, and Prime Editing. (Case I)
- quantify the rate of imprecise editing reads including the reads with donor-reversed integration. (Case II)
- obtain the local alignment of target site and knock-in junction. (Case III)
- analyze the microhmology involved with InDels and imprecise editing. (Case IV)
- find the useful features for the efficiency prediction. (Case V)

If you have the same purpose with the above cases, you can junp to [the section of Case Studies](#CaseStudies)!

# Quick Start

We demonstrate the simple example.
Before that, you have to build the environment according to [the section of Preparation](#Preparation)

#### 1. Download MaChIAto with example files

Click "Download ZIP".

![download_zip.png](https://github.com/KazukiNakamae/temp/blob/temp-images/download_zip.png)

Or, enter the folowing command terminal.

```bash
git clone https://github.com/KazukiNakamae/MaChIAto.git;
```

#### 1. Open terminal and go the directory.

```bash
cd <...>/MaChIAto
```

#### 1. Re-classify the allele frequency table derived from CRISPResso2.

```bash
source activate MaChIAto_env; # *You do not need enter it again once you did it.
```

The re-classification command of MaChIAto (Classifier) 
```bash
python MaChIAto/MaChIAto.py \
-m CRISPResso2 \
-ccf <Alleles_frequency_table.zip of the CRISPResso2 output> \
-o <output directory> \
-a <wt amplicon sequence>\
-g <protospaser sequence>\
<other parameters>;
```

```bash
conda deactivate;
```

The other parameters varies per experiment type {knock-out, homology-based knokc-in, Prime Editing}.

The example is for the knock-out analysis.

##### The case of knock-out
```bash
python MaChIAto/MaChIAto.py -ccf <Alleles_frequency_table.zip of the CRISPResso2 output> \
-o <output directory> \
-a <wt amplicon sequence>\
-g <protospaser sequence>;
```

*The case of the knock-out and knock-in is shown in the section of [MaChIAto_(MaChIAto Classifier)](#MaChIAto_(MaChIAto Classifier)).

#### 1. Check the classification result.

You can get the classification result of the editing. The result is visualized as the pie chart based on the ALL_dataframe,txt.

We provide the detailed description how to read them in the section of [MaChIAtoClassifier result](#MaChIAtoClassifier_Result).

# Quick Start (Optional)

MaChIAto can profile and visualize the MaChIAto result using MaChIAto Aligner, MaChIAto Analyzer, and MaChIAto Reviewer. 

#### 1. Get local alignment and the mutation profiling using MaChIAto Aligner.

```bash
source activate MaChIAto_Aligner_env; # *You do not need enter it again once you did it.
```

The alignment command of MaChIAto Aligner
```bash
Rscript MaChIAto_Aligner/MaChIAtoAligner.R \
<directory of MaChIAto Classifier> \
<output directory>;
```

```bash
conda deactivate;
```

If your donor insert has the extra sequence in the outside of homology arm, you should enter it.
Please refer to the section of [MaChIAtoAligner](#MaChIAtoAligner).

#### 1. Check the alignment result.

You can get the local alignment of the editing. The result is visualized as the map.

We provide the detailed description how to read them in the section of [MaChIAtoAligner result](#MaChIAtoAligner_Result).

#### 1. Aggregate the multiple data.

If you have multiple data (n>3) to profile the characteristics, you should aggregate the data using the following command.

```bash
source activate MaChIAto_Analyzer_env; # *You do not need enter it again once you did it.
python MaChIAto_Analyzer/collect_MaChIAto_data.py \
-i <the prefix of MaChIAto Classifier result> \
-o <output directory> \
-ol <knock-out label> \
-t <calculation target> \
```

```bash
conda deactivate;
```

It is used for the knock-out analysis. The calculation target for MaChIAto Analyzer/Reviewer is protpspacer.
You can target the homology arm and RT template using another command. Please see the section of [MaChIAtoAnalyzer](#MaChIAtoAnalyzer).

#### 1. Investivate the relationship with (epi-)genomic context

You can see the correlation with the >70 (epi-)genomic context using MaChIAto Analyzer.

```bash
Rscript MaChIAto_Analyzer/MaChIAtoAnalyzer.R \
<directory of collect_MaChIAto_data.py> \
<output directory>;
```

```bash
conda deactivate;
```

MaChIAtoAnalyzer basically use > 70 genomic context.
If you want to use epigenomic context or other context, please see the section of [MaChIAtoAnalyzer](#MaChIAtoAnalyzer).

#### 1. Check the correlation.

You can get the correlation between the editing efficacy and the context. The result is visualized as the scatter plot.

We provide the detailed description how to read them in the section of [MaChIAtoAnalyzer result](#MaChIAtoAnalyzer_Result).

#### 1. Profile the mutation tendency.


```bash
source activate MaChIAto_Reviewer_env; # *You do not need enter it again once you did it.
Rscript MaChIAto_Reviewer/MaChIAtoReviewer.R \
<the prefix of MaChIAto Classifier result> \
<the prefix of MaChIAto Aligner result> \
<directory of collect_MaChIAto_data.py> \
<output directory>;
```

```bash
conda deactivate;
```

#### 1. Check the mutation profile.

You can get the local alignment of the editing. The result is visualized as the bar plot and pie chart.

We provide the detailed description how to read them in the section of [MaChIAtoReviewer result](#MaChIAtoReviewer_Result).


# Preparation

Type the following commands.

#### 1. Set bash as default shell
```bash
chsh -s /bin/tcsh;
```

### Install software

#### 1. install miniconda3 from Conda
(https://docs.conda.io/en/latest/miniconda.html)
#### 1. install R (>version R-4.0.1) from the CRAN
(https://cran.ism.ac.jp)
#### 1. install xcode from the App Store
(https://apps.apple.com/jp/app/xcode/id497799835?mt=12)
```bash
sudo xcodebuild -license;
```
#### 1. install XQuartz from the XQuartz project
(xquartz.macosforge.org)

### Build the environment using conda

#### Environment for MaChIAto_(MaChIAto Classifier)
```bash
conda create --name MaChIAto_env;
source activate MaChIAto_env;
conda install -c anaconda python=3.8;
pip install --upgrade pip;
pip install regex tqdm argparse biopython numpy matplotlib GPy gpyopt datetime pandas;
conda deactivate;
```

#### Environment for MaChIAtoAligner
```bash
cd /Volumes/databank2;
conda update -n base -c defaults conda;
conda create --name MaChIAto_Aligner_env;
source activate MaChIAto_Aligner_env;
conda install -c bioconda python=3.8 bwa=0.7.17 samtools=1.9;
conda deactivate;
```

#### Environment for MaChIAtoAnalyzer
```bash
conda create --name MaChIAto_Analyzer_env;
source activate MaChIAto_Analyzer_env;
conda install -c anaconda python=3.8 wget;
conda install -c bioconda emboss;
pip install --upgrade pip;
pip install numpy regex tqdm pandas;
conda deactivate;
```

#### Environment for MaChIAtoReviewer
```bash
conda create -n MaChIAto_Reviewer_env --clone MaChIAto_Aligner_env;
source activate MaChIAto_Reviewer_env;
conda deactivate;
```

### Install OligoArrayAux
```bash
WORKING_DIRECTORY=$PWD;
cd /Volumes/databank2/temp/MaChIAto/MaChIAto_Analyzer;

# For Mac
wget http://unafold.rna.albany.edu/cgi-bin/OligoArrayAux-download.cgi?oligoarrayaux-3.8.tar.gz;
tar zxvf ./OligoArrayAux-download.cgi?oligoarrayaux-3.8.tar.gz;
# Linux user can download the binary version intead: http://unafold.rna.albany.edu/?q=DINAMelt/OligoArrayAux

cd oligoarrayaux-3.8;
./configure --prefix=/Volumes/databank2/temp/MaChIAto/MaChIAto_Analyzer;
make;
make install;
cd $WORKING_DIRECTORY
conda deactivate;
```

# Usage

## Usage of MaChIAto Classifier

### Parameter list

**-m** or **--mode** <str>: This parameter allows for the specification of type of analysis: “CRISPResso” and “CRISPResso2” is allowed in the latest version

**-a** or **--amplicon_seq** <str>: This parameter allows the user to enter the amplicon sequence used for the CRISPResso. The length should be >105bp due to setting for the other parameter.

**-g** or **--guide_seq** <str>: This parameter allows for the specification of the sgRNA sequence used for the CRISPResso. The length of sequence should be 20nt without PAM. The MaChIAto convention is to depict the expected cleavage position using the value of the parameter 3 nt 3' from the end of the guide.", required=True)

**-cf** or **--crispreeso_file** <str> : This parameter allows for the specification of the “Alleles_frequency_table.txt” from CRISPResso. When this parameter is used, “CRISPResso” should be entered as -m parameter.', (default: "./Alleles_frequency_table.txt") (optional)

**-ccf** or **--crispreeso2_file** <str> : This parameter allows for the specification of the “Alleles_frequency_table.zip” from CRISPResso2. When this parameter is used, “CRISPResso2” should be entered as -m parameter.', (default: "./Alleles_frequency_table.zip") (optional)

**-d** or **--donor_seq** <str>: This parameter allows for the specification of the expected HDR amplicon used for the CRISPResso. The length of sequence should be >12bp in knock-out/knock-in analysis. In knock-out analysis, this parameter should not be entered, and then the value will be given fake parameter (“TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT”) and some internal settings is changed for knock-out analysis. However, fake parameter will be poly-C|G|A if amplicon sequence contains poly-T sequence.', (default: "") (optional)

**-e** or **--expected_ki_amplicon_seq** <str>: This parameter allows for the specification of the expected knock-in amplicon sequence which used for the CRISPResso after HDR. The length of sequence should be >12bp. In knock-out analysis, this parameter should not be entered, and then the value will be given fake parameter including fake donor sequence and some internal settings is changed for knock-out analysis.', (default: "") (optional)

**-o** or **--output_folder** <str>: This parameter allows for the specification of the output directory to use for the analysis (default: current directory)', (default: "./") (optional)

**-lh** or **--length_left_homologyarm** <int>: This parameter allows for the specification of the length of 5’ homology arm (default: 20). The length of sequence should be >17bp, and the flanking sequence with the homology arm needs length of >24bp in the expected amplicon sequence.', (default: 20) (optional)

**-rh** or **--length_right_homologyarm** <int>: This parameter allows for the specification of the length of 3’ homology arm (default: 20). The length of sequence should be >17bp, and the flanking sequence with the homology arm needs length of >24bp in the expected amplicon sequence.', (default: 20) (optional)

**-cn** or **--location_comp_nick** <int>: This parameter allows for the specification of the complementary strand nick location [3prime direction is +] (default: 90). This parameter is used in the prime editing and should be over the length of homology arm to which the nickase is adjacent', (default: 90) (optional)

**-n** or **--name**  : This parameter allows for the specification of the name which will be included output directory (default: “untitled-X”). If MaChIAto Analyzer and MaChIAto Reviewer will be used in the following analysis, the value should be “<target_name>-<sample label >” (e.g. DBF4B-C) and underbar "_" should not be used', (default: "untitled-X") (optional)

**--primeediting_analysis** : Re-classify the data as prime editing analysis. This option forces the setting to change for knock-out analysis.', action='store_true (optional)

#### Adcanced parameter (for developers)

**--force_knockout_analysis** : Usually, MaChIAto re-classify the data as knock-in analysis. This option forces the setting to change for knock-out analysis. Under this mode, the length of indicator on knock-in donor is maximum value, and threshold value for alignment of knock-in sequence is 1.0. If this parameter is not entered, MaChIAto can automatically set up this mode by finding some characteristics of knock-out sample. For example, MaChIAto checks that there is no donor sequence or expected knock-in sequence as input, and there is less 3 kinds HDR variants among input data.', action='store_true (optional)

**--skip_optimization** : Usually, MaChIAto run Bayesian optimization for finding the optimized setting. This option allows MaChIAto to skip the process of optimization. The option is made for debugging. So, the option should not be used in usual analysis. However, if the optimization process disturbs an accurate analysis, this option might be useful.', action='store_true (optional)

**--copy_optimization** : Usually, MaChIAto run Bayesian optimization for finding the optimized setting. This option allows MaChIAto to use the provided optimization data instead of the optimization. If you analize NGS data with the setting of previous analized data, this option is useful. Espesially, in substitution editing, we can reccomend to apply the option to negative/positive control sample. However, you should understand that a classification error might frequently occurs.', action='store_true (optional)

**--provided_optimization_file** : This parameter allows for the specification of the “MaChIAto_optimized_param.csv” from the analized MaChIAto folder. When this parameter is used, --copy_optimization parameter is required.', (default: "./MaChIAto_optimized_param.csv") (optional)

### Template Command

##### The case of knock-out
```bash
python MaChIAto/MaChIAto.py \
-m CRISPResso2 \
-ccf <Alleles_frequency_table.zip of the CRISPResso2 output> \
-o <output directory> \
-a <wt amplicon sequence> \
-g <protospaser sequence> \
-n <sample name>-<label name>;
```

##### The case of homology-based knock-in
```bash
python MaChIAto/MaChIAto.py \
-m CRISPResso2 \
-ccf <Alleles_frequency_table.zip of the CRISPResso2 output> \
-o <output directory> \
-a <wt amplicon sequence>\
-g <protospaser sequence> \
-e <expected editing amplicon seqeunce> \
-d <donor insert sequence> \
-lh <length of 5\' homology arm> \
-rh <length of 3\' homology arm> \
-n <sample name>-<label name>;
```

##### The case of Prime Editing (Substitution/Deletion editing)
```bash
python MaChIAto/MaChIAto.py \
-m CRISPResso2 \
-ccf <Alleles_frequency_table.zip of the CRISPResso2 output> \
-o <output directory> \
-a <wt amplicon sequence>\
-g <protospaser sequence> \
-e <expected editing amplicon seqeunce> \
-lh <length of prime binding site> \
-rh <length of RT template> \
-cn <distance between two nick sites> \
--primeediting_analysis \
-n <sample name>-<label name>;
```

##### The case of Prime Editing (Insertion editing)
```bash
python MaChIAto/MaChIAto.py \
-m CRISPResso2 \
-ccf <Alleles_frequency_table.zip of the CRISPResso2 output> \
-o <output directory> \
-a <wt amplicon sequence>\
-g <protospaser sequence> \
-e <expected editing amplicon seqeunce> \
-d <donor insert sequence> \
-lh <length of prime binding site> \
-rh <length of RT template> \
-cn <distance between two nick sites> \
--primeediting_analysis \
-n <sample name>-<label name>;
```

If you want to use the result of CRISPResso version1, the -ccf should be replased into -cf.
Example for CRISPResso version1
```bash
python MaChIAto/MaChIAto.py
-m CRISPResso \
-cf <Alleles_frequency_table.txt of the CRISPResso output> \
-o <output directory> \
-a <wt amplicon sequence>\
-g <protospaser sequence> \
-n <sample name>-<label name>;
```

*The **sample name** and **label name** can be **arbitrary**. If you analyze the multiple experiment (e.g. knock-out and knock-in), the name should be different from others.




# MaChIAtoClassifier_Result

# MaChIAto Aligner

### Parameter list

Input directory:
The directory that MaChIAto Classifier generated.

Output prefix:
The directory into which the output directory is saved.

Left extra sequence (optional):
The outside sequence flanking 5'-homology arm of donor. 

Right extra sequence (optional):
The outside sequence flanking 3'-homology arm of donor. 

### Template Command

```bash
Rscript MaChIAto_Aligner/MaChIAtoAligner.R \
<Input directory> \
<Output prefix> \
<Left extra sequence> \
<Right extra sequence>;
```

*when you use external storage for save output, the process can be too slow. We recommend that you use the tool in the local storage.

# MaChIAtoAligner_Result

# MaChIAto Analyzer

## collect_MaChIAto_data.py

### Parameter list

**-i** or **--indir** <str>: input directory which includes the results of MaChIAto Classifier
**-o** or **--outdir** <str>: output directory
**-ul** or **--untreated_label** <str>: untreated label (This must be one.)' (default: "machiato_dummy_sample") (optional)
**-ol** or **--knock_out_label** <str>: negative control label' (default: []) (optional)
**-il** or **--knock_in_label** <str>: knock-in_label' (default: []) (optional)
**-sc** or **--scaffold_seq** <str>: scaffold sequence of sgRNA (default: "gttttagagctagaaatagcaagttaaaataaggctagtccgttatcaacttgaaaaagtggcaccgagtcggtgc") (optional)
**-t** or **--target_type** <str>: target sequence type {"lmh" | "rmh" | "bmh" | "elmh" | "ermh" | "ebmh" | "protospacer"} (default: "bmh") (optional)
**--ignore_list** <str>: The list of ignore target set contains target names which are not desired to analyze for some reasons. The data (e.g. DBF4B-A, DBF4B-B, DBF4B-C, DBF4B-D) including target name (e.g. DBF4B) shown in the list is skipped through the process of MaChIAto Analyzer. The format should be comma-separated like “TargetA, TargetB, …” “example_data” directory has “ignore_list.csv” as example. (default: "") (optional)

### Template command

##### The case of Double knock-in analysis (*ADVANCE: The analysis includes the comparison between two knock-in methods.)
```bash
python MaChIAto_Analyzer/collect_MaChIAto_data.py \
-i <the prefix of MaChIAto Classifier result> \
-o <output directory> \
-sc <scaffold sequence of sgRNA> \
-ul <unmodified label> \
-ol <knock-out label> \
-il <knock-out label 1> <knock-out label 2> \
-t <calculation target> \
--ignore_list <list of samples ignored>;
```

##### The case of Single knock-in analysis (*STANDARD)
```bash
python MaChIAto_Analyzer/collect_MaChIAto_data.py \
-i <the prefix of MaChIAto Classifier result> \
-o <output directory> \
-sc <scaffold sequence of sgRNA> \
-ul <unmodified label> \
-ol <knock-out label> \
-il <knock-out label> \
-t <calculation target> \
--ignore_list <list of samples ignored>;
```

##### The case of Simple knock-in analysis (*SIMPLE: The analysis can be applied when there is no knock-out sample used as control.)
```bash
python MaChIAto_Analyzer/collect_MaChIAto_data.py \
-i <the prefix of MaChIAto Classifier result> \
-o <output directory> \
-sc <scaffold sequence of sgRNA> \
-il <knock-out label> \
-t <calculation target> \
--ignore_list <list of samples ignored>;
```

##### The case of Double knock-out analysis (*ADVANCE: The analysis includes the comparison between two knock-out methods.)
```bash
python MaChIAto_Analyzer/collect_MaChIAto_data.py \
-i <the prefix of MaChIAto Classifier result> \
-o <output directory> \
-sc <scaffold sequence of sgRNA> \
-ul <unmodified label> \
-ol <knock-out label 1> <knock-out label 2> \
-t <calculation target> \
--ignore_list <list of samples ignored>;
```

##### The case of Single knock-out analysis (*STANDARD)
```bash
python MaChIAto_Analyzer/collect_MaChIAto_data.py \
-i <the prefix of MaChIAto Classifier result> \
-o <output directory> \
-sc <scaffold sequence of sgRNA> \
-ul <unmodified label> \
-ol <knock-out label> \
-t <calculation target> \
--ignore_list <list of samples ignored>;
```

*If you do not enter "-ul <unmodified label>", the process can work. However, some filtering process will be skipped.

### Example of command

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
```

# 3) Simple knock-in analysis : -il PE -t ebmh
```
mkdir ./PE/20200508_MaChIAto2_data_analysis

python3 /Volumes/databank2/temp/MaChIAto/MaChIAto_Analyzer/collect_MaChIAto_data.py \
-i ./PE/for_analyzer \
-o ./PE/20200508_MaChIAto2_data_analysis \
-il PE \
-t ebmh;
```

# 3) Simple knock-in analysis : -il PE -t protospacer
```
mkdir ./PE/20200510_MaChIAto2_data_analysis

python3 /Volumes/databank2/temp/MaChIAto/MaChIAto_Analyzer/collect_MaChIAto_data.py \
-i ./PE/for_analyzer \
-o ./PE/20200510_MaChIAto2_data_analysis \
-il PE \
-t protospacer;
```

50 knock-in | both homology arm
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
```

prime Editing | RT template

```
python3 /Volumes/databank2/temp/MaChIAto/MaChIAto_Analyzer/collect_MaChIAto_data.py \
-i ./PE/MaChIAto_1.8_output \
-o ./PE/20200629_MaChIAto2_data_analysis \
-il C \
-t ermh \
--ignore_list ./PE/20200629_MaChIAto2_data_analysis/ignore_list.csv;
```


## MaChIAto Analyzer

### Parameter list

Summary directory:
The summary directory generated with “collect_MaChIAto_data.py”.

Output prefix:
The directory into which the output directory is saved.

Name of extra data (optional):
The name of feature group that the next argument includes.

Table of extra data (optional):
The pathname of extra data added by the user. The data should be a .csv file.

### Template command

```bash
Rscript MaChIAto_Analyzer/MaChIAtoAnalyzer.R \
<Summary directory> \
<Output prefix> \
<Name of extra data> \
<Table of extra data>;
```

### Example of command

Rscript MaChIAto_Analyzer/MaChIAtoAnalyzer.R \
./20200311_MaChIAto2_data_analysis_mode1 \
./MaChIAtoAnalyzer_beta1.6_output \
InDelphi \
./extra_table/bmh_extra_data_InDelphi.csv \
Accessibility \
./extra_table/bmh_extra_data_accessibility.csv \
Genome_Property \
./extra_table/bmh_extra_data_genomeprop.csv;

# MaChIAto Reviewer

### Parameter list

Summary directory:
The summary directory generated with “collect_MaChIAto_data.py”.

Output prefix:
The directory into which the output directory is saved.

### Template command

Rscript MaChIAto_Reviewer/MaChIAtoReviewer.R \
<Summary directory> \
<Output prefix>;

### Example of command

Rscript /Volumes/databank2/temp/MaChIAto/MaChIAto_Reviewer/MaChIAtoReviewer.R \
./MaChIAto2_v1.7.0_output_200311 \
./MaChIAtoAlignerv14_output_200501_from_MaChIAto2_v1.7.0_output_200311 \
./20200311_MaChIAto2_data_analysis_mode1 \
./MaChIAtoReviewer_beta2.4_output_200610_from_MaChIAto2_v1.7.0_output_200311;

# CaseStudies