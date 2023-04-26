# NeuroVar: A Genetic Expression and Variation data visualization tool for Neurological diseases’ biomarkers

## Table of Contents
1. [Background](#Background)
2. [About NeuroVar](#About NeuroVar)
3. [Implementation and  Operation](#Implementation and  Operation)
4. [Installation](#Installation)
5. [Usage](#Usage)
6. [Case study](#ase study)
7. [Demonstration](#Demonstration)
8. [License](#License)
9. [Contact](#Contact)
10. [Contributors](#Contributors)
<br>

## Background

The expanding availability of large-scale genomic data and the growing interest in uncovering gene-disease associations call for efficient tools to visualize and evaluate gene expression and genetic variation data.  

## About NeuroVar
NeuroVar is a novel tool for visualizing genetic variation (Single nucleotide polymorphisms and insertions/deletions) and gene expression data related to neurological diseases. The tool is available as a desktop application that does not require any computational skills to use.

The current version includes multiple neurological disorders  and diseases with neurological manifestations, including:

★	Epilepsy   ★	amyotrophic lateral sclerosis  ★	autism spectrum disorder   ★	intellectual disability  ★	brain malformation syndrome  ★	Syndromic Disorders   ★	Cerebral Palsy   ★	RASopathy    ★	Amnioacidopathy   ★	Peroxisomal disorders   ★	Hereditary cancer   ★	Mitochondrial disease   ★	Retina related disorders  ★	Hearing Loss Neuropathy   ★	neurodevelopmental disorder   ★	fragile X syndrome  ★	frontotemporal dementia   ★	glycine encephalopathy   ★	Griscelli syndrome   ★	Joubert syndrome    ★	Microcephaly   ★	neurodevelopmental disorder   ★	neurofibromatosis   ★	Parkinson disease   ★	PHARC syndrome    ★	tuberous sclerosis  ★	neuroblastoma.


## Implementation and  Operation

NeuroVar is :

- Desktop application
- Do not require any computational skills to use 
- Platform independent
- All requirement are installed automatically with the tool
- Input :  It requires gene expression CSV files and genetic variants VCF files as inputs.

## Installation

The tool could be installed from this website : https://sites.google.com/view/neurovar 

## Usage

The tool integrates 3 tabs, each for a different functionality,  as explained bellow:

1. First the user select the patient’s disease, next, a list of specific subtypes of the disease with neurological manifestation is provided. More information about each gene is provided including mode of inheritance, description, type, and transcripts. A links for the official online report validating the Gene's association with the disease is also provided. 

![image](https://user-images.githubusercontent.com/73958439/232723944-8e5e658e-bbe5-40e7-92d7-f855ae0400aa.png)


2.Visualize the biomarkers expression profile: 
After importing a CSV file and identifying the key columns, the log2FC value and p-value to define the differential expression profile are requested. The genes’ expression profiles are summarized in a table and represented in a volcano plot.

![image](https://user-images.githubusercontent.com/73958439/232724064-e2803d44-4381-408e-b0b7-0a9553b8a16b.png)


3. Visualize SNPs and Indels data:
The user is requested to define the path to the directory containing the VCF files. The user needs to define the variants type (SNPs or Indels). The VCF files are processed and annotated, and then the variants in the biomarkers are filtered and resumed in a table comparing the reference genome, the controls group and the patients group.

![image](https://user-images.githubusercontent.com/73958439/232724156-3bd91417-89ec-4d1e-a56e-953836a0256b.png)

## Case study

To validate the pipeline, a case study was performed on the public dataset SRP149638 available on the SRA database. The file’s preprocessing, genetic expression analysis, and variant calling were performed using the Exvar R package. 
The dataset corresponds to RNA sequencing data (Expression profiling by high throughput sequencing) from the peripheral blood mononuclear cells (PBMC) from healthy donors and Amyotrophic Lateral Sclerosis (ALS) patients.
The expression analysis shows that only one ALS biomarker is differentially expressed which is TUBA4A gene, a biomarker of ALS type 22. TUBA4A is a protein-coding gene inherited in autosomal dominant mode.

## Demonstration

Demonstration data are provided : https://github.com/omicscodeathon/neurovar/blob/main/demonstration_data.rar

This data will guide you to organize your data correctly and could be used to test the tool.


## License

Artistic license 2.0 

## Contact 

For issues reporting, please send an email to : bioinformaticstool@gmail.com

## Contributors

Hiba Ben Aribi, UTM, Tunisia  (Developer - Writer).

Najla Abassi, UTM, Tunisia  (Developer - Writer).

Olaitan I. Awe, ASBCB, South Africa.
