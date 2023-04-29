# NeuroVar: A Genetic Expression and Variation data visualization tool for Neurological diseases’ biomarkers

## Table of Contents
1. [Background](#Background)
2. [About NeuroVar](#About-NeuroVar)
3. [Implementation and  Operation](#Implementation-and-Operation)
4. [Installation](#Installation)
5. [Usage](#Usage)
6. [Case study](#Case-study)
7. [Demonstration](#Demonstration)
8. [License](#License)
9. [Contact](#Contact)
10. [Contributors](#Contributors)
<br>

## Background

The expanding availability of large-scale genomic data and the growing interest in uncovering gene-disease associations call for efficient tools to visualize and evaluate gene expression and genetic variation data.  

## About NeuroVar
NeuroVar is a novel tool for visualizing genetic variation (Single nucleotide polymorphisms and insertions/deletions) and gene expression data related to neurological diseases. The tool is available as a desktop application that does not require any computational skills to use.

The current version includes :

* 11 Neurological Syndromes are integrated, including:

              ★ Epilepsy 
              ★ Amyotrophic lateral sclerosis 
              ★ Intellectual disability and Autism spectrum disorder 
              ★ Brain malformation syndrome 
              ★ Syndromic Disorders 
              ★ Cerebral Palsy 
              ★ RASopathy 
              ★ Amnioacidopathy  
              ★ CranioFacial Malformations 
              ★ Parkinson disease 
              ★ PHARC syndrome 

* and 7 Non-Neurological diseases with neurological manifestations are integrated, including:

              ★ Peroxisomal disorders  
              ★ Hereditary cancer 
              ★ Mitochondrial disease 
              ★ Retina related disorders 
              ★ General Gene Curation
              ★ Hearing Loss 
              ★ Fatty Acid Oxidation Disordes





## Implementation and  Operation

NeuroVar:

- is a desktop application
- does not require any computational skills to use
- platform-independent
- all requirements are installed automatically with the tool
- input:  NeuroVar requires gene expression CSV files and genetic variants VCF files as inputs.

## Installation

**1 -Download the application**: Download the application file from the [app's website](https://sites.google.com/view/neurovar)                                                           
**2- Extract the application**: Extract the contents of the downloaded file to a directory of your choice                                                                     
**3- Install dependencies**: Install the latest version of the following requirements if your computer does not have them:         

- [Python](https://www.python.org/downloads/windows/)
- [Node.JS](https://nodejs.org/en/download) and npm
- Electron `npm install electron@v23.1.2`
- [R](https://cran.r-project.org/bin/windows/base/)                 


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

Demonstration data are provided: https://github.com/omicscodeathon/neurovar/blob/main/demonstration_data.rar

This data will guide you to organize your data correctly and could be used to test the tool.


## License

Artistic license 2.0

## Contact

For issues reporting, please send an email to : bioinformaticstool@gmail.com

## Contributors

Hiba Ben Aribi, UTM, Tunisia  (Developer - Writer).

Najla Abassi, UTM, Tunisia  (Developer - Writer).

Olaitan I. Awe, African Society for Bioinformatics and Computational Biology (ASBCB), South Africa.
