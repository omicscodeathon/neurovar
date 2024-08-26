# NeuroVar desktop Application User Guide

## Implementation and  Operation

    - does not require any computational skills to use
    - platform: Windows
    - all dependencies are installed automatically with the tool
    - Internet requirement: NO
    - Third-party tools requirement: NO
    - Pre-requirements: Python  installation
    
## Installation


**Step 1 :** Install dependencies

Install the latest version of Python if your computer does not have them:  [Python](https://www.python.org/downloads/windows/)    

**Step 2 :** Download the application
   
Download the application file from the [here](https://drive.google.com/file/d/1UiLqW-ZER6ysMHjTSW033AQwhXlqtuFy/view?usp=sharing)  

**Step 3 :** Extract the application 

Extract the contents of the downloaded file to a directory of your choice  

**Step 4 :** Open The Application

Unzip the folder and open the NeuroVar.exe application file

**Step 5:**

The tool integrates 3 tabs, each for a different functionality.

In the First tab, select the patient’s disease, next, a list of specific subtypes of the disease with neurological manifestation is provided. More information about each gene is provided including mode of inheritance, description, type, and transcripts. A Link to the official online report validating Gene's association with the disease is also provided.


![image](https://github.com/omicscodeathon/neurovar/assets/73958439/b22e51f8-84a2-47b1-9a6b-70d0284d3867)


**Step 6:** Visualize the biomarkers expression profile:

   
In the second tab, import the gene expression data file (as a CSV file) and identify the key columns, the log2FC value and p-value to define the differential expression profile. The genes’ expression profiles are summarized in a table and represented in a volcano plot.


![image](https://github.com/omicscodeathon/neurovar/assets/73958439/c6e976b4-3430-4f54-92d5-7401be22dc76)



**Step 7:** Visualize SNPs and Indels data:

In the third tab,  define the path to the directory containing the genetic variants data files ( as VCF files). Then define the variants type (SNPs or Indels). The VCF files are processed and annotated, and then the variants in the biomarkers are filtered and resumed in a table comparing the reference genome, the control group and the patient group. The data could be downloaded as a  csv file.



#### Demonstration Data 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13375590.svg)](https://doi.org/10.5281/zenodo.13375590)

Demonstration data are provided [here](https://github.com/omicscodeathon/neurovar/blob/main/demonstration_data.rar)

This data will guide you to organize your data correctly and could be used to test the tool.

## License

Artistic license 2.0

## License

Artistic license 2.0
