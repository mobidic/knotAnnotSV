![KnotAnnotSV](https://github.com/mobidic/knotAnnotSV/blob/master/images/logoKNOT.png)

# What is knotAnnotSV?

knotAnnotSV is a simple tool to create a customizable html file (to be displayed on a web browser) from an [AnnotSV](https://lbgi.fr/AnnotSV) output.

The user can customize the order and the number of the annotation columns as well as the visualization mode (direct display or by comment) thanks to a configuration file. He can then visualize, filter and analyze the annotation data thanks to different user friendly available functions (search/filtering box, tooltip, links to public databases, color coded information...).

The available interface is well detailed in the README.knotAnnotSV_latest.pdf file.

At each stage of the analysis process, all the set-up filters are locally stored. The html file can thus be closed at any time by biologists and then re-opened on the same computer to continue the analysis.

Check AnnotSV repo: https://github.com/lgmgeo/AnnotSV


# Installation

Please use git to download the most recent development tree.

The sources can be cloned to any directory:
```bash
cd /path/to/install/
$ git clone https://github.com/mobidic/knotAnnotSV.git
```

# Requirements 

- Linux OS

- Perl library via cpan : YAML::XS, Sort::Key::Natural


# Input

### An AnnotSV output file

Use a classical annotated tsv file from AnnotSV.


### A configuration file

Use an indented yaml file (see config_AnnotSV.yaml for a good start-up) to configure the AnnotSV output file (use space instead of tab for indentation, tabs are not allowed):

- Precise the POSITION (column ordering) of each field you want to display (AnnotSV_ID must be in position 1 and Annotation_mode must be present for good computation) 

- Precise in COMMENTLIST which associated fields you want to display in tooltips (by mouse hovering)

- Precise in RENAME which the column name to display in the output

- Precise in HEADERTIPS some information about this column that you want to display in the header tooltips (by mouse hovering)

- Inactivate fields you don't care either with a starting '#' or 'POSITION: 0' or by deleting line.

```bash
---
AnnotSV_ID:
    POSITION: 1 #mandatory at this position
ACMG_class:
    POSITION: 2
    RENAME: ACMG classification                #change the column name in the output 
    HEADERTIPS: ACMG scoring implementation    #some details
    COMMENTLIST:                               #All the fields I want to knot to the main field 
        - SV length
        - AnnotSV_ranking_score
        - AnnotSV_ranking_criteria
SV_type:
    POSITION: 5
    RENAME: SV type
Annotation_mode:
    POSITION: 4
Gene_name:
    POSITION: 5
    RENAME: Gene Symbol
    HEADERTIPS: Some additionnal explanation about this column
    COMMENTLIST:
        - Gnomad_pLI
        - OMIM_phenotype
SV_chrom:
    POSITION: 0
#SV_start:
#    POSITION: 0
```

# Output

An AnnotSV html file is produced and ready to be displayed on a web browser (Firefox 81.0, Chrome 86.0.4240.75, Edge 83.0.478.54, IE 11, tested so far). 


# USAGE
```
cd /path/to/install/knotAnnotSV

perl ./knotAnnotSV.pl

    --configFile <YAML config file for customizing output>

    --annotSVfile <AnnotSV annotated file> 

    --outDir <output directory (default = current dir)> 

    --outPrefix <output file prefix (default = "")> 
    
    --genomeBuild <Genome Assembly reference (default = hg19)>
    
    --LOEUFcolorRange <Number to define which color gradient to use for LOEUF bin: 1 (red-to-green), 2 (red-shades-only) (default = 1)>

    --datatableDir <Local Path to DataTables directory containing css and js files (default = \"\", requires web connection)>
  
```

# Test knotAnnotSV

To help you get how to make effective use of knotAnnotSV, we have provided an input/output example in the ```example``` folder. 

1. Change to the repo directory, and run the example
```bash

cd /path/to/install/knotAnnotSV

perl ./knotAnnotSV.pl --annotSVfile ./example/example.annotated.tsv --configFile ./config_AnnotSV.yaml --outDir ./example
```
2. Display the html output on a web browser

3. Have fun with the exploring!



--------------------------------------------------------------------------------

**Montpellier Bioinformatique pour le Diagnostic Clinique (MoBiDiC)**

*CHU de Montpellier*

France

![MoBiDiC](https://raw.githubusercontent.com/mobidic/Captain-ACHAB/master/img/logo-mobidic.png)

[Visit our website](https://neuro-2.iurc.montp.inserm.fr/mobidic/)
