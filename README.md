<p align="center">
    <img src="https://github.com/mobidic/knotAnnotSV/blob/master/images/logoKNOT.png" width="600"/>
</p>

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

### Color codes

The full lines are highlighted depending on their SV type. Split lines have a white background. 

<p align="center">
    <img src="https://github.com/mobidic/knotAnnotSV/blob/master/images/SVcolor.png" width="400"/>
</p>


Two different palettes of colors can be used to visualized for each gene its corresponding GnomAD defined LOEUF_bin value (user defined, see the USAGE/OPTIONS section), smaller values reflect loss of function intolerance:

- A red-to-green palette (default): The red color means the lower LOEUF_bin value (0). The green color means the higher LOEUF_bin value (9).
- A sequential palette of red: Darker colors mean lower LOEUF_bin values (e.g 0).

<p align="center">
    <img src="https://github.com/mobidic/knotAnnotSV/blob/master/images/LOEUFbin_ranges.png" width="600"/>
</p>

Depending on the overlapped gene part, the « Gene name » box is fully colored or not:
- Gene totally overlapped: fully colored box
<p align="center">
    <img src="https://github.com/mobidic/knotAnnotSV/blob/master/images/gene_fulloverlap.png" width="200"/>
</p>

- Gene partially overlapped: half colored box depending on whether the 3’ of the gene is overlapped or the 5’ of the gene is overlapped

<p align="center">
    <img src="https://github.com/mobidic/knotAnnotSV/blob/master/images/gene_overlap.png" width="200"/>
</p>


### Main features

Three display modes are available thanks to click action :
- “Compact”: only the “full” AnnotSV lines are displayed (COMPACT button)
- “Expanded”: the “full” and “split” AnnotSV lines are displayed (EXPANDED button)
- “Single SV focus”: only the “full” and “split” AnnotSV lines of a single SV are displayed (double-click on a full line in compact mode)

Additionnal data, links and sorting are available:
- Hover annotation with mouse to display complementary information (tooltips)
- Click on the « AnnotSV ID » to open the SV coordinates in the UCSC Genome browser (the SV region is automatically highlighted in blue and zoomed out by 1.5x)
- Click on the blue hyperlinks to access directly to the corresponding public database (OMIM, genecards)
- By default, the annotation lines are sorted according to these priorization rules: ACMG class > Exomiser Score > OMIM morbid > LOEUF bin (this last one is applied on split lines only)


Column headers have searching and sorting features:
- Searching words or extracting matching records
- Filtering a range of numerical values (with the “>”, “>=”, “<” and “<=” symbols) e.g. to select frequencies smaller than 1%, type “<0.01”
- Filtering out some numerical values (“!= value”) 
- Filtering out some words (“! word”)
- Search is also performed in the tooltips of the OMIM, Pathogenic and Gene Name fields (e.g. you can to match phenotype in the tootip). 
- Column header is clickable to sort values 
- At each stage of the analysis process, all the set-up filters are locally stored
- Reset original sorting by clicking on the header of the first column (AnnotSV ID)



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
