# What is knotAnnotSV?

knotAnnotSV is a simple script to create a customizable html file (to be displayed on a web browser) from an [AnnotSV](https://lbgi.fr/AnnotSV) output.

The user can customize the order and the number of the annotation columns as well as the display mode (direct display or by comment) thanks to a configuration file.

Check AnnotSV repo: https://github.com/lgmgeo/AnnotSV

TODO: explain mouseovering / filter / tab (full-split)


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


### A configuration file

Use an indented yaml file (config_cyto.yaml) to configure the AnnotSV output file (use space instead of tab for indentation, tabs are not allowed):

- Precise the POSITION (column ordering) of each field you want to display

- Precise in COMMENTLIST which associated fields you want to display in tooltips (by mousing over)

- Inactivate fields you don't care either with a starting '#' or 'POSITION: 0' or by deleting line.

```bash
---
AnnotSV ID:
    POSITION: 1
    COMMENTLIST:
        - SV length
        - SV type
AnnotSV ranking:
    POSITION: 2
SV chrom:
    POSITION: 0
#SV start:
#    POSITION: 0
```

# Output

An AnnotSV html file is produced to be displayed on a web browser (Firefox 81.0, Chrome 86.0.4240.75, Edge 83.0.478.54, IE 11, tested so far). 

It should be in the same directory as the ```Datatables``` folder. 

To help you get how to make effective use of knotAnnotSV, we have provided an input/output example in the ```example``` folder. 


# Command line:

1. Change to the repo directory, and run the example
```bash

cd /path/to/install/knotAnnotSV

perl ./knotAnnotSV.pl --annotSVfile ./example/example.annotated.tsv --configFile ./config_cyto.yaml
```
2. Display the html output on a web browser

3. Have fun with the exploring!


# USAGE: arguments
perl knotAnnotSV.pl

```
    --configFile <YAML config file for customizing output>

    --annotSVfile <AnnotSV annotated file> 

    --outDir <output directory (default = current dir)> 

    --outPrefix <output file prefix (default = "")> 

    --datatableDir <directory containing datatables files (default = "")>
  
```


--------------------------------------------------------------------------------

**Montpellier Bioinformatique pour le Diagnostic Clinique (MoBiDiC)**

*CHU de Montpellier*

France

![MoBiDiC](https://raw.githubusercontent.com/mobidic/Captain-ACHAB/master/img/logo-mobidic.png)

[Visit our website](https://neuro-2.iurc.montp.inserm.fr/mobidic/)
