# What is knotAnnotSV?
A simple script to create a customizable html file from an AnnotSV output (So far from a Bed file format).

You can customize column order, number and comments with "mouseovering".


# Installation

To download knotAnnotSV, please use git to download the most recent development tree.

Currently, the tree is hosted on github, and can be obtained via:

```bash
$ git clone https://github.com/thomasguignard/knotAnnotSV.git
```

# Requirements 

- Linux OS
- Perl library via cpan : YAML::XS, Sort::Key::Natural


# Input

## AnnotSV output Files


## configuration of output columns

Use a indented yaml file to configure output (use space instead of tab for indentation, tab are not allowed).

Precise the POSITION field you want to display.

Precise which fields you want to include in COMMENTLIST.

Inactivate fields you don't care with a starting '#' or 'POSITION: 0' or by deleting line.

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

An AnnotSV html File is produced to be display on web browser (firefox 80.0 tested so far). It should be in the same directory as Datatables folder. 

See exemple annotSV.html


# Command:
```bash
#Basic output
./knotAnnotSV.pl --annotSVfile example.annotated.tsv --configFile config_cyto.yaml

#Integrate annotSV anking File
./knotAnnotSV.pl --annotSVfile example.annotated.tsv --configFile config_cyto.yaml --annotSVranking example.ranking.tsv
```


--------------------------------------------------------------------------------

**Montpellier Bioinformatique pour le Diagnostic Clinique (MoBiDiC)**

*CHU de Montpellier*

France

![MoBiDiC](logos/logo-mobidic.png)

[Visit our website](https://neuro-2.iurc.montp.inserm.fr/mobidic/)
