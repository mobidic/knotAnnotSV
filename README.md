# What is parsAnnotSV?
A simple script to create a customizable html file from an AnnotSV output (So far from a Bed file format).

You can customize column order, number and comments with "mouseovering".


# Installation

To download parsAnnot, please use git to download the most recent development tree.

Currently, the tree is hosted on github, and can be obtained via:

```bash
$ git clone https://github.com/thomasguignard/parsAnnot.git
```

# Requirements 

- Linux OS
- Perl library via cpan : Switch, YAML::XS, Sort::Key::Natural


# Output configuration

Use a indented yaml file to configure output (use splace instead of tab for indentation).

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


# Command:
```bash
./parsAnnotSV.pl --annotSVfile CSG202190.xls.annotated.tsv --configFile config_cyto.yaml
```


--------------------------------------------------------------------------------

**Montpellier Bioinformatique pour le Diagnostic Clinique (MoBiDiC)**

*CHU de Montpellier*

France

![MoBiDiC](logos/logo-mobidic.png)

[Visit our website](https://neuro-2.iurc.montp.inserm.fr/mobidic/)
