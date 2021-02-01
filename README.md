# **CroTrait:**<br>A portable tool for *in silico* species identification, serotyping and multilocus sequence typing of *Cronobacter* genus
![icon](assist/icon.jpg)

## Introduction

This is the homepage of **CroTrait**, a software package that performs efficient inference of species, O serotypes and sequence types (STs) of *Cronobacter*. **CroTrait** was developed by *Lu Wang*. **CroTrait** is called from a directory containing one or multiple genome sequences and each genome sequence should be saved to a separate file in FASTA format. Input files can either be complete genomes or draft genomes. It is able to analyze hundreds of genomes data in a matter of hours on a ordinary PC. <br><br>
Moreover, post data analysis and visualization module embedded in **CroTrait** further assist the user in checking and analyzing the data.

## Citation
if you use this software package please cite:<br>
*Lu Wang, Wenxuan Zhu, Gege Lu, et al*. A portable tool for *in silico* species identification, serotyping and multilocus sequence typing of *Cronobacter* genus.
***Journal of Clinical Microbiology*** 2021 Feb (going to be published)

## Environment set up
**CroTrait** is a program written in **python** and the external software **BLAST+** and **MEGA** need to be installed and configured locally:<br>

>**[python](https://www.python.org/)** (version 3.8.0) <br>

>> dependencies <br>
>> **biopython** (version 1.78 or higher) <br>
>> **numpy** <br>
>> **matplotlib** <br>
>> **pandas** <br>

>**[MEGA](https://www.megasoftware.net/)** (version X (64-bit))<br>

>**[BLAST+](https://blast.ncbi.nlm.nih.gov/)** (version 2.9.0)<br>

## User guide
>**1. identify species, O serotypes and STs**<br>
>>**1.1. Assemblies with unknown species**<br>
>>$***CroTrait -d directory [-p prefix]*** <br>
>>-p: represent the prefix of result file <br><br>
>>**1.2. Assemblies with known species (one of the seven species of *Cronobacter*)**<br>
>>$***CroTrait -d directory -s species [-p prefix]*** <br>

>**2. generate O antigen clusters (O-AGCs) pattern**<br>
>>$***CroTrait -t 2 -d genomes -s species***<br>
>>  example figure<br><br>
![icon](assist/icon1.jpg)

>**3. extract O-AGCs sequences in batch**<br>
>>$***CroTrait -t 3 -d directory [-p prefix]*** <br>

>**4. post statistics analysis**<br>
>>$***CroTrait -t 4 -r result_table*** <br>
>> -r: the result created by "1", namely table with identified species, O serotypes and STs.<br>
>> after executing this command, 6 table will generated according to the species, O serotypes and ST.<br>

>**5. post visulization analysis**<br>
>>  example figure<br><br>
![icon](assist/icon2.jpg)

## License
**CroTrait** is a free software package, licensed udner **GPLv3**.

## Feedback/Issues
If you need assistance using CroTrait, you can get in touch by emailing Lu Wang (*wlubio@sina.com*), or by asking on Issues page.







