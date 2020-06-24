tugHall version 2.1
====================

[![License](https://img.shields.io/badge/License-GPLv3-orange.svg)](https://github.com/tugHall/Clone-based/blob/master/LICENSE)


**tugHall** _(**tu**mor **g**ene-**Hall**mark)_ is cancer-cell evolution model simulator, wherein gene mutations are linked to the tumor cell behaviors, which are influenced by the hallmarks of cancer.

This is an _**R**_-based script to simulate the cancer cell evolution in the framework of the model proposed by _**Prof. Mamoru Kato**_,
_Head of Bioinformatics Department, Research Institute, National Cancer Center, Tokyo, JAPAN_.

Authors and contributor list:
---
_**Iurii Nagornov**_

_**Mamoru Kato**_

_Department of Bioinformatics, Research Institute, National Cancer Center Japan, Tokyo, Japan_

All questions and requests can be sent to inagonov@ncc.go.jp

Project source can be downloaded from websites  
--- 
https://github.com/tugHall/Clone-based  -  the developing resource

Short description
---
The wide availability of recent cancer genomic data requires a coherent model that can sort out the relevant findings to systematically explain the clonal evolution and resultant intra-tumor heterogeneity (ITH). Here, we present a new mathematical model designed to computationally simulate the evolution of cancer cells. The model connects well-known cancer hallmarks with the specific mutational states of tumor-related genes. The cell behavior phenotypes are stochastically determined and the hallmarks interfere probabilistically with the phenotypic probabilities. In turn, the hallmark variables depend on the mutational states of tumor-related genes. Thus, it is expected our software can be used to deepen our understanding of cancer-cell evolution and generation of ITH.

How to use **tugHall** and how to analyze data, kindly see user-guides in **Documentation** folder in prefered format *Rmd, pdf, html*.

General changes
---

This version is based on the clone consideration instead cell-based version 1.1.  
Each clone has one or more cells, that allows to accelerate the calculations when number of clones is much less than number of cells.
Definition of clone: the clone is set of cell with same set of genes, which have same mutated / not mutated sites in genes.

Changes in comparison with tugHall v.2.0
---

This version **tugHall 2.1**  has additional changes in the part of binominal distribution to
calculate binominal distribution for big numbers like n=10^14 or more with probability p=0.1 or something. 
For so huge numbers we have to use approximation with normal distribution, because **rbinom()** function has a limitation for numbers **n**x**p**>10^12. 


Content of package
---

* **tugHall_2.1.R** is a R-script of simulation, which uses scripts in **/Code/** folder.
* **/Code/** is the directory with scripts to simulate and to analyze data. 
* **/Input/** and **/Output/** are the directories for input and output data during a simulation. 
* **/Documentation/** is a directory with documentation to use tugHall software and to analyze data of simulation.   
* **/Figures/** is the folder with figures of plots of last simulation, which are used in user-guides also. 
* **/TESTs/** is the directory with tests and it's results. Description of tests is in the folder. 


