# scrna-seq
Single Cell RNA-Seq Workflow

The scripts in this repository perform the basic steps involved in aligning and analyzing single cell RNA-Seq experiemnts.

Single cell RNA-Seq workflows are inherently similar to bulk RNA-Seq workflows. The main difference is each cell in a scRNA-Seq experiment is a unique sequencing library. Also, it should be noted that scRNA-Seq workflows are highly dependent on experimental design and what biological processes are being examined. For example, are you focused of cataloging cell diversity in a heterogenous population, alternative splicing in cell types or detecting differentially expressed genes? Will you be using plate, microfluidic or droplet based methods? Are you using spike-ins, UMIs or neither? Such differences in methodology often result in slightly different analysis workflows.

The steps below cover a semi-generic workflow. I tried to include most common steps, but for clarity the workflow is not comprehensive. As the workflow has numerous steps, I grouped related tasks into stages, which are color coded on the chart below. Relevant scripts and code for each stage is stored in a subdirectory listed to the left of the workflow.

Last, a word about aligners. Current best practices are to use "psuedo-aligners" such as kallisto, salmon or sailfish to quantify transcript abundance prior to making inferences at the gene level. I will focus on kallisto, as that is my personal preference. However, kallisto, salmon or sailfish all perform comprably and are interchangable in most pipelines. Finally, for many biologists learning about bioinformatics the concepts behind pseudo-aligners seem more difficult to grasp than the process of transcript assembly by "true" aligners such as STAR/cufflinks or tophat/cufflinks. To promote the easiest learning of single cell RNA-Seq workflows I have included scripts for both.

![Alt text](https://github.com/ctrhodes/scRNA-seq/blob/master/scRNA-expanded.png?raw=true)

