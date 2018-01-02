# scrna-seq
Single Cell RNA-Seq Workflow

The scripts in this repository perform the basic steps involved in aligning and analyzing single cell RNA-Seq experiemnts.

Single cell RNA-Seq workflows are inherently similar to bulk RNA-Seq workflows. The main difference is each cell in a scRNA-Seq experiment is a unique sequencing library. Also, it should be noted that scRNA-Seq workflows are highly dependent on experimental design and what biological processes are being examined. For example, will you be using plate, microfluidic or droplet based methods? Are you using UMIs? Are you interested in alternative splicing versus differential gene expression? Such differences in methodology often result in different analysis workflows

The steps below cover a semi-generic workflow. The workflow has numerous steps, so the scripts and code for each section is stored in one of 3 subdirectories listed to the left of the workflow.

![Alt text](https://github.com/ctrhodes/scrna-seq/blob/master/scRNSseq_workflow.png?raw=true)

