5/16/25:

Hello! Here is all the information on running scSPRITE



GitHub: https://github.com/caltech-bioinformatics-resource-center/Guttman_Ismagilov_Labs

Paper: Arrastia, M.V., Jachowicz, J.W., Ollikainen, N. et al. Single-cell measurement of higher-order 3D genome organization with scSPRITE. Nat Biotechnol 40, 64–73 (2022). https://doi.org/10.1038/s41587-021-00998-1

https://www.nature.com/articles/s41587-021-00998-1


Data location: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154353

Raw mouse sample 1: https://www.ncbi.nlm.nih.gov/sra/SRX8722756[accn]
Raw mouse sample 2: https://www.ncbi.nlm.nih.gov/sra/SRX8722757[accn]

The raw data used to test functionality of our pipeline was found at the SRA line from data availability.

If you are in the scSPRITE paper go to 
    data availability >
    GSE154353 > 
    SRA	SRP271653 >
    GSM4669509: mESC-scSPRITE-50%-02; Mus musculus; OTHER or GSM4669508: mESC-scSPRITE-50%-01; Mus musculus; OTHER >
    SRR12212044 or SRR12212045
    
    From here you should see tabs for FASTQ download if you ever need to download data.
    
    
Data storage and location: raw fastq files for analysis are stored at raw_fastq. I have kept old or supporting files for these in the folder fastqDir. Paola has said that scSPRITE can only run one sample at a time, so I would leave all other files in fastqDir while you are running analysis on one sample in raw_fastq. There should be two files in raw_fastq: R1 and R2. When it is running, it might create a 3rd file in raw_fastq acknowledging the name of the sample. 

The naming format of your sample reads should match the one I have used, the format set by the creators, or else it will not run. 
i.e. SRR12212045_1.fq.gz
You are more than welcome to edit the code to take different naming formats, but I found it easier to stick to the syntax they wrote the code around. Make sure it is zipped .gz too!


RUNNING THE PIPELINE: 
1) Navigate to /blue/giustirodriguezp/PROJECTS/Sprite/Guttman_Ismagilov_Labs/scSPRITE and open terminal. If you are reading this you are probably already in the right file location.

2) Load snakemake by typing "ml snakemake/9.4.0". Load java by typing "ml java". Load samtools by typing "ml samtools".

3) Do a dry run by typing "snakemake -n". If there are no major errors, you are good to go.

4) First step snakemake run: snakemake --cores 15

5) Second step snakemake run: snakemake -s Snakefile_contact_heatmap --cores 15


People to contact who can help:

Joanna Jachowicz, scSPRITE co-author: joanna.jachowicz@imba.oeaw.ac.at

Andre Chanderbali, HiPerGator Support Ticket: achander@ufl.edu


Notes for why scSPRITE needs edits:

4/10/25:
the code does not make its own new_fastq folder. i made one manually, but mkdir command in the snakefile around line 89 is ideal for future use.

4/18/25:
Only one sample runs at a time. I recommend moving you raw data and other materials to fastqDir until you need it. Move it to raw_fastq when you're ready.

5/19/25:
Andre says theres a mystery python package, possibly custom, for run "Snakefile_contact_heatmap". It is called "contact" and I still haven't gotten to the second step so I'm not really sure what it does.