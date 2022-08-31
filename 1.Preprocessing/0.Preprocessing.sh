perl 2.Cutadapt.pl
gunzip *.gz
pypy Statistics.py ./
perl 1.Statistic_Split.pl
/home/programs/FastQC/fastqc *
mkdir FastQC
mv *fastqc* FastQC/
cd FastQC/
multiqc *
cd ..
gzip *
cd Cutadapt
perl 3.Trimmomatic.pl
gunzip *.gz
pypy Statistics.py ./
perl 1.Statistic_Split.pl
gzip *.fq
cd P/
gunzip *.gz
pypy Statistics.py ./
perl 1.Statistic_Split.pl
/home/programs/FastQC/fastqc *
mkdir FastQC
mv *fastqc* FastQC/
cd FastQC/
multiqc *
cd ..
gzip *.fq

