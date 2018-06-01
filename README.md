# mapandvischainonfiji

This program needs a fastq file as an input and will output quality files, trimed fastq files, sam files, bam files, sorted bam files and read normalized tdfs. 

to download you will run:

git clone https://github.com/maallen3/mapandvischainonfiji.git

To run:

cd mapandvischainonfiji/

#if you are on fiji 
#module load python/2.7.14/pandas

python create_scipts_to_map_on_fiji_argparse.py infile.fastq outdir genome email

For example:

python create_scipts_to_map_on_fiji_argparse.py /scratch/Users/allenma/171025_NB501447_0179_fastq/Allen_Dowell-371/MA_DMSO_Groseq_S7_R1_001.fastq /Users/allenma/proseq/ hg19 allenma@colorado.edu

The output will be these directories in your outdir. 
bams  bedgraphs  cutadapt  e_and_o  qsubscripts  qual  sams  sortedbams  tdfs

By default this script will create 5 slurm scripts and submit them all to the queue.
infilerootname._mom_run.slurm #this script runs all the other scripts
infilerootname._trim.slurm #this script trims the fastq files of adapters
infilerootname._map_.slurm #this script maps with bowtie2
infilerootname._samtobam.slurm #this script changes sams into bam files
infilerootname._bamtobedgraph.slurm #this script creates normalized tdf files



For more options type
python create_scipts_to_map_on_fiji_argparse.py --help

You can turn off the auto submit like this.
python create_scipts_to_map_on_fiji_argparse.py infile outdir genome email --Turnoff_submit

If you turn on --flipreads one other script will be created and run. 



