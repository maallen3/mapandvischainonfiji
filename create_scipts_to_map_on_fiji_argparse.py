import argparse
import sys
import os
import glob 
#this is built for mapping single end nascent data!!! This program expects a infile and a outdir. The outfile must end in a slash.

#What default options do you want to use
bowtieoptions = "--very-sensitive"
processers="1" #this will move to 32 if bowtie is being used
memamount="10"
memtype="gb" #gb, mb, kb
igvdir="/opt/igvtools/2.3.75/"

#-------- these are variables you probably don't need to check
genomesdir="/scratch/Shares/public/genomes/"
def createbowtieindexdic():
	b = {}
	bowtiefiles = [infile for infile in glob.glob(os.path.join(genomesdir, '*/*/*/Sequence/Bowtie2Index/genome*'))]
	for bowtiefile in bowtiefiles:
		dname = os.path.dirname(bowtiefile)
		rootname=bowtiefile.split("/")[7]
		b[rootname]=dname+"/genome"	
	return b
def creategenomefadic():
	g = {}
	genomefastafiles = [infile for infile in glob.glob(os.path.join(genomesdir, '*/*/*/Sequence/WholeGenomeFasta/genome.fa'))]
	for genomefa in genomefastafiles:
		rootname=genomefa.split("/")[7]
		g[rootname] = genomefa
	return g


genomefadic = creategenomefadic()
bowtieindexdic = createbowtieindexdic()
sbatchmem = memamount+memtype
bedtoolsgenomedir="/scratch/Shares/public/genomes/bedtools_genome_files/"



#--------------all scripts below this line



def flipfile(passfile, outdir, wfname):
	wf = open(wfname, "a")
	rootname = passfile.split("/")[-1]
        rootname = rootname.strip(".fastq")
	outfastqdir = outdir+"flipped/"
        ensure_dir(outfastqdir)
	outfastq = outfastqdir+rootname+".flip.fastq"
	wf.write("date\n")
	wf.write("\n")
	wf.write("module load fastx-toolkit/0.0.13\n")
	wf.write("/opt/fastx-toolkit/0.0.13/bin/fastx_reverse_complement -Q33 -i "+passfile+" -o "+outfastq+" \n")	 
	wf.write("echo flipped\n")
	wf.write("date\n")
	wf.close()
	return outfastq

def checkqualityfile(passfile, outdir, wfname):
	wf = open(wfname, "a")
        wf.write("date\n")
	wf.write("module load fastqc/0.11.5\n")
	qualoutdir = outdir+"qual/"
	ensure_dir(qualoutdir)
	wf.write("fastqc "+passfile+" -o "+qualoutdir+"\n")
	wf.write("echo qual\n")
	wf.write("date\n")
        wf.close()

def trimgalorefile(passfile, outdir, wfname):
	rootname = passfile.split("/")[-1]
	rootname = rootname.strip(".fastq")
	cutoutdir = outdir+"cutadapt/"
	ensure_dir(cutoutdir)
	#MA_DMSO_Groseq_S7_R1_001.flip_trimmed.fq
	#MA_DMSO_Groseq_S7_R1_001.flip.fastq  
	trimmedfastq = cutoutdir+rootname+"_trimmed.fq"
	trimmedfastq2 = cutoutdir+rootname+".trimmed.fastq"
	wf = open(wfname, "a")
        wf.write("date\n")
        wf.write("\n")
        wf.write("module load python/2.7.14/cutadapt/1.12\n")
	wf.write("module load trim_galore/0.4.3\n")
        wf.write("echo $PATH\n")
        wf.write("echo $PYTHONPATH\n")
        wf.write("/opt/trim_galore/0.4.3/trim_galore --path_to_cutadapt /opt/cutadapt/python-2.7.14/1.12/bin/cutadapt -o "+cutoutdir+" "+passfile+"\n")
	wf.write("mv "+trimmedfastq+" "+trimmedfastq2+"\n")
	wf.write("echo trimmed\n")
        wf.write("date\n\n\n")
	wf.close()
        return trimmedfastq2


def bowtie2file(passfile, outdir, wfname, genome,processers):
        wf = open(wfname, "a")
        wf.write("date\n")
	wf.write("module load samtools/1.3.1\n")
	wf.write("module load bowtie/2.2.9\n\n")
	samdir = outdir+"sams/"
	ensure_dir(samdir)
	rootname = passfile.split("/")[-1]
	rootname = rootname.strip(".fastq")
	samfile = samdir+rootname+".sam"
	stderrfile = samdir+rootname+".stderr"
	bowtieindex = bowtieindexdic[genome]
	wf.write("bowtie2 -p "+processers+" "+bowtieoptions+" -x "+bowtieindex+" -U "+passfile+" >"+samfile+" 2>"+stderrfile+" \n")
	wf.write("echo mapped\n")
	wf.write("date\n")
	wf.close()
        return samfile


def samtobamfile(passfile, outdir, wfname):
        wf = open(wfname, "a")
        wf.write("date\n")
	wf.write("wc -l "+passfile+" > "+passfile+".wc\n")
	rootname=passfile.split("/")[-1]
	rootname=rootname.strip(".sam")
	bamoutdir = outdir+"bams/"
	ensure_dir(bamoutdir)
	bamfile = bamoutdir+rootname+".bam"
	wf.write("module load samtools/1.3.1\n")
        wf.write("samtools view -S -b -o "+bamfile+" "+passfile+" 2>"+bamfile+".err\n")
	wf.write("samtools flagstat "+bamfile+" > "+bamfile+".flagstat 2>"+bamfile+".flagstat.err\n")
	wf.write("echo bam\n")
	wf.write("date\n")
        wf.close()
        return bamfile


def bamtosortedbamfile(bamfile, outdir, wfname, samtoolsmem):
	wf = open(wfname, "a")
        wf.write("date\n")
	sortedbamfile = bamfile.strip(".bam")+".sorted.bam"
	sortedbamdir = outdir+"sortedbams/"
	ensure_dir(sortedbamdir)
	sortedbamfile = sortedbamdir+sortedbamfile.split("/")[-1]
	wf.write("samtools sort -m "+samtoolsmem+" "+bamfile+" >"+sortedbamfile+"\n")
	wf.write("samtools flagstat "+sortedbamfile+" >"+sortedbamfile+".flagstat 2>"+sortedbamfile+".flagstat.err\n")
	wf.write("echo sorted.bam\n")
        wf.write("date\n")
	return sortedbamfile
                        


def sortedbamtobai(sortedbam, wfname):
	wf = open(wfname, "a")
        wf.write("date\n")
	wf.write("samtools index "+sortedbam+"\n")
        wf.write("echo indexed.bam\n")
        wf.write("date\n")

def sortedbamtobedgraphfile(sortedbam, outdir, wfname, genome):
	wf = open(wfname, "a")
        wf.write("date\n")
        wf.write("module load bedtools/2.25.0\n")
	wf.write("echo bedgraph\n")
	rootname = sortedbam.split("/")[-1]
	rootname = rootname.strip(".sorted.bam")
	Bedgraphoutdir = outdir+"/bedgraphs/"
	ensure_dir(Bedgraphoutdir)
	bedtoolsgenomefile=bedtoolsgenomedir+genome+".chrom.sizes.genome"
	wf.write("genomeCoverageBed -bg -strand + -ibam "+sortedbam+" -g "+bedtoolsgenomefile+" > "+Bedgraphoutdir+rootname+".pos.BedGraph\n")
	wf.write("genomeCoverageBed -bg -strand - -ibam "+sortedbam+" -g "+bedtoolsgenomefile+" | awk -F '\t' -v OFS='\t' '{ $4 = - $4 ; print $0 }' > "+Bedgraphoutdir+rootname+".neg.BedGraph\n")
	wf.write("cat "+Bedgraphoutdir+rootname+".pos.BedGraph "+Bedgraphoutdir+rootname+".neg.BedGraph > "+Bedgraphoutdir+rootname+".unsorted.BedGraph\n")
	wf.write("sortBed -i "+Bedgraphoutdir+rootname+".unsorted.BedGraph >"+Bedgraphoutdir+rootname+".BedGraph\n")
        wf.write("date\n")
	Bedgraphfile = Bedgraphoutdir+rootname+".BedGraph"
	return Bedgraphfile

                        



def readcountcorrectbedgraphfile(Bedgraphfile, sortedbam, outdir, wfname):
	wf = open(wfname, "a")
        wf.write("date\n")
        wf.write("echo readcountcorrectedbedgraph\n")
	wf.write("module load python/2.7.14\n")
	#need to copy readcountcorrectBG.py to out dir
	rootname = Bedgraphfile.strip(".BedGraph")
	flagstatfile = sortedbam+".flagstat"
	readcountcorrBedgraphfile = rootname +".mp.BedGraph" 
	wf.write("python "+outdir+"qsubscripts/readcountcorrectBG.py "+Bedgraphfile+" "+flagstatfile+" "+readcountcorrBedgraphfile+"\n")	
        wf.write("date\n")
	return readcountcorrBedgraphfile
                        

def igvcreatefile(readcountcorrBedgraphfile, outdir, wfname, genome):
	wf = open(wfname, "a")
        wf.write("date\n")
	tdffiledir=outdir+"tdfs/"
	ensure_dir(tdffiledir)
	rootname = readcountcorrBedgraphfile.split("/")[-1]
	rootname = rootname.strip(".mp.BedGraph")
	tdffile = tdffiledir+rootname+".tdf"
	#igvgenomefile = igvdir+"/genomes/"+genome+".chrom.sizes" 
	#igvgenomefile=bedtoolsgenomedir+genome+".chrom.sizes.genome"
	igvgenomefile = genomefadic[genome]
	wf.write(igvdir+"/igvtools toTDF "+readcountcorrBedgraphfile+" "+tdffile+" "+igvgenomefile+"\n")
        wf.write("echo tdf\n")
        wf.write("date\n")



def genebodycov(sortedbam, outdir, wfname, genome):
	wf = open(wfname, "a")
        wf.write("date\n")
	bamqualdir =outdir+"bamqual/"
	ensure_dir(bamqualdir)
	genbodyqualdir = bamqualdir+"genebodycov/"
	ensure_dir(genbodyqualdir)
	bed12 = "/scratch/Shares/public/genomes/allfilesmadefromigenomesfiles/bed12from_genesgtf/"+genome+".bed"
	ensure_dir(genbodyqualdir)
        wf.write("echo "+sortedbam+"\n")
        wf.write("echo "+sortedbamroot+"\n")
        wf.write("module load python/2.7.14/rseqc\n")
	wf.write("cd "+genbodyqualdir+"\n")
	wf.write("geneBody_coverage.py -i "+sortedbam+" -o "+genbodyqualdir+sortedbamroot+".rseqc -r "+bed12+" >"+genbodyqualdir+sortedbamroot+".rseqc.genebodycov.out 2>"+genbodyqualdir+sortedbamroot+".rseqc.genebodycov.err\n")
        wf.write("date\n")
	wf.close()



def createwfile(outdir, rootname,subpartlabel,email,processers,time,queue,sbatchmem):
	wfname = outdir+"qsubscripts/"+rootname+"_"+subpartlabel+".slurm"
        wf = open(wfname, "w")
        wf.write("#!/bin/bash\n")
        wf.write("#SBATCH --job-name="+rootname+"  # Job name\n")
        wf.write("#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)\n")
        wf.write("#SBATCH --mail-user="+email+"# Where to send mail\n")
        wf.write("#SBATCH --nodes=1\n")
        wf.write("#SBATCH --ntasks="+processers+"# Number of CPU (processer cores i.e. tasks) In this example I use 32 for bowtie2. I only need one, since none of the commands I run are parallelized.\n")
        wf.write("#SBATCH --time="+time+" # Time limit hrs:min:sec\n")
        wf.write("#SBATCH -p "+queue+"\n")
        wf.write("#SBATCH --mem="+sbatchmem+" # Memory limit\n")
        wf.write("#SBATCH --output="+outdir+"e_and_o/"+rootname+"_"+subpartlabel+".%j.out\n")
        wf.write("#SBATCH --error="+outdir+"e_and_o/"+rootname+"_"+subpartlabel+".%j.err\n")
        wf.write("\n\n\npwd; hostname; date\n")
        wf.close()
	return wfname


def ensure_dir(file_path):
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        os.makedirs(directory)


def addtosubmit(submitscript, x, runfilename, dependencis):
	wf = open(submitscript, "a")
	line = "jid"+str(x)+"=$(sbatch "
	if len(dependencis)>0:
		line = line+"--dependency=afterany"
		for dependence in dependencis:
			line=line+":$jid"+str(dependence)
	line = line+" "+runfilename+")\n"
	wf.write(line)
	line = "jid"+str(x)+"=$(echo $jid"+str(x)+" | cut -d ' ' -f 4-)\n"
	wf.write(line)
	wf.close()	
	return x

inputfiletype="fastq"

def main(infile, outdir, email, genome,processers=processers):
	#make outdir
	ensure_dir(outdir)
	#make e_and_o
	qsubdir = outdir+"qsubscripts/"
	ensure_dir(qsubdir)
	os.system("cp --no-clobber "+sys.path[0]+"/readcountcorrectBG.py "+qsubdir+"\n")
	ensure_dir(outdir+"e_and_o/")
	ensure_dir(outdir+"qsubscripts/")
	indirfiles=[infile]
	x=0
	for filename in indirfiles:
		passfile = filename
		rootname = filename.split("/")[-1]
		rootname = rootname.strip(inputfiletype)
		#create a master script that controls the chained scripts
		time="95:00:00"
		queue="long"
		subpartlabel="mom_run"
		submitscript = createwfile(outdir, rootname,subpartlabel,email,processers,time,queue,sbatchmem)
		if args.flipreads==True:
			time="5:00:00"
			queue="short"
			#create script header for fliping reads
			subpartlabel="rc"
			filpfilename = createwfile(outdir, rootname,subpartlabel,email,processers,time,queue,sbatchmem)
			#add flip and quality command to script
			passfile = flipfile(passfile, outdir, filpfilename) 
			checkqualityfile(passfile, outdir, filpfilename)
			#add to flip and quality script to submitscript
			filpfilenamejobid = addtosubmit(submitscript, x, filpfilename, dependencis=[])
			x=x+1
		if args.Turnoff_trimgalore==False:
			time="10:00:00"
			queue="short"
			#create script header for triming reads
			subpartlabel="trim"
			trimfilename = createwfile(outdir, rootname,subpartlabel,email,processers,time,queue,sbatchmem)
			#add trim and quality command to script
			passfile = trimgalorefile(passfile, outdir, trimfilename)
			checkqualityfile(passfile, outdir, trimfilename)
			#add to flip and quality script to submitscript
			if  args.flipreads==True:
                        	trimfilenamejobid = addtosubmit(submitscript,x,trimfilename, dependencis=[filpfilenamejobid])
				x=x+1
			else:
				trimfilenamejobid = addtosubmit(submitscript,x,trimfilename, dependencis=[])
				x=x+1
		if args.Turnoff_fastqtosam==False:
			#create script header for mapping reads
			time="48:00:00"
			queue="long"
			processers="32"
			subpartlabel="map_"
                        bowtiescriptfilename = createwfile(outdir, rootname,subpartlabel,email,processers,time,queue,sbatchmem)
			#add mapping command to script
			mappedfile = bowtie2file(passfile,  outdir, bowtiescriptfilename, args.genome,processers)
			#add bowtie script to submit script
			if args.flipreads==True:
				if args.Turnoff_trimgalore==False:
					bowtiescriptjobid=addtosubmit(submitscript,x,bowtiescriptfilename, dependencis=[filpfilenamejobid,trimfilenamejobid])
					x=x+1
				else:
					bowtiescriptjobid=addtosubmit(submitscript,x,bowtiescriptfilename, dependencis=[filpfilenamejobid])
					x=x+1
			else:
				if args.Turnoff_trimgalore==False:
					bowtiescriptjobid=addtosubmit(submitscript,x,bowtiescriptfilename, dependencis=[trimfilenamejobid])
					x=x+1
				else:
					bowtiescriptjobid=addtosubmit(submitscript,x,bowtiescriptfilename, dependencis=[])
					x=x+1
			processers="1"
		#do we need to make sams into something?
		if args.Turnoff_samtobam==False or args.Turnoff_bamtosortedbam==False or args.Turnoff_sortedbamtobai==False:
			#create script header for makingbamuseful reads
			time="24:00:00"
			currmemamount="200"
			queue="short"
			currsbatchmem = currmemamount+memtype
			subpartlabel="samtobam"
			bamtouseful = createwfile(outdir, rootname,subpartlabel,email,processers,time,queue,currsbatchmem)
			if args.Turnoff_samtobam==False:
				mappedfile = samtobamfile(mappedfile, outdir, bamtouseful)
			if args.Turnoff_bamtosortedbam==False:
				memdic = {"gb":"G", "kb":"K", "mb":"M"}
				samtoolsmem = currmemamount+memdic[memtype]				
				sortedbam = bamtosortedbamfile(mappedfile, outdir, bamtouseful,samtoolsmem)
			if args.Turnoff_sortedbamtobai==False:
				sortedbamtobai(sortedbam, bamtouseful)
			#add bedgraphsandtdfs script to submit script
			sortedbamjobid = addtosubmit(submitscript,x,bamtouseful, dependencis=[bowtiescriptjobid])
			x=x+1
		#do we need to make sorted bams into bedgraphs and tdfs?
		if args.Turnoff_sortedbamtobedgraph==False or args.Turnoff_readcountcorrectbedgraph==False or args.Turnoff_igvcreatetdf==False:
			time="10:00:00"
			queue="short"
			currmemamount="200"
			currsbatchmem = currmemamount+memtype
			#create script header for making bedgraphsandtdfs
			subpartlabel="bamtobedgraph"
			bamtotdf = createwfile(outdir, rootname,subpartlabel,email,processers,time,queue,currsbatchmem)
			if args.Turnoff_sortedbamtobedgraph==False:
				Bedgraphfile = sortedbamtobedgraphfile(sortedbam, outdir, bamtotdf, genome)
			if args.Turnoff_readcountcorrectbedgraph==False:
				readcountcorrBedgraphfile = readcountcorrectbedgraphfile(Bedgraphfile, sortedbam, outdir, bamtotdf)
			if args.Turnoff_igvcreatetdf==False:
				igvcreatefile(readcountcorrBedgraphfile, outdir, bamtotdf, genome)
			#add bedgraphsandtdfs script to submit script
			bamtotdfjobid = addtosubmit(submitscript,x,bamtotdf, dependencis=[sortedbamjobid])
			x=x+1
		#create bam quality runs

		wf = open(submitscript, "a")
		wf.write("date\n")
		wf.close()
		if args.Turnoff_submit==False:	
			os.system('sbatch '+submitscript)	
	
	if args.Turnoff_submit==False:
		print "To check if your scripts are still runing, type squeue or qstat on the command line."
		print "If your scripts error the scripts are in "+outdir+"qsubscripts/ and the error files are in "+outdir+"e_and_o/"
		print "Once your jobs have ended check the files in "+outdir+"qual/ for quality reports."
		print "Make sure you rsync what you want backed up to /projects/."
	else:
		print "You now need to check your shell scripts in "+outdir+"qsubscripts/"
		print "Once you decide you like them submit using sbatch"
		




if __name__=="__main__":
	parser = argparse.ArgumentParser(description='This program will make a slurm script for processing GRO-seq data and submit it to the queue on fiji.\n Buy defalut the program will map the reads using bowtie2, then convert the sam to a sorted bam file, then create bedgraphs, then make tdf files which can be used in the program IGV.',usage='python create_scipts_to_map_on_fiji_argparse.py <infile> <outdir> <genome> <email>')
	parser.add_argument("-f","--flipreads", action='store_true', dest="flipreads", default=False)
	parser.add_argument("-q","--Turnoff_checkquality", action='store_true', default=False)
	parser.add_argument("-t","--Turnoff_trimgalore", action='store_true', default=False)
	parser.add_argument("-m","--Turnoff_fastqtosam", action='store_true', default=False)
	parser.add_argument("-s2b","--Turnoff_samtobam", action='store_true', default=False)
	parser.add_argument("-b2sb","--Turnoff_bamtosortedbam", action='store_true', default=False)
	parser.add_argument("-ib","--Turnoff_sortedbamtobai", action='store_true', default=False)
	parser.add_argument("-bg","--Turnoff_sortedbamtobedgraph", action='store_true', default=False)
	parser.add_argument("-c","--Turnoff_readcountcorrectbedgraph", action='store_true', default=False)
	parser.add_argument("-igv","--Turnoff_igvcreatetdf", action='store_true', default=False)
	parser.add_argument("-s","--Turnoff_submit", action='store_true', default=False)
	parser.add_argument("infile", help="Directory with fastq or bam files. Must be located on /scratch/.", type=str)
	parser.add_argument("outdir", help="Directory files will be output into. The program will create this direcoty if it does not exist. This direcory must be located on /scratch/.", type=str)
	genomeschoices = bowtieindexdic.keys()
	parser.add_argument("genome", metavar="genome",help="Which genome do you wish to map too? Allowed values are "+', '.join(genomeschoices)+',', type=str, choices=genomeschoices)
	parser.add_argument("email", help="Email address is required. If you mistype your email then IT gets your emails. Please don't do that. ", type=str)
	args = parser.parse_args()
	print args
##this is built for mapping single end nascent data!!! This program expects a outdir. They must both end in a slash.
	if args.outdir.startswith("/scratch/") and args.infile.startswith("/scratch/"):
			ensure_dir(args.outdir)
                        main(args.infile, args.outdir, args.email, args.genome)
	else:
		print "infile and outdir must be in /scratch/"


 
#test
#python create_scipts_to_map_on_fiji.py /scratch/Users/allenma/171025_NB501447_0179_fastq/Allen_Dowell-371/ /scratch/Users/allenma/pro/ hg19 allenma@colorado.edu 
