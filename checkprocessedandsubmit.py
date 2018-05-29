import glob
import os
import pandas
import subprocess as sub
import time
from datetime import datetime

toprocessesfile="/scratch/Users/allenma/scripts/mapandvischain/GEOandSRA_markkeyword021918.csv"
indir="/scratch/Shares/public/nascentdb/fastqs/"
outdir="/scratch/Shares/public/nascentdb/processedv2.0/"
email="allenma@colorado.edu"

#--------
mappingprefs = {"Homo sapiens":"hg38","Drosophila melanogaster":"dm6", "Mus musculus":"mm10","Saccharomyces cerevisiae":"sacCer3", "Arabidopsis thaliana":"TAIR10","Caenorhabditis elegans":"ce10"}


genomesdir="/scratch/Shares/public/genomes/"
def createbowtieindexdic():
        b = {}
        bowtiefiles = [infile for infile in glob.glob(os.path.join(genomesdir, '*/*/*/Sequence/Bowtie2Index/genome*'))]
        for bowtiefile in bowtiefiles:
                dname = os.path.dirname(bowtiefile)
                rootname=bowtiefile.split("/")[7]
                b[rootname]=dname+"/genome"
        return b

bowtieindexdic = createbowtieindexdic()

print bowtieindexdic

def processthiskeyword(keyword):
	df = pandas.read_csv(toprocessesfile,index_col=0) 
	print df.columns
	print df.index
#	df = df.filter(regex=keyword, axis=0)
	df = df[df['keyword'].str.contains(keyword,na=False)]
	return df




def allkeywords():
        df = pandas.read_csv(toprocessesfile,index_col=0)
        keys = list(set(df['keyword'].tolist()))
	cleanedList = [x for x in keys if str(x) != 'nan']
	return cleanedList

def submitone(keyword, overwrite=False):
	df = processthiskeyword(keyword=keyword)
	bamfilespresent = [infile.split("/")[-1] for infile in glob.glob(os.path.join(outdir+"sortedbams/", '*.bam'))]
	bamfilespresent = [infile.split(".")[0] for infile in bamfilespresent]
	if overwrite==True:
		bamfilespresent=[]
	print bamfilespresent
	for sample in df.index:
		if not sample in bamfilespresent:
			thisrow = df.loc[sample]
			if str(thisrow["ScientificName"])!="nan":
				#print "python create_scipts_to_map_on_fiji_argparse.py "+indir+sample+".fastq "+outdir+" "+mappingprefs[thisrow["ScientificName"]]+" "+email+" --Turnoff_submit"
				os.system("python create_scipts_to_map_on_fiji_argparse.py "+indir+sample+".fastq "+outdir+" "+mappingprefs[thisrow["ScientificName"]]+" "+email)	




def dailylooppernight(numpernight):
	keys = allkeywords()
	print keys
	for i, key in enumerate(keys):
		if i % numpernight == 0:
    			for j in range(numpernight):
				print i, j, keys[j]
		

#dailyloop(5)


def numjobs():
	p = sub.Popen(['qstat',"-u","allenma"],stdout=sub.PIPE,stderr=sub.PIPE)
        output, errors = p.communicate()
        countnow=0
        jobs = []
        for line in output.split("\n"):
                if line.startswith("-"):
                        countnow=1
                elif line=="":
                        countnow=0
                if countnow==1:
                        jobs.append(line)
	return len(jobs)

def numberinqueueloopper(maxnumjobs):
	keys = allkeywords()
	for i, key in enumerate(keys):
		numj = numjobs()
		if numj<maxnumjobs:
			print numj, key
			submitone(key, overwrite=False)	
			

def runnow(maxnumjobs, sleepnow=True):
        now = datetime.now()
        midnight = now.replace(hour=0, minute=0, second=0, microsecond=0)
        seconds = (now - midnight).seconds
        secondsuntilmidnit=86400-seconds
	if sleepnow==True:
		time.sleep(secondsuntilmidnit)#8 hours
	else:
		numberinqueueloopper(maxnumjobs)
	for i in range(25):
		numberinqueueloopper(maxnumjobs)
		now = datetime.now()
		midnight = now.replace(hour=0, minute=0, second=0, microsecond=0)
		seconds = (now - midnight).seconds
		secondsuntilmidnit=86400-seconds
		time.sleep(secondsuntilmidnit)#twenty four hours
		

runnow(500)

#submitone("Compagno2017")
#submitone("Allen2014", overwrite=True)
#python create_scipts_to_map_on_fiji_argparse.py $fastqfile $outdir $genome $email
#this runs without flipping reads and defult is to trim
#python create_scipts_to_map_on_fiji.py /scratch/Users/allenma/171025_NB501447_0179_fastq/Allen_Dowell-371/ /scratch/Users/allenma/pro/ hg19 allenma@colorado.edu
#this filps the reads
#python create_scipts_to_map_on_fiji_argparse.py <infile> <outdir/> hg19 <email> --flipreads
