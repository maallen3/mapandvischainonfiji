#Mary Allen
#May 24, 2018

This package will map GRO-seq data with standard settings on fiji. 
To use

python create_scipts_to_map_on_fiji_argparse.py <indir> <outdir> <genome> <email>


For more options type
python create_scipts_to_map_on_fiji_argparse.py --help

For instance you could turn off the auto submit like this.
python create_scipts_to_map_on_fiji_argparse.py <indir> <outdir> <genome> <email> --Turnoff_submit

By default this script will create 5 slurm scripts and submit them all to the queue. 
<filename>._mom_run.slurm
<filename>._trim.slurm
<filename>._map_.slurm
<filename>._samtobam.slurm
<filename>._bamtobedgraph.slurm

If you turn on --flipreads one other script will be created. The _mom_run runs the other scripts. 

 
