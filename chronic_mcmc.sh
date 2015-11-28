#$ -N CCMCMC
#$ -V -cwd
#$ -M sebastian.funk@lshtm.ac.uk
#$ -m n
#$ -l mem_free=300M,h_vmem=6000M
#$ -t 1-26
Rscript $HOME/code/trypanosomiasis/chronic_run.r -v ${SGE_TASK_ID} -n 1000000
