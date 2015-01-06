#$ -N CCMCMC
#$ -V -cwd
#$ -M sebastian.funk@lshtm.ac.uk
#$ -m n
#$ -l mem_free=300M,h_vmem=6000M
#$ -t 1:26
LD_LIBRARY_PATH=/opt/intel/compiler101/x86_64/lib:$LD_LIBRARY_PATH Rscript $HOME/code/trypanosomiasis/chronic_lhs.R ${SGE_TASK_ID}
