#!/bin/bash
#set -x

JOB=`basename $1 .xyz`

mkdir $JOB

cp $m $JOB/

cd $JOB

WORKDIR=`pwd`

cat << EOF > ${JOB}.run
#!/bin/bash

#SBATCH --time=05:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32 
#SBATCH --partition=thin


RUNDIR="/scratch-shared/${JOB}.\${SLURM_JOBID}"

mkdir -p \$RUNDIR

cd \${RUNDIR}

srun --time=05:00:00 --partition=thin crest \$1 --alpb thf --tstep 2 -origin --verbose -nocross

cp crest_* coord* *.log  ${WORKDIR}

#cd ${WORKDIR}

rm -rf \${RUNDIR}

EOF

chmod 755 ${JOB}.run
sbatch ${JOB}.run $m

cp crest_* coord* *.log .xcontrol ${WORKDIR}
	
