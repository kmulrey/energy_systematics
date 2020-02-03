#! /bin/bash
#SBATCH --time=1-23:50:00
#SBATCH -p normal
#SBATCH -N 1 -n 16
# Take an entire node i.e. 16 cores per instance 
# to get around forced multi-processing in Scipy / LAPack etc. ??
#SBATCH --output /home/acorstanje/slurm_logs/slurm-%j.out
#SBATCH --error /home/acorstanje/slurm_logs/slurm-error-%j.out

echo hostname

export LOFARSOFT=/vol/optcoma/pycrtools
export PYTHONPATH=$LOFARSOFT/release/lib/python:$PYTHONPATH
export LD_LIBRARY_PATH=$LOFARSOFT/release/lib:$LD_LIBRARY_PATH

PYCRTOOLS=$LOFARSOFT/src/PyCRTools

BASE_PATH=/vol/astro7/lofar/sim/pipeline

# BASE_PATH=/vol/astro3/lofar/sim/pipeline
# DB_PATH=$BASE_PATH
# DB_FILE=$DB_PATH/bogus
# DATA_PATH=$BASE_PATH/data
RESULTS_PATH=/vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/sim_tests/results
# LORA_PATH=$BASE_PATH/LORA
DATA_DIR=$BASE_PATH/events
SIMULATION_DIR=$BASE_PATH/run
OUTPUT_DIR=$RESULTS_PATH/production_analysis_Feb2020
OUTPUT_DIR_RADIO_ONLY=$RESULTS_PATH/production_analysis_radio_only_Feb2020
MCVSMC_DIR=$RESULTS_PATH/production_mcvsmc_radio_only_Oct2019 # or _radio_only
LOG_DIR=$RESULTS_PATH/log

#LOG_PATH=$BASE_PATH/log

# Get event id for current task and mark it as NEXT
EVENT_ID=$(awk "NR==$SLURM_ARRAY_TASK_ID" $HOME/cr_physics_new)
#EVENT_ID=$(awk "NR==$SLURM_ARRAY_TASK_ID" $HOME/eventlist_meth_todo.txt)
LOGFILE=$LOG_DIR/cr_xmaxfit-$EVENT_ID.txt

# Set permissions on resulting files
umask 002

echo $EVENT_ID
echo $PYCRTOOLS
# Run the pipeline. Currently set to fetch newest LOFAR (cr-physics pipeline) data from CR Database

#/usr/bin/python -u /vol/astro3/lofar/sim/pipeline/pipeline_job/cr_xmaxfit.py --event=$EVENT_ID --datadir=$DATA_DIR --simulationdir=$SIMULATION_DIR --outputdir=$OUTPUT_DIR --mcvsmcdir=$MCVSMC_DIR --logdir=$LOG_DIR --fetch-lofardata --rewrite-lofardata --mcvsmc --radio-only-mcvsmc &> $LOGFILE
# Adding --force-reprocess 
# temp. removed --mcvsmc --radio-only-mcvsmc --force-reprocess (Feb 26, 2019)
# temp. added --force-reprocess back in (Oct 2019)
/usr/bin/python -u $PYCRTOOLS/xmax_pipeline/cr_xmaxfit.py --event=$EVENT_ID --lorafile-suffix=_GeVfix --datadir=$DATA_DIR --simulationdir=$SIMULATION_DIR --outputdir=$OUTPUT_DIR --outputdir-radio-only=$OUTPUT_DIR_RADIO_ONLY --mcvsmcdir=$MCVSMC_DIR --logdir=$LOG_DIR --fetch-lofardata --rewrite-lofardata --force-reprocess &> $LOGFILE
# Doing --radio-only-mcvsmc now

