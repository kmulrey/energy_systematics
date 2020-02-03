#! /bin/bash
#SBATCH --time=1-23:50:00
#SBATCH -p normal
#SBATCH -N 1 -n 16
# Take an entire node i.e. 16 cores per instance 
# to get around forced multi-processing in Scipy / LAPack etc. ??

echo hostname

export LOFARSOFT=/vol/optcoma/pycrtools
export PYTHONPATH=$LOFARSOFT/release/lib/python:$PYTHONPATH
export LD_LIBRARY_PATH=$LOFARSOFT/release/lib:$LD_LIBRARY_PATH

PYCRTOOLS=$LOFARSOFT/src/PyCRTools

BASE_PATH=/vol/astro3/lofar/sim/pipeline
# BASE_PATH=/vol/astro3/lofar/sim/pipeline
# DB_PATH=$BASE_PATH
# DB_FILE=$DB_PATH/bogus
# DATA_PATH=$BASE_PATH/data
# RESULTS_PATH=$BASE_PATH/results
# LORA_PATH=$BASE_PATH/LORA
DATA_DIR=$BASE_PATH/events
SIMULATION_DIR=$BASE_PATH/run
OUTPUT_DIR=$BASE_PATH/production_analysis
OUTPUT_DIR_RADIO_ONLY=$BASE_PATH/production_analysis_radio_only
MCVSMC_DIR=$BASE_PATH/production_mcvsmc_radio_only # or _radio_only
LOG_DIR=/vol/astro3/lofar/sim/pipeline/log

#LOG_PATH=$BASE_PATH/log

# Get event id for current task and mark it as NEXT
EVENT_ID=$(awk "NR==$SLURM_ARRAY_TASK_ID" $HOME/cr_xmaxfit_lorasim_done_astro3.txt)
#EVENT_ID=$(awk "NR==$SLURM_ARRAY_TASK_ID" $HOME/eventlist_meth_todo.txt)
LOGFILE=$LOG_DIR/cr_xmaxfit-$EVENT_ID.txt

# Set permissions on resulting files
umask 002

echo $EVENT_ID


python cr_xmaxfit.py --event=$EVENT_ID --lorafile-suffix=_GeVfix --datadir=$DATA_DIR --simulationdir=$SIMULATION_DIR --iteration=0 --outputdir=$OUTPUT_DIR --outputdir-radio-only=$OUTPUT_DIR_RADIO_ONLY --mcvsmcdir=$MCVSMC_DIR --logdir=$LOG_DIR --filtdir=

# Run the pipeline. Currently set to fetch newest LOFAR (cr-physics pipeline) data from CR Database

#/usr/bin/python -u /vol/astro3/lofar/sim/pipeline/pipeline_job/cr_xmaxfit.py --event=$EVENT_ID --datadir=$DATA_DIR --simulationdir=$SIMULATION_DIR --outputdir=$OUTPUT_DIR --mcvsmcdir=$MCVSMC_DIR --logdir=$LOG_DIR --fetch-lofardata --rewrite-lofardata --mcvsmc --radio-only-mcvsmc &> $LOGFILE
#/usr/bin/python -u /vol/astro3/lofar/sim/pipeline/pipeline_job/cr_xmaxfit.py --event=$EVENT_ID --force-reprocess --datadir=$DATA_DIR --simulationdir=$SIMULATION_DIR --outputdir=$OUTPUT_DIR --outputdir-radio-only=$OUTPUT_DIR_RADIO_ONLY --mcvsmcdir=$MCVSMC_DIR --logdir=$LOG_DIR --fetch-lofardata --rewrite-lofardata --mcvsmc --radio-only-mcvsmc &> $LOGFILE
# No --lorafile-suffix

## Removed --mcvsmc flag, doing it separately here:
#/usr/bin/python -u /vol/astro3/lofar/sim/pipeline/scripts/mcvsmc/mcvsmc_updated.py --event=$EVENT_ID --iteration=0  --basedir=/vol/astro7/lofar/sim/pipeline/run --recodir=/vol/astro7/lofar/sim/pipeline/test_analysis --outputdir=/vol/astro7/lofar/sim/pipeline/test_mcvsmc >> /vol/astro3/lofar/sim/pipeline/log/cr_xmaxfit-$EVENT_ID-mcvsmc.txt
