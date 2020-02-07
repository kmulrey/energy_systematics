#! /bin/bash
#SBATCH --time=1-23:50:00
#SBATCH -p normal
#SBATCH -N 1 -n 16

echo hostname

export LOFARSOFT=/vol/optcoma/pycrtools
export PYTHONPATH=$LOFARSOFT/release/lib/python:$PYTHONPATH
export LD_LIBRARY_PATH=$LOFARSOFT/release/lib:$LD_LIBRARY_PATH

PYCRTOOLS=$LOFARSOFT/src/PyCRTools
BASE_PATH=/vol/astro7/lofar/sim/pipeline
RESULTS_PATH=/vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/sim_tests/results
NEWSIM_PATH=/vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/sim_tests/
DATA_DIR=$BASE_PATH/events
SIMULATION_DIR=$BASE_PATH/run
#OUTPUT_DIR=$RESULTS_PATH/production_analysis_cal2019
#OUTPUT_DIR_RADIO_ONLY=$RESULTS_PATH/production_analysis_radio_only_Feb2020
OUTPUT_DIR_RADIO_ONLY=$RESULTS_PATH/production_analysis_radio_only_thetaPLUS5
COLLECT_DIR=/vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/sim_tests/filtered_thetaPLUS5
WRITE_FILT=/vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/sim_tests/events_thetaPLUS5
#OUTPUT_DIR_RADIO_ONLY=$RESULTS_PATH/production_analysis_radio_only_cal2019
#MCVSMC_DIR=$RESULTS_PATH/production_mcvsmc_radio_only_Oct2019 # or _radio_only
LOG_DIR=$RESULTS_PATH/log


EVENT_ID=$(awk "NR==$SLURM_ARRAY_TASK_ID" $HOME/cr_physics_new)
#EVENT_ID=158064938
echo $EVENT_ID


LOGFILE=$LOG_DIR/cr_xmaxfit-$EVENT_ID.txt

# Set permissions on resulting files
umask 002

echo $EVENT_ID
echo $PYCRTOOLS

cd /vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/sim_tests/energy_systematics/xmax_pipeline


python cr_xmaxfit.py --event=$EVENT_ID --lorafile-suffix=_GeVfix --datadir=$DATA_DIR --simulationdir=$SIMULATION_DIR --iteration=0  --outputdir-radio-only=$OUTPUT_DIR_RADIO_ONLY --mcvsmcdir=$MCVSMC_DIR --logdir=$LOG_DIR --filtdir=$NEWSIM_PATH --writedir=$WRITE_FILT --collectdir=$COLLECT_DIR
