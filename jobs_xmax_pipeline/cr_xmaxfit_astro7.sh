#! /bin/bash
#SBATCH --time=1-23:50:00
#SBATCH -p normal
#SBATCH -N 1 -n 16

echo hostname
use root
export LOFARSOFT=/vol/optcoma/pycrtools
export PYTHONPATH=$LOFARSOFT/release/lib/python:$PYTHONPATH
export LD_LIBRARY_PATH=$LOFARSOFT/release/lib:$LD_LIBRARY_PATH

#EVENT_ID=$(awk "NR==$SLURM_ARRAY_TASK_ID" $HOME/cr_physics_new)

#EVENT_ID=$SLURM_ARRAY_TASK_ID

echo $SLURM_ARRAY_TASK_ID



# Set permissions on resulting files
umask 002

#echo $EVENT_ID
echo $PYCRTOOLS

cd /vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/sim_tests/energy_systematics/xmax_pipeline


#python cr_xmaxfit.py --event=$EVENT_ID --lorafile-suffix=_GeVfix --datadir=$DATA_DIR --simulationdir=$SIMULATION_DIR --iteration=0  --outputdir-radio-only=$OUTPUT_DIR_RADIO_ONLY --mcvsmcdir=$MCVSMC_DIR --logdir=$LOG_DIR --filtdir=$NEWSIM_PATH --writedir=$WRITE_FILT --collectdir=$COLLECT_DIR

python run_cr_xmax.py -i $SLURM_ARRAY_TASK_ID
