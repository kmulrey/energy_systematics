#! /bin/bash

MAX_SIMULTANEOUS_JOBS=1

export LOFARSOFT=/vol/optcoma/pycrtools
export PYTHONPATH=$LOFARSOFT/release/lib/python:$PYTHONPATH
export LD_LIBRARY_PATH=$LOFARSOFT/release/lib:$LD_LIBRARY_PATH

PYCRTOOLS=$LOFARSOFT/src/PyCRTools

# Remove logfile containing failed event IDs
#rm /vol/astro7/lofar/sim/pipeline/log/failed.txt

# Run cr_physics pipeline
/usr/bin/python -u $PYCRTOOLS/extras/get_event_ids.py --host coma00.science.ru.nl --user crdb --password crdb --dbname crdb --simulation-status=LORASIM_DONE > $HOME/cr_xmaxfit_lorasim_done_astro7.txt
#/usr/bin/python -u /vol/astro2/users/acorstanje/lofarsoft/src/PyCRTools/scripts/acorstanje/eventlist_meth_todo.py > $HOME/cr_xmaxfit_coreasdone_astro7.txt
sleep 5
JOBS=$(sed -n '$=' $HOME/cr_xmaxfit_lorasim_done_astro7.txt)

echo "processing $JOBS astro7 events (LORASIM_DONE)"

/usr/local/slurm/bin/sbatch -p normal --array 1-${JOBS}%${MAX_SIMULTANEOUS_JOBS} cr_xmaxfit_astro7.sh

# Run cr_event pipeline
#/usr/bin/python -u $PYCRTOOLS/extras/get_event_ids.py --host coma00.science.ru.nl --user crdb --password crdb --dbname crdb -a -s NEW > $HOME/cr_event_new
#
#JOBS=$(sed -n '$=' $HOME/cr_event_new)
#
#echo "processing $JOBS NEW events"
#
#/usr/local/slurm/bin/sbatch --array 1-${JOBS}%${MAX_SIMULTANEOUS_JOBS} $PYCRTOOLS/jobs/cr_event.sh

