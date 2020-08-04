#!/bin/bash
#
# prepares runs for mesh dependency study
# the parameter nglob is set to the product of $BASEFACTOR and $run

RUNNAMEBASE="C"
BASEFACTOR=0.01

for run in {01..10}; do
	
	# create copy for individual run
	COMPOSEDFOLDERNAME=${RUNNAMEBASE}${run}
	cp -L -r common $COMPOSEDFOLDERNAME

	# paste individual mesh size factor
	SCALINGFACTOR=$(echo "scale=2; ${BASEFACTOR} * ${run}" | bc)
	sed -i "s/\$MDSNGLOB/${SCALINGFACTOR}/" ${COMPOSEDFOLDERNAME}/in_parameter.tcl

done
