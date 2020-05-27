#!/bin/bash
#
# prepares runs for mesh dependency study
# the parameter nglob is set to the product of $BASEFACTOR and $run

RUNNAMEBASE="M"
BASEFACTOR=0.05

for run in {01..20}; do
	
	# create copy for individual run
	COMPOSEDFOLDERNAME=${RUNNAMEBASE}${run}
	cp -r common $COMPOSEDFOLDERNAME

	# paste individual mesh size factor
	SCALINGFACTOR=$(echo "scale=2; ${BASEFACTOR} * ${run}" | bc)
	sed -i "s/\$MDSNGLOB/${SCALINGFACTOR}/" ${COMPOSEDFOLDERNAME}/in_parameter.tcl

done
