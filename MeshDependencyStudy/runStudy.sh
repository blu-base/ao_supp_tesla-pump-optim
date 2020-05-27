#!/bin/bash
#
# Run cases sequentially on the current machine

RUNNAMEBASE="M"

for run in {01..20}; do
	
	COMPOSEDFOLDERNAME=${RUNNAMEBASE}${run}
	cd $COMPOSEDFOLDERNAME

	./Allrun.sh

	cd ..

done
