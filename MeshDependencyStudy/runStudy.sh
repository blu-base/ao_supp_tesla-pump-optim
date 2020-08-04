#!/bin/bash
#
# Run cases sequentially on the current machine

RUNNAMEBASE="C"

for run in {01..10}; do
	
	COMPOSEDFOLDERNAME=${RUNNAMEBASE}${run}
	cd $COMPOSEDFOLDERNAME

	./Allrun.sh

	cd ..

done
