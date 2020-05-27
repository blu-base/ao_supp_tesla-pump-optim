#!/bin/bash
#
# Resets directory

RUNNAMEBASE="M"
BASEFACTOR=0.05

for run in {01..20}; do
	
	# remove run directories
	COMPOSEDFOLDERNAME=${RUNNAMEBASE}${run}
	rm -r $COMPOSEDFOLDERNAME 2>/dev/null


done
