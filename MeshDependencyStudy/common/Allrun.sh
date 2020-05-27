#!/bin/bash

## Set uo the environment
export LM_LICENSE_FILE=1055@liclux.yourlicenseserver.de
## Additional environment settings might be necessary to
## make your shell aware of the following commands

# run icem to create the mesh
icemcfd -batch -script in_icem.tcl > icemlog.txt 2>&1


# run the simulation
starccm+ -new -batch in_tesrot.java -licpath 1999@flex.cd-adapco.com -podkey XXXXXXXXXXXXXXXXXXXXXX -rsh ssh -np 16 -power -collab > output.txt 2>&1

# Extract values into single value files
./in_extractvalues.sh presDrop
./in_extractvalues.sh torque
./in_extractvalues.sh speed
./in_extractvalues.sh hemolysis
./in_extractvalues.sh eff
