#!/bin/bash

## Set uo the environment
export LM_LICENSE_FILE=1055@liclux.yourlicenseserver.de
## Additional environment settings might be necessary to
## make your shell aware of the following commands

# run icem to create the mesh
/home/software/ansys_inc/v170/icemcfd/linux64_amd/bin/icemcfd -batch -script in_icem.tcl > icemlog.txt 2>&1


# run the simulation
/home/software/CD-adapco/14.04.013-R8/STAR-CCM+14.04.013-R8/star/bin/starccm+ -new -batch in_tesrot.java -licpath 1999@flex.cd-adapco.com -podkey dEqZa/O8SwuZSGcH4MmNHw -rsh ssh -np 6 -power -collab > output.txt 2>&1

# Extract values into single value files
./in_extractvalues.sh presDrop
./in_extractvalues.sh torque
./in_extractvalues.sh speed
./in_extractvalues.sh hemolysis
./in_extractvalues.sh eff


## clean up large files
rm hex.uns
rm slice.msh

#rm slice.sim
