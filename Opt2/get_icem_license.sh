#!/bin/bash
## Checker for free ANSYS Licenses
## Origin by Laszlo Daroczy
## Modified to Check if License for ANSYS ICEM CFD is free


## License Amounts and safty margins for free licenses
CENTRAL_ALL=30
CENTRAL_BORROWED=0
CENTRAL_PROTECTED=0
CENTRAL_AVAILABLE=0
TEACHING_ALL=25
TEACHING_AVAILABLE=0
TEACHING_BORROWED=0
TEACHING_PROTECTED=0
SAFETY_MARGIN=1
SUCCESS=0

sleep 4

export LM_LICENSE_FILE=1055@liclux.yourlicenseserver.de


while [ $SUCCESS -ne 1 ]; do
	# Get system wide values for minimum number of free licenses
	TEACHING_PROTECTED=($(</home/user/licenses/ANSYS/Teaching) )
	CENTRAL_PROTECTED=($(</home/user/licenses/ANSYS/Central) )
	# Get currently used licences
	CENTRAL_BORROWED=$(/ansys_install_dir/fluent/bin/lmstat -f aa_r_cfd | grep "Total of" | awk '{print $11}')
	TEACHING_BORROWED=$(/ansys_install_dir/fluent/bin/lmstat -f aa_t_a | grep "Total of" | awk '{print $11}')

	# compute available
	let CENTRAL_AVAILABLE=$CENTRAL_ALL-$CENTRAL_BORROWED
	let TEACHING_AVAILABLE=$TEACHING_ALL-$TEACHING_BORROWED

	#Preference List: Parallel (1); Chair (2); Teaching (3); CFD (4)


if [ $TEACHING_AVAILABLE -gt $TEACHING_PROTECTED ]; then
	echo 3 > version.txt
				SUCCESS=1
	else

		if [ $CENTRAL_AVAILABLE -gt $CENTRAL_PROTECTED ]; then
							echo 4 > version.txt
							SUCCESS=1
			else
							sleep 4
			fi


fi

done
