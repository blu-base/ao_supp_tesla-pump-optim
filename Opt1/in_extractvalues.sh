#!/bin/bash
# Helper-Script

if [ "$1" = "presDrop" ]
then
	echo "$(tail -n 1  < presDrop.csv | cut -d',' -f2)" > presDrop.val
fi
if [ "$1" = "torque" ]
then
	echo "$(tail -n 1  < torque.csv | cut -d',' -f2)" > torque.val
fi
if [ "$1" = "speed" ]
then
	echo "$(tail -n 1  < rotationSpeed.csv | cut -d',' -f2)" > speed2.val
fi
if [ "$1" = "hemolysis" ]
then
	echo "$(tail -n 1  < HIoutlet.csv | cut -d',' -f2)" > hemolysis.val
fi
if [ "$1" = "eff" ]
then
 	eff=$(tail -n 1  < efficiency.csv | cut -d',' -f2)
	if (( $(bc <<< "$eff >= 1.0") ))
	then
		echo "-99.0" > eff.val
	else
		echo "$eff" > eff.val
	fi
fi
