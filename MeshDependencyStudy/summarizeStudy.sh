#!/bin/bash
## Extracting Data from subdirectories


echo 'case Elements Efficiency Hemolysis PressureDrop rotationSpeed lastZresidual' > dataSummary.csv
for dirs in */
do
	if [ ${dirs%%?} == "common" ]; then
		continue
	fi
	nelem=$(grep "   Number of cells:" ${dirs}output.txt | cut -d':' -f2 | tr -d " ")
        eff=$(tail -n 1 ${dirs}efficiency.csv | cut -d',' -f2)
        hemo=$(tail -n 1 ${dirs}HIoutlet.csv | cut -d',' -f2)
        pres=$(tail -n 1 ${dirs}presDrop.csv | cut -d',' -f2)
        rotspeed=$(tail -n 1 ${dirs}rotationSpeed.csv | cut -d',' -f2)
	lastzresidual=$( tail -n 4 ${dirs}output.txt | head -n 1 | tr -s " " | cut -d' ' -f5)

        echo "${dirs%%?} $nelem $eff $hemo $pres $rotspeed $lastzresidual" >> dataSummary.csv
done

column -t dataSummary.csv
