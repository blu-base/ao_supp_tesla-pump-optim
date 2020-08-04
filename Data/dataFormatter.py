#!/usr/bin/env python3



# Raw Optimization Archives
par3source = "Population_Opt1.csv"
par5source = "Population_Opt2.csv"

headers5par = ["ID", "Validity", "constrainViolation", "efficiency", "hemolysis", "dsp", "rdin", "rdout", "volconstr", "volsp", "bladenum", "dskth", "diffangle", "diffratio", "tonguerd", "initialspeed", "icemversion", "speed", "presdrop", "torque", "effrep", "hemorep" ]
headers3par = ["ID", "Validity", "constrainViolation", "efficiency", "hemolysis", "dsp", "volconstr", "volsp", "rdout", "rdin", "bladenum", "dskth", "diffangle", "diffratio", "tonguerd", "initialspeed", "icemversion", "speed", "presdrop", "torque", "effrep", "hemorep" ]



import numpy as np
import pandas as pd
import math


#gen=pd.read_csv(par5source ,sep=' ',header=None, skiprows=[0])
par5=pd.read_csv(par5source ,sep=' ',header=None)
par5.columns = headers5par
par3=pd.read_csv(par3source ,sep=' ',header=None)
par3.columns = headers3par


par5=par5.drop(columns=["ID", "constrainViolation", "bladenum", "dskth", "diffangle", "diffratio", "tonguerd", "initialspeed", "icemversion", "speed", "presdrop", "torque", "effrep", "hemorep"])
par3=par3.drop(columns=["ID", "constrainViolation", "bladenum", "dskth", "diffangle", "diffratio", "tonguerd", "initialspeed", "icemversion", "speed", "presdrop", "torque", "effrep", "hemorep"])

par3=par3.rename(columns={"efficiency":"hydralic_efficiency(-)", "hemolysis":"Hemolysis_ratio(-)", "dsp":"gap_width(mm)", "rdin":"inner_disk_radius(mm)", "rdout":"outer_disk_radius(mm)", "volconstr":"volute_constructor_parameter(-)"})
par5=par5.rename(columns={"efficiency":"hydralic_efficiency(-)", "hemolysis":"Hemolysis_ratio(-)", "dsp":"gap_width(mm)", "rdin":"inner_disk_radius(mm)", "rdout":"outer_disk_radius(mm)", "volconstr":"volute_constructor_parameter(-)"})

par3.to_csv("Opt1.csv")
par5.to_csv("Opt2.csv")
