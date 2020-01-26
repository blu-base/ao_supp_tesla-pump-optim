#!/usr/bin/env python3



# Raw Optimization Archives
par5source = "../Opt2/Results/FullArchive/Gen.csv"
par3source = "../Opt1/Results/FullArchive/Gen.csv"

headers5par = ["ID", "Validity", "constrainViolation", "efficiency", "hemolysis", "dsp", "rdin", "rdout", "volconstr", "volsp", "bladenum", "dskth", "diffangle", "diffratio", "tonguerd", "initialspeed", "icemversion", "speed", "presdrop", "torque", "effrep", "hemorep" ]
headers3par = ["ID", "Validity", "constrainViolation", "efficiency", "hemolysis", "dsp", "volconstr", "volsp", "rdout", "rdin", "bladenum", "dskth", "diffangle", "diffratio", "tonguerd", "initialspeed", "icemversion", "speed", "presdrop", "torque", "effrep", "hemorep" ]



import numpy as np
import pandas as pd
import math


## Constants
density = 1065 # kg/m^3
gravity = 9.81 # m/s^2
volumeFlowRateLPM = 5 # liter per minute
volumeFlowRate = volumeFlowRateLPM /60 / 1000 # m^3 / s


#gen=pd.read_csv(par5source ,sep=' ',header=None, skiprows=[0])
par5=pd.read_csv(par5source ,sep=' ',header=None)
par5.columns = headers5par
par3=pd.read_csv(par3source ,sep=' ',header=None)
par3.columns = headers3par

# Full collection
gen = pd.concat([par5,par3],sort=True)



## functions for derived variables

## omega to RPM conversion
def getRPMfromOmega(speed):
    # speed in rad/s
    return speed / (2*math.pi) *60

## Pressure head
def head(presdiff):
    return presdiff/(gravity * density)

## specific energy
def specWork(head):
    return gravity * head

## Circumferential velocity at radius rad
def circumVel(omega, rad):
    return omega * rad

## Meridian velocity
def meridianVel(gap, rd, volFlow):
    return volFlow / ( 2 * math.pi * rd * gap)

## Head coefficient
def druckziffer(head, outerRd, omega):
    return head / ( circumVel(omega, outerRd) ** 2 / ( 2 * gravity))

## Specific speed N_qy
def specificSpeedNqy(speed,volumeFlow, specificWork):
    return speed / (2 * math.pi) * pow(volumeFlow, .5) / pow(specificWork, 3/4)

## Specific speed Nq
def specificSpeedNq(rpm, volumeFlow, head):
    return rpm * pow(volumeFlow, .5) / pow(head, 3/4)

## Specific Speed, according to Balje
def specificSpeedNs(speed, volumeFlow, head):
    return speed * pow(volumeFlow, .5) / pow(gravity * head, 3/4)

## Specific Diameter, according to Balje
def specificDiameterds(rd, head, volFlow):
    return 2*rd  *  pow(gravity * head, .25) / pow(volFlow, .5)

## Tip speed ratio
def schnelllaufzahl(druckziffer,durchflussziffer):
    return pow(durchflussziffer, .5) / pow(druckziffer, 3/4)

## Flow coefficient
def durchflussziffer(velMeridian,velCircum):
    return velMeridian / velCircum

def specificDiameter(rd, head, volFlow):
    return 2*rd  *  pow(head, .25) / pow(volFlow, .5)

def durchmesserzahl(druckziffer,durchflussziffer):
    return pow(druckziffer, .25) / pow(durchflussziffer, .5)

def murataRatio(gapwidth, radiusin, omega, volFlow):
    return radiusin * pow(2*math.pi * gapwidth * omega / volFlow, 0.5)

gen['u1'] = gen.apply(lambda x: circumVel(x['speed'], x['rdin']/1000), axis=1)
gen['u2'] = gen.apply(lambda x: circumVel(x['speed'], x['rdout']/1000), axis=1)

gen['um1'] = gen.apply(lambda x: meridianVel(x['volsp']/1000, x['rdin']/1000, volumeFlowRate), axis=1)
gen['um2'] = gen.apply(lambda x: meridianVel(x['volsp']/1000, x['rdout']/1000, volumeFlowRate), axis=1)

gen['rpm'] = gen['speed'].apply(getRPMfromOmega)
gen['head'] = gen['presdrop'].apply(head)
gen['specWork'] = gen['head'].apply(specWork)

gen['druckziffer'] = gen.apply(lambda x: druckziffer(x['head'], x['rdout']/1000, x['speed']), axis=1)
gen['durchflussziffer'] = gen.apply(lambda x: durchflussziffer(x['um2'], x['u2']), axis=1)
gen['schnelllaufzahl'] = gen.apply(lambda x: schnelllaufzahl(x['druckziffer'], x['durchflussziffer']), axis=1)
gen['durchmesserzahl'] = gen.apply(lambda x: durchmesserzahl(x['druckziffer'], x['durchflussziffer']), axis=1)
gen['nqy'] = gen.apply(lambda x: specificSpeedNqy(x['speed'], volumeFlowRate, x['specWork']), axis=1)
gen['nq'] = gen.apply(lambda x: specificSpeedNqy(x['rpm'], volumeFlowRate, x['head']), axis=1)
gen['ns'] = gen.apply(lambda x: specificSpeedNqy(x['speed'], volumeFlowRate, x['head']), axis=1)

gen['specificDia'] = gen.apply(lambda x: specificDiameter(x['rdout']/1000, x['head'], volumeFlowRate), axis=1)
gen['specificDiaBalje'] = gen.apply(lambda x: specificDiameterds(x['rdout']/1000, x['head'], volumeFlowRate), axis=1)

gen['murataRatio'] = gen.apply(lambda x: murataRatio(x['dsp']/1000, x['rdin']/1000, x['speed'], volumeFlowRate), axis=1)


#===============================================================================
# Save the container
#===============================================================================
gen.to_csv('genAll.csv')
