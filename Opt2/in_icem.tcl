## Copyright 2020 Sebastian Engel
## 
## Permission is hereby granted, free of charge, to any person obtaining a copy
## of this software and associated documentation files (the "Software"), to deal
## in the Software without restriction, including without limitation the rights 
## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
## copies of the Software, and to permit persons to whom the Software is 
## furnished to do so, subject to the following conditions:
## 
## The above copyright notice and this permission notice shall be included in 
## all copies or substantial portions of the Software.
## 
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
## IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
## FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
## AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
## LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
## OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN 
## THE SOFTWARE.
## 
## 
## 
## 
## 
################################################################################ 
## Tcl-Script for ANSYS ICEM CFD
## 
## Part of the Work:
## Multi-objective Design Optimization of a Tesla-type Rotary Blood Pump 
## 
## Creates Geometry and Mesh of a 3D Slice of the Tesla pump impeller and the volute
## Outputs a ANSYS Fluent Mesh (V6)
## 
## written by Sebastian Engel
################################################################################ 

# -----------------------------------------------
# Start clean ICEM Session
# -----------------------------------------------#
ic_unload_tetin
ic_hex_unload_blocking 
ic_unload_mesh 
ic_empty_tetin
ic_boco_clear_icons
ic_csystem_display all 0
ic_csystem_set_current global
ic_set_global geo_cad 0.000001 toler
ic_set_meshing_params global 0 gttol 0.000001 gtrel 1
ic_geo_set_units mm

################################################################################
## General settings
################################################################################

# Debug settings
set debug 0

# Output file name
set outputfile "slice"


################################################################################
## Input Parameter settings
################################################################################
# externalize parameter
source in_parameter.tcl
## Disk settings
# Spacing between disks [mm]
#	set dsksp 0.8
# Thickness of a disk [mm]
#	set dskth 0.5
# Outer disk radius [mm]
#	set dskrdout 45.0
# Inner disk radius [mm]
#	set dskrdin 5.0
## Diffusor settings	
# Diffusor inlet/outlet surface ratio [1]
#	set diffratio 3.0
# Diffusor angle [deg]
#	set diffang 8.0


## Volute settings
# Volute spacing around impeller [mm]
#	set volsp 0.4

# Tongue radius [mm]
#	set tonguerd 0.05 

# # Information for volute construction
# # Volumeflow
#	set volumeflowLA 8.3
# # Velocity Moment K_(G,spir) assuming 3m/s 
# # needs to be an optimization parameter
# velCircum between 1 and 12?
# set velCircum 3.
# set cstKGspir [expr ($dskrdout+$volsp)*$velCircum]

# set volConstructor [expr $volumeflowLA/$cstKGspir]
#	set volConstructor 0.2

# How many segments per quarter construct the volute? [1]
set volresobase 8
set volreso [expr $volresobase*4]

## Derived Quantities	
# (radial) Inlet radius [mm]
# should be smaller than dskrdin
set inlrd [expr $dskrdin - 2.0]

# disk tip radius
set dsktiprad [expr $dskth*0.5]
#domain thickness
set domth [expr $dsksp+$dskth]

################################################################################
## Mesh Parameter settings
################################################################################

# Global mesh scale
set nglob 0.3	

# Boundary layer thickness (global)	
# Should be less than a third of $dsksp
if {$volsp > $dsksp} {
	set blth [expr $dsksp/6.0]	
} else {
	set blth [expr $volsp/6.0]
}


# Boundary Layer Nodes
set nbl 15
# Disk surface Nodes [n/mm]
set ndsksrf 10
# Disk space Nodes [1]
set ndsksp 25
# Inlet region Nodes [n/mm]
set ninlreg 10
# Diffusor space Nodes [n/mm]
set ndiffsp 15
# Diffusor Nodes [n/mm]
set ndifflgh 1.
# Inlet length [n/mm]
set ninllgh 8
# Tangential mesh Nodes per quarter
set ntangential 50

# Global factoring
# Disk thickness Nodes [n/mm]
set ndskth [ceil [expr (1+($dskth+2*$blth)/$dsksp)*$ndsksp*$nglob]]
set ndsksrf [ceil [expr $ndsksrf*$nglob]]
set ndsksp [ceil [expr max($ndsksp*$nglob+1,10)]]
set ndiffsp [ceil [expr $ndiffsp*$nglob]]
set ndifflgh [ceil [expr $ndifflgh*(0.25+0.75*$nglob)]]
set ninllgh [ceil [expr $ninllgh*(0.25+0.75*$nglob)]]
set ndsksrf [ceil [expr $ndsksrf*$nglob]]		
set ninlreg [ceil [expr $ninlreg*(0.2+0.8*$nglob)]]	
set ntangential [ceil [expr $ntangential*(0.25+0.75*$nglob)+1]]		



################################################################################
## ICEM initialisation and parameter Checks
################################################################################

# change into working directory - directory of this script
set wdir [pwd]
set staticdir $wdir

ic_chdir $wdir

# Getting ICEM environment (installation dir)
global env
set icemenv $env(ICEM_ACN)

# Get coordinates from point
proc pcrd {point xyz} {
	return [lindex [ic_vset -vec - $point] $xyz]
}

# cylinder coordinates into cartesian
proc cc {radius phi sincos} {
	set degtorad [expr acos(-1)/180.]
	if {$sincos=="s"} {
		return [expr $radius*sin($phi*$degtorad)]
	}
	if {$sincos=="c"} {
		return [expr $radius*cos($phi*$degtorad)]
	}	
	if {$sincos=="cs"} {
		return [expr $radius*cos($phi*$degtorad)],[expr $radius*sin($phi*$degtorad)]
	}
	if {$sincos=="sc"} {
		return [expr $radius*sin($phi*$degtorad)],[expr $radius*cos($phi*$degtorad)]
	}			
}



# Print output to log, with system time option
proc log [list str [list time 1]] {
	set systemTime [clock seconds]
	if {$time} { 
		ic_mess "[clock format $systemTime -format %H:%M:%S]: $str\n"
	} else {
		ic_mess "$str\n"
	}
}

array set clp {1 2 2 3 3 4 4 1}
array set dlp {1 5 2 6 3 7 4 8}
array set llp {0 1 1 2 2 3 3 0}
set pi [expr acos(-1)]
set degtorad [expr $pi/180.]	
# Helper to set Vertices at right radial distance
proc r {ybase yz index} {
	array set siny {1 1. 2 0. 3 -1. 4 0.}
	array set cosz {1 0. 2 1. 3 0. 4 -1.}	
	if {$yz eq "y"} {
		return [expr $ybase*$siny($index)]
	}
	if {$yz eq "z"} {
		return [expr $ybase*$cosz($index)]
	}
}



# Some log describtions	
log "#################################################" 0
log "# Tesla pump Mesh, written by Sebastian Engel   #" 0
log "#################################################" 0
log "Date: [cmd_date]" 0
log "Current User: [cmd_whoami]"
log "System info: [cmd_uname_a]"
log "Free Memory: [cmd_freemem]"
log "ICEM Path: $icemenv"
log "Current Work Directory: $wdir"
if {[ic_batch_mode]} {log "Script executed in batch mode"} else {log "Script executed in Interactive mode"}
log "#################################################" 0	
log "Input Parameter Output"
log "Disks gap \[mm\]: $dsksp"
log "Disk outer radius \[mm\]: $dskrdout"
log "Disk inner radius \[mm\]: $dskrdin"
log "Disk thickness \[mm\]: $dskth"
log "Inlet radius \[mm\]: $inlrd"
log "Diffusor surface rato \[1\]: $diffratio"
log "Diffusor opening angle \[deg\]: $diffang"
log "Volute spacing \[mm\]: $volsp"
log "#################################################" 0	

# -----------------------------------------------
log "Checking some parameters"
# -----------------------------------------------
if {$dskrdout < $dskrdin} { 
	log "Error! Inner disk radius bigger than outer radius."
	break 
}	
if {$inlrd <= 2.0} { 
	log "Error! Radial inlet radius too small."
	break 
}		
if {[expr $dskrdout-$dskth] < $dskrdin} { 
	log "Error! Check disk radii. Outer disk radius is too small in relation to inner radius."
	break 
}				
if {[expr $dsksp*3.0] < $blth} { 
	log "Error! Check boundary layer thickness. It's very thick. Probably causing mesh errors."
	log "Setting boundary layer thickness to [expr $dsksp/4.0]"
	set blth [expr $dsksp/4.0]
}	
if {$tonguerd < 0.05} { 
	log "Error! Tongue radius is very small. Posssibly causing mesh errors."
	break 
}	
if {$volsp < 0.1} { 
	log "Error! Disk-Tongue spacing is very small. Posssibly causing mesh errors."
	break 
}				

log "Checks completed"




################################################################################
log "Start Building of Geometry"
################################################################################

# Center axis
ic_point {} GEOM p.axis.b "0,0,0"
ic_point {} GEOM p.axis.t "0,0,[expr $dsksp+$dskth]"	
ic_curve point GEOM c.axis "p.axis.b p.axis.t"	

# Creating Slice geometry
ic_point {} INLET p.inl.b "$inlrd,0,0"
ic_point {} INLET p.inl.t "$inlrd,0,[expr $dsksp+$dskth]"


# Disk constructors
ic_point {} DISKS p.dsk.in.b "$dskrdin,0,0"	
ic_point {} DISKS p.dsk.in.t "$dskrdin,0,$domth"	
ic_point {} GEOM p.dsk.in.bl.b "[expr $dskrdin-$blth],0,0"	
ic_point {} GEOM p.dsk.in.bl.t "[expr $dskrdin-$blth],0,$domth"	
ic_point {} GEOM p.dsk.in.archelp.b "[expr $dskrdin+$dsktiprad],0,0"	
ic_point {} GEOM p.dsk.in.archelp.t "[expr $dskrdin+$dsktiprad],0,$domth"	
ic_point {} DISKS p.dsk.in.arc.b "[expr $dskrdin+$dsktiprad],0,$dsktiprad"	
ic_point {} DISKS p.dsk.in.arc.t "[expr $dskrdin+$dsktiprad],0,[expr $dsksp+$dsktiprad]"	

ic_point {} DISKS p.dsk.out.b "$dskrdout,0,0"	
ic_point {} DISKS p.dsk.out.t "$dskrdout,0,[expr $dsksp+$dskth]"	
ic_point {} GEOM p.dsk.out.bl.b "[expr $dskrdout+$blth],0,0"	
ic_point {} GEOM p.dsk.out.bl.t "[expr $dskrdout+$blth],0,[expr $dsksp+$dskth]"		
ic_point {} GEOM p.dsk.out.archelp.b "[expr $dskrdout-$dsktiprad],0,0"	
ic_point {} GEOM p.dsk.out.archelp.t "[expr $dskrdout-$dsktiprad],0,$domth"		
ic_point {} DISKS p.dsk.out.arc.b "[expr $dskrdout-$dsktiprad],0,$dsktiprad"	
ic_point {} DISKS p.dsk.out.arc.t "[expr $dskrdout-$dsktiprad],0,[expr $dsksp+$dsktiprad]"	

# Volute spacing
ic_point {} GEOM p.vol.sp.b "[expr $dskrdout+$volsp],0,0"	
ic_point {} GEOM p.vol.sp.t "[expr $dskrdout+$volsp],0,$domth"

ic_curve point SYM2 c.sym.volsp.t "p.vol.sp.t p.dsk.out.t"
ic_curve point SYM1 c.sym.volsp.b "p.vol.sp.b p.dsk.out.b"	

# Inlet region lines
ic_curve point SYM2 c.sym.inl.t "p.inl.t p.dsk.in.t"
ic_curve point SYM1 c.sym.inl.b "p.inl.b p.dsk.in.b"
ic_curve point INLET c.inl "p.inl.b p.inl.t"	

# Disk lines - Arcs and Rest
ic_curve arc_ctr_rad DISKS c.disk.in.rd.b "p.dsk.in.archelp.b p.dsk.in.b p.dsk.in.arc.b $dsktiprad {} {} 0"
ic_curve arc_ctr_rad DISKS c.disk.in.rd.t "p.dsk.in.archelp.t p.dsk.in.t p.dsk.in.arc.t $dsktiprad {} {} 0"
ic_curve arc_ctr_rad DISKS c.disk.out.rd.b "p.dsk.out.archelp.b p.dsk.out.b p.dsk.out.arc.b $dsktiprad {} {} 0"
ic_curve arc_ctr_rad DISKS c.disk.out.rd.t "p.dsk.out.archelp.t p.dsk.out.t p.dsk.out.arc.t $dsktiprad {} {} 0"

ic_curve point DISKS c.dsk.b "p.dsk.in.arc.b p.dsk.out.arc.b"
ic_curve point DISKS c.dsk.t "p.dsk.in.arc.t p.dsk.out.arc.t"


# Disk Boundary Layer Helper
foreach pnt {p.dsk.out.b p.dsk.out.t p.dsk.in.b p.dsk.in.t p.dsk.out.bl.b p.dsk.out.bl.t p.dsk.in.bl.b p.dsk.in.bl.t} {
	ic_geo_duplicate_set_fam_and_data point $pnt $pnt.help {} _0	
	ic_geo_duplicate_set_fam_and_data point $pnt $pnt.perp {} _0			
	ic_geo_set_part point $pnt.help GEOM 0		
	ic_geo_set_part point $pnt.perp GEOM 0		
}

# Helper points on surface
ic_move_geometry point names p.dsk.out.b.help rotate -45 rotate_axis {0 1 0} cent p.dsk.out.archelp.b
ic_move_geometry point names p.dsk.out.t.help rotate 45 rotate_axis {0 1 0} cent p.dsk.out.archelp.t	
ic_move_geometry point names p.dsk.in.b.help rotate 45 rotate_axis {0 1 0} cent p.dsk.in.archelp.b
ic_move_geometry point names p.dsk.in.t.help rotate -45 rotate_axis {0 1 0} cent p.dsk.in.archelp.t	
# Helper points in BL
ic_move_geometry point names p.dsk.out.bl.b.help rotate -45 rotate_axis {0 1 0} cent p.dsk.out.archelp.b
ic_move_geometry point names p.dsk.out.bl.t.help rotate 45 rotate_axis {0 1 0} cent p.dsk.out.archelp.t		
ic_move_geometry point names p.dsk.in.bl.b.help rotate 45 rotate_axis {0 1 0} cent p.dsk.in.archelp.b
ic_move_geometry point names p.dsk.in.bl.t.help rotate -45 rotate_axis {0 1 0} cent p.dsk.in.archelp.t
ic_move_geometry point names p.dsk.out.bl.b.perp rotate -90 rotate_axis {0 1 0} cent p.dsk.out.archelp.b
ic_move_geometry point names p.dsk.out.bl.t.perp rotate 90 rotate_axis {0 1 0} cent p.dsk.out.archelp.t		
ic_move_geometry point names p.dsk.in.bl.b.perp rotate 90 rotate_axis {0 1 0} cent p.dsk.in.archelp.b
ic_move_geometry point names p.dsk.in.bl.t.perp rotate -90 rotate_axis {0 1 0} cent p.dsk.in.archelp.t	

# -----------------------------------------------
log "Revole Geometry"
# -----------------------------------------------

# Catch all entities
set pointlist [ic_geo_get_objects point]
set crvlist [ic_geo_get_objects curve]

# Define Points which don't need to be revolved 
set ptsToRemove "p.dsk.out.archelp.b p.dsk.out.archelp.t p.dsk.in.archelp.b p.dsk.in.archelp.t"

foreach pnt $pointlist {
	if {[pcrd $pnt 0] == 0} {
		lappend ptsToRemove $pnt
	}
}

# Removing $ptsToRemove from $pointlist
foreach pt $ptsToRemove {
	set id [lsearch $pointlist $pt]
	set pointlist [lreplace $pointlist $id $id]
}

# Revolving Points and creating Arcs, every 90°
foreach pnt $pointlist {
	for {set x 1} {$x <= 4} {incr x} {
		ic_geo_duplicate_set_fam_and_data point $pnt $pnt.$clp($x) {} 
		ic_move_geometry point names $pnt.$clp($x) rotate [expr $x*90] rotate_axis {0 0 1} cent {0 0 0}
	}
	for {set x 1} {$x <= 4} {incr x} {
		ic_curve arc_ctr_rad [ic_geo_get_part point $pnt] c.$pnt.$x "\{0 0 [pcrd $pnt 2] \} $pnt.$x $pnt.$clp($x) {} {} {} 0"
	}
	ic_delete_geometry point $pnt
}

# Revolving Curves without "build topologie", > faster and entity names won't change
foreach crv $crvlist {
	if {$crv != "c.axis"} {
		for {set x 1} {$x <= 4} {incr x} {
			ic_geo_duplicate_set_fam_and_data curve $crv $crv.$clp($x) {} 
			ic_move_geometry curve names $crv.$clp($x) rotate [expr $x*90] rotate_axis {0 0 1} cent {0 0 0} _0
			ic_geo_cre_srf_rev [ic_geo_get_part curve $crv] s.$crv.$x $crv {0 0 0} {0 0 1} [expr ($x-1)*90] [expr $x*90] c 0
		}	
		ic_delete_geometry curve $crv 
	}
}

# -----------------------------------------------
log "Creating Volute geometry"
# -----------------------------------------------
# Determine volute segment division
# degree


set volsegdeg [expr 360./$volreso]
# radian
set volseglgh [expr $volsegdeg*$degtorad]


# # Information for volute construction
# # Volumeflow
# set volumeflowLA 8.3
# # Velocity Moment K_(G,spir) assuming 3m/s 
# # needs to be an optimization parameter
# set cstKGspir [expr ($dskrdout+$volsp)*3.]

# Constructing Volute points
for {set x 1} {$x <= $volreso} {incr x} {
	set volsegrad [expr ($dskrdout+$volsp)*exp($volConstructor/($domth)*$x*$volsegdeg/360.0)]
	set volx [expr cos($volseglgh*$x)*$volsegrad]
	set voly [expr sin($volseglgh*$x)*$volsegrad]
	if {$volConstructor < 0.02} {
		set blvolx [expr cos($volseglgh*$x)*($volsegrad-$blth)]
		set blvoly [expr sin($volseglgh*$x)*($volsegrad-$blth)]
	} else {
		set blvolx [expr cos($volseglgh*$x)*(min(($volsegrad-$dskrdout)/2.,$volsp*3.)+$dskrdout)]
		set blvoly [expr sin($volseglgh*$x)*(min(($volsegrad-$dskrdout)/2.,$volsp*3.)+$dskrdout)]	
	}
	ic_point {} VOLUTE p.vol.wll.b.$x "$volx,$voly,0"			
	ic_point {} VOLUTE p.vol.wll.t.$x "$volx,$voly,$domth"		
	ic_point {} GEOM p.vol.wll.bl.b.$x "$blvolx,$blvoly,0"			
	ic_point {} GEOM p.vol.wll.bl.t.$x "$blvolx,$blvoly,$domth"				
}

# Tongue geometry
ic_curve arc_ctr_rad DIFFUSOR c.vol.tongue.arc.diff.b "p.vol.sp.b.1 \{1. 0 0\} \{0 1. 0\} $tonguerd 185. 270."
ic_curve arc_ctr_rad DIFFUSOR c.vol.tongue.arc.diff.t "p.vol.sp.t.1 \{1. 0 $domth\} \{0 1. $domth\} $tonguerd 185. 270."

ic_curve arc_ctr_rad VOLUTE c.vol.tongue.arc.vol.b "p.vol.sp.b.1 \{1. 0 0\} \{0 1. 0\} $tonguerd 270. 360."
ic_curve arc_ctr_rad VOLUTE c.vol.tongue.arc.vol.t "p.vol.sp.t.1 \{1. 0 $domth\} \{0 1. $domth\} $tonguerd 270. 360."	

ic_point {} VOLUTE p.vol.wll.b.0 "[expr [pcrd p.vol.sp.b.1 0]-$tonguerd],0,0"		
ic_point {} VOLUTE p.vol.wll.t.0 "[expr [pcrd p.vol.sp.t.1 0]-$tonguerd],0,$domth"	
ic_point {} GEOM p.vol.wll.bl.b.0 "[expr [pcrd p.vol.sp.b.1 0]-$tonguerd-$blth],0,0"		
ic_point {} GEOM p.vol.wll.bl.t.0 "[expr [pcrd p.vol.sp.t.1 0]-$tonguerd-$blth],0,$domth"	

ic_point {} DIFFUSOR p.diff.tongue.b "[expr [pcrd p.vol.sp.b.1 0]+[cc $tonguerd 95 s]],[cc $tonguerd 95 c],0"		
ic_point {} DIFFUSOR p.diff.tongue.t "[expr [pcrd p.vol.sp.t.1 0]+[cc $tonguerd 95 s]],[cc $tonguerd 95 c],$domth"	
ic_point {} GEOM p.diff.tongue.bl.b "[expr [pcrd p.vol.sp.b.1 0]+[cc [expr $tonguerd+$blth] 92 s]],[cc [expr $tonguerd+$blth] 92 c],0"		
ic_point {} GEOM p.diff.tongue.bl.t "[expr [pcrd p.vol.sp.t.1 0]+[cc [expr $tonguerd+$blth] 92 s]],[cc [expr $tonguerd+$blth] 92 c],$domth"	

ic_point {} DIFFUSOR p.diff.tongue.tip.b "[pcrd p.vol.sp.b.1 0],-$tonguerd,0"		
ic_point {} DIFFUSOR p.diff.tongue.tip.t "[pcrd p.vol.sp.t.1 0],-$tonguerd,$domth"		
ic_point {} GEOM p.diff.tongue.tip.bl.b.1 "[pcrd p.vol.sp.b.1 0],[expr -$tonguerd-$blth],0"		
ic_point {} GEOM p.diff.tongue.tip.bl.t.1 "[pcrd p.vol.sp.t.1 0],[expr -$tonguerd-$blth],$domth"	


ic_curve point DIFFUSOR c.diff.tongue "p.diff.tongue.b p.diff.tongue.t"	
ic_geo_cre_srf_rev DIFFUSOR s.diff.tongue.1 c.diff.tongue p.vol.sp.b.1 {0 0 1} 5 90 c 0
ic_geo_cre_srf_rev VOLUTE s.diff.tongue.2 c.diff.tongue p.vol.sp.b.1 {0 0 1} 90 180 c 0



# Creating Volute Curves - Connecting points
for {set i 1} {$i <= 4} {incr i} {
	set baselist_t ""
	set baselist_b ""		
	for {set j 0} {$j <= $volresobase} {incr j} {
		lappend baselist_t "p.vol.wll.t.[expr ($i-1)*$volresobase+$j]"		
		lappend baselist_b "p.vol.wll.b.[expr ($i-1)*$volresobase+$j]"
		lappend baselist_bl_t "p.vol.wll.bl.t.[expr ($i-1)*$volresobase+$j]"		
		lappend baselist_bl_b "p.vol.wll.bl.b.[expr ($i-1)*$volresobase+$j]"			
	}
	ic_curve point VOLUTE c.vol.wll.b.$i "$baselist_b"	
	ic_curve point VOLUTE c.vol.wll.t.$i "$baselist_t"	
	ic_curve point GEOM c.vol.wll.bl.b.$i "$baselist_bl_b"	
	ic_curve point GEOM c.vol.wll.bl.t.$i "$baselist_bl_t"			
}

for {set i 1} {$i <= 4} {incr i} {

	ic_curve point VOLUTE c.vol.wll.con.t.$i "p.vol.sp.t.$i p.vol.wll.t.[expr ($i-1)*$volresobase]"	
	ic_curve point VOLUTE c.vol.wll.con.b.$i "p.vol.sp.b.$i p.vol.wll.b.[expr ($i-1)*$volresobase]"	
	ic_curve point VOLUTE c.vol.wll.con.$i "p.vol.wll.t.[expr ($i-1)*$volresobase] p.vol.wll.b.[expr ($i-1)*$volresobase]"			
}	
#		ic_curve point DIFFUSOR c.diff.in.b "p.vol.wll.b.0 p.vol.wll.b.$volreso"	


# Creating Volute Surfaces
for {set i 1} {$i <= 4} {incr i} {
	ic_geo_create_surface_from_curves VOLUTE s.vol.wll.$i 0.001 "c.vol.wll.t.$i c.vol.wll.b.$i" 0		
	ic_geo_create_surface_from_curves SYM1 s.vol.wll.b.$i 0.001 "c.p.vol.sp.b.$i c.vol.wll.b.$i" 0	
	ic_geo_create_surface_from_curves SYM2 s.vol.wll.t.$i 0.001 "c.p.vol.sp.t.$i c.vol.wll.t.$i" 0		
}	


# -----------------------------------------------
log "Creating Diffusor geometry"
# -----------------------------------------------	

set diffradout [expr [pcrd p.vol.wll.b.$volreso 0]]
set diffradbase [expr [pcrd p.vol.wll.b.0 0]]	
#log $diffradbase
set difflgh [expr $diffratio*$diffradbase]
#log "[expr $diffradbase+$diffradbase*sin($diffang/2)]"


# Points and Lines
# ic_point {} DIFFUSOR p.diff.wll.out.b "[expr $diffradout+$difflgh*tan($diffang/2.*$degtorad)],$difflgh,0"	
# ic_point {} DIFFUSOR p.diff.wll.out.t "[expr $diffradout+$difflgh*tan($diffang/2.*$degtorad)],$difflgh,$domth"	
ic_point {} DIFFUSOR p.diff.wll.out.b "[expr $diffradout+$difflgh*tan($diffang*$degtorad)],$difflgh,0"	
ic_point {} DIFFUSOR p.diff.wll.out.t "[expr $diffradout+$difflgh*tan($diffang*$degtorad)],$difflgh,$domth"	
ic_curve point DIFFUSOR c.diff.wll.out.b "p.vol.wll.b.$volreso p.diff.wll.out.b"	
ic_curve point DIFFUSOR c.diff.wll.out.t "p.vol.wll.t.$volreso p.diff.wll.out.t"		


# Diffusor with defined opening angle(diffang) and length (difflgh) including smooth transition helper
# ic_point {} DIFFUSOR p.diff.wll.in.b "[expr $diffradbase+$tonguerd-$difflgh*sin($diffang/2*$degtorad)],$difflgh,0"	
# ic_point {} DIFFUSOR p.diff.wll.in.t "[expr $diffradbase+$tonguerd-$difflgh*sin($diffang/2*$degtorad)],$difflgh,$domth"		
# ic_point {} DIFFUSOR p.diff.wll.in.b.help "[expr $diffradbase+$tonguerd-$difflgh/2.0*sin($diffang/2.*$degtorad)],[expr $difflgh/2.0],0"	
# ic_point {} DIFFUSOR p.diff.wll.in.t.help "[expr $diffradbase+$tonguerd-$difflgh/2.0*sin($diffang/2.*$degtorad)],[expr $difflgh/2.0],$domth"	

ic_point {} DIFFUSOR p.diff.wll.in.b "[expr $diffradbase+$tonguerd],$difflgh,0"	
ic_point {} DIFFUSOR p.diff.wll.in.t "[expr $diffradbase+$tonguerd],$difflgh,$domth"		
ic_point {} DIFFUSOR p.diff.wll.in.b.help "[expr $diffradbase+$tonguerd],[expr $difflgh/2.0],0"	
ic_point {} DIFFUSOR p.diff.wll.in.t.help "[expr $diffradbase+$tonguerd],[expr $difflgh/2.0],$domth"	
ic_curve point DIFFUSOR c.diff.wll.in.b.help "p.diff.wll.in.b.help  p.diff.wll.in.b"	
ic_curve point DIFFUSOR c.diff.wll.in.t.help "p.diff.wll.in.t.help  p.diff.wll.in.t"	


# Make smooth transition from tongue to diffusor by bridging
foreach level {t b} {	
	ic_geo_cre_bridge_crv DIFFUSOR c.diff.tongue.bdg.$level c.vol.tongue.arc.diff.$level c.diff.wll.in.$level.help 1 1 .2 .5
}



# Combine Diffusor edge curve
foreach level {t b} {
	ic_curve concat DIFFUSOR c.diff.wll.in.$level "c.diff.tongue.bdg.$level c.diff.wll.in.$level.help"		
}


# Diffusor BL helper points
# ic_point {} GEOM p.diff.wll.out.bl.b "[expr [pcrd p.diff.wll.out.b 0]-$difflgh/2.0*sin($diffang/2.*$degtorad)*0.1],[pcrd p.diff.wll.out.b 1],0"	
# ic_point {} GEOM p.diff.wll.out.bl.t "[expr [pcrd p.diff.wll.out.t 0]-$difflgh/2.0*sin($diffang/2.*$degtorad)*0.1],[pcrd p.diff.wll.out.t 1],$domth"	
ic_point {} GEOM p.diff.wll.out.bl.b "[expr [pcrd p.diff.wll.out.b 0]-([pcrd p.diff.wll.out.b 0]-[pcrd p.diff.wll.in.b 0])*0.6],[pcrd p.diff.wll.out.b 1],0"	
ic_point {} GEOM p.diff.wll.out.bl.t "[expr [pcrd p.diff.wll.out.t 0]-([pcrd p.diff.wll.out.t 0]-[pcrd p.diff.wll.in.t 0])*0.6],[pcrd p.diff.wll.out.t 1],$domth"	
ic_point {} GEOM p.diff.wll.in.bl.b "[expr [pcrd p.diff.wll.in.b 0]+$difflgh/2.0*sin($diffang/2.*$degtorad)*0.1],[pcrd p.diff.wll.in.b 1],0"	
ic_point {} GEOM p.diff.wll.in.bl.t "[expr [pcrd p.diff.wll.in.t 0]+$difflgh/2.0*sin($diffang/2.*$degtorad)*0.1],[pcrd p.diff.wll.in.t 1],$domth"		





# Outlet Curves and Surfaces
ic_curve point OUTLET c.outlet.wll.out {p.diff.wll.out.b p.diff.wll.out.t}
ic_curve point OUTLET c.outlet.wll.in {p.diff.wll.in.b p.diff.wll.in.t}
ic_curve point OUTLET c.outlet.wll.t {p.diff.wll.out.t p.diff.wll.in.t}
ic_curve point OUTLET c.outlet.wll.b {p.diff.wll.out.b p.diff.wll.in.b}
ic_geo_create_surface_from_curves OUTLET s.outlet 0 "c.outlet.wll.out c.outlet.wll.in c.outlet.wll.t c.outlet.wll.b" 0	





# Get some helper points for tongue blocking
foreach level {t b} {
	set j 0; foreach angle "115 [expr 145-2./$volsp] [expr 155-2./$volsp] 105" {incr j
		ic_point {} GEOM p.block.tip.help.$level.$j p.vol.sp.$level.1+vector([cc [expr 20./$volConstructor] $angle sc],0)
		ic_curve point GEOM c.block.tip.$level.$j "p.diff.tongue.tip.bl.$level.1 p.block.tip.help.$level.$j"
		ic_point intersect GEOM p.block.tip.bl.$level.$j "c.block.tip.$level.$j c.vol.wll.bl.$level.4" tol 0.0001	
		ic_point projcurv GEOM p.block.tip.wll.$level.$j "p.block.tip.bl.$level.$j c.vol.wll.$level.4"
	}
}
# Single tip block setting
foreach level {t b} {
	set j 1; foreach angle {100 270} {incr j	
		ic_point {} GEOM p.diff.tongue.tip.bl.$level.$j p.diff.tongue.tip.$level+vector([cc [expr $blth/2.] $angle sc],0)+vector(0,[expr -$blth/2.],0)
	}
}

# Diffusor BL helper curves	
foreach level {t b} {
	#	ic_curve point GEOM c.diff.wll.in.bl.$level "p.diff.wll.in.bl.$level  p.diff.wll.in.$level.help.bl"		
	ic_curve point GEOM c.diff.wll.out.bl.$level "p.diff.wll.out.bl.$level  p.vol.wll.bl.$level.$volreso"		
}

# Tongue Boundary Layer Helper

ic_curve point GEOM c.vol.tongue.arc.vol.bl.t {p.vol.wll.bl.t.0 p.diff.tongue.tip.bl.t.3 p.diff.tongue.tip.bl.t.2 p.diff.tongue.bl.t}
ic_point projcurv GEOM p.diff.tongue.tip.bl.help.t.1 {p.diff.tongue.tip.bl.t.1 c.vol.tongue.arc.vol.bl.t}
ic_curve split GEOM c.vol.tongue.arc.diff.bl.t {c.vol.tongue.arc.vol.bl.t p.diff.tongue.tip.bl.help.t.1}

ic_curve point GEOM c.vol.tongue.arc.vol.bl.b {p.vol.wll.bl.b.0 p.diff.tongue.tip.bl.b.3 p.diff.tongue.tip.bl.b.2 p.diff.tongue.bl.b}
ic_point projcurv GEOM p.diff.tongue.tip.bl.help.b.1 {p.diff.tongue.tip.bl.b.1 c.vol.tongue.arc.vol.bl.b}
ic_curve split GEOM c.vol.tongue.arc.diff.bl.b {c.vol.tongue.arc.vol.bl.b p.diff.tongue.tip.bl.help.b.1}


# Helper Geometry for Diffusor Boundary layer
foreach level {t b} {	

	#	ic_geo_duplicate_set_fam_and_data curve c.diff.tongue.bdg.$level c.diff.tongue.bdg.bl.$level {} _0
	ic_geo_duplicate_set_fam_and_data point p.diff.wll.in.$level.help p.diff.wll.in.$level.help.bl {} _0		
	ic_move_geometry point names p.diff.wll.in.$level.help.bl translate "[expr $difflgh/2.0*sin($diffang/2.*$degtorad)*0.05] 0 0"	
	ic_curve point GEOM c.diff.wll.in.$level.help.bl "p.diff.wll.in.$level.help.bl  p.diff.wll.in.bl.$level"		
	ic_geo_cre_bridge_crv GEOM c.diff.tongue.bdg.bl.$level c.vol.tongue.arc.diff.bl.$level c.diff.wll.in.$level.help.bl 0 1 .15 .5	
	# ic_geo_set_part curve c.diff.tongue.bdg.bl.$level GEOM 0
	# ic_geo_set_part point p.diff.wll.in.$level.help.bl GEOM 0	
}	

# Diffusor BL helper curves	
foreach level {t b} {
	ic_curve point GEOM c.diff.wll.in.bl.$level "p.diff.wll.in.bl.$level  p.diff.wll.in.$level.help.bl"		
	#	ic_curve point GEOM c.diff.wll.out.bl.$level "p.diff.wll.out.bl.$level  p.vol.wll.bl.$level.$volreso"		
}	

# Combine Diffusor edge curve
foreach level {t b} {
	ic_curve concat GEOM c.diff.wll.in.$level.bl "c.diff.tongue.bdg.bl.$level c.diff.wll.in.$level.help.bl"		

}

# Surfaces of Diffusor
ic_geo_create_surface_from_curves DIFFUSOR s.diff.wll.out 0 "c.diff.wll.out.b c.diff.wll.out.t" 0
ic_geo_create_surface_from_curves DIFFUSOR s.diff.wll.in 0 "c.diff.wll.in.b c.diff.wll.in.t" 0
ic_geo_create_surface_from_curves SYM2 s.diff.wll.t 0 "c.diff.wll.out.t c.diff.wll.in.t" 0
ic_geo_create_surface_from_curves SYM1 s.diff.wll.b 0 "c.diff.wll.out.b c.diff.wll.in.b" 0	


# -----------------------------------------------	
log "Geometry creation completed"	
# -----------------------------------------------	

################################################################################
log "Start Blocking"
################################################################################



# Inits (not all needed, untested)
ic_geo_new_family FLUID
ic_hex_initialize_blocking {} FLUID 0 101
ic_hex_switch_blocking root
ic_hex_unblank_blocks 
ic_hex_multi_grid_level 0
ic_hex_projection_limit 0	
ic_hex_default_bunching_law default 1.2
ic_hex_floating_grid off
ic_hex_transfinite_degree 1
ic_hex_set_mesh_params GEOM INLET DISKS VOLUTE SYM2 SYM1 DIFFUSOR OUTLET FLUID -version 110
ic_hex_error_messages off_minor
ic_hex_switch_blocking root	


ic_hex_mark_blocks unmark
ic_hex_mark_blocks superblock 13
ic_hex_mark_blocks face_neighbors corners { 21 37 25 41 } { 22 38 26 42 }	
ic_hex_ogrid 1 m GEOM INLET DISKS VOLUTE SYM2 SYM1 DIFFUSOR OUTLET FLUID -version 50	


# Building basic Volute Block Structure
ic_hex_mark_blocks unmark
ic_hex_mark_blocks superblock 34
ic_hex_split_grid 25 41 0.7 m GEOM INLET DISKS VOLUTE SYM2 SYM1 DIFFUSOR OUTLET FLUID marked
ic_hex_mark_blocks unmark
ic_hex_delete_blocks numbers 13 39	
ic_hex_extrude_faces 1 { 41 42 77 78 } FLUID 4.0 -nsub 1 -type original FLUID GEOM INLET OUTLET SYM2 SYM1 DIFFUSOR DISKS VOLUTE -version 50
ic_hex_collapse_edge 78 96 join keep_first -version 101 fix_first
ic_hex_collapse_edge 77 90 join keep_first -version 101 fix_first
ic_hex_collapse_edge 133 95 join keep_first -version 101 fix_first
ic_hex_collapse_edge 129 89 join keep_first -version 101 fix_first
ic_hex_extrude_faces 1 {128 132 129 133} FLUID 4.0 -nsub 1 -type original FLUID GEOM INLET OUTLET SYM2 SYM1 DIFFUSOR DISKS VOLUTE -version 50


# Recreating center
ic_hex_create_block FLUID 73 65 77 69 74 66 78 70 GEOM INLET DISKS VOLUTE SYM2 SYM1 DIFFUSOR OUTLET FLUID
ic_hex_mark_blocks unmark	
ic_hex_mark_blocks superblock 112
ic_hex_mark_blocks face_neighbors corners { 66 74 70 78 } { 65 73 69 77 }
ic_hex_ogrid 1 m GEOM INLET DISKS VOLUTE SYM2 SYM1 DIFFUSOR OUTLET FLUID -version 50
#ic_hex_delete_blocks numbers 112
ic_hex_delete_blocks numbers 113



# Moving some nodes to get an easy overview
ic_hex_move_node 21 p.vol.wll.b.[expr $volresobase*2]
ic_hex_move_node 22 p.vol.wll.t.[expr $volresobase*2]
ic_hex_move_node 25 p.vol.wll.b.[expr $volresobase*1]
ic_hex_move_node 26 p.vol.wll.t.[expr $volresobase*1]
ic_hex_move_node 37 p.vol.wll.b.[expr $volresobase*3]
ic_hex_move_node 38 p.vol.wll.t.[expr $volresobase*3]



ic_hex_move_node 77 p.dsk.out.b.1
ic_hex_move_node 78 p.dsk.out.t.1	
ic_hex_move_node 69 p.dsk.out.b.2
ic_hex_move_node 70 p.dsk.out.t.2	
ic_hex_move_node 65 p.dsk.out.b.3
ic_hex_move_node 66 p.dsk.out.t.3
ic_hex_move_node 73 p.dsk.out.b.4
ic_hex_move_node 74 p.dsk.out.t.4







# Preparing boundary layer blocks
ic_hex_split_grid 78 133 0.5 m GEOM INLET DISKS VOLUTE SYM2 SYM1 DIFFUSOR OUTLET FLUID
ic_hex_split_grid 292 300 0.25 m GEOM INLET DISKS VOLUTE SYM2 SYM1 DIFFUSOR OUTLET FLUID




# Diffusor Outlet + boundary layer
ic_hex_move_node 164 p.diff.wll.out.b
ic_hex_move_node 166 p.diff.wll.out.t
ic_hex_move_node 165 p.diff.wll.in.b
ic_hex_move_node 167 p.diff.wll.in.t

ic_hex_move_node 274 p.diff.wll.out.bl.b
ic_hex_move_node 290 p.diff.wll.out.bl.t	
ic_hex_move_node 310 p.diff.wll.in.bl.b
ic_hex_move_node 330 p.diff.wll.in.bl.t





# moving tip block vertices to tip.bl points
ic_hex_move_node 341 p.diff.tongue.tip.bl.t.1
ic_hex_move_node 321 p.diff.tongue.tip.bl.b.1	
ic_hex_move_node 329 p.diff.tongue.tip.bl.t.2
ic_hex_move_node 309 p.diff.tongue.tip.bl.b.2
ic_hex_move_node 277 p.diff.tongue.tip.bl.b.3
ic_hex_move_node 292 p.diff.tongue.tip.bl.t.3	
ic_hex_move_node 133 p.diff.tongue.tip.t
ic_hex_move_node 129 p.diff.tongue.tip.b



# moving tip blocking vertices to volute wall
ic_hex_move_node 300 p.block.tip.bl.t.1	
ic_hex_move_node 340 p.block.tip.bl.t.2
ic_hex_move_node 291 p.block.tip.bl.t.3	
ic_hex_move_node 289 p.block.tip.bl.t.4	

ic_hex_set_node_projection 300 c.vol.wll.bl.t.4	
ic_hex_set_node_projection 340 c.vol.wll.bl.t.4
ic_hex_set_node_projection 291 c.vol.wll.bl.t.4	
ic_hex_set_node_projection 289 c.vol.wll.bl.t.4	

ic_hex_move_node 303 p.block.tip.bl.b.1	
ic_hex_move_node 320 p.block.tip.bl.b.2	
ic_hex_move_node 276 p.block.tip.bl.b.3	
ic_hex_move_node 273 p.block.tip.bl.b.4

ic_hex_set_node_projection 303 c.vol.wll.bl.b.4		
ic_hex_set_node_projection 320 c.vol.wll.bl.b.4		
ic_hex_set_node_projection 276 c.vol.wll.bl.b.4		
ic_hex_set_node_projection 273 c.vol.wll.bl.b.4	


# volute tongue wall boundary layer
ic_hex_move_node 275 p.block.tip.wll.b.1
ic_hex_move_node 293 p.block.tip.wll.t.1	
ic_hex_move_node 311 p.block.tip.wll.b.2
ic_hex_move_node 333 p.block.tip.wll.t.2
ic_hex_move_node 42  p.block.tip.wll.t.3
ic_hex_move_node 41  p.block.tip.wll.b.3	
ic_hex_move_node 132 p.block.tip.wll.t.4
ic_hex_move_node 128 p.block.tip.wll.b.4

ic_hex_move_node 268 p.vol.wll.bl.b.[expr $volresobase*3]
ic_hex_move_node 284 p.vol.wll.bl.t.[expr $volresobase*3]
ic_hex_move_node 286 p.vol.wll.bl.t.[expr $volresobase*2]
ic_hex_move_node 270 p.vol.wll.bl.b.[expr $volresobase*2]
ic_hex_move_node 288 p.vol.wll.bl.t.[expr $volresobase*1]
ic_hex_move_node 272 p.vol.wll.bl.b.[expr $volresobase*1]	





foreach level {t b} {
	ic_point {} GEOM p.block.belowtip.help.$level.1 p.diff.tongue.tip.bl.$level.3+vector([cc [expr 20./$volConstructor] [expr 250+min(6./$volsp,5.)] sc],0)
	ic_point {} GEOM p.block.belowtip.help.$level.2 p.diff.tongue.tip.bl.$level.1+vector([cc [expr 20./$volConstructor] [expr 240-min(6./$volsp,30.)] sc],0)	
	ic_point {} GEOM p.block.belowtip.help.$level.3 p.diff.tongue.tip.bl.$level.3+vector([cc [expr 20./$volConstructor] [expr 250-min(2./$volsp,20.)] sc],0)		

	ic_curve point GEOM c.block.belowtip.help.$level.1 "p.diff.tongue.tip.bl.$level.3 p.block.belowtip.help.$level.1"	
	ic_curve point GEOM c.block.belowtip.help.$level.2 "p.diff.tongue.tip.bl.$level.1 p.block.belowtip.help.$level.2"
	ic_curve point GEOM c.block.belowtip.help.$level.3 "p.diff.tongue.tip.bl.$level.3 p.block.belowtip.help.$level.3"	

	for {set j 1} {$j <=3} {incr j} {			
		ic_point intersect GEOM p.block.belowtip.disk.bl.$level.$j "c.block.belowtip.help.$level.$j c.p.dsk.out.bl.$level.4" tol 0.0001	
		ic_point intersect GEOM p.block.belowtip.disk.$level.$j "c.block.belowtip.help.$level.$j c.p.dsk.out.$level.4" tol 0.0001	
	}
}	

ic_curve point GEOM c.yblock.center.help.t {p.dsk.out.t.help p.block.belowtip.disk.bl.t.2}
ic_curve point GEOM c.yblock.center.help.b {p.dsk.out.b.help p.block.belowtip.disk.bl.b.2}	
ic_point crv_par GEOM p.yblock.center.help.t {c.yblock.center.help.t 0.5}
ic_point crv_par GEOM p.yblock.center.help.b {c.yblock.center.help.b 0.5}	
ic_point projsurf DISKS p.yblock.center.t {p.yblock.center.help.t s.c.disk.out.rd.t.1}
ic_point projsurf DISKS p.yblock.center.b {p.yblock.center.help.b s.c.disk.out.rd.b.1}

ic_point projcurv GEOM p.block.belowtip.disk.bl.proj.t.1 {p.block.belowtip.disk.bl.t.1 c.p.dsk.out.t.perp.4}		
ic_point projcurv GEOM p.block.belowtip.disk.bl.proj.t.2 {p.block.belowtip.disk.bl.t.2 c.p.dsk.out.t.perp.4}
ic_point projcurv GEOM p.block.belowtip.disk.bl.proj.t.3 {p.block.belowtip.disk.bl.t.3 c.p.dsk.out.t.perp.4}

ic_point projcurv GEOM p.block.belowtip.disk.bl.proj.b.1 {p.block.belowtip.disk.bl.b.1 c.p.dsk.out.b.perp.4}		
ic_point projcurv GEOM p.block.belowtip.disk.bl.proj.b.2 {p.block.belowtip.disk.bl.b.2 c.p.dsk.out.b.perp.4}
ic_point projcurv GEOM p.block.belowtip.disk.bl.proj.b.3 {p.block.belowtip.disk.bl.b.3 c.p.dsk.out.b.perp.4}



ic_hex_move_node 312 p.block.belowtip.disk.bl.b.2
ic_hex_move_node 331 p.block.belowtip.disk.bl.t.2
ic_hex_move_node 347 p.block.belowtip.disk.bl.b.3
ic_hex_move_node 344 p.block.belowtip.disk.bl.t.3

ic_hex_move_node 313 p.block.belowtip.disk.bl.b.1
ic_hex_move_node 332 p.block.belowtip.disk.bl.t.1
ic_hex_move_node 77 p.block.belowtip.disk.b.1
ic_hex_move_node 78 p.block.belowtip.disk.t.1	



# separating inlet blockings
ic_hex_split_grid 70 226 0.851687 m GEOM INLET DISKS VOLUTE SYM2 SYM1 DIFFUSOR OUTLET FLUID

ic_hex_move_node 361 p.dsk.in.bl.t.1
ic_hex_move_node 360 p.dsk.in.bl.t.2
ic_hex_move_node 358 p.dsk.in.bl.t.3
ic_hex_move_node 356 p.dsk.in.bl.t.4

ic_hex_move_node 353 p.dsk.in.bl.b.1		
ic_hex_move_node 352 p.dsk.in.bl.b.2
ic_hex_move_node 350 p.dsk.in.bl.b.3
ic_hex_move_node 348 p.dsk.in.bl.b.4








# Disk Boundary layer mesh

# preparation block splits
ic_hex_mark_blocks unmark
ic_hex_mark_blocks superblock 31	
ic_hex_mark_blocks superblock 32
ic_hex_mark_blocks superblock 33
ic_hex_mark_blocks superblock 34	
ic_hex_mark_blocks superblock 62	
ic_hex_mark_blocks superblock 158
ic_hex_mark_blocks superblock 159	
ic_hex_mark_blocks superblock 160	
ic_hex_mark_blocks superblock 161	

ic_hex_split_grid 348 356 0.8 3:0,1 6:0,1 m GEOM INLET DISKS VOLUTE SYM2 SYM1 DIFFUSOR OUTLET FLUID marked


ic_hex_mark_blocks unmark
ic_hex_mark_blocks superblock 31	
ic_hex_mark_blocks superblock 32
ic_hex_mark_blocks superblock 33
ic_hex_mark_blocks superblock 34	
ic_hex_mark_blocks superblock 62	
ic_hex_mark_blocks superblock 158
ic_hex_mark_blocks superblock 159	
ic_hex_mark_blocks superblock 160	
ic_hex_mark_blocks superblock 161	
ic_hex_split_grid 350 379 0.3 3:0,1 6:0,1 m GEOM INLET DISKS VOLUTE SYM2 SYM1 DIFFUSOR OUTLET FLUID marked



# Disk Outlet boundary layer

ic_hex_create_composite {c.p.dsk.out.bl.t.1 c.p.dsk.out.bl.t.2 c.p.dsk.out.bl.t.3 c.p.dsk.out.bl.t.4}	
ic_hex_set_edge_projection 328 332 0 1 c.p.dsk.out.bl.t.1
ic_hex_set_edge_projection 326 328 0 1 c.p.dsk.out.bl.t.1
ic_hex_set_edge_projection 324 326 0 1 c.p.dsk.out.bl.t.1
ic_hex_set_edge_projection 331 324 0 1 c.p.dsk.out.bl.t.1	
ic_hex_set_edge_projection 331 344 0 1 c.p.dsk.out.bl.t.1		

ic_hex_create_composite {c.p.dsk.out.bl.b.1 c.p.dsk.out.bl.b.2 c.p.dsk.out.bl.b.3 c.p.dsk.out.bl.b.4}		
ic_hex_set_edge_projection 313 308 0 1 c.p.dsk.out.bl.b.1
ic_hex_set_edge_projection 308 306 0 1 c.p.dsk.out.bl.b.1
ic_hex_set_edge_projection 306 304 0 1 c.p.dsk.out.bl.b.1
ic_hex_set_edge_projection 304 312 0 1 c.p.dsk.out.bl.b.1	



#	ic_hex_project_to_surface FLUID INLET GEOM OUTLET SYM DIFFUSOR DISKS VOLUTE

# # making disk ogrid top
ic_hex_mark_blocks unmark
ic_hex_mark_blocks face_neighbors corners { 74 356 78 361 } { 66 74 358 356 } { 66 358 70 360 } { 70 78 360 361 } { 70 78 328 332 } { 66 326 70 328 } { 66 74 326 324 } { 78 331 332 344 } { 74 324 78 331 }	
ic_hex_ogrid 1 m GEOM INLET DISKS VOLUTE SYM2 SYM1 DIFFUSOR OUTLET FLUID -version 50
ic_hex_mark_blocks unmark
ic_hex_delete_blocks numbers 216 217 218 219 220 221 222 223 224

# # making disk ogrid bottom
ic_hex_mark_blocks unmark
ic_hex_mark_blocks face_neighbors corners { 65 350 69 352 } { 65 306 69 308 } { 65 73 350 348 } { 65 73 306 304 } { 69 77 352 353 } { 69 77 308 313 } { 73 348 77 353 } { 73 304 77 312 } { 77 312 313 347 }
ic_hex_ogrid 1 m GEOM INLET DISKS VOLUTE SYM2 SYM1 DIFFUSOR OUTLET FLUID -version 50
ic_hex_delete_blocks numbers 160 158 33 32 31 159 34 161 62

# Adding blocking at disk inlet for radius resolution
ic_hex_split_grid 490 496 p.dsk.in.arc.t.2 m GEOM INLET DISKS SYM2 SYM1 VOLUTE DIFFUSOR OUTLET FLUID




#Inlet associations
set vrt_t {234 226 222 230} 
set vrt_b {232 224 220 228} 		
for  {set z 0} {$z < 4} {incr z} {
	ic_hex_set_edge_projection [lindex $vrt_t $z] [lindex $vrt_t $llp($z)] 0 1 c.p.inl.t.[expr $z+1]
	ic_hex_move_node [lindex $vrt_t $z]  "p.inl.t.[expr $z+1]"
	ic_hex_set_edge_projection [lindex $vrt_b $z] [lindex $vrt_b $llp($z)] 0 1 c.p.inl.b.[expr $z+1]
	ic_hex_move_node [lindex $vrt_b $z]  p.inl.b.[expr $z+1]		
}		



# Disk Inlet boundary layer
set vrt_t {419 418 379 371} 
set vrt_b {475 474 435 427} 	
for  {set z 0} {$z < 4} {incr z} {
	ic_hex_set_edge_projection [lindex $vrt_t $z] [lindex $vrt_t $llp($z)] 0 1 c.p.dsk.in.bl.t.help.[expr $z+1]	
	ic_hex_move_node [lindex $vrt_t $z]  p.dsk.in.bl.t.help.[expr $z+1]	
	ic_hex_set_edge_projection [lindex $vrt_b $z] [lindex $vrt_b $llp($z)] 0 1 c.p.dsk.in.bl.b.help.[expr $z+1]				
	ic_hex_move_node [lindex $vrt_b $z]  p.dsk.in.bl.b.help.[expr $z+1]	
}		



# Disk Inlet boundary layer edge in symmetry
set vrt_t {361 360 358 356} 
set vrt_b {353 352 350 348} 	
for  {set z 0} {$z < 4} {incr z} {
	ic_hex_set_edge_projection [lindex $vrt_t $z] [lindex $vrt_t $llp($z)] 0 1 c.p.dsk.in.bl.t.[expr $z+1]	
	ic_hex_move_node [lindex $vrt_t $z]  p.dsk.in.bl.t.[expr $z+1]				
	ic_hex_set_edge_projection [lindex $vrt_b $z] [lindex $vrt_b $llp($z)] 0 1 c.p.dsk.in.bl.b.[expr $z+1]		
	ic_hex_move_node [lindex $vrt_b $z]  p.dsk.in.bl.b.[expr $z+1]				
}



#Disk wall at Inlet
set vrt_t {507 506 481 488} 
set vrt_b {661 660 643 650} 
for  {set z 0} {$z < 4} {incr z} {
	ic_hex_set_edge_projection [lindex $vrt_t $z] [lindex $vrt_t $llp($z)] 0 1 c.p.dsk.in.t.[expr $z+1]	
	ic_hex_move_node [lindex $vrt_t $z]  p.dsk.in.t.[expr $z+1]					
	ic_hex_set_edge_projection [lindex $vrt_b $z] [lindex $vrt_b $llp($z)] 0 1 c.p.dsk.in.b.[expr $z+1]	
	ic_hex_move_node [lindex $vrt_b $z]  p.dsk.in.b.[expr $z+1]				
}	

#Disk radius ogrid at Inlet with helper points
set vrt_t {497 496 478 485} 
set vrt_b {675 674 648 655} 
for  {set z 0} {$z < 4} {incr z} {
	ic_hex_set_edge_projection [lindex $vrt_t $z] [lindex $vrt_t $llp($z)] 0 1 c.p.dsk.in.t.help.[expr $z+1]	
	ic_hex_move_node [lindex $vrt_t $z]  p.dsk.in.t.help.[expr $z+1]			
	ic_hex_set_edge_projection [lindex $vrt_b $z] [lindex $vrt_b $llp($z)] 0 1 c.p.dsk.in.b.help.[expr $z+1]		
	ic_hex_move_node [lindex $vrt_b $z]  p.dsk.in.b.help.[expr $z+1]				
}		



#Disk radius ogrid at Inlet with helper points	boundary layer	
set vrt_t {847 846 843 840} 
set vrt_b {835 834 831 828} 
for  {set z 0} {$z < 4} {incr z} {
	ic_hex_set_edge_projection [lindex $vrt_t $z] [lindex $vrt_t $llp($z)] 0 1 c.p.dsk.in.bl.t.perp.[expr $z+1]	
	ic_hex_move_node [lindex $vrt_t $z]  p.dsk.in.bl.t.perp.[expr $z+1]			
	ic_hex_set_edge_projection [lindex $vrt_b $z] [lindex $vrt_b $llp($z)] 0 1 c.p.dsk.in.bl.b.perp.[expr $z+1]		
	ic_hex_move_node [lindex $vrt_b $z]  p.dsk.in.bl.b.perp.[expr $z+1]				
}	

#Disk radius ogrid at Inlet with helper points		
set vrt_t {849 848 844 841} 
set vrt_b {837 836 832 829} 
for  {set z 0} {$z < 4} {incr z} {
	ic_hex_set_edge_projection [lindex $vrt_t $z] [lindex $vrt_t $llp($z)] 0 1 c.p.dsk.in.arc.t.[expr $z+1]	
	ic_hex_move_node [lindex $vrt_t $z]  p.dsk.in.arc.t.[expr $z+1]			
	ic_hex_set_edge_projection [lindex $vrt_b $z] [lindex $vrt_b $llp($z)] 0 1 c.p.dsk.in.arc.b.[expr $z+1]		
	ic_hex_move_node [lindex $vrt_b $z]  p.dsk.in.arc.b.[expr $z+1]				
}		


# Volute Wall
ic_hex_create_composite {c.vol.wll.t.1 c.vol.wll.t.2 c.vol.wll.t.3 c.vol.wll.t.4 c.vol.tongue.arc.vol.t c.diff.wll.out.t}
ic_hex_set_edge_projection 293 132 0 1 c.vol.wll.t.1
ic_hex_set_edge_projection 333 293 0 1 c.vol.wll.t.1
ic_hex_set_edge_projection 42  333 0 1 c.vol.wll.t.1
ic_hex_set_edge_projection 38  42  0 1 c.vol.wll.t.1
ic_hex_set_edge_projection 22  38  0 1 c.vol.wll.t.1
ic_hex_set_edge_projection 22  26  0 1 c.vol.wll.t.1
ic_hex_set_edge_projection 26  133 0 1 c.vol.wll.t.1

ic_hex_create_composite {c.vol.wll.b.1 c.vol.wll.b.2 c.vol.wll.b.3 c.vol.wll.b.4 c.vol.tongue.arc.vol.b c.diff.wll.out.b}
ic_hex_set_edge_projection 25  129 0 1 c.vol.wll.b.1
ic_hex_set_edge_projection 21  25  0 1 c.vol.wll.b.1
ic_hex_set_edge_projection 21  37  0 1 c.vol.wll.b.1
ic_hex_set_edge_projection 37  41  0 1 c.vol.wll.b.1
ic_hex_set_edge_projection 41  311 0 1 c.vol.wll.b.1
ic_hex_set_edge_projection 311 275 0 1 c.vol.wll.b.1
ic_hex_set_edge_projection 275 128 0 1 c.vol.wll.b.1






log "....."

# Outlet Associations
ic_hex_set_edge_projection 167 330 0 1 c.outlet.wll.t
ic_hex_set_edge_projection 330 290 0 1 c.outlet.wll.t
ic_hex_set_edge_projection 290 166 0 1 c.outlet.wll.t
ic_hex_set_edge_projection 165 310 0 1 c.outlet.wll.b
ic_hex_set_edge_projection 310 274 0 1 c.outlet.wll.b
ic_hex_set_edge_projection 274 164 0 1 c.outlet.wll.b

#Diffusor boundary layer	
ic_hex_create_composite {c.vol.tongue.arc.diff.bl.b c.diff.wll.in.b.bl}
ic_hex_set_edge_projection 309 310 0 1 c.vol.tongue.arc.diff.bl.b
ic_hex_create_composite {c.vol.tongue.arc.diff.bl.t c.diff.wll.in.t.bl}
ic_hex_set_edge_projection 329 330 0 1 c.vol.tongue.arc.diff.bl.t


ic_hex_create_composite {c.vol.wll.bl.t.4 c.diff.wll.out.bl.t}
ic_hex_set_edge_projection 284 291 0 1 c.vol.wll.bl.t.4
ic_hex_set_edge_projection 291 340 0 1 c.vol.wll.bl.t.4
ic_hex_set_edge_projection 340 300 0 1 c.vol.wll.bl.t.4
ic_hex_set_edge_projection 300 289 0 1 c.vol.wll.bl.t.4
ic_hex_set_edge_projection 289 290 0 1 c.vol.wll.bl.t.4

ic_hex_create_composite {c.vol.wll.bl.b.4 c.diff.wll.out.bl.b}
ic_hex_set_edge_projection 268 276 0 1 c.vol.wll.bl.b.4
ic_hex_set_edge_projection 276 320 0 1 c.vol.wll.bl.b.4
ic_hex_set_edge_projection 320 303 0 1 c.vol.wll.bl.b.4
ic_hex_set_edge_projection 303 273 0 1 c.vol.wll.bl.b.4
ic_hex_set_edge_projection 273 274 0 1 c.vol.wll.bl.b.4	
ic_hex_set_edge_projection 286 284 0 1 c.vol.wll.bl.t.4
ic_hex_set_edge_projection 270 268 0 1 c.vol.wll.bl.b.4

# Disk outlet radius associations
ic_hex_create_composite {c.p.dsk.out.t.1 c.p.dsk.out.t.2 c.p.dsk.out.t.3 c.p.dsk.out.t.4 }
ic_hex_set_edge_projection 500 504 0 1 c.p.dsk.out.t.1
ic_hex_set_edge_projection 499 504 0 1 c.p.dsk.out.t.1
ic_hex_set_edge_projection 486 499 0 1 c.p.dsk.out.t.1
ic_hex_set_edge_projection 498 500 0 1 c.p.dsk.out.t.1
ic_hex_set_edge_projection 479 498 0 1 c.p.dsk.out.t.1
ic_hex_set_edge_projection 479 486 0 1 c.p.dsk.out.t.1	

ic_hex_create_composite {c.p.dsk.out.b.1 c.p.dsk.out.b.4 c.p.dsk.out.b.2 c.p.dsk.out.b.3}
ic_hex_set_edge_projection 658 659 0 1 c.p.dsk.out.b.1
ic_hex_set_edge_projection 657 659 0 1 c.p.dsk.out.b.1
ic_hex_set_edge_projection 649 657 0 1 c.p.dsk.out.b.1
ic_hex_set_edge_projection 656 658 0 1 c.p.dsk.out.b.1
ic_hex_set_edge_projection 642 656 0 1 c.p.dsk.out.b.1
ic_hex_set_edge_projection 642 649 0 1 c.p.dsk.out.b.1	


ic_hex_create_composite {c.p.dsk.out.t.help.1 c.p.dsk.out.t.help.2 c.p.dsk.out.t.help.3 c.p.dsk.out.t.help.4 }
ic_hex_set_edge_projection 492 494 0 1 c.p.dsk.out.t.help.1
ic_hex_set_edge_projection 484 493 0 1 c.p.dsk.out.t.help.1
ic_hex_set_edge_projection 477 492 0 1 c.p.dsk.out.t.help.1
ic_hex_set_edge_projection 477 484 0 1 c.p.dsk.out.t.help.1

ic_hex_create_composite {c.p.dsk.out.b.help.1 c.p.dsk.out.b.help.4 c.p.dsk.out.b.help.2 c.p.dsk.out.b.help.3}
ic_hex_set_edge_projection 670 672 0 1 c.p.dsk.out.b.help.1
ic_hex_set_edge_projection 654 671 0 1 c.p.dsk.out.b.help.1
ic_hex_set_edge_projection 647 670 0 1 c.p.dsk.out.b.help.1
ic_hex_set_edge_projection 647 654 0 1 c.p.dsk.out.b.help.1	


ic_hex_create_composite {c.p.dsk.out.arc.t.1 c.p.dsk.out.arc.t.2 c.p.dsk.out.arc.t.3 c.p.dsk.out.arc.t.4}
ic_hex_set_edge_projection 490 491 0 1 c.p.dsk.out.arc.t.1
ic_hex_set_edge_projection 491 483 0 1 c.p.dsk.out.arc.t.1
ic_hex_set_edge_projection 476 483 0 1 c.p.dsk.out.arc.t.1
ic_hex_set_edge_projection 476 490 0 1 c.p.dsk.out.arc.t.1	

ic_hex_create_composite {c.p.dsk.out.arc.b.1 c.p.dsk.out.arc.b.4 c.p.dsk.out.arc.b.2 c.p.dsk.out.arc.b.3}
ic_hex_set_edge_projection 668 669 0 1 c.p.dsk.out.arc.b.1
ic_hex_set_edge_projection 653 669 0 1 c.p.dsk.out.arc.b.1
ic_hex_set_edge_projection 646 668 0 1 c.p.dsk.out.arc.b.1
ic_hex_set_edge_projection 646 653 0 1 c.p.dsk.out.arc.b.1	


# Disk outlet boundary layer
ic_hex_create_composite {c.p.dsk.out.bl.b.help.1 c.p.dsk.out.bl.b.help.2 c.p.dsk.out.bl.b.help.3 c.p.dsk.out.bl.b.help.4 }
ic_hex_set_edge_projection 426 467 0 1 c.p.dsk.out.bl.b.help.1
ic_hex_set_edge_projection 463 468 0 1 c.p.dsk.out.bl.b.help.1
ic_hex_set_edge_projection 434 426 0 1 c.p.dsk.out.bl.b.help.1
ic_hex_set_edge_projection 434 463 0 1 c.p.dsk.out.bl.b.help.1	

ic_hex_create_composite {c.p.dsk.out.bl.t.help.4 c.p.dsk.out.bl.t.help.1 c.p.dsk.out.bl.t.help.2 c.p.dsk.out.bl.t.help.3}
ic_hex_set_edge_projection 378 370 0 1 c.p.dsk.out.bl.t.help.4
ic_hex_set_edge_projection 378 407 0 1 c.p.dsk.out.bl.t.help.4
ic_hex_set_edge_projection 407 412 0 1 c.p.dsk.out.bl.t.help.4
ic_hex_set_edge_projection 370 411 0 1 c.p.dsk.out.bl.t.help.4


ic_hex_create_composite {c.p.dsk.out.bl.t.perp.1 c.p.dsk.out.bl.t.perp.2 c.p.dsk.out.bl.t.perp.3 c.p.dsk.out.bl.t.perp.4}
ic_hex_set_edge_projection 381 388 0 1 c.p.dsk.out.bl.t.perp.1
ic_hex_set_edge_projection 365 388 0 1 c.p.dsk.out.bl.t.perp.1
ic_hex_set_edge_projection 373 365 0 1 c.p.dsk.out.bl.t.perp.1
ic_hex_set_edge_projection 373 381 0 1 c.p.dsk.out.bl.t.perp.1	

ic_hex_create_composite {c.p.dsk.out.bl.b.perp.1 c.p.dsk.out.bl.b.perp.2 c.p.dsk.out.bl.b.perp.3 c.p.dsk.out.bl.b.perp.4}
ic_hex_set_edge_projection 437 444 0 1 c.p.dsk.out.bl.b.perp.1
ic_hex_set_edge_projection 429 437 0 1 c.p.dsk.out.bl.b.perp.1
ic_hex_set_edge_projection 429 421 0 1 c.p.dsk.out.bl.b.perp.1
ic_hex_set_edge_projection 421 444 0 1 c.p.dsk.out.bl.b.perp.1


ic_hex_create_composite {c.vol.wll.bl.t.4 c.vol.tongue.arc.vol.bl.t}	
ic_hex_create_composite {c.vol.wll.bl.b.3 c.vol.tongue.arc.vol.bl.b}
ic_hex_set_edge_projection 288 292 0 1 c.vol.wll.bl.t.4
ic_hex_set_edge_projection 272 277 0 1 c.vol.wll.bl.b.3
ic_hex_set_edge_projection 270 272 0 1 c.vol.wll.bl.b.2
ic_hex_set_edge_projection 286 288 0 1 c.vol.wll.bl.t.2


ic_hex_move_node 498 p.dsk.out.t.2
ic_hex_move_node 328 p.dsk.out.bl.t.2
ic_hex_move_node 492 p.dsk.out.t.help.2
ic_hex_move_node 407 p.dsk.out.bl.t.help.2
ic_hex_move_node 490 p.dsk.out.arc.t.2
ic_hex_move_node 381 p.dsk.out.bl.t.perp.2
ic_hex_move_node 656 p.dsk.out.b.2
ic_hex_move_node 308 p.dsk.out.bl.b.2
ic_hex_move_node 670 p.dsk.out.b.help.2
ic_hex_move_node 668 p.dsk.out.arc.b.2
ic_hex_move_node 463 p.dsk.out.bl.b.help.2
ic_hex_move_node 437 p.dsk.out.bl.b.perp.2
ic_hex_move_node 642 p.dsk.out.b.3
ic_hex_move_node 306 p.dsk.out.bl.b.3
ic_hex_move_node 647 p.dsk.out.b.help.3
ic_hex_move_node 434 p.dsk.out.bl.b.help.3
ic_hex_move_node 646 p.dsk.out.arc.b.3
ic_hex_move_node 429 p.dsk.out.bl.b.perp.3
ic_hex_move_node 479 p.dsk.out.t.perp.3
ic_hex_move_node 326 p.dsk.out.bl.t.3
ic_hex_move_node 477 p.dsk.out.t.help.3
ic_hex_move_node 476 p.dsk.out.arc.t.3
ic_hex_move_node 373 p.dsk.out.bl.t.perp.3
ic_hex_move_node 378 p.dsk.out.bl.t.help.3
ic_hex_move_node 649 p.dsk.out.b.perp.4
ic_hex_move_node 304 p.dsk.out.bl.b.4
ic_hex_move_node 654 p.dsk.out.b.help.4
ic_hex_move_node 426 p.dsk.out.bl.b.help.4
ic_hex_move_node 653 p.dsk.out.arc.b.4
ic_hex_move_node 421 p.dsk.out.bl.b.perp.4
ic_hex_move_node 486 p.dsk.out.t.perp.4
ic_hex_move_node 484 p.dsk.out.t.help.4
ic_hex_move_node 483 p.dsk.out.arc.t.4
ic_hex_move_node 365 p.dsk.out.bl.t.perp.4
ic_hex_move_node 370 p.dsk.out.bl.t.help.4
ic_hex_move_node 324 p.dsk.out.bl.t.4

log "............"	

ic_point projcurv GEOM p.yblock.south.t {p.yblock.center.t c.p.dsk.out.arc.t.4}
ic_point projcurv GEOM p.yblock.south.b {p.yblock.center.b c.p.dsk.out.arc.b.4}


# Y block center  vertices	
ic_point {} GEOM p.yblock.center.center.t p.yblock.south.t+(p.block.belowtip.disk.bl.proj.t.3-p.yblock.south.t)*0.55
ic_point {} GEOM p.yblock.center.center.b p.yblock.south.b+(p.block.belowtip.disk.bl.proj.b.3-p.yblock.south.b)*0.55

ic_point projsurf GEOM p.yblock.center.center.proj.t {p.yblock.center.center.t s.c.disk.out.rd.t.1}	
ic_point projsurf GEOM p.yblock.center.center.proj.b {p.yblock.center.center.b s.c.disk.out.rd.b.1}		




# Yblock south corner vertices

ic_point projcurv GEOM p.yblock.south.center.t {p.yblock.center.center.proj.t c.p.dsk.out.arc.t.4}
ic_point projcurv GEOM p.yblock.south.center.b {p.yblock.center.center.proj.b c.p.dsk.out.arc.b.4}	
ic_hex_move_node 491 p.yblock.south.center.t
ic_hex_move_node 669 p.yblock.south.center.b		

# Yblock south corner boundary layer vertices	
ic_point {} GEOM p.yblock.south.center.bl.t p.yblock.south.center.t+vector(0,0,-$blth)
ic_point {} GEOM p.yblock.south.center.bl.b p.yblock.south.center.b+vector(0,0,$blth)
ic_hex_move_node 388 p.yblock.south.center.bl.t
ic_hex_move_node 444 p.yblock.south.center.bl.b	






ic_hex_move_node 499 p.block.belowtip.disk.bl.proj.t.2
ic_hex_move_node 504 p.block.belowtip.disk.bl.proj.t.3
ic_hex_move_node 500 p.block.belowtip.disk.bl.proj.t.1

ic_hex_move_node 657 p.block.belowtip.disk.bl.proj.b.2
ic_hex_move_node 659 p.block.belowtip.disk.bl.proj.b.3
ic_hex_move_node 658 p.block.belowtip.disk.bl.proj.b.1	


# ic_hex_move_node 671 p.block.belowtip.disk.bl.b.2
# ic_hex_move_node 673 p.block.belowtip.disk.bl.b.3
# ic_hex_move_node 672 p.block.belowtip.disk.bl.b.1	

# Yblock

log "....................................."	



ic_point projcurv GEOM p.yblock.belowtip.disk.bl.proj.help.t.1 {p.block.belowtip.disk.bl.proj.t.1 c.p.dsk.out.t.help.4}
ic_point projcurv GEOM p.yblock.belowtip.disk.bl.proj.help.b.1 {p.block.belowtip.disk.bl.proj.b.1 c.p.dsk.out.b.help.4}

ic_point {} GEOM p.yblock.belowtip.disk.bl.proj.t.1 p.yblock.belowtip.disk.bl.proj.help.t.1+vector(0,[expr -$blth/10.],0)
ic_point {} GEOM p.yblock.belowtip.disk.bl.proj.b.1 p.yblock.belowtip.disk.bl.proj.help.b.1+vector(0,[expr -$blth/10.],0)


ic_hex_move_node 494 p.yblock.belowtip.disk.bl.proj.t.1
ic_hex_move_node 672 p.yblock.belowtip.disk.bl.proj.b.1

# Y block center Boundary layer vertices		
ic_point {} GEOM p.yblock.center.center.proj.bl.t p.yblock.belowtip.disk.bl.proj.t.1+vector([expr $blth/sqrt(2.0)],[expr -$blth/10.],[expr -$blth/sqrt(2.0)])
ic_point {} GEOM p.yblock.center.center.proj.bl.b p.yblock.belowtip.disk.bl.proj.b.1+vector([expr $blth/sqrt(2.0)],[expr -$blth/10.],[expr $blth/sqrt(2.0)])
ic_hex_move_node 412 p.yblock.center.center.proj.bl.t
ic_hex_move_node 468 p.yblock.center.center.proj.bl.b
ic_hex_set_node_projection 412 c.p.dsk.out.bl.t.help.1		
ic_hex_set_node_projection 468 c.p.dsk.out.bl.b.help.1		

# Y block east edge center vertices	
ic_point {} GEOM p.yblock.east.center.help.t p.yblock.south.center.t+(p.block.belowtip.disk.bl.proj.t.2-p.yblock.south.center.t)*0.5
ic_point {} GEOM p.yblock.east.center.help.b p.yblock.south.center.b+(p.block.belowtip.disk.bl.proj.b.2-p.yblock.south.center.b)*0.5
ic_point projsurf GEOM p.yblock.east.center.t {p.yblock.east.center.help.t s.c.disk.out.rd.t.1}
ic_point projsurf GEOM p.yblock.east.center.b {p.yblock.east.center.help.b s.c.disk.out.rd.b.1}
ic_hex_move_node 493 p.yblock.east.center.t	
ic_hex_move_node 671 p.yblock.east.center.b	


# Y block east edge center Boundary layer vertices
ic_point {} GEOM p.yblock.east.center.bl.t p.yblock.east.center.t+(p.yblock.east.center.t-p.yblock.east.center.help.t)*[expr $blth/[ic_vcalc distance p.yblock.east.center.help.t p.yblock.east.center.t]]
ic_point {} GEOM p.yblock.east.center.bl.b p.yblock.east.center.b+(p.yblock.east.center.b-p.yblock.east.center.help.b)*[expr $blth/[ic_vcalc distance p.yblock.east.center.help.b p.yblock.east.center.b]]	
ic_hex_move_node 411 p.yblock.east.center.bl.t	
ic_hex_move_node 467 p.yblock.east.center.bl.b	
ic_hex_set_node_projection 411 c.p.dsk.out.bl.t.help.1	
ic_hex_set_node_projection 467 c.p.dsk.out.bl.b.help.1

ic_point {} GEOM p.yblock.centerfin.help.t p.yblock.south.center.t+(p.block.belowtip.disk.bl.t.3-p.yblock.south.center.t)*0.55
ic_point {} GEOM p.yblock.centerfin.help.b p.yblock.south.center.b+(p.block.belowtip.disk.bl.b.3-p.yblock.south.center.b)*0.55

ic_point projsurf GEOM p.yblock.centerfin.t {p.yblock.centerfin.help.t s.c.disk.out.rd.t.1}
ic_point projsurf GEOM p.yblock.centerfin.b {p.yblock.centerfin.help.b s.c.disk.out.rd.b.1}

ic_hex_move_node 495 p.yblock.centerfin.t
ic_hex_move_node 673 p.yblock.centerfin.b	

ic_point {} GEOM p.yblock.centerfin.bl.t p.yblock.centerfin.t+vector([expr $blth],0,[expr -$blth/sqrt(2.0)])
ic_point {} GEOM p.yblock.centerfin.bl.b p.yblock.centerfin.b+vector([expr $blth],0,[expr $blth/sqrt(2.0)])

ic_hex_move_node 417 p.yblock.centerfin.bl.t
ic_hex_move_node 473 p.yblock.centerfin.bl.b

# Diffusor Associations
#outer wall
ic_hex_set_edge_projection 132 166 0 1 c.vol.wll.t.1	
ic_hex_set_edge_projection 128 164 0 1 c.vol.wll.b.1	
# inner wall
ic_hex_create_composite { c.diff.wll.in.t c.vol.tongue.arc.diff.t}	
ic_hex_set_edge_projection 133 167 0 1 c.diff.wll.in.t
ic_hex_create_composite {c.diff.wll.in.b c.vol.tongue.arc.diff.b }	
ic_hex_set_edge_projection 129 165 0 1 c.diff.wll.in.b


################################################################################
log "Mesh settings"
################################################################################	


# General settings

# Axial mesh count
ic_hex_set_mesh 427 371 n $ndsksp h1rel 0.0 h2rel 0.0 r1 1.2 r2 1.2 lmax 0 default copy_to_parallel unlocked

# Disk Boundary layer settings
ic_hex_set_mesh 353 661 n $nbl h1rel 0.0 h2rel 0.03 r1 1.2 r2 1.05 lmax 0 default copy_to_parallel unlocked
ic_hex_set_mesh 361 507 n $nbl h1rel 0.0 h2rel 0.03 r1 1.2 r2 1.05 lmax 0 default copy_to_parallel unlocked

# Volute/Diffusor Boundary layer settings
ic_hex_set_mesh 292 133 n $nbl h1rel 0.0 h2rel 0.03 r1 1.2 r2 1.2 lmax 0 default copy_to_parallel unlocked
ic_hex_set_mesh 133 329 n $nbl h1rel 0.03 h2rel 0.0 r1 1.2 r2 1.2 lmax 0 default copy_to_parallel unlocked


# Diffusor length + Tongue
ic_hex_set_mesh 329 330 n [ceil [expr $difflgh*$ndifflgh]] h1rel linked 329 341 h2rel 0.0 r1 1.05 r2 1.2 lmax 0 default copy_to_parallel unlocked

# Tongue mesh + Disk
ic_hex_set_mesh 26 133 n [ceil [expr 100*($nglob*0.25+0.75)]] h1rel 0 h2rel linked 292 341 r1 1.2 r2 1.055 lmax 1e+010 default copy_to_parallel unlocked

# ic_hex_set_mesh 26 133 n 100 h1rel 0.0410849763966 h2rel linked 292 341 r1 1.2 r2 1.08 lmax 9.96142e+009 default unlocked
# ic_hex_set_mesh 288 292 n 100 h1rel 0.0410849763966 h2rel linked 292 341 r1 1.2 r2 1.08 lmax 9.96142e+009 default unlocked
# ic_hex_set_mesh 25 129 n 100 h1rel 0.0410849763966 h2rel linked 292 341 r1 1.2 r2 1.08 lmax 9.96142e+009 default unlocked
# ic_hex_set_mesh 272 277 n 100 h1rel 0.0410849763966 h2rel linked 292 341 r1 1.2 r2 1.08 lmax 9.96142e+009 default unlocked	

# Volute before tongue
ic_hex_set_mesh 324 331 n [ceil [expr 100*($nglob*0.25+0.75)]] h1rel 0.0 h2rel 0.01 r1 1.2 r2 1.2 lmax 0 default copy_to_parallel unlocked

# Volute edges above and before tongue in outlet direction
ic_hex_set_mesh 344 341 n [ceil [expr 45*($nglob*0.25+0.75)]] h2rel linked 277 129 r2 1.04
ic_hex_set_mesh 313 277 h2rel linked 277 129 r2 1.04
ic_hex_set_mesh 340 300 h2rel linked 277 129 r2 1.04
ic_hex_set_mesh 333 293 h2rel linked 277 129 r2 1.04
ic_hex_set_mesh 320 303 h2rel linked 277 129 r2 1.04
ic_hex_set_mesh 311 275 h2rel linked 277 129 r2 1.04
ic_hex_set_mesh 332 292 h2rel linked 277 129 r2 1.04
ic_hex_set_mesh 347 321 h2rel linked 277 129 r2 1.04


# Volute width mesh
ic_hex_set_mesh 331 291 n [ceil [expr 40*($nglob*0.25+0.75)]] h1rel linked 331 499 h2rel linked 291 42 r1 1.2 r2 1.2 lmax 0 default copy_to_parallel unlocked
ic_hex_renew
ic_hex_set_mesh 329 289 n [ceil [expr 40*($nglob*0.25+0.75)]] h1rel linked 331 499 h2rel linked 291 42 r1 1.2 r2 1.2 lmax 0 default unlocked
ic_hex_set_mesh 310 274 n [ceil [expr 40*($nglob*0.25+0.75)]] h1rel linked 331 499 h2rel linked 291 42 r1 1.2 r2 1.2 lmax 0 default unlocked
ic_hex_set_mesh 330 290 n [ceil [expr 40*($nglob*0.25+0.75)]] h1rel linked 331 499 h2rel linked 291 42 r1 1.2 r2 1.2 lmax 0 default unlocked

# Disk radial mesh
ic_hex_set_mesh 491 849 n [ceil [expr 45*($nglob*0.25+0.75)]] h1rel linked 444 468 h2rel linked 835 475 r1 1.2 r2 1.2 lmax 0 default copy_to_parallel unlocked


# Inlet mesh count
ic_hex_set_mesh 361 234 n [ceil [expr 15*($nglob*0.25+0.75)]] h1rel 0.0 h2rel 0.0 r1 1.2 r2 1.2 lmax 0 default copy_to_parallel unlocked

#Volute tangential mesh
ic_hex_set_mesh 483 491 n [ceil [expr 100*($nglob*0.25+0.75)]] h1rel 0.0 h2rel linked 491 490 r1 1.2 r2 1.2 lmax 0 default copy_to_parallel unlocked

ic_hex_set_mesh 476 483 h2rel linked 483 491 default copy_to_parallel unlocked
ic_hex_set_mesh 476 490 h2rel linked 490 491 default copy_to_parallel unlocked	



# ic_hex_set_mesh_params INLET DISKS SYM2 VOLUTE DIFFUSOR OUTLET FLUID -version 110	
ic_hex_set_mesh_params INLET DISKS SYM2 VOLUTE DIFFUSOR OUTLET FLUID fix_counts -version 110	
ic_hex_compute_mesh_size GEOM INLET DISKS SYM2 VOLUTE DIFFUSOR OUTLET FLUID



# Volute tongue mesh correction
ic_hex_set_mesh 25 129 h2rel linked 292 341 r1 1.2 r2 1.05  default copy_to_parallel unlocked





# Disk tangential mesh off quarters
ic_hex_set_mesh 476 490 n [ceil [expr 45*($nglob*0.25+0.75)]] h1rel 0.0 h2rel 0.0 r1 1.2 r2 1.2 lmax 0 default copy_to_parallel unlocked
ic_hex_set_mesh 476 483 n [ceil [expr 45*($nglob*0.25+0.75)]] h1rel 0.0 h2rel 0.0 r1 1.2 r2 1.2 lmax 0 default copy_to_parallel unlocked

# Disk inlet radius
ic_hex_set_mesh 427 348 n [ceil [expr $ndsksp/2.+3]] h1rel 0.0 h2rel 0.0 r1 1.2 r2 1.2 lmax 0 default copy_to_parallel unlocked	
ic_hex_set_mesh 485 488 n [ceil [expr $ndsksp/2.+3]] h1rel 0.0 h2rel 0.0 r1 1.2 r2 1.2 lmax 0 default copy_to_parallel unlocked
ic_hex_set_mesh 841 485 n [ceil [expr $ndsksp/2.+3]] h1rel 0.0 h2rel 0.0 r1 1.2 r2 1.2 lmax 0 default copy_to_parallel unlocked

# Disk outlet radius
ic_hex_set_mesh 378 326 n [ceil [expr $ndsksp/2.+3]] h1rel 0.0 h2rel 0.0 r1 1.2 r2 1.2 lmax 0 default copy_to_parallel unlocked
ic_hex_set_mesh 642 647 n [ceil [expr $ndsksp/2.+3]] h1rel 0.0 h2rel 0.0 r1 1.2 r2 1.2 lmax 0 default copy_to_parallel unlocked


# Yblock fixings
ic_hex_set_mesh 500 504 h1rel 0.0 h2rel 0.0 
ic_hex_set_mesh 412 417 h1rel 0.0 h2rel 0.0 
ic_hex_set_mesh 332 344 h1rel 0.0 h2rel 0.0 
ic_hex_set_mesh 468 473 h1rel 0.0 h2rel 0.0 
ic_hex_set_mesh 672 673 h1rel 0.0 h2rel 0.0 
ic_hex_set_mesh 658 659 h1rel 0.0 h2rel 0.0 
ic_hex_set_mesh 313 347 h1rel 0.0 h2rel 0.0 	

ic_hex_set_mesh 656 658 h2rel linked 313 347 r2 1.04
ic_hex_set_mesh 670 672 h2rel linked 313 347 r2 1.04
ic_hex_set_mesh 463 468 h2rel linked 313 347 r2 1.04
ic_hex_set_mesh 328 332 h2rel linked 313 347 r2 1.04
ic_hex_set_mesh 407 412 h2rel linked 313 347 r2 1.04
ic_hex_set_mesh 492 494 h2rel linked 313 347 r2 1.04
ic_hex_set_mesh 498 500 h2rel linked 313 347 r2 1.04
ic_hex_set_mesh 668 669 h2rel linked 313 347 r2 1.04
ic_hex_set_mesh 437 444 h2rel linked 313 347 r2 1.04
ic_hex_set_mesh 381 388 h2rel linked 313 347 r2 1.04
ic_hex_set_mesh 490 491 h2rel linked 313 347 r2 1.04
ic_hex_set_mesh 308 313 h2rel linked 313 347 r2 1.04	



#Inlet Sym Correction
ic_hex_project_face node_numbers { 228 232 348 353 } { 348 353 650 661 } -to_surface s.c.sym.inl.b.1
ic_hex_project_face node_numbers { 220 350 228 348 } { 348 650 350 643 } -to_surface s.c.sym.inl.b.2
ic_hex_project_face node_numbers { 220 224 350 352 } { 350 352 643 660 } -to_surface s.c.sym.inl.b.3
ic_hex_project_face node_numbers { 224 352 232 353 } { 352 660 353 661 } -to_surface s.c.sym.inl.b.4	

# Fixing inverted blocks
ic_hex_mark_blocks unmark
ic_hex_mark_blocks inverted
ic_hex_invert_super_block marked	
ic_hex_renew

# Fixing inverted blocks #2
ic_hex_mark_blocks unmark
ic_hex_mark_blocks inverted
ic_hex_invert_super_block marked	

ic_hex_renew


# Fixing inverted blocks #3
ic_hex_mark_blocks unmark
ic_hex_mark_blocks inverted
ic_hex_invert_super_block marked	

if {$volConstructor < 0.1} {
	ic_hex_set_edge_projection 273 274 0 2 0
	ic_hex_set_edge_projection 289 290 0 2 0
}

if {$volConstructor < 0.02 } {
	ic_hex_set_edge_projection 309 310 0 2 0
	ic_hex_set_edge_projection 329 330 0 2 0
}

ic_hex_set_mesh_params GEOM INLET DISKS SYM2 SYM1 VOLUTE DIFFUSOR OUTLET FLUID fix_counts -version 110
ic_hex_compute_mesh_size GEOM INLET DISKS SYM2 SYM1 VOLUTE DIFFUSOR OUTLET FLUID	


## Hexa Mesh
ic_hex_create_mesh INLET DISKS SYM2 SYM1 VOLUTE DIFFUSOR OUTLET FLUID proj 2 dim_to_mesh 3 nproc 4


log "Meshing complete"	

if {$debug == 0} {
	################################################################################
	log "Try Mesh Smoothing"
	################################################################################

	set determinante [split [ic_hex_rh 20 FLUID proj 2 minval 0 -type determinant_27 maxval 1 new_format] " "]
	set skewness [split [ic_hex_rh 20 FLUID proj 2 minval 0 -type eriksson maxval 1 new_format] " "]
	log "Worst Determinant: [lindex $determinante 0]"
	log "Worst Skewness: [lindex $skewness 0]"

	# Try Smoothing

	ic_hex_smooth 5 GEOM INLET DISKS SYM2 SYM1 VOLUTE DIFFUSOR OUTLET FLUID elliptic iter_srf 15 iter_vol 5 exp_srf 0.0 exp_vol 0.0 niter_post 3 limit_post 0.2 smooth_type 202 nfix_layers -1 rebunch_edges 0 treat_unstruct 2 stabilize_srf 1.0 stabilize_vol 2.0 ortho_distance_srf 1 ortho_distance_vol 0 surface_fitting 1 keep_per_geom 1

	# Smoothing succeded?
	set determinante [split [ic_hex_rh 20 FLUID proj 2 minval 0 -type determinant_27 maxval 1 new_format] " "]
	set skewness [split [ic_hex_rh 20 FLUID proj 2 minval 0 -type eriksson maxval 1 new_format] " "]
	if {[lindex $determinante 0]<0.4} {
		ic_hex_create_mesh INLET DISKS SYM2 SYM1 VOLUTE DIFFUSOR OUTLET FLUID proj 2 dim_to_mesh 3 nproc 4
		log "No smoothing applied. Determinante too low, [lindex $determinante 0]"
	} else {
		log "After Smoothing Worst Determinant: [lindex $determinante 0]"
		if {[lindex $skewness 0]<0.3} {
			ic_hex_create_mesh INLET DISKS SYM2 SYM1 VOLUTE DIFFUSOR OUTLET FLUID proj 2 dim_to_mesh 3 nproc 4
			log "No smoothing applied. Skewness not good enough, [lindex $skewness 0]"
		} else {
			log "After Smooting Worst Skewness: [lindex $skewness 0]"  
		}			
	}


} 
# debug encasing

if {$debug == 1} {
	################################################################################
	log "Export Final Mesh"
	################################################################################
	# Boundary Conditions
	ic_boco_solver Fluent_V6
	ic_solution_set_solver Fluent_V6 1
	ic_boco_solver {ANSYS Fluent}
	ic_solver_mesh_info {ANSYS Fluent}

	ic_boco_set VOLUTE { { 1  {WALL}  0  } }
	ic_boco_set DIFFUSOR { { 1  {WALL}  0  } }		
	ic_boco_set INLET { { 1  {MASFI}  0  } }
	ic_boco_set FLUID { { 1  {FLUID}  0  } }
	ic_boco_set SYM1 {{1 PER 0}}
	ic_boco_set SYM2 {{1 PER 0}}		
	ic_boco_set OUTLET { { 1  {PRESO}  0  } }
	ic_boco_set DISKS { { 1  {WALL}  0  } }
	ic_boco_save $outputfile.fbc
	ic_boco_save_atr $outputfile.atr	

	#Convert Final to Unstruct	
	cmd_rm "hex.uns"	
	ic_hex_write_file hex.uns INLET DISKS SYM2 SYM1 VOLUTE DIFFUSOR OUTLET FLUID proj 2 dim_to_mesh 3 -family_boco $outputfile.fbc	

	log "Export Command: $icemenv/icemcfd/output-interfaces/fluent6 -dom hex.uns -b $outputfile.fbc -scale 0.001,0.001,0.001 $outputfile.msh"
	ic_exec "$icemenv/icemcfd/output-interfaces/fluent6" -dom hex.uns -b $outputfile.fbc -scale 0.001,0.001,0.001 -per "type trans base \{0 0 0\} axis \{0 0 $domth\} angle 0" $outputfile.msh
	log "Mesh exported to $outputfile.msh"
	# debug encasing	
} 

