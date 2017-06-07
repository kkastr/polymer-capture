# Chain size
# Unitless
set N            [lindex $argv 0]

# Applied voltage drop
# V/m to eps/sigma
set dVexp        [lindex $argv 1]

# Geometry Parameters (radius and full thickness)
# sigma
set reff         [lindex $argv 2]
set leff         [lindex $argv 3]

# Electrolyte Choice
set muflag      [lindex $argv 4]

# Run Parameters
# Unitless
set nevents      [lindex $argv 5]
set rseed        [lindex $argv 6]
t_random seed    $rseed

# Visualization on/off
set visflag     [lindex $argv 7]


########################################################################################
#### Input Parameters

# Domain size
set box_l 400

# Experimental quantities
# Currently tuning to Ntune=100, but can set Ntune=N if we want to try that
set Ntune    100
set sigmaexp [expr {5.0*pow(10.,-9.)}]
set Dexp     [expr {2.38*pow(10.,-12.)/pow($Ntune*$sigmaexp*1e6,0.608)}]
if {$muflag == 1} {
    # NaCl high-salt limit
    set muexp [expr {3.14*pow(10.,-8.)}]
} elseif {$muflag == 2} {
    # LiCl 40 M
    set muexp [expr {0.6*pow(10.,-8.)}]
}

# Scales
set sigma 1.0
set eps   1.0
set kT    1.0
set gamma 1.0

# Peclet tuning
set dVtot  [expr {$kT*$muexp*$dVexp/($Ntune*$Dexp)}]
set tauexp [expr {$kT*$sigmaexp*$sigmaexp/($Ntune*$gamma*$Dexp)}]

# Time Parameters
set eqltime    [expr {1e2*$N}]
set cutofftime 1e6

# Distance Parameters
# Failure distance is propto sqrt(N)
set eqldist     [expr {sqrt($N) * 10}]
set captdist    [expr {sqrt($N) * 2} ]
set cutoffdist  [expr {sqrt($N) * 8}]

# FENE Potential Parameters
set kap [expr {30.0*$eps/($sigma*$sigma)}]
set lam [expr {1.5*$sigma}]

# LJ Potential Parameters
set cut   [expr {pow(2.0,1.0/6.0)*$sigma}]
set shift [expr {0.25*$eps}]

# Angular Potential Parameters
# Lp ~~ angular spring constant
set Lp [expr {10.0*$eps/($sigma*$sigma)}]


########################################################################################
#### Initialization (only runs once)

# Spatial domain creation
setmd box_l $box_l $box_l $box_l

# Temporal domain creation
setmd time_step 0.01
setmd skin 0.4

# Interaction creations
inter 0 fene $kap $lam
inter 7 angle $Lp 3.14159
inter 0 0  lennard-jones $eps $sigma $cut $shift 0.
inter 0 76 lennard-jones $eps $sigma $cut $shift 0.

# Derived Geometric Quantities (do not edit)
set rnom [expr {$reff + 0.5*$sigma}]
set lnom [expr {$leff - 1.0*$sigma + 1e-3}]

# Pore creation
constraint pore center [expr {$box_l/2.}] [expr {$box_l/2.}] [expr {$box_l/2.}] axis 0 0 1 radius $rnom length [expr $lnom/2.] type 76

# Polymer creation
set x [expr {$box_l/2.}]
set y [expr {$box_l/2.}]
set z [expr {$box_l/2. - $eqldist}]
part 0 pos $x $y $z type 0 fix 0 0 0 ext_force 0. 0. 0.
# Build the chain (without angular potential)
for { set i 1 } { $i < $N } { incr i } {
    set z [expr {$z - $sigma}]
    part $i pos $x $y $z type 0 fix 0 0 0 ext_force 0. 0. 0.
    part $i bond 0 [expr {$i - 1}] 
}
# After all particles are placed, apply the angular potential
for { set i 1 } { $i < [expr {$N-1}] } { incr i } {
    part $i bond 7 [expr {$i - 1}] [expr {$i + 1}]
}

# Initialize Visualization
if { $visflag == 1 } {
    prepare_vmd_connection vmdout
    imd listen 100
    imd positions
}

# Output files
set Fout  [open "data/${muflag}_${dVexp}_${reff}_${leff}_${N}_${rseed}.dat"  "a+"]
puts $Fout "muflag,dVexp,reff,leff,N,rseed,t,nretract,ntimeout"
close $Fout


########################################################################################
########################################################################################

# Main Loop 
set casenum 1
set nretract 0
set ntimeout 0
while { $casenum <= $nevents } {

    ##############################################################################
    # Equilibrate far from the pore

    # Randomly sample hemisphere
    # https://tinyurl.com/zkmpszk
    set rejectY 1
    while {$rejectY == 1} {
	# Sample [-1,1] three times
	set Y1 [expr {2*rand()-1}]
	set Y2 [expr {2*rand()-1}]
	set Y3 [expr {2*rand()-1}]
	# Reflect to negative hemisphere
	if {$Y3>0} {set Y3 [expr {-$Y3}]}
	# Reject when Ynorm>1
	set Ynorm [expr {sqrt($Y1*$Y1+$Y2*$Y2+$Y3*$Y3)}]
	if {$Ynorm <= 1.0} {
	    set Y1 [expr {$eqldist*$Y1/$Ynorm}]
	    set Y2 [expr {$eqldist*$Y2/$Ynorm}]
	    set Y3 [expr {$eqldist*$Y3/$Ynorm}]
	    # Reject when Y3 too close to wall (overlaps)
	    if {$Y3 <= -$sigma/2.} {set rejectY 0}
	}
    }
    
    # Place the polymer 
    set x [expr {$box_l/2. + $Y1}]
    set y [expr {$box_l/2. + $Y2}]
    set z [expr {$box_l/2. - $leff/2. + $Y3}]
    part 0 pos $x $y $z
    for { set i 1 } { $i < $N } { incr i } {
	set x [expr {$x - $sigma}]
	part $i pos $x $y $z 
    }

    # Equilibrate (fix middle monomer)
    puts "Equilibrating..."
    part [expr {$N/2}] fix
    thermostat langevin $kT 0.1
    for { set i 0 } { $i < $N } { incr i } {
	# Remove any electric field via hacked forces.c
	part $i ext_force 0 0 0
    }    
    for {set i 0} {$i < $eqltime} {incr i 1} {
	integrate 100
	if { $visflag == 1 } {imd positions}
    }
    puts "Equilibrated"
    part [expr {$N/2}] unfix

    ##############################################################################
    # Translate to near the pore

    # Calculate closest monomer
    set rclosest $box_l
    set iclosest 0
    for { set i 0 } { $i < $N } { incr i } {
    	# Monomer position relative to pore face
    	set x [expr {[lindex [part $i print pos] 0] - $box_l/2.}]
    	set y [expr {[lindex [part $i print pos] 1] - $box_l/2.}] 
    	set z [expr {[lindex [part $i print pos] 2] - ($box_l/2.-$leff/2.)}]
    	# Count closest monomer
    	set r [expr {sqrt($x*$x+$y*$y+$z*$z)}]
    	if { $r < $rclosest } {
    	    set rclosest $r
    	    set iclosest $i
    	}
    }

    # Translate closest monomer to capture radius from pore face
    set oldx [expr {[lindex [part $iclosest print pos] 0] - $box_l/2.}]
    set oldy [expr {[lindex [part $iclosest print pos] 1] - $box_l/2.}]
    set oldz [expr {[lindex [part $iclosest print pos] 2] - ($box_l/2.-$leff/2.)}]
    set newx [expr {$captdist*$oldx/$rclosest + $box_l/2.}]
    set newy [expr {$captdist*$oldy/$rclosest + $box_l/2.}]
    set newz [expr {$captdist*$oldz/$rclosest + ($box_l/2.-$leff/2.)}]
    part $iclosest pos $newx $newy $newz

    # Move the rest of the chain accordingly
    set shiftx [expr {$newx - ($oldx + $box_l/2.)}]
    set shifty [expr {$newy - ($oldy + $box_l/2.)}]
    set shiftz [expr {$newz - ($oldz + ($box_l/2.-$leff/2.))}]
    for { set i 0 } { $i < $N } { incr i } {
	if {$i != $iclosest} {
	    set oldx [lindex [part $i print pos] 0]
	    set oldy [lindex [part $i print pos] 1]
	    set oldz [lindex [part $i print pos] 2]
	    set newx [expr {$oldx + $shiftx}]
	    set newy [expr {$oldy + $shifty}]
	    set newz [expr {$oldz + $shiftz}]
	    part $i pos $newx $newy $newz 
	}
    }

    # If overlaps occur, kill event
    set zmax [expr {-$box_l}]
    for { set i 0 } { $i < $N } { incr i } {
    	# Monomer position relative to pore face
    	set z [expr {[lindex [part $i print pos] 2] - ($box_l/2.-$leff/2.)}]
    	# Count closest monomer
    	if { $z > $zmax } {
    	    set zmax $z
    	}
    }
    if {$zmax > -$sigma/2.} {
	puts "Bad overlap after transport."
	continue
    }


    ##############################################################################
    # Conduct translocation

    # Simulate
    thermostat langevin $kT $gamma
    for { set i 0 } { $i < $N } { incr i } {
	# Apply electric field via hacked forces.c
	part $i ext_force $dVtot $reff $leff
    }
    set ntrans 0
    set t 0
    while {1} {

	##################################################################
	# Compute Metrics

	# Analyse Position
	set ntrans 0
	set rmin $box_l
	for { set i 0 } { $i < $N } { incr i } {
	    # Monomer position relative to pore
	    set x [expr {[lindex [part $i print pos] 0] - $box_l/2.}]
	    set y [expr {[lindex [part $i print pos] 1] - $box_l/2.}]
	    set z [expr {[lindex [part $i print pos] 2] - ($box_l/2.-$leff/2.)}]
	    # Count Ntrans
	    if { $z >= $leff/2. } {incr ntrans}
	    # Count closest monomer
	    set r [expr {sqrt($x*$x+$y*$y+$z*$z)}]
	    if { $r < $rmin } {set rmin $r}
	}

	##################################################################
	# Interpret Metrics

	# If event succeeds, increment casenum, track tau, reset nfails, and stop
	if { $ntrans == $N } {
	    puts "Event succeeded."
	    set Fout  [open "data/${muflag}_${dVexp}_${reff}_${leff}_${N}_${rseed}.dat"  "a+"]
	    puts $Fout "$muflag,$dVexp,$reff,$leff,$N,$rseed,$t,$nretract,$ntimeout"
	    close $Fout
	    incr casenum 1
	    set nretract 0
	    set ntimeout 0
	    break
	}

	# If chain is too far, consider it a failure
	if { $rmin > $cutoffdist } {
	    puts "Chain too far from pore."
	    incr nretract 1
	    break
	}

	# If time limit exceed, consider event a failure
	if { $t >= $cutofftime} {
	    puts "Time limit exceeded."
	    incr ntimeout 1
	    break
	}

	# Integrate
	integrate 100
	incr t 1
	if { $visflag == 1 } {imd positions}
    }
}






