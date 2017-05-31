# set N      [lindex $argv 0]
# set rseed  [lindex $argv 1]
# #set filename [lindex $argv 2]
# set vis [lindex $argv 3]

#t_random seed    $rseed


#puts $N


set N 10

set boxx 100
set boxy 100
set boxz 100
set cx [expr $boxx/2.0]
set cy [expr $boxy/2.0]
set cz [expr $boxz/2.0]

set temp 1.0; set gamma 1.0; set gamma_equilibration 0.1


setmd box_l $boxx $boxy $boxz
setmd periodic 1 1 1
setmd time_step 0.01; setmd skin 0.4
thermostat langevin $temp $gamma



#FENE potential
set k_fene 30.0
set r_fene 1.5

#Angular Potential - Harmonic
set k_angle 6.0
set pi 3.14159

#Shifted Lennard-Jones
set eps 1.0
set sigma 1.0
set lj_cutoff 1.12246
set lj_shift 0.25
set lj_offset 0 

inter 0 fene $k_fene $r_fene
inter 1 angle $k_angle $pi
inter 0 0 lennard-jones $eps $sigma $lj_cutoff $lj_shift $lj_offset
inter 1 0 lennard-jones $eps $sigma $lj_cutoff $lj_shift $lj_offset


proc rand_range {min max} { 
	return [expr rand() * ($max - $min) + $min] 
}

proc randgen {} {
	set sumsq 1.1
	
	while {$sumsq >= 1.0} {
		set x1 [rand_range -1 1]
		set x2 [rand_range -1 1]
		set sumsq [expr $x1*$x1 + $x2*$x2]
	}
	set x [expr 2.0 * $x1 *sqrt( 1.0 - $sumsq )]
	set y [expr 2.0 * $x2 *sqrt( 1.0 - $sumsq )]
	set z [expr abs(1.0 - 2.0 * ($sumsq))]
	set r [expr $x*$x + $y*$y + $z*$z]
	set randlist "$x $y $z"
	return $randlist
}

proc fbuild {} {
	global N cx cy cz

	set rlist [randgen]
	set mx [expr $cx + 20*[lindex $rlist 0]]
	set my [expr $cy + 20*[lindex $rlist 1]]
	set mz [expr $cz + 20*[lindex $rlist 2]]

	part [expr $N/2] pos $mx $my $mz type 0

	set kx [lindex [part [expr $N/2] print pos] 0]
	set ky [lindex [part [expr $N/2] print pos] 1]
	set kz [lindex [part [expr $N/2] print pos] 2]
	
	for { set i [expr ($N/2) + 1] } { $i < $N } { incr i } {
		set x [expr $kx + ( $i - $N/2 )]
		set y [expr $ky]
		set z [expr $kz]
		part $i pos $x $y $z type 0 fix 0 0 0
	}
	for { set i [expr ($N/2) - 1] } { $i > 0 } { set i [expr $i - 1] } {
		set x [expr $kx + ( $i - $N/2 )]
		set y [expr $ky]
		set z [expr $kz]
		part $i pos $x $y $z type 0 fix 0 0 0
	}
	# if { $i > 0 } {
	# 	part $i bond 0 [expr $i - 1]
	# }
	# if { $i > 1} {
	# 	part [expr $i - 1] bond 1 $i [expr $i - 2]
	# }
}


fbuild


part [expr $N + 1] pos $cx $cy $cz type 1

constraint pore center [expr $cx] [expr $cy] [expr $cz] axis 0 0 1 radius 1.8 length 1.0 type 1






set visu 1
if [expr $visu == 1] {
	prepare_vmd_connection nano 11000
	#exec sleep 3
	imd positions
	imd listen 500
}	



# proc fmd{} {
# 	while { $flag == 0 } {
		
# 	}
# }


