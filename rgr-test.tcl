set N [lindex $argv 0]
set rseed  [lindex $argv 1]
set filename [lindex $argv 2]
#set transportdist [lindex $argv 3]
set vis [lindex $argv 3]


t_random seed $rseed



set boxx 500
set boxy 500
set boxz 500 

set cx [expr $boxx/2.0]
set cy [expr $boxy/2.0]
set cz [expr $boxz/2.0]
set nmax 10
set temp 1; set gamma 1.0; set gamma_equilibration 0.01



setmd box_l $boxx $boxy $boxz
setmd periodic 1 1 1
setmd time_step 0.01; setmd skin 0.4
thermostat langevin $temp $gamma


#FENE potential
set k_fene 30.0
set r_fene 1.5

#Angular Potential - Harmonic
set k_angle [expr 10.0 ]
set pi 3.14159

#Shifted Lennard-Jones
set eps 1.0
set sigma 1.0
set lj_cutoff 1.12246
set lj_shift 0.25
set lj_offset 0 


#mesh 
set Nr 200
set Nz 200
set dr 1
set dz 1



inter 0 fene $k_fene $r_fene
inter 1 angle $k_angle $pi
inter 0 0 lennard-jones $eps $sigma $lj_cutoff $lj_shift $lj_offset
inter 1 0 lennard-jones $eps $sigma $lj_cutoff $lj_shift $lj_offset


set rgwlc2 [expr (1.0/3.0) * $k_angle * $N - pow($k_angle,2) + 2.0 * (pow($k_angle,3)/$N) * (1 - ($k_angle/$N)*(1- exp(-$N/$k_angle)) )  ]
set rgwlc [expr sqrt($rgwlc2)]
# puts $rgwlc

set t_trans 0
set trans_flag 0
set fixed_N [expr $N/2]
set equil_time [expr 50.0 * $N]
set tpore 0.5
set z_line [expr $cz - $tpore/2]
set force [expr -3.0]
set n_attempt 0
set illegal_mov 0 
set rpore 1.8
set reff [expr $rpore - 0.5]
set cutofftime 1e6
set transportdist [expr 40]
set cutoffdist [expr $transportdist + $rgwlc]
set contactdist [expr 3]

puts "$transportdist, $rpore, $force"

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

for {set k 0} { $k < $Nr} {incr k} {
	#puts ""
	for {set l 0} {$l < $Nz} {incr l} {
		set rarray($k,$l) 0
		# puts -nonewline "[set array($k,$l)] \ "
	}
}
# puts ""

for {set k 0} { $k < $Nr} {incr k} {
	#puts ""
	for {set l 0} {$l < $Nz} {incr l} {
		set rgarray($k,$l) 0
		# puts -nonewline "[set array($k,$l)] \ "
	}
}
# puts ""

for {set k 0} { $k < $Nr} {incr k} {
	#puts ""
	for {set l 0} {$l < $Nz} {incr l} {
		set rgparray($k,$l) 0
		# puts -nonewline "[set array($k,$l)] \ "
	}
}




part [expr $N] pos $cx $cy $cz type 99 fix
#part [expr $N+1] pos [expr 250 + 80/sqrt(2)] [ expr 250 + 80/sqrt(2)] 250 type 99 fix
constraint pore center [expr $cx] [expr $cy] [expr $cz] axis 0 0 1 radius $rpore length $tpore type 1

for { set i 0 } { $i < $N } { incr i } {
	set x [expr $cx - $N/2 + $i]
	set y [expr $cy]
	set z [expr $cz  + $transportdist]
	part $i pos $x $y $z type 0 fix 


	if { $i > 0 } {
		part $i bond 0 [expr $i - 1]
	}
	if { $i > 1} {
		part [expr $i - 1] bond 1 $i [expr $i - 2]
	}
}

if { $vis == 1 } {
    prepare_vmd_connection vmdout
    imd listen 100
    imd positions
}

# set part_pos_contact [open "data/${filename}_$N/part_pos_contact-$N-$rseed.xyz" "a"]
# set part_pos_trans [open "data/${filename}_$N/part_pos_trans-$N-$rseed.xyz" "a"]
#set part_pos_z [open "data/${filename}_$N/part_pos_z-$N-$rseed.xyz" "a"]
#set metric_csv [open "data/${filename}_$N/metric-${filename}-$N-$rseed.csv" "a"]



set flag 0
set t 0
set n 0
set rg_flag 0
set fail 0
set success 0
set trans_flag 0
set stuck 0
set drawrand 0


set position_flag 0
set contactflag 0
set ncontact 0
set st 0


#set positions_csv [open "data/${filename}_$N/traj-${filename}-$N-$rseed.csv" "a"]
while {$flag == 0} {

	while {1} {

		set rg_test [analyze rg 0 1 $N]

		puts "$rg_test"

		integrate 100
		imd positions
		
		incr t
	}
}