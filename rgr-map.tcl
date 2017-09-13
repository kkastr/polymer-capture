set N [lindex $argv 0]
set rseed  [lindex $argv 1]
set filename [lindex $argv 2]
set transportdist [lindex $argv 3]
set vis [lindex $argv 4]


t_random seed $rseed



set boxx 500
set boxy 500
set boxz 500 

set cx [expr $boxx/2.0]
set cy [expr $boxy/2.0]
set cz [expr $boxz/2.0]
set nmax 10
set temp 1.0; set gamma 1.0; set gamma_equilibration 0.01



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
# set transportdist [expr 60]
set cutoffdist [expr $transportdist + 60]
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




part [expr $N] pos $cx $cy $cz type 99 fix
constraint pore center [expr $cx] [expr $cy] [expr $cz] axis 0 0 1 radius $rpore length $tpore type 1

for { set i 0 } { $i < $N } { incr i } {
	set x [expr $cx - $N/2 + $i]
	set y [expr $cy]
	set z [expr $cz  + $transportdist]
	part $i pos $x $y $z type 0 


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

set translation_part [expr int($N*0.5)]
set rgr_csv [open "data/${filename}_$N/rgr-${filename}-$N-$rseed-$transportdist.csv" "a"]
while {$flag == 0} {

	set overlap 0
	if {$position_flag == 1} {
		for {set i 0} {$i < $N} {incr i} {
			part $i ext_force 0 0 0
		}
	for { set i 0 } { $i < $N } { incr i } {
		set x [expr $cx - $N/2 + $i]
		set y [expr $cy]
		set z [expr $cz  + $transportdist]
		part $i pos $x $y $z type 0 
		}
	set position_flag 0
	}	


	# #puts "equilibrated."
	
	
	
	# #puts $drawrand

	# set rlist {}

	# set randlist [randgen]


	# set px [lindex [part $translation_part print pos] 0]
	# set py [lindex [part $translation_part print pos] 1]
	# set pz [lindex [part $translation_part print pos] 2]


	# set mx [expr $cx + $lj_cutoff + $transportdist*[lindex $randlist 0]]
	# set my [expr $cy + $lj_cutoff + $transportdist*[lindex $randlist 1]]
	# set mz [expr $cz + $lj_cutoff + $transportdist*[lindex $randlist 2]]

	# #puts "$mx $my $mz"
	# set mr2  [expr ($mx - $cx)*($mx -$cx) + ($my - $cy)*($my -$cy) + ($mz - $cz)*($mz -$cz)]
	# set mr [expr sqrt($mr2)]
	# #puts $mr
	
	# #puts $min_part
	# part $translation_part pos $mx $my $mz type 0

	# #set translation_part [expr $N*0.5]
	# set tx [lindex [part $translation_part print pos] 0]
	# set ty [lindex [part $translation_part print pos] 1]
	# set tz [lindex [part $translation_part print pos] 2]
	

	# #part $min_part fix
	# # for {set i 0} {$i < $N} {incr i} {
	# # 	part $i bond delete
	# # }
	# set tr [expr sqrt($tx*$tx + $ty*$ty + $tz*$tz)]

	# for { set i [expr 0] } { $i < $N } { incr i } {
	# 	if {$i != $translation_part} {
	# 		set ix [expr [lindex [part $i print pos] 0] - $px + $tx] 
	# 		set iy [expr [lindex [part $i print pos] 1] - $py + $ty] 
	# 		set iz [expr [lindex [part $i print pos] 2] - $pz + $tz] 
	# 		set ir [expr sqrt($ix*$ix + $iy*$iy + $iz*$iz)]
	# 		#puts "$ix $iy $iz"
	# 		#puts $ir
	# 		part $i pos $ix $iy $iz type 0 
	# 	}
	# }

	
	# set zlist {}
	# set zmincheck 32768
	# for {set i 0} {$i < $N} {incr i} {
	# 	set z [lindex [part $i print pos] 2]
	# 	if {$z < $zmincheck} {
	# 		set zmincheck $z
	# 	}
	# }

	if {$zmincheck < [expr $cz + $lj_cutoff + 0.5 + 0.5] } {
		puts "Bad overlap with pore"
		set position_flag 1
		incr overlap
	
		#puts  "overlap $overlap"
		continue
	}
	puts "overlap $overlap"

	puts "positionflag $position_flag"
	
	#set drawrand 0

	for {set i 0} {$i < $N} {incr i} {
		part $i ext_force $force 1.6 1.6
	}

	part $fixed_N fix

	thermostat langevin $temp $gamma_equilibration
	for {set i 0} {$i < $equil_time} {incr i} {
	    
	    integrate 100
	    imd positions
	}

	thermostat langevin $temp $gamma

	part $fixed_N unfix

	set rg_at_equil [analyze rg 0 1 $N]

	puts $rgr_csv "$rg_at_equil"
	#puts $rg_at_equil
	#puts "I'm back in the first while"
	set tstart $t

	incr t
	if {$t > 1000} {
		set flag 1
	}
}
# close $part_pos_contact
# close $part_pos_trans
#close $part_pos_z
# close $metric_csv
close $rgr_csv
