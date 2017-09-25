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
set nmax 5
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
set force [expr -5.0]
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



set rmatrix_csv [open "data/${filename}_$N/rmatrix-${filename}-$N-$rseed.csv" "a"]
set rgrmatrix_csv [open "data/${filename}_$N/rgrmatrix-${filename}-$N-$rseed.csv" "a"]
set rgpmatrix_csv [open "data/${filename}_$N/rgpmatrix-${filename}-$N-$rseed.csv" "a"]

#set positions_csv [open "data/${filename}_$N/traj-${filename}-$N-$rseed.csv" "a"]
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
			set rgrarray($k,$l) 0
			# puts -nonewline "[set array($k,$l)] \ "
		}
	}
	for {set k 0} { $k < $Nr} {incr k} {
		#puts ""
		for {set l 0} {$l < $Nz} {incr l} {
			set rgparray($k,$l) 0
			# puts -nonewline "[set array($k,$l)] \ "
		}
	}


	set position_flag 0
	}	
	part $fixed_N fix

	thermostat langevin $temp $gamma_equilibration
	for {set i 0} {$i < $equil_time} {incr i} {
	    
	    integrate 100
	    imd positions
	}

	thermostat langevin $temp $gamma

	part $fixed_N unfix

	#puts "equilibrated."
	
	set rg_at_equil [analyze rg 0 1 $N]
	
	#puts $drawrand

	# set rlist {}

	# for {set j 0} {$j < $N} {incr j} {
	# 	set x [lindex [part $j print pos] 0]
	# 	set y [lindex [part $j print pos] 1]
	#     set z [lindex [part $j print pos] 2]
	#     set r [expr sqrt($x*$x + $y*$y + $z*$z)]
	#     lappend rlist $r
	# }

	# set rmin [::tcl::mathfunc::min {*}$rlist]
	# set rmax [::tcl::mathfunc::max {*}$rlist]


	# for {set i 0} {$i < $N} {incr i} {
	#     set z [lindex [part $i print pos] 2]
	#     set y [lindex [part $i print pos] 1]
	#     set x [lindex [part $i print pos] 0]
	#     set r [expr sqrt($x*$x + $y*$y + $z*$z)]
	#     if {$r == $rmin} {
	#       set min_part $i
	#     } 
	# }


	# set randlist [randgen]


	# set px [lindex [part $min_part print pos] 0]
	# set py [lindex [part $min_part print pos] 1]
	# set pz [lindex [part $min_part print pos] 2]


	# set mx [expr $cx + $lj_cutoff + $transportdist*[lindex $randlist 0]]
	# set my [expr $cy + $lj_cutoff + $transportdist*[lindex $randlist 1]]
	# set mz [expr $cz + $lj_cutoff + $transportdist*[lindex $randlist 2]]

	# #puts "$mx $my $mz"
	# set mr2  [expr ($mx - $cx)*($mx -$cx) + ($my - $cy)*($my -$cy) + ($mz - $cz)*($mz -$cz)]
	# set mr [expr sqrt($mr2)]
	# #puts $mr

	# #puts $min_part
	# part $min_part pos $mx $my $mz type 0

	# set tx [lindex [part $min_part print pos] 0]
	# set ty [lindex [part $min_part print pos] 1]
	# set tz [lindex [part $min_part print pos] 2]
	

	# #part $min_part fix
	# # for {set i 0} {$i < $N} {incr i} {
	# # 	part $i bond delete
	# # }
	# set tr [expr sqrt($tx*$tx + $ty*$ty + $tz*$tz)]

	# for { set i [expr 0] } { $i < $N } { incr i } {
	# 	if {$i != $min_part} {
	# 		set ix [expr [lindex [part $i print pos] 0] - $px + $tx] 
	# 		set iy [expr [lindex [part $i print pos] 1] - $py + $ty] 
	# 		set iz [expr [lindex [part $i print pos] 2] - $pz + $tz] 
	# 		set ir [expr sqrt($ix*$ix + $iy*$iy + $iz*$iz)]
	# 		#puts "$ix $iy $iz"
	# 		#puts $ir
	# 		part $i pos $ix $iy $iz type 0 
	# 		}
	# 	}

	
	set zlist {}
	set zmincheck 32768
	for {set i 0} {$i < $N} {incr i} {
		set z [lindex [part $i print pos] 2]
		if {$z < $zmincheck} {
			set zmincheck $z
		}
	}

	if {$zmincheck < [expr $cz + $lj_cutoff + 0.5 + 0.5] } {
		puts "Bad overlap with pore"
		set position_flag 1
		incr overlap
	
		#puts  "overlap $overlap"
		continue
	}
	puts "overlap $overlap"

	puts "positionflag $position_flag"
	
	set drawrand 0


	#puts "I'm back in the first while"
	set tstart $t
	


	while {1} {

		for {set i 0} {$i < $N} {incr i} {
			part $i ext_force $force $reff $reff
		}

		set z_list {}
		set r_list {}		
		if { $n > $nmax } {
			puts "n > nmax"
			set flag 1
			puts "$fail $success"
			break
		}
		set rmincheck 32768
		set imin 0
		for {set i 0} { $i < $N } {incr i} {
			set x [lindex [part $i print pos] 0]
			set y [lindex [part $i print pos] 1]
			set z [lindex [part $i print pos] 2]
			set r [expr sqrt(($x-$cx)*($x-$cx) + ($y-$cy)*($y-$cy) + ($z-$cz)*($z-$cz))]
			if {$r < $rmincheck} {
				set imin $i
				#puts "$imin"
				set rmincheck $r
			}
			lappend z_list $z
			lappend r_list $r			
		}
		
		
		
		set z_min [::tcl::mathfunc::min {*}$z_list]
		set z_max [::tcl::mathfunc::max {*}$z_list]
		set r_min [::tcl::mathfunc::min {*}$r_list]
		set r_max [::tcl::mathfunc::max {*}$r_list]

		set iminx [lindex [part $imin print pos] 0]
		set iminy [lindex [part $imin print pos] 1]
		set iminz [lindex [part $imin print pos] 2]
		set iminr [expr sqrt(pow([expr $iminx - $cx],2) + pow([expr $iminy - $cy],2) + pow([expr $iminz - $cz],2))]




		set x_2 0
		set y_2 0
		set z_2 0
		set xs 0
		set ys 0
		set zs 0

		for { set i 0 } { $i < $N } { incr i } {
			set x [lindex [part $i print pos] 0]
			set y [lindex [part $i print pos] 1]
			set z [lindex [part $i print pos] 2]
			
			set x_2 [expr $x_2 + pow($x,2)] 
			set y_2 [expr $y_2 + pow($y,2)]
			set z_2 [expr $z_2 + pow($z,2)]
			set xs [expr $xs + $x] 
			set ys [expr $ys + $y]
			set zs [expr $zs + $z]
		}
		
		set xcom [expr $xs/$N]
		set ycom [expr $ys/$N]
		set zcom [expr $zs/$N]
		set rcom [expr sqrt(pow($xcom - (250 + 80/sqrt(2)),2) + pow($ycom - (250 + 80/sqrt(2)),2))]
		set rg_per_tu [analyze rg 0 1 $N]


		# part [expr $N + 1] pos $xcom $ycom $zcom type 99


		# set p1x [lindex [part 51 print pos] 0]
		# set p1y [lindex [part 51 print pos] 1]
		# set p1z [lindex [part 51 print pos] 2]
		
		# puts "$p1x $p1y $p1z"
		set idx1 [expr int( $rcom )]
		set idx2 [expr abs(int( $zcom - 250 ))] 

		# puts $idx1
		# puts $idx2

		set rgr [lindex $rg_per_tu 0]
		set rgp [expr 0.5*([lindex $rg_per_tu 1] + [lindex $rg_per_tu 2]) ]
		set rarray($idx1,$idx2) [expr [set rarray($idx1,$idx2)] + 1]
		set rgrarray($idx1,$idx2) [expr [set rarray($idx1,$idx2)] + $rgr]
		set rgparray($idx1,$idx2) [expr [set rarray($idx1,$idx2)] + $rgp]



		# puts "$idx1 $idx2 [set rarray($idx1,$idx2)]"

		# if {$t%10 == 0} { 	
		# 	puts $positions_csv "$imin,$r_min,$iminx,$iminy,$iminz"
		# }

		if {$r_min > $cutoffdist} {
			puts "Dist greater than r = $cutoffdist from pore"
			incr fail
			set position_flag 1
			break
		}

		if {$r_min > [expr 0.5 * $rgwlc] } {
			set tstart $t
		}

		if {[expr $t - $tstart ] > $cutofftime} {
			puts "cutoff time exceeded - stuck event"
			incr stuck
			set position_flag 1 
			break
		}



		if {$r_min <= $contactdist && $contactflag == 0} {
			puts "Contact with pore"
			for {set i 0} {$i < $N} {incr i} {
				set x [lindex [part $i print pos] 0]
				set y [lindex [part $i print pos] 1]
				set z [lindex [part $i print pos] 2]
				set r [expr sqrt(($x-$cx)*($x-$cx) + ($y-$cy)*($y-$cy) + ($z-$cz)*($z-$cz))]
				if {$r == $r_min} {
					set ncontact $i
				}
			}
			set contactflag 1
		}

		if {$z_min > $z_line && $trans_flag == 1} {
			puts "zmin greater than zline"
			set trans_flag 0
		}

		if {$z_min < $z_line && $trans_flag == 0} {
			
		
			#puts "inside translocation if"
			set t_last_thread $t
			
			if { $n_attempt == 0 } {
				puts "n_attempt $n_attempt"
				set t_first_thread $t
			
  				set n_attempt [expr $n_attempt + 1]
			}
			set rg_calc_trans [analyze rg 0 1 $N]

			for {set i 0} {$i < $N} {incr i} {
				set x [lindex [part $i print pos] 0]
				set y [lindex [part $i print pos] 1]
				set z [lindex [part $i print pos] 2]
				set r [expr sqrt(($x-$cx)*($x-$cx) + ($y-$cy)*($y-$cy) + ($z-$cz)*($z-$cz))]
				if {$z < $z_line} {
					set st $i
				}
			}


			set trans_flag [expr $trans_flag + 1 ]
		}
		#puts $z_max
		if {$z_max < $z_line} {
			puts "zmax less than zline"
			set t_thread $t
			puts $t_last_thread
			set t_trans [expr $t_thread - $t_last_thread]
			
			set metric_csv [open "data/${filename}_$N/metric-${filename}-$N-$rseed.csv" "a"]

			puts $metric_csv "$N,$t_trans,$rg_calc_trans,$rg_at_equil,$t_first_thread,$t_thread,$t_last_thread,$fail,$stuck,$ncontact,$st"
			close $metric_csv

	      	set n_attempt 0
	      	set rg_flag 0  
			set n [expr $n + 1.0]
			set position_flag 1
			set contactflag	0
			break
		}




		if {$t_trans != 0} {
			set t_trans 0

			for {set k 0} { $k < $Nr} {incr k} {
				puts $rmatrix_csv ""
				for {set l 0} {$l < $Nz} {incr l} {
					puts -nonewline  $rmatrix_csv "[set rarray($k,$l)] \ "
				}
			}
			puts $rmatrix_csv " Next Event "


			for {set k 0} { $k < $Nr} {incr k} {
				puts $rgrmatrix_csv ""
				for {set l 0} {$l < $Nz} {incr l} {
					puts -nonewline $rgrmatrix_csv "[set rgrarray($k,$l)] \ "
				}
			}
			puts $rgrmatrix_csv " Next Event "

			for {set k 0} { $k < $Nr} {incr k} {
				puts $rgpmatrix_csv ""
				for {set l 0} {$l < $Nz} {incr l} {
					puts -nonewline $rgpmatrix_csv "[set rgparray($k,$l)] \ "
				}
			}
			puts $rgpmatrix_csv " Next Event "

		}
		integrate 100
		imd positions
		
		incr t
	}
}
# close $part_pos_contact
# close $part_pos_trans
#close $part_pos_z
# close $metric_csv
#close $positions_csv
close rmatrix_csv
close rgrmatrix_csv 
close rgpmatrix_csv 
