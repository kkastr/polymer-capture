set N [lindex $argv 0]
set rseed  [lindex $argv 1]
set filename [lindex $argv 2]
set vis [lindex $argv 3]


t_random seed $rseed



set boxx 400
set boxy 400
set boxz 400 

set cx [expr $boxx/8.0]
set cy [expr $boxy/8.0]
set cz [expr $boxz/8.0]
set tmax 10000
set nmax 10
set temp 1; set gamma 1.0; set gamma_equilibration 0.1



setmd box_l $boxx $boxy $boxz
setmd periodic 1 1 1
setmd time_step 0.01; setmd skin 0.4
thermostat langevin $temp $gamma


#FENE potential
set k_fene 30.0
set r_fene 1.5

#Angular Potential - Harmonic
set k_angle [expr 10.0 / $temp]
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

set t_trans 0
set trans_flag 0
set fixed_N [expr 0]
set equil_time [expr 10.0 * $N]
set t_pore 1
set z_line [expr $cz - $t_pore/2]
set force [expr -5.0]
set n_attempt 0
set illegal_mov 0 



part [expr $N] pos $cx $cy $cz type 99 fix
constraint pore center [expr $cx] [expr $cy] [expr $cz] axis 0 0 1 radius 1.5 length 1 type 1

for { set i 0 } { $i < $N } { incr i } {
	set x [expr $cx]
	set y [expr $cy]
	set z [expr $cz + 0.97*$i]
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

set part_pos_contact [open "data/${filename}_$N/part_pos_contact-$N-$rseed.xyz" "a"]
set part_pos_trans [open "data/${filename}_$N/part_pos_trans-$N-$rseed.xyz" "a"]
set part_pos_z [open "data/${filename}_$N/part_pos_z-$N-$rseed.xyz" "a"]

set flag 0
set t 0	
set n 0
set rg_flag 0
set fail 0
set success 0
set trans_flag 0

while {$flag == 0} {

	for { set i 0 } { $i < $N } { incr i } {
		set x [expr $cx]
		set y [expr $cy]
		set z [expr $cz + 0.97*$i -0.97]
		part $i pos $x $y $z type 0 
	}

	part $fixed_N fix

	thermostat langevin $temp $gamma_equilibration

	for {set i 0} {$i < $equil_time} {incr i} {
	
	    integrate 100
	    imd positions
	}
	thermostat langevin $temp $gamma

	part $fixed_N unfix


	for {set i 0} {$i < $N} {incr i} {
		part $i ext_force $force 1 1
	}

	#puts "I'm back in the first while"
	while {1} {
		set z_list {}
		set r_list {}		
		if { $n > $nmax } {
			puts "n > nmax"
			set flag 1
			puts "$fail $success"
			break
		}

		for {set i 0} { $i < $N } {incr i} {
			set x [lindex [part $i print pos] 0]
			set y [lindex [part $i print pos] 1]
			set z [lindex [part $i print pos] 2]
			set r [expr sqrt(($x-$cx)*($x-$cx) + ($y-$cy)*($y-$cy) + ($z-$cz)*($z-$cz))]
			#puts $r
			lappend z_list $z
			lappend r_list $r
			#puts "hi, I'm getting particle positions"
		}
		set z_min [::tcl::mathfunc::min {*}$z_list]
		set z_max [::tcl::mathfunc::max {*}$z_list]
		set r_min [::tcl::mathfunc::min {*}$r_list]
		set r_max [::tcl::mathfunc::max {*}$r_list]


		if {$r_min > 30.0} {
			puts "Dist greater than r = 30 from pore"
			incr fail
			break
		}


		if {$z_min > $z_line && $trans_flag == 1} {
			puts "zmin greater than zline"
			set trans_flag 0
		}

		if {$z_min < $z_line} {
			
			if {$trans_flag == 0 } {
				#puts "inside translocation if"
				set t_thread $t
				
				if { $n_attempt == 0 } {
					puts "n_attempt $n_attempt"
					set t_first_thread $t
				
	        
					set n_attempt [expr $n_attempt + 1]
				}
				set rg_calc_trans [analyze rg 0 1 $N]

				
				puts $part_pos_trans "$N"
	    		puts $part_pos_trans "Position trans starting $t_thread"
	    		set n_cis 0
				for {set l 0} {$l < $N} {incr l} {
					puts $part_pos_trans "a$l [part $l print pos]"
					set z [lindex [part $l print pos] 2]
					if { $z <= $z_min } {
						set thread_index $l
					}
				}
				set trans_flag [expr $trans_flag + 1 ]
			}
		}
		#puts $z_max
		if {$z_max < $z_line} {
			puts "zmax less than zline"
			set t_last_thread $t
			puts $t_last_thread
			set t_trans [expr $t_last_thread - $t_thread]
			set trans_time [open "data/${filename}_$N/trans_time_$N-$rseed.dat" "a"]
			puts $trans_time "$t_trans $N $t_first_thread $t_last_thread $t_thread"
			close $trans_time
			set rg_trans [open "data/${filename}_$N/rg_trans-$N-$rseed.dat" "a"]
			puts $rg_trans "$rg_calc_trans $N $t_thread"
			close $rg_trans
			set thread_indexing [open "data/${filename}_$N/thread_indexing-$N-$rseed.dat" "a"]
			puts $thread_indexing "$thread_index"
			close $thread_indexing
	      	set n_attempt 0
	      	set rg_flag 0  
			set n [expr $n + 1.0]
			break
		}
		if {$t_trans != 0} {
			set t_trans 0 
		}
		integrate 100
		imd positions
		
		incr t
	}
}
close $part_pos_contact
close $part_pos_trans
close $part_pos_z

