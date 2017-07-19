set N [lindex $argv 0]
set rseed  [lindex $argv 1]
set filename [lindex $argv 2]
set vis [lindex $argv 3]


t_random seed $rseed



set boxx 500
set boxy 500
set boxz 500 

set cx [expr $boxx/2.0]
set cy [expr $boxy/2.0]
set cz [expr $boxz/2.0]

set tmax 10000
set nmax 10
set temp 0.1; set gamma 1.0; set gamma_equilibration 0.1



setmd box_l $boxx $boxy $boxz
setmd periodic 1 1 1
setmd time_step 0.01; setmd skin 0.4
thermostat langevin $temp $gamma


#FENE potential
set k_fene 30.0
set r_fene 1.5

#Angular Potential - Harmonic
set k_angle [expr $temp * 10.0 ]
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
set fixed_N [expr $N/2]
set equil_time [expr $N * $N]
set tpore 1
set z_line [expr $cz - $tpore/2]
set force [expr -1.0]
set n_attempt 0
set illegal_mov 0 
set rpore 1.5
set cutofftime 1e6
set cutoffdist 2
set transportdist 20



part [expr $N] pos $cx $cy $cz type 99 fix
constraint pore center [expr $cx] [expr $cy] [expr $cz] axis 0 0 1 radius $rpore length $tpore type 1

for { set i 0 } { $i < $N } { incr i } {
	set x [expr $cx]
	set y [expr $cy]
	set z [expr $cz  + $i]
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
#set part_pos_z [open "data/${filename}_$N/part_pos_z-$N-$rseed.xyz" "a"]
set metric_csv [open "data/${filename}_$N/metric-${filename}-$N-$rseed.csv" "a"]

#puts $metric_csv "N,tau,rgxtrans,rgytrans,rgztrans,rgxycorrtrans,rgxequil,rgyequil,rgzequil,rgxycorrequil,tfirstthread,tthread,tlastthread"


set flag 0
set t 0
set n 0
set rg_flag 0
set fail 0
set success 0
set trans_flag 0
set stuck 0


set position_flag 0
while {$flag == 0} {


	if {$position_flag == 1} {
		for {set i 0} {$i < $N} {incr i} {
			part $i ext_force 0 0 0
		}

	for { set i 0 } { $i < $N } { incr i } {
		set x [expr $cx]
		set y [expr $cy]
		set z [expr $cz  + $i]
		part $i pos $x $y $z type 0 
	}
	set position_flag 0
	}
	part 0 fix
	#part 1 fix

	thermostat langevin $temp $gamma_equilibration
	for {set i 0} {$i < $equil_time} {incr i} {
	    
	    integrate 100
	    imd positions
	
	}
	thermostat langevin $temp $gamma

	part 0 unfix
	#part 1 unfix

	puts "equilibrated."
	
	set rg_at_equil [analyze rg 0 1 $N]

	
	#puts "I'm back in the first while"

	while {1} {
		for {set i 0} {$i < $N} {incr i} {
			part $i ext_force $force 1 1
		}

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
		
			lappend z_list $z
			lappend r_list $r
			
		}
		set z_min [::tcl::mathfunc::min {*}$z_list]
		set z_max [::tcl::mathfunc::max {*}$z_list]
		set r_min [::tcl::mathfunc::min {*}$r_list]
		set r_max [::tcl::mathfunc::max {*}$r_list]


		if {$r_min > $cutoffdist} {
			puts "Dist greater than r = $cutoffdist from pore"
			incr fail
			set position_flag 1
			break
		}

		if {$t > $cutofftime} {
			puts "cut off time exceeded - stuck event"
			incr stuck
			set position_flag 1 
			break
		}


		if {$z_min > $z_line && $trans_flag == 1} {
			puts "zmin greater than zline"
			set trans_flag 0
		}

		if {$z_min < $z_line} {
			
			if {$trans_flag == 0 } {
				#puts "inside translocation if"
				set t_last_thread $t
				
				if { $n_attempt == 0 } {
					puts "n_attempt $n_attempt"
					set t_first_thread $t
				
	        
					set n_attempt [expr $n_attempt + 1]
				}
				set rg_calc_trans [analyze rg 0 1 $N]

				
				puts $part_pos_trans "$N"
	    		puts $part_pos_trans "Position trans starting $t_last_thread"
	         
				for {set l 0} {$l < $N} {incr l} {
					puts $part_pos_trans "a$l [part $l print pos]"
				}
				set trans_flag [expr $trans_flag + 1 ]
			}
		}
		#puts $z_max
		if {$z_max < $z_line} {
			puts "zmax less than zline"
			set t_thread $t
			puts $t_last_thread
			set t_trans [expr $t_thread - $t_last_thread]

			set metric_csv [open "data/${filename}_$N/metric-${filename}-$N-$rseed.csv" "a"]
			puts $metric_csv "$N,$t_trans,$rg_calc_trans,$rg_at_equil,$t_first_thread,$t_thread,$t_last_thread,$fail,$stuck"
			close $metric_csv

	      	set n_attempt 0
	      	set rg_flag 0  
			set n [expr $n + 1.0]
			set position_flag 1
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
#close $part_pos_z
#close $metric_csv
