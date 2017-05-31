set N      [lindex $argv 0]
set rseed  [lindex $argv 1]
set filename [lindex $argv 2]
set vis [lindex $argv 3]

t_random seed    $rseed


puts $N

set box_l 1000

setmd box_l $box_l $box_l $box_l

setmd periodic 1 1 1

#puts $box_l

set box_center 100.0

puts $box_center

setmd time_step 0.01; setmd skin 0.4
set temp 1.0; set gamma 1.0; set gamma_equilibration 0.1
thermostat langevin $temp $gamma

#potential parameters.

#FENE potential
set k_fene 20.0
set r_fene [expr 2.0 * 0.97]

#Angular Potential - Harmonic
set k_angle 6.0
set pi 3.14159

#Shifted Lennard-Jones
set eps 1.0
set sigma 1.0
set lj_cutoff 1.12246
set lj_shift 0.25
set lj_offset 0.0

#nanopore sizes.
set r_ex [expr 1.8 * $sigma]
set r_ent [expr 1.8 * $sigma]
set r_gap [expr 100 * $sigma]
set l_ex [expr  1.5 * $sigma]
set l_ent [expr 4.5 * $sigma]
set l_gap [expr 19.5 * $sigma]

#particle parameters.
#number of fixed particles, 5% of polymer
set fixed_part [expr floor(( 0.05 * $N )) + 1] 

set part_type 0
set pore_type 1

#integration parameters.
set int_steps 100
set n 0
set n_max 10.0
set trans_flag 0

#force parameters.
set force [expr -3.848]

#translocation time parameters
set t 0
set t_trans 0
set n_attempt 0 
set t_thread 0
set t_first_thread 0
set t_last_thread 0 
set t_contact 0
set z_line 198
set z_top 242
set rg_flag 0


#equilibration time
set equil_time [expr 100 * $N]


#pore offset
set xy_offset [expr 30.0/sqrt(2.0)]

#increased seeding into cavity
set shift [expr 5]
#30 / sqrt(2) for off-axis
#interactions between types of particles, fene and shifted LJ.
inter 0 fene $k_fene $r_fene
inter 1 angle $k_angle $pi
inter 0 0 lennard-jones $eps $sigma $lj_cutoff $lj_shift $lj_offset
inter 1 0 lennard-jones $eps $sigma $lj_cutoff $lj_shift $lj_offset

#particle initialization.

proc particle_positions {N box_center xy_offset l_gap l_ex shift part_type force} {
	for {set i 0} {$i < [expr $N]} {incr i} {
		set posx [expr $box_center + $xy_offset]
		set posy [expr $box_center + $xy_offset]
		set posz [expr ( 0.97 * $i) + 2.0 * $box_center + 2.0 * $l_gap + $l_ex - $shift ]
		part [expr $i] pos $posx $posy $posz type $part_type ext_force $force 1 1
	}
}

particle_positions $N $box_center $xy_offset $l_gap $l_ex $shift $part_type $force

#bond initialization.

for {set k 1} {$k < $N} {incr k} {
    part [expr $k-1] bond 0 $k
   }
for {set k 2} {$k < $N} {incr k} {
    part [expr $k-1] bond 1 [expr $k-2] $k
}



#dummy particles

part [expr $N] pos $box_center $box_center [expr 2.0 * $box_center] type 99 fix
#part [expr $N + 1] pos $box_center $box_center [expr 2.0 * $box_center - $l_ex] type 98 fix
#part [expr $N + 2] pos $box_center $box_center [expr 2.0 * $box_center + 2.0 *$l_gap + $l_ex] type 99 fix
#part [expr $N + 2] pos $box_center $box_center [expr 2.0 * $box_center + $l_ex] type 98 fix
#part [expr $N + 4] pos $box_center $box_center [expr 2.0 * $box_center +  2.0 * $l_gap + $l_ex + 2.0 * $l_ent] type 99 fixi


#constraint initialization, 3 nanopores. 
constraint pore center [expr $box_center] [expr $box_center] [expr 2.0 * $box_center] axis 0 0 1 radius $r_ex length $l_ex type $pore_type
constraint pore center [expr $box_center] [expr $box_center] [expr 2.0 * $box_center + $l_gap + $l_ex] axis 0 0 1 radius $r_gap length $l_gap type $pore_type
constraint pore center [expr $box_center + $xy_offset] [expr $box_center + $xy_offset] [expr 2.0 * $box_center + 2.0 * $l_gap + $l_ex + $l_ent] axis 0 0 1 radius $r_ent length $l_ent type $pore_type




if {$vis==1} {
	set vmd "yes"

	if {$vmd == "yes"} {
		prepare_vmd_connection cavity 3000
		exec sleep 4
		imd positions
	}
}



#Equilibration. 

proc equilibration {fixed_part temp gamma_equilibration gamma equil_time} {
	
	for {set i 0} { $i < [expr $fixed_part]} {incr i} {
		part $i fix
	}

	thermostat langevin $temp $gamma_equilibration
	for {set i 0} {$i < $equil_time} {incr i} {
	    integrate 100
	    imd positions
	}
	thermostat langevin $temp $gamma


	for {set i 0} { $i < [expr $fixed_part]} {incr i} {
		part $i unfix
	}

} 


equilibration $fixed_part $temp $gamma_equilibration $gamma $equil_time

puts "equilibrated."
set rg_equil [open "data/N_${filename}_$N/rg_equil-$N-$rseed" "a"]
set rg_calc [analyze rg 0 1 $N]
puts $rg_equil "$rg_calc $N"
close $rg_equil



#force init after equilibration

# for {set i 0} {$i < $N} {incr i} {
# 	part $i ext_force $force 1 1
# }

set part_pos_contact [open "data/N_${filename}_$N/part_pos_contact-$N-$rseed.xyz" "a"]
set part_pos_trans [open "data/N_${filename}_$N/part_pos_trans-$N-$rseed.xyz" "a"]
set part_pos_z [open "data/N_${filename}_$N/part_pos_z-$N-$rseed.xyz" "a"]

set top_thread 0
set top_thread_flag 1

while { $n < $n_max } {
	
	integrate $int_steps
	imd positions
	set z_list {}


	for {set j 0} {$j < $N} {incr j} {
		set z [lindex [part $j print pos] 2]
		lappend z_list $z
		if {$z <= $z_top && $top_thread_flag==1} {
			
			set top_thread $j
			puts $top_thread
      		set top_thread_flag 0
        }	
	}
	
	set z_min [::tcl::mathfunc::min {*}$z_list]
	set z_max [::tcl::mathfunc::max {*}$z_list]
	set rmin [analyze mindist 0 99]

	if { $z_min > 252.0 } {
		puts "Above z = 252"
		
		particle_positions $N $box_center $xy_offset $l_gap $l_ex $shift $part_type $force
		
		equilibration $fixed_part $temp $gamma_equilibration $gamma $equil_time

		puts "Failed event, re - equilibrated."
		set rg_equil [open "data/N_${filename}_$N/rg_equil-$N-$rseed" "a"]
		set rg_nalc [analyze rg 0 1 $N]
		puts $rg_equil "$rg_calc $N"
		close $rg_equil
	}

  
	if { $rmin < 4.0 && $rg_flag == 0 } {
		set t_contact $t
		set rg_calc [analyze rg 0 1 $N]
		set radius_of_gyration [open "data/N_${filename}_$N/radius_of_gyration_$N-$rseed" "a"]
		puts $radius_of_gyration "$rg_calc $N $t_contact"
		close $radius_of_gyration
		
		puts $part_pos_contact "$N"
    	puts $part_pos_contact "Positions After rmin less than four"	
		for {set l 0} {$l < $N} {incr l} {
			puts $part_pos_contact "a$l [part $l print pos]"
		}
		
    	set rg_flag 1
	}
  	if {$z_min < 205 && [expr int($t) % 100] == 0} {
  		puts $part_pos_z "$N"
    	puts $part_pos_z "Positions after z less than 205"
    	for {set j 0} {$j < $N} {incr j} {
      		puts $part_pos_z "a$j [part $j print pos]"
    	}
	}
	if {$z_min > $z_line && $trans_flag == 1} {
		puts "zmin greater than zline"
		set trans_flag 0
	}	
	if {$z_min < $z_line} {
		if {$trans_flag == 0 } {
			set t_thread $t
			
			if { $n_attempt == 0 } {
				puts "n_attempt $n_attempt"
				set t_first_thread $t
			
        
				set n_attempt [expr $n_attempt + 1.0]
			}
			set rg_calc_trans [analyze rg 0 1 $N]

			
			puts $part_pos_trans "$N"
    		puts $part_pos_trans "Position trans starting $t_thread"
    		set n_cis 0
			for {set l 0} {$l < $N} {incr l} {
				puts $part_pos_trans "a$l [part $l print pos]"
				set z [lindex [part $l print pos] 2]
				if { $z >= 252 } {
					set n_cis [expr $n_cis + 1]
				}
				if { $z <= $z_min } {
					set thread_index $l
				}
			}
			set trans_flag [expr $trans_flag + 1 ]
		}
	}
	
	if {$z_max < $z_line} {
		puts "zmax less than zline"
		
		set t_last_thread $t
		puts $t_last_thread
		set t_trans [expr $t_last_thread - $t_thread]
		set n_cis_print [open "data/N_${filename}_$N/n_cis-$N-$rseed.txt" "a"]
		puts $n_cis_print "$n_cis $N $t_thread"
		close $n_cis_print
		set trans_time [open "data/N_${filename}_$N/trans_time_$N-$rseed.dat" "a"]
		puts $trans_time "$t_trans $N $t_first_thread $t_last_thread $t_thread"
    	close $trans_time
    	set rg_trans [open "data/N_${filename}_$N/rg_trans-$N-$rseed.dat" "a"]
      	puts $rg_trans "$rg_calc_trans $N $t_thread"
      	close $rg_trans
      	set thread_indexing [open "data/N_${filename}_$N/thread_indexing-$N-$rseed.dat" "a"]
      	puts $thread_indexing "$thread_index"
      	close $thread_indexing
      	set top_thread_index [open "data/N_${filename}_$N/top_thread_index-$N-$rseed.dat" "a"]
      	puts $top_thread_index "$top_thread"
      	close $top_thread_index
		set n_attempt 0
    	set top_thread_flag 1
		set rg_flag 0  
		set n [expr $n + 1.0]
	}
  
	if {$t_trans != 0} { 
		
		#re-init
		particle_positions $N $box_center $xy_offset $l_gap $l_ex $shift $part_type $force

		equilibration $fixed_part $temp $gamma_equilibration $gamma $equil_time

		puts "Translocated - re - equilibrated."
		set rg_equil [open "data/N_${filename}_$N/rg_equil-$N-$rseed" "a"]
		set rg_calc [analyze rg 0 1 $N]
		puts $rg_equil "$rg_calc $N"
		close $rg_equil
		
		set t_trans 0
	}
	set t [expr $t + 1.0]
	
	#puts $n
	#puts $t
}	

close $part_pos_contact
close $part_pos_trans
close $part_pos_z
