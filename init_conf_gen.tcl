set N      [lindex $argv 0]
set rseed  [lindex $argv 1]

t_random seed    $rseed

#md initialization.#35211123456 issue seed	
#set rseed 3021123456
#set N 100

 
puts $N

set box_l 1500

setmd box_l $box_l $box_l $box_l

setmd periodic 1 1 1

puts $box_l

set box_center 100.0

puts $box_center




setmd time_step 0.01; setmd skin 0.4
set temp 1.0; set gamma 1; set gamma_equilibration 0.1
thermostat langevin $temp $gamma

#potential parameters.

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
set lj_offset 0.0

#nanopore sizes.
set r_ex [expr 1.8 * $sigma]
set r_ent [expr 5.5 * $sigma]
set r_gap [expr 99.5 * $sigma]
set l_ex [expr  1.5 * $sigma]
set l_ent [expr 4.5 * $sigma]
set l_gap [expr 19.5 * $sigma]

#particle parameters.
#set fixed_part [expr int(floor($N/2))]
#puts "$fixed_part"
#number of fixed particles, 10% of polymer length
set fixed_part [expr int ($N / 2.0)] 
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
set rg_flag 0

set equil_time [expr 100*$N]


#interactions between types of particles, fene and shifted LJ.
inter 0 fene $k_fene $r_fene
inter 1 angle $k_angle $pi
inter 0 0 lennard-jones $eps $sigma $lj_cutoff $lj_shift $lj_offset
inter 1 0 lennard-jones $eps $sigma $lj_cutoff $lj_shift $lj_offset

#particle initialization.
#force on particles - enumerated in field.dat defined in forces.c format magnitude stuff stuff (keep stuff=1).


set z_top [expr 2.0 * $box_center + 2.0 * $l_gap + $l_ex + 2.0 * $l_ent + 2.5 ]

set xy_offset [expr 30/sqrt(2)]
for {set i 0} {$i < [expr $N]} {incr i} {
  set posx [expr   $box_center  + ( 0.97/sqrt(2) * $i ) + $xy_offset - (($N)/(2.)-1)/sqrt(2) ]
  set posy [expr   $box_center  + ( 0.97/sqrt(2) * $i ) + $xy_offset - (($N)/(2.)-1)/sqrt(2) ]
	set posz [expr  2.0 * $box_center +  2.0 * $l_gap + $l_ex + 2.0 * $l_ent + $lj_cutoff + 200]
	part [expr $i] pos $posx $posy $posz type $part_type ext_force 0 0 0
}


#bond initialization.
for {set k 1} {$k < $N} {incr k} {
    part [expr $k-1] bond 0 $k
   }
for {set k 2} {$k < $N} {incr k} {
    part [expr $k-1] bond 1 [expr $k-2] $k
}

part [expr $N] pos [expr $box_center + $xy_offset] [expr $xy_offset + $box_center] [expr 2.0 * $box_center + 2.0 * $l_gap + $l_ex + $l_ent + 2.0] type 98 fix

constraint pore center [expr $box_center + $xy_offset] [expr $xy_offset + $box_center] [expr 2.0 * $box_center] axis 0 0 1 radius $r_ex length $l_ex type $pore_type
constraint pore center [expr $box_center + $xy_offset] [expr $xy_offset + $box_center] [expr 2.0 * $box_center + 2.0 * $l_gap + $l_ex + $l_ent] axis 0 0 1 radius $r_ent length $l_ent type $pore_type





# set vmd "yes"

# if {$vmd == "yes"} {
#    prepare_vmd_connection init_config_gen 3000
#    exec sleep 4
#     imd positions
# }


#Equilibration. 

part $fixed_part fix

thermostat langevin $temp $gamma_equilibration
for {set i 0} {$i < $equil_time} {incr i} {
    integrate 100
    imd positions
}
thermostat langevin $temp $gamma

puts "equilibrated $n"

part $fixed_part unfix


set z_list {}
set x_list {}
set y_list {}

for {set j 0} {$j < $N} {incr j} {
    set z [lindex [part $j print pos] 2]
    lappend z_list $z
}

set z_min [::tcl::mathfunc::min {*}$z_list]
set z_max [::tcl::mathfunc::max {*}$z_list]

for {set i 0} {$i < $N} {incr i} {
    set z [lindex [part $i print pos] 2]
    set y [lindex [part $i print pos] 1]
    set x [lindex [part $i print pos] 0]
    if {$z == $z_min} {
      set min_part $i
    } 
}

set min_partx [lindex [part $min_part print pos] 0]
set min_party [lindex [part $min_part print pos] 1]

set rmin [analyze mindist 0 98]

set rg2 [expr  ( 1.0 / 3.0 ) * $k_angle * ($N*0.97 ) - pow($k_angle, 2) + ( 2.0 * ( pow($k_angle, 3)/($N * 0.97 ) ) ) * ( 1.0 - ($k_angle/( $N * 0.97 )) *( 1.0 - exp(-($N * 0.97)/$k_angle))) ]

set rg [expr sqrt($rg2)]
puts [expr 0.1* $rg]
for {set i 0} {$i < [expr $N]} {incr i} {
  set x_transport [expr [lindex [part $i print pos] 0] - ( $min_partx - 100 - $xy_offset ) ]
  set y_transport [expr [lindex [part $i print pos] 1] - ( $min_party - 100 - $xy_offset ) ]
  set z_transport [expr [lindex [part $i print pos] 2] - ($z_min - 0.0 ) + $z_top]
  part [expr $i] pos $x_transport $y_transport $z_transport ext_force $force 1 1
} 

puts [part $min_part print pos]


set flag 0

while { $n < $n_max } {

	imd positions
	integrate $int_steps

  if {$flag == 1} {
    for {set i 0} {$i < [expr $N]} {incr i} {
      set posx [expr   $box_center  + ( 0.97/sqrt(2) * $i ) + $xy_offset - (($N)/(2.)-1)/sqrt(2) ]
      set posy [expr   $box_center  + ( 0.97/sqrt(2) * $i ) + $xy_offset - (($N)/(2.)-1)/sqrt(2) ]
      set posz [expr  2.0 * $box_center +  2.0 * $l_gap + $l_ex + 2.0 * $l_ent + $lj_cutoff + 200]
      part [expr $i] pos $posx $posy $posz type $part_type ext_force 0 0 0
    }


  part $fixed_part fix

  thermostat langevin $temp $gamma
  for {set i 0} {$i < [expr $equil_time/100.]} {incr i} {
      integrate 10
      imd positions
  }
  thermostat langevin $temp $gamma_equilibration
  for {set i 0} {$i < $equil_time} {incr i} {
      integrate 100
      imd positions
  }
  thermostat langevin $temp $gamma

  puts "equilibrated $n"

  part $fixed_part unfix


  set z_list {}

  for {set j 0} {$j < $N} {incr j} {
      set z [lindex [part $j print pos] 2]
      lappend z_list $z
  }

  set z_min [::tcl::mathfunc::min {*}$z_list]
  set z_max [::tcl::mathfunc::max {*}$z_list]
  

  for {set i 0} {$i < $N} {incr i} {
      set z [lindex [part $i print pos] 2]
      set y [lindex [part $i print pos] 1]
      set x [lindex [part $i print pos] 0]
      if {$z == $z_min} {
        set min_part $i
      } 
  }

  set min_partx [lindex [part $min_part print pos] 0]
  set min_party [lindex [part $min_part print pos] 1]

  set rmin [analyze mindist 0 98]
  
  for {set i 0} {$i < [expr $N]} {incr i} {
    set x_transport [expr [lindex [part $i print pos] 0] - ( $min_partx - 100 -$xy_offset ) ]
    set y_transport [expr [lindex [part $i print pos] 1] - ( $min_party - 100 -$xy_offset) ]
    set z_transport [expr [lindex [part $i print pos] 2] - ($z_min - 0.0 ) + $z_top]
    part [expr $i] pos $x_transport $y_transport $z_transport ext_force $force 1 1
  } 
  set flag 0
  }

  set z_list {}

  for {set j 0} {$j < $N} {incr j} {
    set z [lindex [part $j print pos] 2]
    lappend z_list $z
  }
  

  set z_min [::tcl::mathfunc::min {*}$z_list]
  set z_max [::tcl::mathfunc::max {*}$z_list]
  set rmin [analyze mindist 0 98]

  if { $z_min > [expr $z_top + 10.0]} {
    puts "z_min > 260"
    for {set i 0} {$i < [expr $N]} {incr i} {
      set posx [expr   $box_center  + ( 0.97/sqrt(2) * $i ) + $xy_offset - (($N)/(2.)-1)/sqrt(2) ]
      set posy [expr   $box_center  + ( 0.97/sqrt(2) * $i ) + $xy_offset - (($N)/(2.)-1)/sqrt(2) ]
      set posz [expr  2.0 * $box_center +  2.0 * $l_gap + $l_ex + 2.0 * $l_ent + $lj_cutoff + 200]
      part [expr $i] pos $posx $posy $posz type $part_type ext_force 0 0 0
    }

  part $fixed_part fix

  thermostat langevin $temp $gamma
  for {set i 0} {$i < [expr $equil_time/100.]} {incr i} {
      integrate 10
      imd positions
  }
  thermostat langevin $temp $gamma_equilibration
  for {set i 0} {$i < $equil_time} {incr i} {
      integrate 100
      imd positions
  }
  thermostat langevin $temp $gamma

  puts "equilibrated $n"

  part $fixed_part unfix


  set z_list {}

  for {set j 0} {$j < $N} {incr j} {
      set z [lindex [part $j print pos] 2]
      lappend z_list $z
  }

  set z_min [::tcl::mathfunc::min {*}$z_list]
  set z_max [::tcl::mathfunc::max {*}$z_list]

  for {set i 0} {$i < $N} {incr i} {
      set z [lindex [part $i print pos] 2]
      set y [lindex [part $i print pos] 1]
      set x [lindex [part $i print pos] 0]
      if {$z == $z_min} {
        set min_part $i
      } 
  }

  set min_partx [lindex [part $min_part print pos] 0]
  set min_party [lindex [part $min_part print pos] 1]

  set rmin [analyze mindist 0 98]
  
  for {set i 0} {$i < [expr $N]} {incr i} {
    set x_transport [expr [lindex [part $i print pos] 0] - ( $min_partx - 100 -$xy_offset) ]
    set y_transport [expr [lindex [part $i print pos] 1] - ( $min_party - 100 -$xy_offset) ]
    set z_transport [expr [lindex [part $i print pos] 2] - ($z_min - 0.0 ) + $z_top]
    part [expr $i] pos $x_transport $y_transport $z_transport ext_force $force 1 1
  } 
  }

  if { $z_min < [expr 2.0 * $box_center + 2.0 * $l_gap + $l_ex + $l_ent + 2.0] } {
   
    set n [expr $n + 1.0]
    set flag 1

  }
 
  set out [open "|gzip -c - > data/N_gen_fixfene_$N/checkpoint_N_[expr $N]_seed[expr $rseed]_num[expr int($n)].block.gz" "w"]
  blockfile $out write variable all
  blockfile $out write interactions
  blockfile $out write random
  blockfile $out write bitrandom
  blockfile $out write particles "id pos type" all
  blockfile $out write bonds all
  blockfile $out write configs
  close $out
}	


