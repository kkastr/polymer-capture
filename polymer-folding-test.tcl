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
set temp 0.00; set gamma 1.0; set gamma_equilibration 0.01



setmd box_l $boxx $boxy $boxz
setmd periodic 1 1 1
setmd time_step 0.01; setmd skin 0.4
thermostat langevin $temp $gamma


#FENE potential
set k_fene 30.0
set r_fene 1.5

#Angular Potential - Harmonic
set k_angle [expr 10.0]
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
set equil_time [expr 20.0 * $N]
set tpore 1
set z_line [expr $cz - $tpore/2]
set force [expr -3.0]
set n_attempt 0
set illegal_mov 0 
set rpore 2.1
set cutofftime 1e6
set transportdist [expr 60]
set cutoffdist [expr $transportdist + 60]
set contactdist [expr 3]


part [expr $N] pos $cx $cy $cz type 99 fix
constraint pore center [expr $cx] [expr $cy] [expr $cz] axis 0 0 1 radius $rpore length $tpore type 1



for { set i 0 } { $i < $N } { incr i } {
	set x [expr $cx - $N/2 + $i]
	set y [expr $cy]
	set z [expr $cz  + 2.5]
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

for {set i 0} {$i < $N} {incr i} {
	part $i ext_force $force 1.6 1.6
}

while {1} {
	integrate 100
	imd positions

	incr t
}

