set N      [lindex $argv 0]
set rseed  [lindex $argv 1]
set filename [lindex $argv 2]
set vis [lindex $argv 3]

#set rseed 18
t_random seed    $rseed


#puts $N


#set N 20

set boxx 400
set boxy 400
set boxz 400
set cx [expr $boxx/8.0]
set cy [expr $boxy/8.0]
set cz [expr $boxz/8.0]
set tmax 10000
set nmax 5
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
set fixed_N [expr $N/2]
set equil_time [expr 10.0 * $N]
set t_pore 1
set z_line [expr $cz - $t_pore/2]
set force [expr -10.0]
set n_attempt 0
set illegal_mov 0 


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

for { set i 0 } { $i < $N } { incr i } {
	set x [expr $cx - $N/2 + $i]
		set y [expr $cy]
	set z [expr $cz  + 100]
	part $i pos $x $y $z type 0 


	if { $i > 0 } {
		part $i bond 0 [expr $i - 1]
	}
	if { $i > 1} {
		part [expr $i - 1] bond 1 $i [expr $i - 2]
	}
}


set rlist {}

for {set j 0} {$j < $N} {incr j} {
	set x [lindex [part $j print pos] 0]
	set y [lindex [part $j print pos] 1]
    set z [lindex [part $j print pos] 2]
    set r [expr sqrt($x*$x + $y*$y + $z*$z)]
    lappend rlist $r
}

set rmin [::tcl::mathfunc::min {*}$rlist]
set rmax [::tcl::mathfunc::max {*}$rlist]


for {set i 0} {$i < $N} {incr i} {
    set z [lindex [part $i print pos] 2]
    set y [lindex [part $i print pos] 1]
    set x [lindex [part $i print pos] 0]
    set r [expr sqrt($x*$x + $y*$y + $z*$z)]
    if {$r == $rmin} {
      set min_part $i
    } 
}


set randlist [randgen]


set px [lindex [part $min_part print pos] 0]
set py [lindex [part $min_part print pos] 1]
set pz [lindex [part $min_part print pos] 2]


set mx [expr $cx + 40*[lindex $randlist 0]]
set my [expr $cy + 40*[lindex $randlist 1]]
set mz [expr $cz + 40*[lindex $randlist 2]]

puts "$mx $my $mz"
set mr2  [expr ($mx - $cx)*($mx -$cx) + ($my - $cy)*($my -$cy) + ($mz - $cz)*($mz -$cz)]
set mr [expr sqrt($mr2)]
puts $mr
#set mindist [analyze mindist 0 99]
#puts $mindist
puts $min_part
part $min_part pos $mx $my $mz type 0

set tx [lindex [part $min_part print pos] 0]
set ty [lindex [part $min_part print pos] 1]
set tz [lindex [part $min_part print pos] 2]


#part $min_part fix
# for {set i 0} {$i < $N} {incr i} {
# 	part $i bond delete
# }
set tr [expr sqrt($tx*$tx + $ty*$ty + $tz*$tz)]

for { set i [expr 0] } { $i < $N } { incr i } {
	if {$i != $min_part} {
		set ix [expr [lindex [part $i print pos] 0] - $px + $tx] 
		set iy [expr [lindex [part $i print pos] 1] - $py + $ty] 
		set iz [expr [lindex [part $i print pos] 2] - $pz + $tz] 
		set ir [expr sqrt($ix*$ix + $iy*$iy + $iz*$iz)]
		#puts "$ix $iy $iz"
		puts $ir
		part $i pos $ix $iy $iz type 0 
	}
}

set zlist {}

for {set i 0} {$i < $N} { incr i} {
	set z [lindex [part $i print pos] 2]
	if {$z < [expr $cz + $lj_cutoff]} {
		puts "Nooooooo - Overlap"
		set illegal_mov 1
	}
}





