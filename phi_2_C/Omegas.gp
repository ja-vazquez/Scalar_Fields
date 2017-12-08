#set term post eps enhanced "Times Roman" 20
set term post eps enhanced color "Times Roman" 18
set output "Omegas.eps"
#set output "n21.ps"

#set terminal vdi 13 0 fff f00 f0 f ff f0f
#set logscale x
#set mxtics
set key top
set xlabel "{/symbol=24 a}"
set ylabel "{/Symbol=24 W}" 
#set xzeroaxis
#set grid
set xtics 
set grid
set mxtics 9
show mxtics
#set ytics 0.01
#set with lines
set logscale x

plot [:][:]'Ode_phi2.dat' title "{/Symbol W}_{DE}" w l lt 2,'Odm_phi2.dat' title "{/Symbol W}_{DM}" w l lt 3,'Ob_phi2.dat' title "{/Symbol W}_{b}" w l lt 1,'Or_phi2.dat' title "{/Symbol W}_{r}" w l lt 5,'On_phi2.dat' title "{/Symbol W}_{/Symbol nu}" w l lt 4 

