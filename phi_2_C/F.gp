#set term post eps enhanced "Times Roman" 20
set term post eps enhanced color "Times Roman" 18
set output "F.eps"
#set output "n21.ps"

#set terminal vdi 13 0 fff f00 f0 f ff f0f
#set logscale x
#set mxtics
set key top
set xlabel "{/symbol=24 a}"
set ylabel "{/symbol=24 F}" 
#set xzeroaxis
#set grid
set xtics 
set grid
set mxtics 9
show mxtics
#set ytics 0.01
#set with lines
set logscale x

plot [:][:]'F.dat' title "{/symbol F}" w l lt 2 

