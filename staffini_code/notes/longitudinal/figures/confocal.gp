unset key

set term fig
set out "waist.fig"

set xrange [-1 to 1]
set border 0
unset xtics
unset ytics

plot sqrt(1+x*x) lt 1 lc 1 lw 5 , -sqrt(1+x*x) lt 1 lc 1 lw 5

######################################################################
set term postscript enhanced eps "Times-Roman,44"  
epsname="composite.eps"
set out epsname


set xrange [0 to 1]
ymax=6
set yrange [-ymax to ymax]

set param
set trange [-ymax to ymax]
set zeroaxis lw 4

set multiplot
set size 0.25, 1
# Spacing between panels
mpdx=0.18;
# Offset of first panel
mpx0=0.24


f1(x)=exp(-x*x/(2*5))
f0(x)=exp(-x*x*5/2)
o=3

set label 2 at graph -1.5, 1.05 "b)   Object"

set origin 0, 0
plot f0(t-o), t lt 1 lc rgb("#0000d0") lw 11

set label 1 at graph -2, 0.51 "{/Symbol \336}"
set label 2 at graph 0.3, 1.05 " (1)"

set origin mpx0+0*mpdx, 0
plot f1(t+o), t lt 1 lc rgb("#0000d0") lw 11

set label 1 at graph -1.3, 0.51 "+"
set label 2 " (2)"

set origin mpx0+1*mpdx, 0
plot f0(t+o), t lt 1 lc rgb("#0000d0") lw 11

set label 2 " (3)"

set origin mpx0+2*mpdx, 0
plot f1(t-o), t lt 1 lc rgb("#0000d0") lw 11

set label 2 " (4)"

set origin mpx0+3*mpdx, 0
plot f0(t-o), t lt 1 lc rgb("#0000d0") lw 11

unset multiplot
set macros
!epstopdf @epsname
