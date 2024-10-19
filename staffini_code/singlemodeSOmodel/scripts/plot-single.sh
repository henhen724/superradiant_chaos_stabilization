#!/bin/bash

  
P=0.400
OmegaTil=1.520
#plotkind="chaos"
plotkind="evolution"


if [ "$plotkind" == "chaos" ]; then

############# CHAOS PLOT #############

gnuplot <<EOF
dirname="Esslinger_single_OmTil=${OmegaTil}_P=${P}"
rt=dirname."/Evo-Pump_${P}-OmegaTil_${OmegaTil}.dat"

######################################################################
tmin=0
tmax=5000
range=8
fname="Omtil=${OmegaTil}_P=${P}_chaos-map"
xrange=range
yrange=range
######################################################################

labstr="P=${P}, {/Symbol w}=${OmegaTil} E_0"

set out fname."-orig.eps"
set term postscript enhanced color eps "Times-Roman,24" solid

rl="Re[{/Symbol l}]"
il="Im[{/Symbol l}]"

set xlabel rl
set ylabel il

set xrange [-xrange to xrange]
set yrange [-yrange to yrange]

set pal mod HSV
set pal def (0 0 1 1, 1 1 1 1)
set cbrange [tmin to tmax]

set label at graph 0.05, 0.95 labstr front
set grid front; unset grid
plot x w l lc pal, x w l lc rgb("magenta")
unset label

load "~/Code/utilities/png.gp"
set out "a.png"

sp=tmin*500
ep=tmax*500


plot rt every ::sp::ep  u 2:3:1 w p lc pal lw 3

set macros # allows to use gnuplot variable names in system commands
!insert-graphs-magenta.pl @fname a.png

epsname=fname.".eps"
!epstopdf @epsname
EOF

echo "Created eps/pdf file Omtil=${OmegaTil}_P=${P}_chaos-map"


else

############# TIME EVO PLOT #############

gnuplot <<EOF

dirname="Esslinger_single_OmTil=${OmegaTil}_P=${P}"
rt=dirname."/Evo-Pump_${P}-OmegaTil_${OmegaTil}.dat"

lmax=8
tmax=5E3

fname="Omtil=${OmegaTil}_P=${P}_time-evo"

labstr='~{/Symbol w}{0.7\~}=${P}, P=${OmegaTil}'

set term postscript enhanced color eps "Times-Roman,24" solid
set out fname."-orig.eps"

xlab="Time, t (ms)"

unset key

set xrange [0 to tmax]
set xtics 1E3

##################################################
set size 1,0.9
set origin 0,0

unset label 1

set xlabel xlab
set format x "%g"

rl="Re[{/Symbol l}]"
il="Im[{/Symbol l}]"

set label at graph 0.05, 0.95 labstr front
set grid front; unset grid

set key top right 

set ylabel rl.", ".il offset 0.8, 0
set yrange [-lmax to lmax]
set ytics 2
plot x w l lc rgb("magenta"), \
     rt every 1::1::2 u (0):(0) w l lw 2 lt 1 t "Re[{/Symbol l}]", \
     rt every 1::1::2 u (0):(0) w l lw 2 lt 3 t "Im[{/Symbol l}]"

######################################################################
load "~/bin/png.gp"
set term png size 4096, 1024 crop truecolor 

set out "a.png"
set xrange [0 to tmax]
set yrange [-lmax to lmax]
plot rt u 1:3 w l lw 4 lt 3 t il, \
     rt u 1:2 w l lw 4 lt 1 t rl


set macros
!insert-graphs-magenta.pl  @fname a.png

epsname=fname.".eps"
!epstopdf @epsname
EOF

echo "Created eps/pdf file Omtil=${OmegaTil}_P=${P}_time-evo"

fi


