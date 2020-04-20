set terminal eps enhanced color
set output "nd_t.eps"

L=10
dt=0.1
nt=200

set multiplot layout 2,1 title 'Hubbard, L=10 sites, Sz=0. Changing from U=-4 to U=-2'


set xlabel 'time'
set ylabel '<n_{up} n_{dw} >'
#set xrange[0:100]
p "nd_t.dat" u ($0*dt):1 w l t ''
set xlabel 'frequency'
set ylabel 'log(Pw)'
set xrange[0:11]
p "nd_w.dat" u ($0*2*pi/(nt*dt)):(log($1)) w lp t 'Power spectral',0 dt 2
set xrange[*:*]

unset multiplot


set multiplot layout 2,1
set xlabel 'time'
set ylabel '<c_i^+ c_1>'
#set xrange[0:100]
p for[i=1:L:2] "cic0_t.dat" u ($0*dt):(column(i)) w l t 'i='.i
p for[i=2:L:2] "cic0_t.dat" u ($0*dt):(column(i)) w l t 'i='.i

unset multiplot

set multiplot layout 2,1

set ylabel 'nk'
p for[i=0:(L/2)] "nk_t.dat" u ($0*dt):(column(i+1)) w l t 'k='.i,0 dt 2 t ''
set xlabel 'frequency'
set ylabel 'log( nk(w) )+k'
set xrange[0:11]
p for[i=0:(L/2)] "nwk.dat" u ($0*2*pi/(nt*dt)):(log(column(i+1))+i) w l t 'k='.i
set xrange[*:*]

unset multiplot
