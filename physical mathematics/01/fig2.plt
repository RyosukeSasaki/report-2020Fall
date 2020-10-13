set terminal postscript eps enhanced color "GothicBBB-Medium-UniJIS-UTF8-H,20"
set output 'fig2.eps'
#set terminal wxt
set size square
set samples 10000
set xrange [-20:20]
set yrange [0:1.2]

set label 1 at -13,0.04 "-√x"
set arrow from  -sqrt(100),0.0 to -sqrt(100),1.2 nohead lw 1 dt (10,5) lc rgb "red"
g(x,t)=exp(-(x**2)/2+(x**3)/(3*sqrt(t))-(x**4)/(4*t))
t=100
plot g(x, t) title "x=100"