set terminal postscript eps enhanced color "GothicBBB-Medium-UniJIS-UTF8-H,20"
set output 'fig1.eps'
#set terminal wxt
set size square
set ylabel "e^{-xh(τ)}"
set xlabel "τ"
set samples 10000
set xrange [0:2]

g(x,t)=exp(-t*h(x))
h(x)=x-log(x)

plot g(x,100) title "x=100"