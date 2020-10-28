set terminal postscript eps enhanced color "GothicBBB-Medium-UniJIS-UTF8-H,20"
set output 'fig2.eps'
set samples 10000
set size square
set xrange [-20:20]
set zeroaxis dt (5,5)
set xtics axis

A=1.0
T=1.0/2.0
d=1.0/20.0
o=T/(2.0*pi)

p=d*A/T
set yrange [-0.1:p+0.1]
set label left at 1,p "dA/T_0"
set label left at 1,p+0.07 "|C_n|"
set label right at 22.,0.01 "ω"

f(x)=abs((d*A/T)*(sin((x*pi*d)/T))/((x*pi*d)/T))

plot f(x)