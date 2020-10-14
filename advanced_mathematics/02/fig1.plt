set terminal postscript eps enhanced color "GothicBBB-Medium-UniJIS-UTF8-H,20"
set output 'fig1.eps'
set samples 10000
set xrange [-pi:pi]
set key bottom right
set size square

series(x, n)=(n>1 ? series(x, n-1)+func(x,n) : func(x,1))
func(x, n)=((-1)**(n-1))*(2*sin(n*x))/n
f(x)=x

plot f(x) title "f(x)",series(x, 1) title "n=1", series(x, 5) title "n=5", series(x, 9) title "n=9"