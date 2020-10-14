set terminal postscript eps enhanced color "GothicBBB-Medium-UniJIS-UTF8-H,20"
set output 'fig2.eps'
set samples 10000
set xrange [-pi:pi]
set key bottom right
set size square

series(x,n)=(n>1 ? series(x, n-1)+func(x,n) : func(x,1)+pi/4)
func(x,n)=((-1)**(n-1))*(sin(n*x))/n - (n % 2)*(2*cos((n)*x))/((n)**2*pi)
f(x)=(x>0 ? x : 0)

plot f(x), series(x,1) title "n=1", series(x,5) title "n=5", series(x,9) title "n=9"