#HyperVolume
reset
set term post eps m enhanced
set output 'Mut_1d.eps'
set size 0.5,0.5
set pointsize 0.2 
set grid
set title ''
set xr [0:1]
#set yr [0:4]
set xlabel 'X'
set ylabel 'P'
#set parametric
#et samples 10000
sigma=0.1
sigma2=sigma*(1-0.9)**0.2
mu=0.4
f(x)=1/(sigma*sqrt(2*3.14))*exp(-(x-mu)**2/(2*sigma**2))
f2(x)=1/(sigma2*sqrt(2*3.14))*exp(-(x-mu)**2/(2*sigma2**2))
p f(x) t'g/gmax=.0',f2(x) t'g/gmax=.9' 

#HyperVolume
reset
set isosample 500,500
set term post eps color
set output 'Mut_2d.eps'
set pm3d map
set size 0.5,0.5
set pointsize 0.2 
set grid
set title ''
set xr [0:1]
set yr [0:1]
set ylabel 'P'
set dgrid3d 500,500,40
set cntrparam levels 20 
#set parametric
#et samples 10000
sigma  = 0.12
sigma2 = 0.033
mu = 0.4
mu2= 0.9
f(x)=1/(sigma*sqrt(2*3.14))*exp(-(x-mu)**2/(2*sigma**2))
fy(x)=1/(sigma2*sqrt(2*3.14))*exp(-(x-mu2)**2/(2*sigma2**2))

sp f(x)*fy(y) with pm3d t'' 
