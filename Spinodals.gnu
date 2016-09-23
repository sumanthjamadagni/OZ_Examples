eps_crit_50_6 = 0.45 
eps_crit_50_8 = 0.70
eps_crit_50_10 = 0.85 
eps_crit_50_12 = 1.05 
eps_crit_50_18 = 1.65 

eps_crit_100_6 =  0.45
eps_crit_100_12 =  0.85
eps_crit_100_24 =  1.5

set key outside

set term postscript enhanced solid portrait 20 

set size 1.2, 0.60
set output 'T-rho.eps'

set xlabel '{/Symbol r}' 
set ylabel 'T'

plot [0:0.78] \
'm50-n6/LiqSpinodal.dat' u 4:1 w p lc 1 pt 6 ps 1 t '50/6',\
'm50-n6/VapSpinodal.dat' u 4:1 w p lc 1 pt 5 ps 1 noti ,\
'm50-n8/LiqSpinodal.dat' u 4:1 w p lc 2 pt 6 ps 1 t '50/8',\
'm50-n8/VapSpinodal.dat' u 4:1 w p lc 2 pt 5 ps 1 noti ,\
'm50-n10/LiqSpinodal.dat' u 4:1 w p lc 3 pt 6 ps 1 t '50/10',\
'm50-n10/VapSpinodal.dat' u 4:1 w p lc 3 pt 5 ps 1 noti ,\
'm50-n12/LiqSpinodal.dat' u 4:1 w p lc 4 pt 6 ps 1 t '50/12',\
'm50-n12/VapSpinodal.dat' u 4:1 w p lc 4 pt 5 ps 1 noti ,\
'm50-n18/LiqSpinodal.dat' u 4:1 w p lc 5 pt 6 ps 1 t '50/18',\
'm50-n18/VapSpinodal.dat' u 4:1 w p lc 5 pt 5 ps 1 noti 


set output 'Tr-rho.eps'
set xlabel '{/Symbol r}' 
set ylabel 'T_r = T/T@_c^{est}' 

plot [0:0.78] [:1.1]\
'm50-n6/LiqSpinodal.dat' u 4:($1 * eps_crit_50_6) w p lc 1 pt 6 ps 1 t '50/6',\
'm50-n6/VapSpinodal.dat' u 4:($1 * eps_crit_50_6) w p lc 1 pt 5 ps 1 noti ,\
'm50-n8/LiqSpinodal.dat' u 4:($1 * eps_crit_50_8) w p lc 2 pt 6 ps 1 t '50/8',\
'm50-n8/VapSpinodal.dat' u 4:($1 * eps_crit_50_8) w p lc 2 pt 5 ps 1 noti ,\
'm50-n10/LiqSpinodal.dat' u 4:($1 * eps_crit_50_10) w p lc 3 pt 6 ps 1 t '50/10',\
'm50-n10/VapSpinodal.dat' u 4:($1 * eps_crit_50_10) w p lc 3 pt 5 ps 1 noti ,\
'm50-n12/LiqSpinodal.dat' u 4:($1 * eps_crit_50_12) w p lc 4 pt 6 ps 1 t '50/12',\
'm50-n12/VapSpinodal.dat' u 4:($1 * eps_crit_50_12) w p lc 4 pt 5 ps 1 noti ,\
'm50-n18/LiqSpinodal.dat' u 4:($1 * eps_crit_50_18) w p lc 5 pt 6 ps 1 t '50/18',\
'm50-n18/VapSpinodal.dat' u 4:($1 * eps_crit_50_18) w p lc 5 pt 5 ps 1 noti 

#---------
set output 'eps-rho.eps'

set ylabel '{/Symbol e}' 
set xlabel '{/Symbol r}' 

plot [0:0.78] \
'm50-n6/LiqSpinodal.dat' u 4:2 w p lc 1 pt 6 ps 1 t '50/6',\
'm50-n6/VapSpinodal.dat' u 4:2 w p lc 1 pt 5 ps 1 noti ,\
'm50-n8/LiqSpinodal.dat' u 4:2 w p lc 2 pt 6 ps 1 t '50/8',\
'm50-n8/VapSpinodal.dat' u 4:2 w p lc 2 pt 5 ps 1 noti ,\
'm50-n10/LiqSpinodal.dat' u 4:2 w p lc 3 pt 6 ps 1 t '50/10',\
'm50-n10/VapSpinodal.dat' u 4:2 w p lc 3 pt 5 ps 1 noti ,\
'm50-n12/LiqSpinodal.dat' u 4:2 w p lc 4 pt 6 ps 1 t '50/12',\
'm50-n12/VapSpinodal.dat' u 4:2 w p lc 4 pt 5 ps 1 noti ,\
'm50-n18/LiqSpinodal.dat' u 4:2 w p lc 5 pt 6 ps 1 t '50/18',\
'm50-n18/VapSpinodal.dat' u 4:2 w p lc 5 pt 5 ps 1 noti 

!display T-rho.eps & 
!display Tr-rho.eps & 
!display eps-rho.eps & 

