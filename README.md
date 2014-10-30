CUBISM-MPCF
===========
Please cite:

1- Rossinelli D., Hejazialhosseini B. , Hadjidoukas P., Bekas C., Curioni A., Bertsch A., Futral S., Schmidt S.J., Adams N.A. and Koumoutsakos P., 11 PFLOP/s simulations of cloud cavitation collapse, Proceedings of the International Conference on High Performance Computing, Networking, Storage and Analysis (SC '13), 2013.

2- Hejazialhosseini B., Rossinelli D., Conti C., Koumoutsakos P., High throughput software for direct numerical simulations of compressible two-phase flows. In Proceedings of the International Conference on High Performance Computing, Networking, Storage and Analysis (SC '12). IEEE Computer Society Press, Los Alamitos, CA, USA, Article 16, 2012.


==========================
# AIM week 

# 1 Initial Conditions parameters :
  * Box size 
  * Initial T = 300 K
  * Initial Pressure (p) = 1 bar
  * T_s_m
  * T_s_e
  * I - current = 1250000 A
  * U = 240 V
  * V_arc + arc box
  * gamma = 1.4
  * flap area
    * rho_flap
    * s_flap
    * L_flap
  * zeta_grid = 0.219
  * threshold_pressure = 6000 Pa
  

# 1.2  Parameters to executable
```
app -sim flap -tend 2 -dumpperiod 10 -saveperiod 10 -cfl 0.1 -bpdx 100 -bpdz 1 -pOutside 100000 -pCrit 6000 -tInit 300 -tSM 0.015 -tSE 0.015 -iI 1250000 -iU 240 -arcX 0.1 -arcWidth 0.1 -zetaGrid 0.219 -threshP 6000  -flRho 3500 -flS 0.022 -flL 0.850 -g1 1.4 -verb 0 -nsteps 100
```
  
