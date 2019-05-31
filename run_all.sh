
module purge
module use /g/data3/hh5/public/modules
module load conda
module load parallel

parallel python run_test.py $* {} <<EOF
arm_sgp_summer_1997_A_2017_GFS
arm_sgp_summer_1997_A
arm_sgp_summer_1997_B
arm_sgp_summer_1997_C
arm_sgp_summer_1997_R
arm_sgp_summer_1997_S
arm_sgp_summer_1997_T
arm_sgp_summer_1997_U
arm_sgp_summer_1997_X
astex
bomex
default
LASSO_2016051812_MSDA
LASSO_2016051812
LASSO_2016051812_VARA
twpice_2017_GFS
twpice
EOF
