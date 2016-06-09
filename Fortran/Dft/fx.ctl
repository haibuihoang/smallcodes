dset ^fx.dat
undef -999000000.0
title  OUTPUT FROM WRF V3.4 MODEL
xdef  128 linear 0 0.0001
ydef  1 linear 0 0.0001
zdef  1 levels  0.00000
tdef  1 linear 00Z10APR2013    1440MN      
VARS   5
fx    1  0  orginal function
fr    1  0  Re Fourier comp.
fi    1  0  Im Fourier comp.
fx1    1  0  restore function
fx2    1  0  restore function based on first 10 coef
ENDVARS
