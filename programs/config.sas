libname INDATA "[esim data]" ;
libname OUTDATA "[simulation folder]" ;

%let outdat = [outdat] ;

/*log normal parameters */
%let mean = 11.00255 ;
%let sd = 0.75275;

/* relevant percentiles */
%let pctiles= 25 50 75 ;


options fullstimer mlogic  sasautos=("../macros"); 
