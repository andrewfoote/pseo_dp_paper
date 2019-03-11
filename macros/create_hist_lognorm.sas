%macro create_hist_lognorm(nbin,celldat) ;



data histogram_lognorm (keep=bincounter lower upper);
   retain bincounter step lower upper width cdf_check;
   lower=10000;
   step=1/&nbin;
   do bincounter=1 to %eval(&nbin-1);
      upper=round(exp(&mean + &sd*probit(step*bincounter)));
      cdf_check=cdf('NORMAL',log(upper),&mean,&sd);
      width=upper-lower;
      output;
      lower=upper;
   end;
   /* Can censor data at a reasonable large value and it 
      will not affect percentiles */
   upper=round(exp(&mean + &sd*probit(0.975))) ;
   width = upper-lower ; 
   output ;
   bincounter + 1;
   lower=upper ;
       
   upper=round(exp(&mean +  &sd*probit(0.999)));
   width=upper-lower;
   cdf_check=cdf('NORMAL',log(upper),&mean,&sd);
   output;
run;


data _NULL_ ;
    set histogram_lognorm end=eof;
    if eof then do ;
        call symput('maxbin_lognorm',bincounter);
    end;
run;

data shell_lognorm (keep= cellid  bincounter) ;
    set &celldat.;

    do ii = 1 to %eval(&nbin.+2);
        bincounter = ii ;
        output;
    end;
run;

proc sort data=shell_lognorm;
    by bincounter;
run;

data shell_hist_lognorm (keep=cellid bincounter lower upper) ;
    merge shell_lognorm histogram_lognorm ;
    by bincounter ;
run;

proc sort data=shell_hist_lognorm;
     by cellid bincounter ;
run;



%mend create_hist_lognorm; 
