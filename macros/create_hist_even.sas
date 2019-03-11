/* macro creates an evenly spaced histogram */

%macro create_hist_even(numbins,celldat) ;

%let min = 10000 ;
%let max = 300000;
%let totalmax = 500000 ;
    
%let binwidth = (&max-&min)/&numbins;

data histogram_even (keep=bincounter lower upper) ;
    bincounter = 0 ;
    do ii = &min to &max by &binwidth ;
        bincounter + 1 ;
        lower = ii  ;
        upper = ii + &binwidth ;
        output ;
    end;

    /* Making a final bin to catch all extra observations -
        bin is &max and above */

    lower = upper ;
    upper = &totalmax. ;
    bincounter = bincounter + 1 ;
    output;
    lower = &totalmax.;
    upper = &totalmax.+1;
    bincounter = bincounter + 1 ;
    output ;
    call symput('maxbin',bincounter) ;
run;

data _NULL_ ;
    set histogram_even end=eof;
    if eof then do ;
        call symput('maxbin_even',bincounter);
    end;
run;

data shell_even (keep= cellid  bincounter) ;
    set &celldat.;

    do ii = 1 to %eval(&numbins.+2);
        bincounter = ii ;
        output;
    end;
run;

proc sort data=shell_even;
    by bincounter;
run;

data shell_hist_even (keep=cellid bincounter lower upper) ;
    merge shell_even histogram_even ;
    by bincounter ;
run;

proc sort data=shell_hist_even;
     by cellid bincounter ;
run;



%mend create_hist_even;
