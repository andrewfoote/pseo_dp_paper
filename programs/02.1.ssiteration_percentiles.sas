/**************************
This program computes the smooth
sensitivity for each cell, and then draws
20 draws for each cell from gamma
distribution
 ***************************/
%macro printlog ;


proc printto
    log="02.iteration_percentiles.log" 
    print="02.iteration_percentiles.lst" ;
run;


%mend printlog ;

%macro iteration ;
*%printlog;
    %include "config.sas" ;

data synth ;
    set INDATA.synth ;
run;

proc sort data=synth out=temp2 nodupkey ;
    by cellid;
run;


data _NULL_ ;
    set temp2 end=eof;
    if eof then call symput('maxcell',cellid) ;
run;
    
proc print data= synth (obs = 100 ) ;
    title 'temp2' ;
run;

%let testing = yes ;

/* looping over cellnumber */
%do cellnum = 1 %to  &maxcell. ;
    
    data temp3 ;
        set synth (where=(cellid = &cellnum.)) ;
    run;

    %let nobs = 0 ;
    data _NULL_ ;
        set temp3 end=eof;
        if eof then call symput('nobs',_N_) ;
    run;

    %if %eval(&nobs ne 0) %then %do ;

    %dp_iteration(temp3 , ann_earn, 25 50 75, 10, 300, 10,iteration_step,noprint=NO);


    data iteration_step ;
        set iteration_step;
        cellid = &cellnum. ;
    run;
    
    %if %eval(&cellnum.>1) %then %do ;
        proc append base=smoothsens_errors
                    data=iteration_step ;
        run;
    %end;
    %else %do ;
        data smoothsens_errors;
            set iteration_step;
        run;
    %end;

    %end;
%end;


proc sort data=smoothsens_errors ;
    by cellid eps ;
run;

proc print data=smoothsens_errors;
    title 'smoothsens errors' ;
run;

data errors ;
    set smoothsens_errors;

    abs_error_p25 = abs(mean_error_p25) ;
    abs_error_p50 = abs(mean_error_p50) ;
    abs_error_p75 = abs(mean_error_p75) ;

    rel_error_p25 = 1- (abs_error_p25/p25_value) ;
    rel_error_p50 = 1- (abs_error_p50/p50_value) ;
    rel_error_p75 = 1- (abs_error_p75/p75_value) ;

    if cellcount<30 then cellbin=1 ;
        else if 30 <= cellcount < 50 then cellbin =2 ;
        else if 50 <= cellcount < 80 then cellbin =3 ;
        else if 80 <= cellcount < 100 then cellbin = 4 ;
        else if 100 <= cellcount < 200 then cellbin = 5 ;
        else if 200 <= cellcount < 300 then cellbin = 6 ;
        else if 300 <= cellcount then cellbin = 7 ;
     epsilon = eps/100 ;
run;

proc print data=errors;
    title 'errors' ;
run;

/* summarizing all the errors */
proc summary data=errors;
    class  epsilon cellbin ;
    var abs_: rel_: ;
    output out=OUTDATA.mean_errors_ss

        mean(abs_error_p25
             abs_error_p50
             abs_error_p75
             rel_error_p25
             rel_error_p50
             rel_error_p75) =

        meanabs_error_p25_ss
        meanabs_error_p50_ss
        meanabs_error_p75_ss

        meanrel_error_p25_ss
        meanrel_error_p50_ss
        meanrel_error_p75_ss
        ;
run;


%mend iteration;

%iteration;
