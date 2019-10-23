/**************************
    This program iterates over three loops for two different types of
    histograms:

    log normal
    evenly spaced bins

    We iterate each combination of numbins and epsilon 50 times.

    Each time, we save the results, and store them for later.

    At end, we merge in the true percentiles
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


%let flag  = 0 ;
%do numbin = 10 %to 30 ;
    %create_hist_lognorm(&numbin,temp2) ;
    %create_hist_even(&numbin.,temp2) ;
    data _NULL_ ;
        set histogram_lognorm end=eof;
        if eof then do ;
            call symput('maxbin_lognorm',bincounter);
            end;
    run;

    data _NULL_ ;
        set histogram_even end=eof;
        if eof then do ;
            call symput('maxbin_even',bincounter);
        end;
    run;
    
    %do epsilon = 1 %to 30 ;
        %do iter = 1 %to 50 ;
            %let eps = %sysevalf(&epsilon./10) ;
            %let runtime = %sysfunc(datetime()) ;
            %put "******************************************************" ;
            %put "EPSILON: &eps. --- NUMBIN: &numbin. --- ITERATION: &iter. " ;
            %put "THIS LOOP STARTED: &runtime." ;
            %put "******************************************************" ;

            %printjunk;
            
            %histogram_lognormal(&numbin.,&eps,synth) ;

            %histogram_evenspaced(&numbin,&eps,synth) ;

            %printlog;
            data pctiles_&iter. ;
                 merge pctiles_lognorm pctiles_even;
                 by cellid ;
                 numbin = &numbin;
                 epsilon = &epsilon;
                 iteration = &iter;
            run;

            proc print data=pctiles_&iter. (obs=20) ;
                title "pctiles &iter." ;
            run;            

        %end; /* end of iteration loop */
        data pctiles_app ;
            set
                %do jj= 1 %to %eval(&iter.-1) ;
                 pctiles_&jj.
                %end;
                ;
        run;
            %if %eval(&flag ne 0) %then %do ;
                proc append base=pctiles_stored
                            data=pctiles_app;
                run;
            %end;
            %else %do ;
                data pctiles_stored ;
                    set pctiles_app;
                run;
            %end;

        %let flag = 1 ;
    %end;
%end;


proc sort data=pctiles_stored ;
    by cellid ;
run;

data OUTDATA.errors ;
    merge pctiles_stored OUTDATA.simulated_percentiles;
    by cellid ;

    /* L1 error */
    L1_p25_lognorm = p25_earn - p25_lognorm;
    L1_p50_lognorm = p50_earn - p50_lognorm;
    L1_p75_lognorm = p75_earn - p75_lognorm;

    L1_p25_even = p25_earn - p25_even ;
    L1_p50_even = p50_earn - p50_even;
    L1_p75_even = p75_earn - p75_even ;

    absL1_p25_lognorm =  abs(L1_p25_lognorm) ;
    absL1_p50_lognorm = abs(L1_p50_lognorm) ;
    absL1_p75_lognorm = abs(L1_p75_lognorm) ;

    absL1_p25_even = abs(L1_p25_even) ;
    absL1_p50_even = abs(L1_p50_even) ;
    absL1_p75_even = abs(L1_p75_even) ;

    /* Relative error */
    rel_p25_lognorm = 1- (absL1_p25_lognorm/p25_earn) ;
    rel_p50_lognorm = 1- (absL1_p50_lognorm/p50_earn) ;
    rel_p75_lognorm = 1- (absL1_p75_lognorm/p75_earn) ;

    rel_p25_even = 1- (absL1_p25_even/p25_earn) ;
    rel_p50_even = 1- (absL1_p50_even/p50_earn) ;
    rel_p75_even = 1- (absL1_p75_even/p75_earn) ;

    
    tot_error_lognorm = (L1_p25_lognorm/p25_earn) + (L1_p50_lognorm/p50_earn) + (L1_p75_lognorm/p75_earn);
    tot_error_even = (L1_p25_even/p25_earn) + (L1_p50_even/p50_earn) + (L1_p75_even/p75_earn) ;

    tot_abserror_lognorm = absL1_p25_lognorm + absL1_p50_lognorm + absL1_p75_lognorm ;
    tot_abserror_even = absL1_p25_even + absL1_p50_lognorm + absL1_p75_lognorm ;

    if cellcount_lognorm<30 then cellbin=1 ;
        else if 30 <= cellcount_lognorm < 50 then cellbin =2 ;
        else if 50 <= cellcount_lognorm < 80 then cellbin =3 ;
        else if 80 <= cellcount_lognorm < 100 then cellbin = 4 ;
        else if 100 <= cellcount_lognorm < 200 then cellbin = 5 ;
        else if 200 <= cellcount_lognorm < 300 then cellbin = 6 ;
        else if 300 <= cellcount_lognorm then cellbin = 7 ;

run;

/* adding ranks to the data */

proc sort data=OUTDATA.errors out=error_1 ;
    by iteration epsilon numbin ;
run;
    
proc rank data=error_1 out=error_rank ;
     by iteration epsilon numbin ;

     var p25_earn p50_earn p75_earn p25_lognorm p50_lognorm p75_lognorm p25_even p50_even p75_even ;
     ranks p25_truerank p50_truerank p75_truerank
           p25_logrank p50_logrank p75_logrank
           p25_evenrank p50_evenrank p75_evenrank;
run;

proc corr data=error_rank spearman outp=pearson_measures ;
    var p25_truerank p50_truerank p75_truerank
           p25_logrank p50_logrank p75_logrank
           p25_evenrank p50_evenrank p75_evenrank ;
    by iteration epsilon numbin ; 
run;

/* Going to measure the spearman correlation value here */

proc export data=pearson_measures outfile="&outdat./pearson_measures.dta" replace;
run;    


proc summary data=OUTDATA.errors;
    class iteration numbin epsilon cellbin ;
    var L1_: rel_: ;
    output out=OUTDATA.mean_errors_protection

        mean(L1_p25_lognorm
             L1_p50_lognorm
             L1_p75_lognorm
             L1_p25_even
             L1_p50_even
             L1_p75_even
             rel_p25_lognorm
             rel_p50_lognorm
             rel_p75_lognorm
             rel_p25_even
             rel_p50_even
             rel_p75_even

             tot_error_lognorm
             tot_error_even
             tot_abserror_lognorm
             tot_abserror_even

                          ) =

        meanL1_p25_lognorm
        meanL1_p50_lognorm
        meanL1_p75_lognorm
        meanL1_p25_even
        meanL1_p50_even
        meanL1_p75_even

        meanrel_p25_lognorm
        meanrel_p50_lognorm
        meanrel_p75_lognorm
        meanrel_p25_even
        meanrel_p50_even
        meanrel_p75_even

        meantot_err_lognorm
        meantot_err_even
        meantot_abserr_lognorm
        meantot_abserr_even

        ;
run;


%mend iteration;

%iteration;
