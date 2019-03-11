/**********************************
    This macro takes the following inputs:
    - dataset
    - variable to calculate things from
    - list of percentiles
    - value of epsilon (privacy loss parameter) for percentiles
    - value of epsilon for N
    - SEED
    - output dataset
    

    First, it calculates the true values for the percentiles
    Second, it calculates the smooth sensitivity for each of these
    percentiles

    Third, it draws from the eligible distribution
    Finally,  it adds TRUE+ (SS*eta) for each value of percentile
    Outputs that dataset

 ************************************/

%macro dp_iteration(indset,calcvar,ptiles,mineps,maxeps,epsstep,outdset , noprint =YES) ;

%create_cdf(WORK);

%let epsilon_cell = 1.5 ; 

%if "&noprint."="YES" %then %do ;
    %put "ENTERING DIFF PRIVACY MACRO" ;
    %printjunk;
%end;

proc summary data=&indset. ;
    var &calcvar. ;
    output out= ptile_stats
           %do ii = 1 %to %sysfunc(countw(&ptiles.,' ')) ;
                p%scan(&ptiles,&ii.)(&calcvar.)=p%scan(&ptiles.,&ii.)_value  
           %end;
               N(&calcvar.)= cellcount ;
run;

%if %eval("&testing."="yes") %then %do ;

proc print data=&indset (obs=1) ;
    title 'crossvars';
run;
proc print data=ptile_stats ;
    title 'true statistics' ;
run;

%end;

/******************************************/
/*  Step 2: calculating smooth sensitivity - for all values of epsilon*/
/******************************************/
    


%do ii = 1 %to %sysfunc(countw(&ptiles.,' ')) ;
    %compute_ss_iteration(&indset.,&calcvar.,%scan(&ptiles,&ii.),&mineps.,&maxeps.,&epsstep.) ;

    proc print data=smoothsens;
        title 'smoothsensitivity' ;
    run;

    data smoothsens_%scan(&ptiles,&ii.) (drop = ss epsilon epsval ) ;
        set smoothsens (rename=(col1 = ss col2=epsilon)) end=eof ;

        array smoothsens{&mineps.:&maxeps.} ss_%scan(&ptiles.,&ii.)_&mineps.-ss_%scan(&ptiles.,&ii.)_&maxeps. ;
        retain ss_%scan(&ptiles.,&ii.)_&mineps.-ss_%scan(&ptiles.,&ii.)_&maxeps. ;
        
        epsval = &mineps.+ (_N_*&epsstep.) ;
        if ss = 0 then ss = . ; 
        smoothsens{epsval} = ss ;
        if eof then output ;
    run;
        

        
%end; /* end of loop through percentiles */

    data ss_all;
        merge smoothsens_25 smoothsens_50 smoothsens_75 ;
    run;

    proc print data=ss_all;
           title 'smoothsensitivity';
    run;
 

/*************************************/    
/*  Step 3: add these all together   */
/*************************************/


data _NULL_ ;
    set WORK.gamma_cdf end = eof ;
    if eof then call symput('maxbin',bincounter) ;
run;


data errors_bydraw (keep=error_p25 error_p50 error_p75 eps p25_value p50_value p75_value draw cellcount cellcount_protected cellcount_iter);
    if _n_ eq 0 then set WORK.gamma_cdf ;
    if _n_ eq 1 then do ;
        declare hash h_cdf(dataset: "WORK.gamma_cdf", ordered: 'Y') ;
        rc = h_cdf.defineKey('bincounter', 'eta') ;
        rc = h_cdf.defineData('eta','cdf','cdf_lag','bincounter') ;
        rc = h_cdf.defineDone() ;
        declare hiter i_dist('h_cdf') ;
    end;

    
    merge ptile_stats ss_all  ;


    array ss_25{&mineps.:&maxeps.} ss_25_&mineps.-ss_25_&maxeps. ;
    array ss_50{&mineps.:&maxeps.} ss_50_&mineps.-ss_50_&maxeps. ;
    array ss_75{&mineps.:&maxeps.} ss_75_&mineps.-ss_75_&maxeps. ;
    eta_cellcount = ranexp(0)*(1/&epsilon_cell.)-ranexp(0)*(1/&epsilon_cell.) ;
    cellcount_protected = round(cellcount+ eta_cellcount,1)  ;
    if cellcount_protected < 0 then do ;
        cellcount_protected = 0 ;
        flag_negative = 0 ;
    end;
    /* calculating protected values */
        call streaminit(1) ;
    do eps = &mineps. to &maxeps. by &epsstep. ;
        p = 1-(1/exp(eps/100)) ;
        /* doing 20 draws of the gamma distribution to get values of eta */
        do draw = 1 to 30 ;
        
        %do ii = 1 %to %sysfunc(countw(&ptiles., ' ' )) ;
                %let ptile = %scan(&ptiles.,&ii.) ;
                r = ranuni(-1) ;
                bincounter =1 ;
                rc=h_cdf.find() ;
                binsteps = int(log2(&maxcell.)) ;
                do while (binsteps>0 ) ;
                    if cdf<r then bincounter = bincounter + int(2**binsteps) ;
                    else bincounter = bincounter - int(2**binsteps) ;
                    rc = h_cdf.find() ;
                    binsteps=binsteps-1 ; 
                 end;
                 rc = i_dist.setcur() ;
                 do until (cdf_lag le r le cdf) ;
                     if cdf<r then rc=i_dist.next() ;
                     else if cdf_lag > r then rc=i_dist.prev() ;               
                 end;
                
                 error_p&ptile. =  (ss_&ptile.{eps}/(eps/1600))*eta;
        %end;
    
        x = rand('geometric',p) ;
        y = rand('geometric',p) ;
        noise = x-y  ;
        
        cellcount_iter = cellcount+noise ;
        cellcount_iter = max(cellcount_iter,0) ;
        output;
        end;
    end;
    /* NOTE: Rounding values to reasonable values, given what they are;
    earnings rounded to nearest cent; cellcount rounded to nearest integer. */
run;

proc summary data=errors_bydraw nway;
    class eps;
    var error_p25 error_p50 error_p75 cellcount_protected
        p25_value p50_value p75_value ;
    output out=&outdset.  mean(error_p25 error_p50 error_p75 cellcount_protected cellcount_iter) =
                          mean_error_p25 mean_error_p50 mean_error_p75 mean_protectedcount
                          mean_count_iter
                          mean(p25_value p50_value p75_value cellcount) =
                                      p25_value p50_value p75_value cellcount;
run;

/**************************************
    one thing to add here -
    pull out eps = 6,
    merge to all rows as "pXX_value_eps6"
    (a "protected" but very close value)
    which I will use as my baseline for
    relative error calculations
****************************************/
/*
 data eps6 (keep=p25_value_eps6 p50_value_eps6 p75_value_eps6);
     set &outdset. (where=(eps=600)) ;
     p25_value_eps6 = p25_value + mean_error_p25  ;
     p50_value_eps6 = p50_value + mean_error_p50  ;
     p75_value_eps6 = p75_value + mean_error_p75  ;
     err_cellcount_eps6 = mean_count_iter-cellcount ; 
     call symput('p25_6',p25_value_eps6) ;
     call symput('p50_6',p50_value_eps6) ;
     call symput('p75_6',p75_value_eps6) ;
     call symput('p25_6_err',mean_error_p25) ;
     call symput('p50_6_err',mean_error_p50) ;
     call symput('p75_6_err',mean_error_p75) ;
     call symput('err_cellcount_6',err_cellcount_eps6) ;
 run;

 data &outdset.;
     set &outdset. ;
     p25_value_eps6 = &p25_6. ;
     p50_value_eps6 = &p50_6. ;
     p75_value_eps6 = &p75_6. ;
     p25_err_eps6 = &p25_6_err. ;
     p50_err_eps6 = &p50_6_err. ;
     p75_err_eps6 = &p75_6_err. ;
     err_cellcount_eps6 = &err_cellcount_6. ; 
 run;
*/
 /* merging that into everything */

%if "&noprint."="YES" %then %do ;
%printlog ;
%put "LEAVING OPTIMAL SPECTRAL" ;
        %put "Best clusters are: &optclus." ;
%end;         

%mend dp_iteration ;    
