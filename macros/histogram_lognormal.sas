/******************************
    This macro takes a few inputs:

    1. number of bins
    2. input dataset
    

    And then puts things in histogram bins,
    adds noise, and then calculates the percentiles

    The output dataset (pctiles_lognorm) is cellid with the three calculated
    percentiles and a protected cell count
*****************************/

%macro histogram_lognormal(numbins,eps,indata) ;


proc sort data=&indata out=earnings ;
    by cellid;
run;

proc sort data=earnings out=temp2  nodupkey;
    by cellid ;
run;

/******************************
    create empirical histogram
******************************/
 data histogram_constructed ;

    /*creating hashtable lookup */
    if _n_ eq 0 then set histogram_lognorm ;
    if _n_ eq 1 then do ;
        declare hash h_hist(dataset: "histogram_lognorm",ordered:'Y');
        rc = h_hist.defineKey('bincounter') ;
        rc = h_hist.defineData('bincounter','lower','upper');
        rc = h_hist.defineDone() ;
        declare hiter i_hist('h_hist') ;
     end;
     set earnings ;

     /* creating a binary search tree */
     
         
     bincounter = 1 ;
     rc = h_hist.find();
     binsteps = int(log2(&maxbin_lognorm.)) ; /* max steps to get close */
        
     do while (binsteps>0) ;        
         if ann_earn > upper then bincounter = bincounter + int(2**binsteps) ;
         else bincounter = bincounter - int(2**binsteps) ;
         if bincounter > &maxbin_lognorm. then bincounter = &maxbin_lognorm. ;
         else if bincounter < 1 then bincounter = 1 ;
         rc = h_hist.find() ;
         binsteps = binsteps-1 ; 
     end;
     rc = i_hist.setcur() ;
     do until (lower le ann_earn < upper) ;
        if upper le ann_earn then rc=i_hist.next() ;
        else if lower > ann_earn then rc=i_hist.prev() ;
        if bincounter = &maxbin_lognorm. and ann_earn > upper then ann_earn = upper-1 ; 
     end;
     output ;
run;   


proc summary data=histogram_constructed nway ;
    class cellid bincounter;
    var ann_earn ;
    output out=hist_counts N(ann_earn) = count;
run;

data hist_counts2 (drop= x y );
    merge shell_hist_lognorm (in=a) hist_counts (in=b) ;
    by cellid bincounter ;

    if a and not b then count = 0  ;

    /* Add noise */
    p = 1- (1/exp(&eps.)) ;
    x = rand('geometric',p);
    y = rand('geometric',p) ;
    noisy_count = count+(x-y) ;

run;

proc summary data=hist_counts2;
    class cellid ;
    var noisy_count ;
    output out=totals sum(noisy_count) = total_count_noisy ;
run;

data hist_noisy  ;
    merge hist_counts2 totals ;
    by cellid ;
run;

proc sort data=hist_noisy ;
    by cellid bincounter ;
run;


/******************************************
    Now extract the percentiles
******************************************/
/*
data ptiles_fuzz (keep=cellid total_count_noisy p25_fuzz p50_fuzz p75_fuzz) ;
    array bhist{&maxbin_lognorm.} _TEMPORARY_ ;
    array thist{&maxbin_lognorm.} _TEMPORARY_ ;
    array ebc{&maxbin_lognorm.} _TEMPORARY_ ;

    if _N_ = 1 then do;
        set histogram_lognorm (keep=bincounter  bottom top);
        do i = 1 to &maxbin_lognorm. ;
            bhist{bincounter}=bottom;
            thist{bincounter}=top ;
        end;
    end;    
*/
%do ii = 1 %to %sysfunc(countw(&pctiles.,' ')) ;
    %let pctile = %scan(&pctiles.,&ii.);
    
    data temp_&pctile. (drop = flag lag_count total_count );

        set hist_noisy ;
        by cellid;

        retain total_count lag_count flag pctile_&pctile._protected ;

        if first.cellid then do ;
            total_count = 0 ;
            flag = 0 ;
            pctile_&pctile._protected = . ;
            lag_count = 0 ;
            if _N_ < 21 then put flag= ;
        end;
        if total_count_noisy > 0  then do ;           
            lag_count = total_count ;
            total_count = total_count + noisy_count ;
            /*  if current bin includes the pctile */
                if (total_count/total_count_noisy >=  &pctile./100 and flag = 0) then do ;
                    flag = 1 ;
                    x = (&pctile./100)*total_count_noisy - lag_count ;
            /* here, we assume that the earnings are uniformly distributed inside a bin */
                    pctile_&pctile._protected = lower+(x/noisy_count)*(upper-lower)  ;
                end; 
        end;                         
        if last.cellid then output ; 

    run;
%end; 

data pctiles_lognorm
       (keep=cellid p25_lognorm p50_lognorm p75_lognorm cellcount_lognorm) ;
    merge
        %do ii=1 %to %sysfunc(countw(&pctiles.,' ')) ;
           %let pctile = %scan(&pctiles.,&ii.) ;
           temp_&pctile.
        %end;
        ;
        ;
     by cellid ;     
     cellcount_lognorm = total_count_noisy  ;
     

     p25_lognorm = round(pctile_25_protected,0.01) ;
     p50_lognorm = round(pctile_50_protected,0.01) ;
     p75_lognorm = round(pctile_75_protected,0.01) ;

run;

proc sort data=pctiles_lognorm;
    by cellid;
run;


%mend histogram_lognormal ;

    
