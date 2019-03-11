/******************************
    This macro takes a few inputs:

    1. number of bins
    2. input dataset
    

    And then puts things in histogram bins,
    adds noise, and then calculates the percentiles

    The output dataset (pctiles_even) is cellid with the three calculated
    percentiles and a protected cell count
*****************************/

%macro histogram_evenspaced(numbins,eps,indata) ;


proc sort data=&indata out=temp ;
    by cellid;
run;


/******************************
    create empirical histogram
******************************/
 data histogram_constructed ;

    /*creating hashtable lookup */
    if _n_ eq 0 then set histogram_even ;
    if _n_ eq 1 then do ;
        declare hash h_hist(dataset: "histogram_even",ordered:'Y');
        rc = h_hist.defineKey('bincounter') ;
        rc = h_hist.defineData('bincounter','lower','upper');
        rc = h_hist.defineDone() ;
        declare hiter i_hist('h_hist') ;
     end;
     set earnings ;

     /* creating a binary search tree */
     
         
     bincounter = 1 ;
     rc = h_hist.find();
     binsteps = int(log2(&maxbin_even.+1)) ; /* max steps to get close */
        
     do while (binsteps>0) ;        
         if ann_earn > upper then bincounter = bincounter + int(2**binsteps) ;
         else bincounter = bincounter - int(2**binsteps) ;
         if bincounter > &maxbin_even. then bincounter = &maxbin_even. ;
         else if bincounter < 1 then bincounter = 1 ;
         rc = h_hist.find() ;
         binsteps = binsteps-1 ; 
     end;
     rc = i_hist.setcur() ;
     do until (lower le ann_earn < upper) ;
        if upper le ann_earn then rc=i_hist.next() ;
        else if lower > ann_earn then rc=i_hist.prev() ;
        if bincounter = &maxbin_even. and ann_earn > upper then ann_earn = upper-1 ; 
     end;
     output ;
run;   


proc summary data=histogram_constructed nway ;
    class cellid bincounter;
    var ann_earn ;
    output out=hist_counts N(ann_earn) = count;
run;

proc sort data=hist_counts ;
    by cellid bincounter;
run;

data hist_counts2 (drop= x y );
    merge shell_hist_even (in=a) hist_counts (in=b) ;
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
proc print data=hist_noisy (obs=20) ;
    title 'hist noisy' ;
run;

/******************************************
    Now extract the percentiles
******************************************/

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

proc print data=temp_25 (obs=20) ;
    title 'temp 25 - obs 20 ' ;
run;

data pctiles_even
       (keep=cellid p25_even p50_even p75_even cellcount_even) ;
    merge
        %do ii=1 %to %sysfunc(countw(&pctiles.,' ')) ;
           %let pctile = %scan(&pctiles.,&ii.) ;
           temp_&pctile.
        %end;
        ;
        ;
     by cellid ;     
     cellcount_even = total_count_noisy  ;
     

     p25_even = round(pctile_25_protected,0.01) ;
     p50_even = round(pctile_50_protected,0.01) ;
     p75_even = round(pctile_75_protected,0.01) ;

run;

proc sort data=pctiles_even ;
    by cellid ;
run;

proc print data=pctiles_even (obs=20);
    title 'pctiles even';
run;

%mend histogram_evenspaced ;

    
