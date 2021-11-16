/*******************
    This program extracts the true percentiles by cell ID
  ******************/

%macro percentiles;

%include "config.sas" ;

proc summary data=INDATA.synthetic_data nway ;
    class cellid ;
    var ann_earn ;
    outputs out=OUTDATA.simulated_percentiles
        p25(ann_earn) = p25_earn
        p50(ann_earn) = p50_earn
        p75(ann_earn) = p75_earn;
run;

%mend percentiles ;

%percentiles;
