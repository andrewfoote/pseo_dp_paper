/*************************
    This code creates the CDF
    of 1/(1+|x|^4) (gamma distribution for SS)
   *************************/


%macro create_cdf(outlib) ;

%let bottom = -1000 ;
%let top = 1000 ;
%let step = 0.001;

%let norm = pi/sqrt(2) ;

data &outlib..gamma_cdf  (keep=eta cdf cdf_lag bincounter ) ;

    format eta BEST12.5;
    format cdf cdf_lag pdf pdf_lag BEST30.25  ;
    retain cdf cdf_lag pdf pdf_lag counter norm;
    norm = sqrt(2)/constant("pi") ;
    counter = 1 ; 
    cdf = 0 ;
    cdf_lag = 0 ;
    pdf_lag = 0 ;
    do eta = &bottom. to &top. by &step.  ;
        pdf = (1/(1+(abs(eta)**4)))*norm  ;
        counter + 1 ;
        bincounter = int(counter/20) + 1 ;        
        cdf = cdf_lag + ((pdf+pdf_lag)*&step./2) ;
        
        output;
        cdf_lag = cdf ;
        pdf_lag = pdf ;
     end;
run;
    


proc print data=&outlib..gamma_cdf (obs = 10) ;
run;

proc means data=&outlib..gamma_cdf N mean min max p50 ;
    var cdf cdf_lag bincounter  ;
run;

%mend create_cdf ;


