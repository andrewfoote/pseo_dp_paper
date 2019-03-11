%include "config.sas" ;


proc export data=OUTDATA.mean_errors_protection
    outfile="&outdat./mean_errors_protection.dta"  replace;
run;

proc export data=OUTDATA.errors
    outfile="&outdat./simulation_errors.dta" replace;
run;


data eps_ss_errors eps_bin_ss_errors;
    set OUTDATA.mean_errors_ss ;

    if _TYPE_ = 2 then output eps_ss_errors;
    if _TYPE_ = 3 then output eps_bin_ss_errors;
run;

proc export data=eps_ss_errors
    outfile="&outdat./sim_eps_ss_errors.dta" replace;
run;

proc export data=eps_bin_ss_errors
    outfile="&outdat./ss_bin_eps_ss_errors.dta" replace;
run;
