/***************************************
This macro computes the smooth sensitivity
    for a given percentile

    It requires four inputs:
    - dataset from which to compute SS
    - X, the variable for which SS is computed
    - percentile, the percentile of the SS
              NOTE: percentile is expressed as a whole number
    - epsilon, the privacy loss parameter
******************************************/

%macro compute_ss_iteration(dset,X,pctile,mineps,maxeps,epsstep) ;

*begin proc iml ;
proc iml ; 

*read in dataset ;
use &dset.  ;
read all var {&X.} into Y;

*variable of interest that we want to calculate SS over ;



*here I sort in ascending order ; 
call sort(Y,1) ;
X = t(Y);

*number of elements in X ;
len = NCOL(X) ;

*percentile position;
* does this need to be floored? ;
ptile = len*&pctile./100 ;

*ptile_comp is the opposite pctile in dist;
ptile_comp=len*(100-&pctile.)/100 ;



/*
n is the number of observations between the given percentile
the end of the distribution (top or bottom)

To that end, we need to do a min() function to properly assign n
(since we need the minimum gap between pctile and end of list)

in the case where pctile = 50, n is just what it was in Yan's code (len-med-2)
*/
n1 = len-ptile-2 ;
n2 = len-ptile_comp-2 ;
n3 = 50 ; 

n = min(n1,n2,n3) ; /* only looping over data sets differing by 50 - it becomes irrelevant
    and small after that */


/************* CHECKING CONDITION - ENOUGH OBSERVATIONS? **********/
       col = (&maxeps.-&mineps.)/&epsstep. ;
    SS = repeat(.,col,2) ;  
if n > 1 then do ;
  
*outputs of SS procedure ;
A = repeat(0,1,n+1) ;


e = constant("e") ;

do jj = 1 to  col ;
    eps = &mineps. + (jj*&epsstep.) ;
    beta= (eps/100)/4; 
do k = 1 to n ;
	cur = 0 ;
      
	do t = 1 to k+2 ;
           
                x1 = X[1,ptile+t];
                x2 = X[1,ptile+t-k-1] ;
            
            /****** THERE ARE A BUNCH OF PROBLEMS HERE *****/
                    if x1-x2 > cur then  cur=x1-x2 ;
            
	end;
        *putting the max local scaled sensitivity in matrix A ;
        int = k*beta; 
	A[1,k] = cur/(e**(int));

        /*
        In cases where int is VERY
        small, SAS doesn't have the
        machine precision , and considers
        it missing.
        
        This step is a fix for the problem.        
        */
        if A[1,k]=. then A[1,k]=0 ; 
end; 
reg = max(A) ; 
SS[jj,1] = reg ;
SS[jj,2] = eps ;
end;



end;  

/*outputting SS to a dataset, titled "smoothsens" */ 
    create smoothsens from SS;
%put "HERE CREATE" ;
append from SS ;
%put "HERE APPEND" ;
close smoothsens;

quit;

%mend ;
