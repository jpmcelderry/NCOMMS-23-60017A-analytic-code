/* MD Dec 2020 - Test Hierarchical modeling on smaller DS */

/* creating a generalized version of macro
   Feb 2021 - added the ability to pass in fixed and random effects
    May 2021 - added the ability to pass in class and variable list
 */

/* just testing something */
%macro test_covars(UserFixedE=,UserRandomE=/*,numRandomEffects=,numFixedEffects=*/);


%LET none=0;
%if (&UserRandomE ^= %str()) %then %let numRandomEffects=%SYSFUNC(countw(&UserRandomE)); %else %let numRandomEffects=0;
%put The number of random effects=&numRandomEffects;


%if (&UserFixedE ^= %str()) %then %let numFixedEffects=%SYSFUNC(countw(&UserFixedE)); %else %let numFixedeEffects=0;
%put The number of fixed effects=&numFixedEffects;


/*create n variables of n fixed effects */
%DO n=1 %TO &numFixedEffects; 
    %let fixed&n = %SYSFUNC(STRIP(%SCAN(&UserFixedE,&n)));
    %put The variable fixed&n=&&fixed&n;
%END;

/*create n variables of n random effects */
%DO n=1 %TO &numRandomEffects; 
    %let random&n = %SYSFUNC(STRIP(%SCAN(&UserRandomE,&n)));
    %put The variable random&n=&&random&n;
%END;


%mend;

/*************
   Main macro: run_proc_glimmix_general
    this macro takes in: statistical parameters, 
    a list of fixed effects (UserFixedE),
     a list of random effects (UserRandomE)
     a list of class variables (classlist)
     a list of all variables (varlist)
     the name of a file to which esimtates will be written (outfile).  The file will be writen in excel format
     
**********/ 
title '***************   GENERAL MACRO: *********************';
%macro run_proc_glimmix_general(maxiter=500, method=rspl, technique=quanew, update=dbfgs,
                        data=fam_hist.swdk3_both_all,show_output=1,UserFixedE=,UserRandomE=,outfile=,varlist=,classlist=,cholesky=, scoring=,
                        inititer=, maxtime=, initglm=, linesearch=, lsprecision=,
                        restart=);


%if (&linesearch ^= %str()) %then %do;
  %if (&technique ^= congra) & (&technique ^= quanew) & (&technique ^= newrap) %then
    %let linesearch = ;
%end; 
%if (&technique ^= congra) & (&technique ^= quanew) %then
  %let restart = ;


%put METHOD=&method  TECHNIQUE=&technique  UPDATE=&update;
%put CHOLESKY=&cholesky  SCORING=&scoring  INITITER=&inititer  MAXTIME=&maxtime;
%put INITGLM=&initglm  LINESEARCH=&linesearch  LSPRECISION=&lsprecision;
%put RESTART=&restart SHOW_OUTPUT=&show_output;
%put FixedEffects=&UserFixedE;

ods exclude all;

%if &show_output = 1 %then %do;
  ods exclude ClassLevels SolutionR ParameterEstimates;
  ods output classlevels=_class;
  ods output solutionr=randomEffects;
  ods output parameterestimates=fixedEffects;
  ods output estimates=Estimates; /* this statement is necessary to print out estimates of the model; */
%end;
 
%let none=0;

/* get the number of fixed and random effects */
%if (&UserFixedE ^= %str()) %then %let numFixedEffects=%SYSFUNC(countw(&UserFixedE)); %else %let numFixedeEffects=0;
%put The number of fixed effects=&numFixedEffects;

%if (&UserRandomE ^= %str()) %then %let numRandomEffects=%SYSFUNC(countw(&UserRandomE)); %else %let numRandomEffects=0;
%put The number of random effects=&numRandomEffects;


/*create n variables of n fixed effects	 This is to automate printing out the "estimates"  */

%DO n=1 %TO &numFixedEffects; 
    %let fixed&n = %SYSFUNC(STRIP(%SCAN(&UserFixedE,&n)));
    %put The variable fixed&n=&&fixed&n;
%END;
 

%DO j=1 %TO &numRandomEffects; 
    %let random&j = %SYSFUNC(STRIP(%SCAN(&UserRandomE,&j)));
    %put The variable random&j=&&random&j;
%END;
 

ods trace on;
proc glimmix data=&data
             method=&method 
			 maxopt=&maxiter 
			 &cholesky
			 &inititer
             &initglm
			 &scoring
			 ;

class %IF (&classlist ^=%str()) %THEN &classlist; %else REGION GENDER;
      ;
 /* full variable list passed to model below*/
 model CASE_STATUS (event='1') = %IF (&varlist ^= %str()) %THEN &varlist; %else indexage REGION GENDER REGISTRATIONYEAR INDEXYEAR; &UserFixedE  / solution dist=binary link=logit;
	
		%IF &numRandomEffects > &none %THEN %DO;
		  random &UserRandomE / solution; %END; 
	  
 
 nloptions maxfunc=&maxiter maxiter=&maxiter technique=&technique
           maxfunc=5000
           &maxtime
		   &linesearch
		   &lsprecision
		   &restart
           %if %str(&update) ^= %str() %then %do;
             update=&update
		   %end;
           ;


	 /* print out the expnentiated point estimate of the n fixed effects*/  
	   %do n=1 %to &numFixedEffects;
	       	   %put printing estimate of fixed effect &&fixed&n;
		   estimate "&&fixed&n" &&fixed&n 1  / exp cl;
          %end; 

	  /* print out the j exponentiated random effect estimates for the single fixed effects passed into the macro*/
	   %do j=1 %to &numRandomEffects;
	       	   %put printing estimate of random effect &&random&j;
		   estimate "&&random&j" &&fixed1 1 | &&random&j 1  / exp cl;
          %end; 
   run;
ods trace off;

/* writing out the exponentiated esimates to the file (outfile) */
ods excel file="&outfile..xlsx";
proc print data=Estimates;
run;
ods excel close;
   
*%exit: ;

%mend;
