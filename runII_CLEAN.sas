/* M D'Arcy 
   This file will use the un_proc_glimmix_general.sas
   The file containing the relevant data for this example is: .txt 

   For any "statistical" questions (e.g. - why are those parameters set the way they are) ask RP

   Jan 2021 - added in the ability to write estimates to a file.
 */


/* Below is a more generalized version of proc glimmix */
%include "C:\DIRECTORY/run_proc_glimmix_general_CLEAN.sas";
libname t "C:\DIRECTORY";

*options MEMSIZE=MAX;
options ps=32000 source source2 nomprint;



/* read in Inflammations and infections data formatted for these analyses in CleanHierarchy.sas */

/* make sure that text file columns look correct */
 
   proc import datafile='C:\DIRECTORY/HierarchyInfectInflam_SHORT.csv'
     out= t.testdata 
     dbms=csv 
     replace;
     getnames=yes ;
	 guessingrows = 15900; /*IMPORTANT*/
  run;   


data tmp; set t.testdata; myid = patid;  run;
 

title '******************  Running glimmix general passing in fixed and random effects';
/* For information about the main statistical parameters, talk to RM
   The user only needs to pass in the data, UserFixedE and UserRandomE variables.
  The random and fixed effects are all columns in the dataset
  The outfile is the name of the file (minus the .xlsx) you want to write estimates to.  It currently must be included in the macro call
  */
 

title '********************************* Infections and Inflammation hierarchical analyses';

/*NOT BMI adjusted */
/* running the macro for 1-10 years primary analysis; the same code was used for the 10-32 years prior to index date */

option spool; 
 %run_proc_glimmix_general(maxiter=400, method=rspl, technique=newrap, update=,
                  cholesky=%str(), initglm=%str(), scoring=%str(),
                  maxtime=%str(maxtime=9999), inititer=%str(),
                  linesearch=%str(), lsprecision=%str(),
                  data=tmp,show_output=1,UserFixedE=InfectInflam_lt10,
		  UserRandomE=InfectInflam_Pneum_lt10 InfectInflam_Heart_lt10 InfectInflam_HN_lt10 InfectInflam_Menin_lt10 InfectInflam_Hep_lt10 InfectInflam_TB_lt10 InfectInflam_Herpes_lt10 InfectInflam_UT_lt10 InfectInflam_GI_Diver_lt10 InfectInflam_GI_Haem_lt10  InfectInflam_GI_Dys_lt10 InfectInflam_GI_Ulc_lt10 InfectInflam_GI_Other_lt10 InfectInflam_GI_Ose_lt10 InfectInflam_GI_Inf_lt10 InfectInflam_GI_Chol_lt10 InfectInflam_GI_IChol_lt10 InfectInflam_GI_GERD_lt10 InfectInflam_GI_HP_lt10 InfectInflam_GI_GastrNIIGC_lt10 InfectInflam_GI_Aden_lt10 InfectInflam_GI_Calc_lt10 InfectInflam_GI_NoCAT_lt10 InfectInflam_OBGYN_lt10 InfectInflam_Limb_lt10 InfectInflam_Skin_lt10 InfectInflam_Bone_lt10 InfectInflam_Septi_lt10 InfectInflam_Enceph_lt10 InfectInflam_Flu_lt10 InfectInflam_Malaria_lt10 InfectInflam_Shingles_lt10 InfectInflam_Dental_lt10 InfectInflam_Resp_lt10 InfectInflam_UResp_lt10 InfectInflam_Arthrit_lt10 InfectInflam_Gout_lt10 InfectInflam_AI_IBD_lt10 InfectInflam_AI_RA_lt10 InfectInflam_AI_Psor_lt10 InfectInflam_AI_LThy_lt10 InfectInflam_AI_DMT1_lt10 InfectInflam_AI_Lupus_lt10 InfectInflam_AI_IDP_lt10 InfectInflam_AI_oth2_lt10 InfectInflam_AI_Alo_lt10 InfectInflam_AI_PolyR_lt10 InfectInflam_AI_Rose_lt10 InfectInflam_AI_NoCAT_lt10  InfectInflam_COPD_lt10  InfectInflam_Other_lt10, outfile=InfectInflamFullFeb2022, varlist=indexage REGION GENDER REGISTRATIONYEAR INDEXYEAR,classlist=REGION GENDER);

 
