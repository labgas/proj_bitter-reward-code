* AUTHOR

Lukas Van Oudenhove @LaBGAS, KU Leuven

https://gbiomed.kuleuven.be/english/research/50000625/50000628/labgas

September 2023;


* RESEARCH QUESTION AND HYPOTHESES



* IMPORT AND LIBRARY DEF
------------------------;

%let path=\\gbw-s-labgas01.luna.kuleuven.be\data\proj_bitter-reward\BIDS\phenotype;
* NOTE: path to the location on the server through the Windows filesystem, SAS is only available for RedHat Linux systems, not Ubuntu (or at least so it seems);

proc import
	datafile="&path\food_images_wanting_liking.tsv"
	out=work.ratings_all
	dbms=tab
	replace;
run;

* GET RID OF NANS AND CONVERT LIKING VAR TO NUMERIC;
data work.ratings_all;
set work.ratings_all;
 if liking = 'NaN' then liking = .;
run;

data work.ratings_all;
set work.ratings_all;
 liking_numeric = input(liking,3.);
 drop liking;
run;

data work.ratings_all;
set work.ratings_all;
 liking = liking_numeric;
 drop liking_numeric;
run;

* CREATE COMBINED TRIAL_TYPE & CONDITION VARIABLE;
data work.ratings_all;
set work.ratings_all;
	if trial_type = 'high calorie' AND condition = 'bitter' then cond_trial_type = 'bit_high_calorie';
	if trial_type = 'low calorie' AND condition = 'bitter' then cond_trial_type = 'bit_low_calorie';
	if trial_type = 'neutral' AND condition = 'bitter' then cond_trial_type = 'bit_neutral';
	if trial_type = 'high calorie' AND condition = 'saline' then cond_trial_type = 'pla_high_calorie';
	if trial_type = 'low calorie' AND condition = 'saline' then cond_trial_type = 'pla_low_calorie';
	if trial_type = 'neutral' AND condition = 'saline' then cond_trial_type = 'pla_neutral';
run;


* PRIMARY OUTCOME
-----------------;

* CHECK DISTRIBUTION

* across subjects & conditions;
proc univariate data=work.ratings_all;
var rating;
histogram rating / kernel normal (mu=est sigma=est) lognormal (sigma=est theta=est zeta=est);
run;
* NOTE: distribution looks ok, even though not formally normal;

* by condition;
proc sort data=work.ratings_all;
by trial_type;
run;

proc univariate data=work.ratings_all;
by trial_type;
var rating;
histogram rating / kernel normal (mu=est sigma=est) lognormal (sigma=est theta=est zeta=est);
run;

proc univariate data=work.ratings_all;
class trial_type;
var rating;
histogram rating / kernel overlay;
run;

proc univariate data=work.ratings_all noprint;
by trial_type;
class participant_id;
var rating;
histogram rating / kernel overlay;
run;

* by subject;
proc sort data=work.ratings_all;
by participant_id;
run;

proc univariate data=work.ratings_all;
by participant_id;
var rating;
histogram rating / kernel normal (mu=est sigma=est) lognormal (sigma=est theta=est zeta=est);
run;

proc univariate data=work.ratings_all;
class participant_id;
var rating;
histogram rating / kernel overlay;
run;

proc univariate data=work.ratings_all noprint;
by participant_id;
class trial_type;
var rating;
histogram rating / kernel overlay;
run;

* summary multi-panel histogram;
proc univariate data=work.ratings_all noprint;
class trial_type participant_id;
histogram rating / ncols=4 nrows=30 kernel overlay;
run;
* too many tiles ;


* MIXED MODELS

* model including trial id and its interaction with trial type;
proc mixed data=work.ratings_all;
class trial_type participant_id;
model rating = trial_type | trial_id_cond_concat / solution residual influence ddfm=kenwardroger2;
random intercept trial_id_cond_concat / subject=participant_id type=un g gcorr;
lsmeans trial_type;
lsmestimate trial_type 'erythritol versus sucrose' 1 0 -1,
	'erythritol versus sucralose' 1 -1 / adjust=bon stepdown adjdfe=row;
run;

proc mixed data=work.ratings_all;
class trial_type participant_id;
model rating = trial_type | trial_id_cond_concat / solution residual influence ddfm=kenwardroger2;
random intercept trial_id_cond_concat / subject=participant_id group=trial_type type=un g gcorr;
lsmeans trial_type;
lsmestimate trial_type 'erythritol versus sucrose' 1 0 -1,
	'erythritol versus sucralose' 1 -1 / adjust=bon stepdown adjdfe=row;
run;
* NOTES: 
1. the second model, with a trial- rather than subject-level covariance structure, fits a lot better
2. the interaction effect is far from significant, so we omit it;

* model including trial id but omitting its interaction;
proc mixed data=work.ratings_all;
class trial_type participant_id;
model rating = trial_type trial_id_cond_concat / solution residual influence ddfm=kenwardroger2;
random intercept trial_id_cond_concat / subject=participant_id group=trial_type type=un g gcorr;
lsmeans trial_type;
lsmestimate trial_type 'erythritol versus sucrose' 1 0 -1,
	'erythritol versus sucralose' 1 -1 / adjust=bon stepdown adjdfe=row;
run;
* NOTES: 
1. the fit of this model without the interaction is slightly better than the model with the interaction, supporting its omission
2. the trial_id main effect is far from significant, so we omit it;

* model only including main effect of trial type;
proc mixed data=work.ratings_all;
class trial_type participant_id;
model rating = trial_type / solution residual influence ddfm=kenwardroger2;
random intercept / subject=participant_id group=trial_type type=un g gcorr;
lsmeans trial_type;
lsmestimate trial_type 'erythritol versus sucrose' 1 0 -1,
	'erythritol versus sucralose' 1 -1 / adjust=bon stepdown adjdfe=row;
run;
* NOTE: 
1. the fit is slightly worse than the model with the main effect of trial_id, but we do want to omit it for parsimony purposes, and the results are virtually identical
2. this is the final model to be reported;


* COVARIATES
------------;

* CHECK DISTRIBUTIONS;
* across conditions;
proc univariate data=work.ratings_all;
where trial_id_cond_concat = 1;
var hunger intensity liking concentration;
histogram hunger intensity liking concentration / normal (mu=est sigma=est) lognormal (sigma=est theta=est zeta=est);
run;
* distributions quite funky - check residuals of models;

* by condition;
proc univariate data=work.ratings_all;
where trial_id_cond_concat = 1;
by trial_type;
var intensity liking concentration;
histogram intensity liking concentration / kernel normal (mu=est sigma=est) lognormal (sigma=est theta=est zeta=est);
run;

proc univariate data=work.ratings_all;
where trial_id_cond_concat = 1;
class trial_type;
var intensity liking concentration;
histogram intensity liking concentration / kernel overlay;
run;


* MIXED MODELS;
proc mixed data=work.ratings_all;
where trial_id_cond_concat = 1;
class trial_type participant_id;
model intensity = trial_type / solution residual influence ddfm=kenwardroger2;
random participant_id / type=cs g gcorr;
lsmeans trial_type;
lsmestimate trial_type 'erythritol versus sucrose' 1 0 -1,
	'erythritol versus sucralose' 1 -1 / adjust=bon stepdown adjdfe=row;
run;

proc mixed data=work.ratings_all;
class trial_type participant_id;
where trial_id_cond_concat = 1;
model liking = trial_type / solution residual influence ddfm=kenwardroger2;
random participant_id / type=cs g gcorr;
lsmeans trial_type;
lsmestimate trial_type 'erythritol versus sucrose' 1 0 -1,
	'erythritol versus sucralose' 1 -1 / adjust=bon stepdown adjdfe=row;
run;

* adjust final model for covariates;
proc mixed data=work.ratings_all;
class trial_type participant_id;
model rating = trial_type hunger intensity / solution residual influence ddfm=kenwardroger2;
random intercept / subject=participant_id group=trial_type type=un g gcorr;
lsmeans trial_type;
lsmestimate trial_type 'erythritol versus sucrose' 1 0 -1,
	'erythritol versus sucralose' 1 -1 / adjust=bon stepdown adjdfe=row;
run;
* NOTE: these models excludes the water condition since it does not have intensity ratings;

* test interaction between intensity and trial type;
proc mixed data=work.ratings_all;
class trial_type participant_id;
model rating = trial_type | intensity hunger / solution residual influence ddfm=kenwardroger2;
random intercept / subject=participant_id group=trial_type type=un g gcorr;
lsmeans trial_type;
lsmestimate trial_type 'erythritol versus sucrose' 1 0 -1,
	'erythritol versus sucralose' 1 -1 / adjust=bon stepdown adjdfe=row;
run;
* effect of intensity does not differ between trial types!

* concentration (excluding sucrose, as it always has the same concentration);
proc mixed data=work.ratings_all;
where trial_type NE 'sucrose';
class trial_type participant_id;
model rating = trial_type | concentration / solution residual influence ddfm=kenwardroger2;
random intercept / subject=participant_id group=trial_type type=un g gcorr;
lsmeans trial_type;
estimate 'slope of concentration in erythritol' concentration 1 concentration*trial_type 1;
estimate 'slope of concentration in sucralose' concentration 1 concentration*trial_type 0 1;
run;
* no effect of concentration, and it does not differ between ery and sucra!;


* CREATE FILE WITH AVERAGES PER CONDITION PER SUBJECT
-----------------------------------------------------;

proc means data=work.ratings_all noprint;
var wanting liking;
output out=work.ratings_means_long mean= ;
class cond_trial_type participant_id;
types () cond_trial_type*participant_id;
run;

data work.ratings_means_long;
set work.ratings_means_long;
if _TYPE_ = 0 then delete;
drop _TYPE_ _FREQ_;
run;

proc sort data=work.ratings_means_long;
by participant_id cond_trial_type;
run;

proc transpose data=work.ratings_means_long out=work.means_wide_wanting (drop=_NAME_ rename=(bit_high_calorie=want_bit_high_calorie bit_low_calorie=want_bit_low_calorie bit_neutral=want_bit_neutral pla_high_calorie=want_pla_high_calorie pla_low_calorie=want_pla_low_calorie pla_neutral=want_pla_neutral));
var wanting;
by participant_id; 
id cond_trial_type;
run;

proc transpose data=work.ratings_means_long out=work.means_wide_liking (drop=_NAME_ rename=(bit_high_calorie=like_bit_high_calorie bit_low_calorie=like_bit_low_calorie bit_neutral=like_bit_neutral pla_high_calorie=like_pla_high_calorie pla_low_calorie=like_pla_low_calorie pla_neutral=like_pla_neutral));
var liking;
by participant_id; 
id cond_trial_type;
run;

data work.ratings_means;
merge work.means_wide_liking work.means_wide_wanting;
by participant_id;
run;

* sanity check;
proc univariate data=work.ratings_means;
var concentration_erythritol concentration_sucralose;
run;

* simple correlations;
proc corr data=work.ratings_means plots=all;
var concentration_erythritol intensity_erythritol liking_erythritol rating_erythritol hunger;
run;

proc corr data=work.ratings_means plots=all;
var concentration_sucralose intensity_sucralose liking_sucralose rating_sucralose hunger;
run;

proc corr data=work.ratings_means plots=all;
var intensity_sucrose liking_sucrose rating_sucrose hunger;
run;

proc corr data=work.ratings_means plots=all;
var liking_sucrose liking_erythritol liking_sucralose;
run;

proc corr data=work.ratings_means plots=all;
var rating_sucrose rating_erythritol rating_sucralose rating_water;
run;

proc corr data=work.ratings_means plots=all;
var intensity_sucrose intensity_erythritol intensity_sucralose;
run;


* tidy up library;
*proc datasets library=ery_4a nodetails nolist;
*delete ratings_means_long means_wide_hunger means_wide_intensity means_wide_liking means_wide_rating;
*run;

* export to excel/tsv for use in Matlab;
proc export data=work.ratings_means outfile="&path\ratings_means.tsv" DBMS=TAB;
PUTNAMES=YES;
run;
