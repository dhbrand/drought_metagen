proc import datafile="Z:\GitHub\drought_metagen\SAS\hf.sas.csv"
     out=hf
     dbms=csv
     replace;
     getnames=no;
run;

proc print data=hf (obs=10);
run;

proc import datafile="Z:\GitHub\drought_metagen\SAS\Y.hf.csv"
     out=Y
     dbms=csv
     replace;
     getnames=yes;
run;

proc print data=Y (obs=10);
run;

proc transpose data=X out=X_t;
run;
proc print data=hf;
var V1-V10;
run;

ods graphics on;

proc glmselect data=hf plots=all;
model V4 = V5-V20 /
        selection=elasticnet (choose=bic);
run;
