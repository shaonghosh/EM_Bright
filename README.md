# EM_Bright
This repository contains the tools required for the "EM-Bright analysis" of LIGO-Virgo gravitational wave triggers. You will need LValert tools and pycbc installation to run the EM-Bright analysis.


*** Subscribe to the cbc_gstlal node ***
lvalert_admin -a <user.name> -b <your password> --subscribe --node cbc_gstlal

*** Start listening ***
lvalert_listen -a <user.name> -b <your password> -c /home/shaon/analysis/o2/EM_Bright_classification/auto_EM_Bright/myLVAlertListen.ini -r alertInstance & 


Important: Need PyCBC installation 
