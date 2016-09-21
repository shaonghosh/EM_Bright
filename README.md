# EM_Bright
This repository contains the tools required for the "EM-Bright analysis" of LIGO-Virgo gravitational wave triggers. You will need LValert tools and pycbc installation to run the EM-Bright analysis.

Run >> python initiate.py --user user.name --password <ligo.org password>
to create the LV alert ini file and start the jobs.

*** Subscribe to the cbc_gstlal node ***
lvalert_admin --user <user.name> --password <your password> --subscribe --node cbc_gstlal

*** Start listening ***
lvalert_listen --user <user.name> --password <your password> -c /home/shaon/analysis/o2/EM_Bright_classification/auto_EM_Bright/myLVAlertListen.ini -r alertInstance & 

Notes: 

Need PyCBC installation 

Make sure your lalsuite installation is later than 23rd June.


Sample run command:

python initiate.py --user shaon.ghosh --netrcpath /home/shaon/.netrc



