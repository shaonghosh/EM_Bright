# EM_Bright
This repository contains the tools required for the "EM-Bright analysis" of LIGO-Virgo gravitational wave triggers. 


1. Sample run command: If running in older version of LALsuite (before Nov 2016) or if running on lalinference_o2 branch then you can directly run em_bright_classifier. 

em_bright_classifier --configfile etc/emBright.ini --graceid M262774

2. If running on latest version of LALsuite, then run the driver script as follows:

./embright_driver.sh etc/emBright.ini <graceid>

You might have to chmod 777 the embright_driver.sh script after you pull the repository.


