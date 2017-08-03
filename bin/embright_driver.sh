#!/bin/bash
source /home/shaon/.profile
source /home/gracedb.processor/users/shaon/setup.sh
/home/gracedb.processor/users/shaon/EM_Bright/bin/em_bright_classifier --configfile $2 --graceid $4
