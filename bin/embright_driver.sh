#!/bin/bash
. ~cbc/pe/local/etc/lalinference-user-env.sh
/home/gracedb.processor/opt/bin/em_bright_classifier --configfile $2 --graceid $4
