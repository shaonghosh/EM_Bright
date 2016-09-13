#!/usr/bin/env python
import os
import json
import numpy as np
from sys import stdin
from sys import exit
import time as Time

import ligo.gracedb.rest
from pylal import SnglInspiralUtils

import genDiskMassProbability
from getEllipsoidSamples import getSamples



def getCoinc(graceid, coinc_path, psd_path):
    '''
    Attempts to fetch a coinc and psd files of a given graceID from graceDB.
    Saves the files in the supplied coinc and psd paths.
    Returns zero if the fetching is successful.
    Returns one upon failure.
    '''
    gracedb = ligo.gracedb.rest.GraceDb(os.environ.get('GRACEDB_SERVICE_URL', ligo.gracedb.rest.DEFAULT_SERVICE_URL))

    try:
        coinc_object = gracedb.files(graceid, "coinc.xml")
        psd_object = gracedb.files(graceid, "psd.xml.gz")

        coinc_file = open(coinc_path + '/coinc_' + graceid + '.xml', 'w')
        psd_file = open(psd_path + '/psd_' + graceid + '.xml.gz', 'w')

        coinc_file.writelines(coinc_object.read())
        psd_file.writelines(psd_object.read())

        coinc_file.close()
        psd_file.close()
        return 0
    except:
        return 1


def readCoinc(CoincFile):
    '''
    Reads the coinc file and returns the masses, chi1 and gps times as a list.
    '''
    coinc = SnglInspiralUtils.ReadSnglInspiralFromFiles(CoincFile)
    for coinc_index, coinc_row in enumerate(coinc):
        mass1 = coinc_row.mass1
        mass2 = coinc_row.mass2
        chi1 = coinc_row.spin1z
        time = coinc_row.end_time
    return [mass1, mass2, chi1, time]



#########################################################################################
'''
Receives alerts from graceDB, obtains the required coinc and psd files and then launches
the EM-Bright classification jobs.
'''
# Load the LVAlert message contents into a dictionary
streamdata = json.loads(stdin.read())
gdb = ligo.gracedb.rest.GraceDb()

#print streamdata
# Do something with new events having FAR below threshold

alert_type = 'None'
if streamdata['alert_type']:
    alert_type = streamdata['alert_type']

#if alert_type == 'new':
    # The object is a serialized event. Get the FAR

graceid = str(streamdata['uid'])
coinc_path = 'all_coincs' ### Currently hardcoded
psd_path = 'all_psds' ### Currently hardcoded
source_class_path = 'all_source_classifications' ### Currently hardcoded
os.system('mkdir -p ' + coinc_path)
os.system('mkdir -p ' + psd_path)
os.system('mkdir -p ' + source_class_path)

x = 1
###CHECK! This can cause an infinite loop. Discuss this with the reviewers
while x == 1:
    x = getCoinc(graceid, coinc_path, psd_path)

### Check if this event  has been analyzed ###
if os.path.isfile(source_class_path + '/Source_Classification_' + graceid + '_.dat'):
        print 'Event already analyzed... skipping'
        exit(0)


### Check if the required coinc file and the psd file is present ###
if os.path.isfile(coinc_path + '/coinc_' + graceid + '.xml') and os.path.isfile(psd_path + '/psd_' + graceid + '.xml.gz'):
    start = Time.time()
    coincFileName = [coinc_path + '/coinc_' + graceid + '.xml']
    [mass1, mass2, chi1, time] = readCoinc(coincFileName)

else:
    print 'Did not find the coinc and psd files... quitting...'
    exit(1)


File = open(coinc_path + '/masses_chi1_' + graceid + '_.dat', 'w')
File.writelines(graceid + '\t' + str(mass1) + '\t' +  str(mass2) + '\t' + str(chi1) + '\n')

File.close()

samples_sngl = getSamples(graceid, mass1, mass2, chi1, 1000, {'H1=' + psd_path + '/psd_' + graceid + '.xml.gz'}, {'L1=' + psd_path + '/psd_' + graceid + '.xml.gz'}, saveData=True)
if ~np.any( np.isnan(samples_sngl[0]) ):
    diskMassObject_sngl = genDiskMassProbability.genDiskMass(samples_sngl, 'test', 0.03)
    [NS_prob_1_sngl, NS_prob_2_sngl, diskMass_sngl] = diskMassObject_sngl.fromEllipsoidSample()
    em_bright_prob_sngl = np.sum((diskMass_sngl > 0.)*100./len(diskMass_sngl))

else:
    print 'Return was NaNs'
    [NS_prob_2_sngl, em_bright_prob_sngl] = [0., 0.]

end = Time.time()
print '*** Time taken in computing probabilities = ' + str(end - start)

source_classification = open(source_class_path + '/Source_Classification_' + graceid + '_.dat', 'w')
source_classification.writelines('The probability of second object being a neutron star for the trigger ' + graceid + ' = ' + str(NS_prob_2_sngl) + '\n')
source_classification.writelines('The probability of remnant mass in the disk in excess of 0.03 M_sun for the trigger ' + graceid + ' = '  + str(NS_prob_2_sngl) + '\n')

source_classification.close()

message = 'Computed from detection pipeline: The probability of second object being a neutron star  = ' + str(NS_prob_2_sngl) + '\n The probability of remnant mass in the disk in excess of 0.03 M_sun = '  + str(NS_prob_2_sngl) + '\n'

gdb.writeLog(graceid, message)
gdb.writeLog(graceid, message, tagname='em_follow')






