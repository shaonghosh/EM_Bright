#!/usr/bin/env python
import os
import json
import numpy as np
from sys import stdin
from sys import exit
import time as Time
import datetime
import ConfigParser

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
    ### Reed: "reading in an environmental variable here seems fine, but also overkill. 
    ### I don't see where this variable is set in initiate.py, so this could be fragile. 
    ### I'd suggest specifying it as a keyword argument instead." - Clarify this from Reed.
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
    mass1 = coinc.get_column('mass1')[0]
    mass2 = coinc.get_column('mass2')[0]
    chi1 = coinc.get_column('spin1z')[0]
    time = coinc.get_column('end_time')[0]
    snr1 = coinc.get_column('snr')[0]
    snr2 = coinc.get_column('snr')[1]
    snr = np.max([snr1, snr2]) ## Only returning the largest snr

    return [float(mass1), float(mass2), float(chi1), int(time), float(snr)]



#########################################################################################

### Reading information from config file ###
configParser = ConfigParser.ConfigParser()
configFile = 'configFile.ini'
configParser.read(configFile)
coinc_path = configParser.get('Paths', 'coincPath')
psd_path = configParser.get('Paths', 'psdPath')
source_class_path = configParser.get('Paths', 'results')
log_path = configParser.get('Paths', 'logs')

ellipsoidSample = int( configParser.get('EMBright', 'elipsoidSample') )
diskMassThreshold = float( configParser.get('EMBright', 'diskMassThreshold') )
f_low = float( configParser.get('EMBright', 'fmin') )

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

if alert_type == 'new':
    graceid = str(streamdata['uid'])
    try:
        os.system('mkdir -p ' + log_path)
    except:
        print 'Could not create logs directory...'
        exit(1)
    log = open(log_path + '/log_' + graceid + '.txt', 'a')
    log.writelines('\n' + str(datetime.datetime.today()) + '\t' + 'Analyzing event: ' + graceid + '\n')
    try:
        os.system('mkdir -p ' + coinc_path)
        log.writelines(str(datetime.datetime.today()) + '\t' + 'Successfully created coinc directory\n')

        os.system('mkdir -p ' + psd_path)
        log.writelines(str(datetime.datetime.today()) + '\t' + 'Successfully created psd directory\n')

        os.system('mkdir -p ' + source_class_path)
        log.writelines(str(datetime.datetime.today()) + '\t' + 'Successfully created results directory\n')

    except:
        log.writelines(str(datetime.datetime.today()) + '\t' + 'Failed to creat coinc and/or psd and/or results directory, Check write privilege to the given path\n')
        exit(1)


    x = 1
    countTrials = 0
    while x == 1:
        log.writelines(str(datetime.datetime.today()) + '\t' + 'Fetching coinc and psd file. Trial number: ' +  str(countTrials+1) + '\n')
        x = getCoinc(graceid, coinc_path, psd_path)
        if countTrials >= 5:
            log.writelines(str(datetime.datetime.today()) + '\t' + 'Could not fetch coinc and/or psd files\n')
            exit(1)
        if x == 1: Time.sleep(5) ### Wait for five seconds if getCoinc is unsuccessful
        countTrials += 1

    log.writelines(str(datetime.datetime.today()) + '\t' + 'Successfully fetched coinc and/or psd files\n')


    ### Check if this event  has been analyzed ### 
    ### This is not required anymoure since I am now only analyzing new events. Remove this ###
    if os.path.isfile(source_class_path + '/Source_Classification_' + graceid + '_.dat'):
        log.writelines(str(datetime.datetime.today()) + '\t' + 'Event already analyzed... skipping\n')
        exit(0)


    start = Time.time()
    coincFileName = [coinc_path + '/coinc_' + graceid + '.xml']
    [mass1, mass2, chi1, time, snr] = readCoinc(coincFileName)

    File = open(coinc_path + '/masses_chi1_' + graceid + '_.dat', 'w')
    File.writelines(graceid + '\t' + str(mass1) + '\t' +  str(mass2) + '\t' + str(chi1) + '\t' + str(snr) + '\n')

    File.close()

    samples_sngl = getSamples(graceid, mass1, mass2, chi1, snr, ellipsoidSample, {'H1=' + psd_path + '/psd_' + graceid + '.xml.gz'}, {'L1=' + psd_path + '/psd_' + graceid + '.xml.gz'}, fmin=f_low, saveData=True)
    log.writelines(str(datetime.datetime.today()) + '\t' + 'Created ambiguity ellipsoid samples\n')

    ### Currently NaNs are generated when the ellipsoid generation failed. This will be changed in subsequent version. ###
    if ~np.any( np.isnan(samples_sngl[0]) ): 
        diskMassObject_sngl = genDiskMassProbability.genDiskMass(samples_sngl, 'test', diskMassThreshold)
        [NS_prob_1_sngl, NS_prob_2_sngl, diskMass_sngl] = diskMassObject_sngl.fromEllipsoidSample()
        em_bright_prob_sngl = np.sum((diskMass_sngl > 0.)*100./len(diskMass_sngl))

    else:
        log.writelines(str(datetime.datetime.today()) + '\t' + 'Return was NaNs\n')
        [NS_prob_2_sngl, em_bright_prob_sngl] = [0., 0.]
        message = 'EM-Bright probabilities computation failed for trigger + ' + graceid + '\n'

        gdb.writeLog(graceid, message, tagname='em_follow')
        end = Time.time()
        log.writelines(str(datetime.datetime.today()) + '\t' + 'Time taken in computing EM-Bright probabilities = ' + str(end - start) + '\n')        
        exit(0)
        

    end = Time.time()
    log.writelines(str(datetime.datetime.today()) + '\t' + 'Time taken in computing EM-Bright probabilities = ' + str(end - start) + '\n')

    
    ### Find an appropriate use of this file (e.g. uploading this as JSON) else get rid of it ###
    source_classification = open(source_class_path + '/Source_Classification_' + graceid + '_.dat', 'w')
    source_classification.writelines('The probability of second object being a neutron star for the trigger ' + graceid + ' = ' + str(NS_prob_2_sngl) + '% \n')
    source_classification.writelines('The probability of remnant mass outside the black hole in excess of ' + str(diskMassThreshold) + ' M_sun for the trigger ' + graceid + ' = '  + str(NS_prob_2_sngl) + '% \n')

    source_classification.close()

    message = 'Computed from detection pipeline: The probability of second object being a neutron star  = ' + str(NS_prob_2_sngl) + '% \n The probability of remnant mass outside the black hole in excess of ' + str(diskMassThreshold) + ' M_sun = '  + str(NS_prob_2_sngl) + '% \n'

    gdb.writeLog(graceid, message, tagname='em_follow')






