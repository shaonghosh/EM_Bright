#!/usr/bin/env python
import os
import json
import numpy as np
import sys
from sys import stdin
from sys import exit
import time as Time
import datetime
import ConfigParser

import ligo.gracedb.rest
from pylal import SnglInspiralUtils

import genDiskMassProbability
from getEllipsoidSamples import getSamples



def getCoinc(graceid, gracedb_url, coinc_path, psd_path):
    '''
    Attempts to fetch a coinc and psd files of a given graceID from graceDB.
    Saves the files in the supplied coinc and psd paths.
    Returns zero if the fetching is successful.
    Returns one upon failure.
    '''
    gracedb = ligo.gracedb.rest.GraceDb( gracedb_url )

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
    Reads the coinc file. Finds the highest SNR IFO and returns the point estimates of the same.
    '''
    coinc = SnglInspiralUtils.ReadSnglInspiralFromFiles(CoincFile)
    ifo=[]; mass1=[]; mass2=[]; chi1=[]; snr= np.array([])
    for row in coinc:
        ifo.append(row.ifo)
        mass1.append(row.mass1)
        mass2.append(row.mass2)
        chi1.append(row.spin1z)
        snr = np.append(snr, row.snr)        
    index = np.argmax(snr) ## To return highest SNR IFO point estimates
    return [mass1[index], mass2[index], chi1[index], snr[index], str(ifo[index])]



#########################################################################################

### Reading information from config file ###
configParser = ConfigParser.ConfigParser()
configParser.read( sys.argv[1] )
gracedb_url = configParser.get('gracedb', 'gracedb_url')
coinc_path = configParser.get('Paths', 'coincPath')
psd_path = configParser.get('Paths', 'psdPath')
source_class_path = configParser.get('Paths', 'results')
log_path = configParser.get('Paths', 'logs')
numTrials = int( configParser.get('Paths', 'numTrials') )
wait = float( configParser.get('Paths', 'wait') )

ellipsoidSample = int( configParser.get('EMBright', 'elipsoidSample') )
diskMassThreshold = float( configParser.get('EMBright', 'diskMassThreshold') )
forced = configParser.getboolean('EMBright', 'Forced')
f_low = float( configParser.get('EMBright', 'fmin') )
mass1_cut = float( configParser.get('EMBright', 'mass1_cut') )
chi1_cut = float( configParser.get('EMBright', 'chi1_cut') )
lowMass_approx = configParser.get('EMBright', 'lowMass_approx')
highMass_approx = configParser.get('EMBright', 'highMass_approx')
tagnames = configParser.get('gracedb', 'tagnames').split() 

'''
Receives alerts from graceDB, obtains the required coinc and psd files and then launches
the EM-Bright classification jobs.
'''
# Load the LVAlert message contents into a dictionary
streamdata = json.loads(stdin.read())
gdb = ligo.gracedb.rest.GraceDb( gracedb_url )

#print streamdata
# Do something with new events having FAR below threshold

# alert_type = 'None'
# if streamdata['alert_type']:
#     alert_type = streamdata['alert_type']

if streamdata['alert_type'] == 'new':
    graceid = str(streamdata['uid'])
    try:
        os.system('mkdir -p ' + log_path)
    except:
        print 'Could not create logs directory...'
        exit(1)
    
    logFileName = log_path + '/log_' + graceid + '.txt'
    log = open(logFileName, 'a')
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
    
    for countTrials in xrange(numTrials): ### iterate a maximum of numTrials times
        log.writelines(str(datetime.datetime.today()) + '\t' + 'Fetching coinc and psd file. Trial number: ' +  str(countTrials+1) + '\n')
        if getCoinc(graceid, gracedb_url, coinc_path, psd_path): ### this failed, so we sleep
            Time.sleep(wait)
        else: ### success! so we exit the loop
            break
    else: ### we did not break from the loop, so we must have timed out
        log.writelines(str(datetime.datetime.today()) + '\t' + 'Could not fetch coinc and/or psd files\n')       
        exit(1)

    log.writelines(str(datetime.datetime.today()) + '\t' + 'Successfully fetched coinc and/or psd files\n')



    start = Time.time()
    coincFileName = [coinc_path + '/coinc_' + graceid + '.xml']
    [mass1, mass2, chi1, snr, ifo] = readCoinc(coincFileName)

    File = open(coinc_path + '/masses_chi1_' + graceid + '_.dat', 'w')
    File.writelines(graceid + '\t' + str(mass1) + '\t' +  str(mass2) + '\t' + str(chi1) + '\t' + str(snr) + '\n')

    File.close()

    samples_sngl = getSamples(graceid, mass1, mass2, chi1, snr, ellipsoidSample, {ifo + '=' + psd_path + '/psd_' + graceid + '.xml.gz'}, fmin=f_low, NMcs=10, NEtas=10, NChis=10, mass1_cut=mass1_cut, chi1_cut=chi1_cut, lowMass_approx=lowMass_approx, highMass_approx=highMass_approx, Forced=forced, logFile=logFileName, saveData=True)
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
    source_classification.writelines('The probability of remnant mass outside the black hole in excess of ' + str(diskMassThreshold) + ' M_sun for the trigger ' + graceid + ' = '  + str(em_bright_prob_sngl) + '% \n')

    source_classification.close()

    message = 'EM-Bright probabilities computed from detection pipeline: The probability of second object being a neutron star  = ' + str(NS_prob_2_sngl) + '% \n The probability of remnant mass outside the black hole in excess of ' + str(diskMassThreshold) + ' M_sun = '  + str(em_bright_prob_sngl) + '% \n'

    gdb.writeLog(graceid, message, tagname=tagnames)

