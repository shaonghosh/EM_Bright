import os
import ligo.gracedb.logging
import ligo.gracedb.rest


def getCoinc(graceid, coinc_path, psd_path):
    gracedb = ligo.gracedb.rest.GraceDb(os.environ.get('GRACEDB_SERVICE_URL',ligo.gracedb.rest.DEFAULT_SERVICE_URL))

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





