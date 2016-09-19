import os
import sys
import optparse

parser = optparse.OptionParser()
parser.add_option("-U", "--user", action="store", type="string", metavar="NAME", help="username albert.einstein")
parser.add_option("-P", "--password", action="store", type="string", metavar="NAME", help="lv alert password")

(opts,args) = parser.parse_args()


work_dir = os.getcwd()

text = '''[cbc_gstlal_mdc]
nodes = test_gstlal
executable = ''' + work_dir + '''/sourceClassDriver.py

[cbc_gstlal_highmass]
nodes = test_gstlal
executable = ''' + work_dir + '''/sourceClassDriver.py

'''

iniFile = open('myLVAlertListen.ini', 'w')

iniFile.writelines(text)
iniFile.close()

os.system('chmod 777 sourceClassDriver.py')
lvalert_cmd = 'lvalert_listen -a ' + opts.user + ' -b ' + opts.password + ' -c ' + work_dir + '/myLVAlertListen.ini -r alertInstance &'

print 'Launching LValert job...'
os.system(lvalert_cmd)
print 'Done...'
