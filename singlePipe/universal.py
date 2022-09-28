from termcolor import colored
import subprocess
from datetime import datetime

def run_cmd (c,outputdir):
    print(colored(c, 'green'))
    with open(outputdir+'/'+'log.txt','a') as logfile:
        logfile.write(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        process = subprocess.Popen(c,
                                   stdout=subprocess.PIPE, 
                                   stderr=logfile)
        stdout, stderr = process.communicate()
        stdout, stderr
    if process.returncode != 0:
        err_msg = ["error in the code, check log file:"]
        raise Exception(err_msg)
        
def run_cmd_file (mycmd,f,outputdir):
    print(colored(mycmd, 'green'))
    with open(outputdir+'/'+'log.txt','a') as logfile:
        logfile.write(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        process = subprocess.Popen(mycmd,
                                   stdout=f, 
                                   stderr=logfile)
        stdout, stderr = process.communicate()
        stdout, stderr 
    if process.returncode != 0:
        err_msg = ["error in the code, check log file:"]
        raise Exception(err_msg)
