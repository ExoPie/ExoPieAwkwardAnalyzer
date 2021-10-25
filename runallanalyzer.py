import os 
import subprocess


def submitprocess(filename):
    command="python bbDMAnalyzer_onezip_iterate.py "+filename+" &"
    os.system(command)
    return True

def getNprocesses():
    n_running_processes = subprocess.check_output("ps -aux | grep \"python bbDMAnalyzer_onezip_iterate.py\" | wc -l",shell=True)
    return n_running_processes

counter=0
n_running_processes=0

for ifile in open('list.txt'):
    issubmitted=False
    print ("\n   --------------- ----- ")
    print ("runnning process for: ", ifile.rstrip())
    print ("number of running processes", int(n_running_processes))

    canbesubmitted=False
    while (canbesubmitted==False) & (issubmitted==False) :
        n_running_processes=getNprocesses()
        #print ("waiting inside while loop untill processes finish")
        if (int(n_running_processes)>12):
            os.system("sleep 5")
            print ("sleeping")
        if (int(n_running_processes)<12):
            print ("processes under 10, will submit now")
            issubmitted=submitprocess(ifile.rstrip())
            canbesubmitted=True
    if issubmitted:
        print ("job number %f submitted",counter)
    counter=counter+1
    
        
        
        

