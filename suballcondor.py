import os 

eospath="/eos/cms/store/group/phys_exotica/bbMET/2017_SkimmedFiles/skim_v17_12-01_type1met_taupt18_pct"
tosubmit="condor_test"
inputtextfile="skimfileslist/samplespplit/"
os.system("cp pyfiles.tar "+tosubmit)
os.system("cp data.tar "+tosubmit)
os.system("cp runanalysis.sh "+tosubmit)


def filetolist(textfile):
    flist=[]
    for iline in open(textfile,'r'):
        flist.append(iline.rstrip())
    return flist

def createsubmitfile(filename):
    TEMP_SUB_FILE="""
    universe = vanilla
    request_memory = 8192
    Proxy_filename = x509up
    Proxy_path = /afs/cern.ch/user/k/khurana/private/$(Proxy_filename)
    request_cpus = 4
    +JobFlavour = "nextweek"
    executable = runanalysis.sh
    should_transfer_files = YES
    output = output_c/condor.$(Cluster).$(Process).out
    error = error_c/condor.$(Cluster).$(Process).err
    log = log_c/condor.$(Cluster).$(Process).log
    transfer_input_files = ./pyfiles.tar, ./data.tar, requirements.txt
    on_exit_remove = (ExitBySignal == False) && (ExitCode == 0)
    on_exit_hold = ( (ExitBySignal == True) || (ExitCode != 0) )
    on_exit_hold_reason = strcat("Job held by ON_EXIT_HOLD due to ",ifThenElse((ExitBySignal == True), "exit by signal",strcat("exit code ",ExitCode)), ".")
    periodic_release =  (NumJobStarts < 5) && ((CurrentTime - EnteredCurrentStatus) > (60*60))
    arguments =  """+eospath+"""/$(jobid) 2017 $(jobid) $(Proxy_path)
    queue jobid from """+filename +"""
    """
    
    fsub_out =  tosubmit+"/"+filename.split("/")[-1].replace(".txt",".sub") 
    fsub=open(fsub_out,'w')
    fsub.write(TEMP_SUB_FILE)
    fsub.close()
    
    os.system("condor_submit "+fsub_out)

    return 0 



os.system("ls -1 "+inputtextfile+" >tmp.txt")
sampleliist=filetolist("tmp.txt")

for isample  in  sampleliist:
    fullsamplepath=inputtextfile+isample
    print ("submitting jobs for: ",fullsamplepath)
    createsubmitfile(fullsamplepath)#"skimfileslist/samplespplit/QCD_bEnriched_HT1000to1500_TuneCP5_13TeV-madgraph-pythia8.txt")

