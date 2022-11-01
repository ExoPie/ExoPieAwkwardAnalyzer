#!/bin/sh

export SCRAM_ARCH=slc7_amd64_gcc820

inputfilename=$1
year=$2

echo "inputfilename is:" $inputfilename
echo "running the code for: "$year

#cp /eos/cms/store/group/phys_exotica/bbMET/2017testanalyzer/mypack.tar.bz2 .
#echo "untarring the virtual box now"
#tar -xf mypack.tar.bz2
#scl enable rh-python36 bash
#echo "enabled python 3"
#ls -ltr /opt/rh/rh-python36/

basepath=$PWD

##########################################
## setup the environment 
export PYTHONPATH=/usr/bin/python3.6
python3 -m venv venv2
ls -ltr 
source venv2/bin/activate
pip install -U pip requests
pip install uproot4
pip install -r requirements.txt
ls -ltr 
source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.22.02/x86_64-centos7-gcc48-opt/bin/thisroot.sh
root -l -b -q 
echo "root worked"
source venv2/bin/activate

#########
## packaged installed

##########################################

#locate python3

#virtualenv --system-site-packages mypack
#cd mypack 
#source bin/activate
#echo "installing requirements.txt"

#python test.py 


#echo "untarrng done, now activating"
#source mypack/bin/activate 

echo "virtual env activated, now untarring the python files and yaml files, needed to run the code"
cd $basepath
tar -xf pyfiles.tar 
tar -xf data.tar
echo "getting the proxy files"
export X509_USER_PROXY=$4
voms-proxy-info -all
voms-proxy-info -all -file $4
ls -ltr 
echo "time for action-----------"

echo "present path" $PWD
echo "base path" $basepath
mkdir output 

ls -ltr 
python bbDMAnalyzer_onezip_iterate.py $inputfilename $year

echo "listing after running "
ls -ltr 
echo "checking for output"
ls -ltr output

xrdcp -f output/"$3" root://eoscms.cern.ch//eos/cms/store/group/phys_exotica/bbMET/2017testanalyzer/"$3"
#if [ -e "$3" ]; then
#  until xrdcp -f output/"$3" root://eoscms.cern.ch//eos/cms/store/group/phys_exotica/bbMET/2017testanalyzer/"$3"; do 
#    sleep 5
#    echo "Retrying"
#  done

#fi

#exitcode=$?

#if [ ! -e "$3" ]; then
#  echo "Error: The python script failed, could not create the output file."
#  
#fi
#exit $exitcode
