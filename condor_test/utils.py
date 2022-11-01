import sys 
import uproot4 
import numpy 
import math 
import time
import awkward1 as ak 
import numpy as np
import pandas as pd
from ROOT import TH1F, TFile
import  matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from array import array


def FileToList(filename):
    return ([i.rstrip() for i in open(filename)])
def SetHist(HISTNAME, binning):
    h = TH1F()
    if len(binning) == 3:
        h = TH1F(HISTNAME, HISTNAME, binning[0], binning[1], binning[2])
    else:
        nBins = len(binning) - 1
        #h = TH1F(HISTNAME, HISTNAME, binning[0], binning[1], binning[2])  ## make it variable binning histogram                                                                 \
                                                                                                                                                                                  
        h = TH1F(HISTNAME, HISTNAME, nBins, array('d', binning))
    return h



import copy
def VarToHist(df_var, df_weight, HISTNAME, binning):
    binning_ = copy.deepcopy(binning)
    df_var = pd.Series(df_var)

    h_var = SetHist(HISTNAME, binning_)
    weight = df_weight

    if len(binning_) >3:
        binning_.append(10000) ## to take care of overflow                                                                                                                         
        n, bins, patches = plt.hist(df_var, binning_, histtype='step', weights=weight)


    if len(binning_)==3:
        binning_.append(binning_[-1]*3) ## to take care of overflow                                                                                                                 
        n, bins, patches = plt.hist(df_var, binning_[0], range=(binning_[1], binning_[2]), histtype='step', weights=weight)
    ## this is outside if                                                                                                                                                        \
                                                                                                                                                                                  
    n=list(n)
    n_last = n[-1]
    n.remove(n_last)
    n[-1]  = n[-1]  + n_last
    for ibin in range(len((n))):
        h_var.SetBinContent(ibin+1, n[ibin])

    return h_var



'''
year='2016'
year_file = open("Year.py", "w")
if year == '2016':
    print('code is running for 2016')
    year_file.write('era="2016"')
elif year == '2017':
    print('code is running for 2017')
    year_file.write('era="2017"')
elif year == '2018':
    print('code is running for 2018')
    year_file.write('era="2018"')
else:
    print('please provide year')
    sys.exit()
year_file.close()

import ana_weight as wgt



def weight_( ep_WmunuRecoil, nEle, nMu, ep_muPt, ep_muEta):#, ep_ZmumuRecoil, ep_WmunuRecoil, nEle, ep_elePt, ep_eleEta, ep_eleIsPTight, nMu, ep_muPt, ep_muEta, ep_isTightMuon):
    total_weight=1.0
    if (nEle==0) & (nMu==1) :
        total_weight = weight_W1mu_(ep_WmunuRecoil, nEle, nMu, ep_muPt, ep_muEta)
        
    
    return total_weight

def weight_W1mu_(ep_WmunuRecoil, nEle, nMu, ep_muPt, ep_muEta):
    weightMET, weightMET_up, weightMET_down     = wgt.getMETtrig_First( ep_WmunuRecoil, 'R')
    weightMu, weightMu_up, weightMu_down        = wgt.mu_weight(ep_muPt, ep_muEta, 'T')
    return ( weightMET * weightMu)

'''

def getpt_eta_phi(mupx, mupy,mupz):
    mupt = numpy.sqrt(mupx**2 + mupy**2)
    mup = numpy.sqrt(mupx**2 + mupy**2 + mupz**2)
    mueta = numpy.log((mup + mupz)/(mup - mupz))/2
    muphi = numpy.arctan2(mupy, mupx)
    return (mupt, mueta, muphi)

def geteta(mupx, mupy,mupz):
    mup = numpy.sqrt(mupx**2 + mupy**2 + mupz**2)
    mueta = numpy.log((mup + mupz)/(mup - mupz))/2
    return (mueta)



def getphi(mupx, mupy):
    muphi = numpy.arctan2(mupy, mupx)
    return (muphi)



def getpt(mupx, mupy):
    mupt = numpy.sqrt(mupx**2 + mupy**2)
    return (mupt)

##def Phi_mpi_pi(x):
##    kPI=numpy.array(3.14159265)
##    kPI = kPI.repeat(len(x))
##    kTWOPI = 2 * kPI
##
##    #while ((x.any() >= kPI).any()): x = x - kTWOPI;
##    #while ((x.any() < -kPI).any()): x = x + kTWOPI;
##    while (ak.any(x >= kPI)): x = x - kTWOPI;
##    while (ak.any(x < -kPI)): x = x + kTWOPI;
##    return x;
##

## solution from 
## https://stackoverflow.com/questions/64282188/while-loop-equivalent-function-snippet-in-awkward-array/64285047#64285047

def Phi_mpi_pi(x):
    y = numpy.add(x, numpy.pi)
    y = numpy.mod(y, 2*numpy.pi)
    y = numpy.subtract(y, numpy.pi)
    return y

def DeltaPhi(phi1,phi2):
    #print (phi1-phi2)
    phi = Phi_mpi_pi(phi1-phi2)
    return abs(phi)

def getrecoil(nEle,elept,elephi,elepx_,elepy_,met_,metphi_):
    WenuRecoilPx = -( met_*numpy.cos(metphi_) + elepx_)
    WenuRecoilPy = -( met_*numpy.sin(metphi_) + elepy_)
    WenuRecoilPt = (numpy.sqrt(WenuRecoilPx**2  +  WenuRecoilPy**2))
    return WenuRecoilPt


def getrecoil1(elepx_,elepy_,met_,metphi_):
    WenuRecoilPx = -( met_*numpy.cos(metphi_) + elepx_)
    WenuRecoilPy = -( met_*numpy.sin(metphi_) + elepy_)
    WenuRecoilPt = (numpy.sqrt(WenuRecoilPx**2  +  WenuRecoilPy**2))
    return WenuRecoilPt



def getMT(nEle,elept,elephi,elepx_,elepy_,met_,metphi_):
    dphi = DeltaPhi(elephi,metphi_)
    MT = numpy.sqrt( 2 * elept * met_ * (1.0 - numpy.cos(dphi)) )
    return MT
    


def getRecoilPhi(nEle,elept,elephi,elepx_,elepy_,met_,metphi_):
    WenuRecoilPx = -( met_*numpy.cos(metphi_) + elepx_)
    WenuRecoilPy = -( met_*numpy.sin(metphi_) + elepy_)
    WenurecoilPhi = numpy.arctan2(WenuRecoilPx,WenuRecoilPy)
    return WenurecoilPhi
    
def Delta_R(eta1, eta2, phi1,phi2):
    deltaeta = eta1-eta2
    deltaphi = DeltaPhi(phi1,phi2)
    DR = numpy.sqrt ( deltaeta**2 + deltaphi**2 )
    return DR 

def jetcleaning(ak4eta, lepeta, ak4phi, lepphi, DRCut):
    ## usage: (obj_to_clean, obj_cleaned_against, so on
    dr_ = Delta_R(ak4eta, lepeta, ak4phi, lepphi)
    
    return (dr_ > DRCut)



def getFirstElement(x):
    if len(x)>0: return x[0]

def getSecondElement(x):
    if len(x)>1: return x[1]


def getTwoElement(x):
    if len(x)==1: return (x[0], x[0])
    if len(x)>1: return [x[0], x[1]]
    
def getNthElement(x,n):
    if len(x)>n: return x[n]

def getMinimum(x):
    if len(x)>0: return min(x)
        
def countTrue(x):
    if len(x)>0: return numpy.sum(x)

def deltaR(phoeta, jeteta, phophi, jetphi, cut_=0.4):
    phoeta_unzip, jeteta_unzip = ak.unzip(ak.cartesian([phoeta,jeteta], nested=True))
    phophi_unzip, jetphi_unzip = ak.unzip(ak.cartesian([phophi,jetphi], nested=True))
    deta_unzip = phoeta_unzip - jeteta_unzip
    dphi_unzip = Phi_mpi_pi(phophi_unzip - jetphi_unzip)
    dr_unzip = numpy.sqrt(deta_unzip**2 + dphi_unzip**2)
    dr_pho_jet_status = ak.any(dr_unzip<=cut_,axis=-1)  ## use axis in new version of awkward
    return dr_pho_jet_status



def getN(var_, i):
    return ak.mask(var_, ak.num(var_, axis=1)>i, highlevel=False)[:,i]
