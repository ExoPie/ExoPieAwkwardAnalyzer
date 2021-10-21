import copy
binning={}
binning["SR_1b"]={"metpt":[40,0.,1000.],
                  "metphi":[40,-3.14,3.14],
                  #"jetpt0":[50,0,1000.],
                  #"jetpt1":[50,0,1000.],
                  #"jeteta0":[50,-2.5,2.5],
                  #"jeteta1":[50,-2.5,2.5],
                  #"csv0": [20,0,1.],
                  "nTrueInt": [50,0,100.],
                  "nJetLoose":[5,0,5],
                  "nEleLoose":[3,0,3],
                  "min_dphi_jet_met":[50,-5,5]
                  
              }

binning["SR_2b"] =copy.deepcopy(binning["SR_1b"])
binning["SR_2b"]["jetpt0"] = [25,0,1000.]

binning["ZeeCR_2b"]=copy.deepcopy(binning["SR_1b"])

binning["ZeeCR_1b"]=copy.deepcopy(binning["SR_1b"])

binning["ZmumuCR_2b"]=copy.deepcopy(binning["SR_1b"])

binning["ZmumuCR_1b"]=copy.deepcopy(binning["SR_1b"])

binning["TopenuCR_2b"]=copy.deepcopy(binning["SR_1b"])

binning["TopenuCR_1b"]=copy.deepcopy(binning["SR_1b"])

binning["TopmunuCR_2b"]=copy.deepcopy(binning["SR_1b"])

binning["TopmunuCR_1b"]=copy.deepcopy(binning["SR_1b"])

binning["WenuCR_2b"]=copy.deepcopy(binning["SR_1b"])

binning["WenuCR_1b"]=copy.deepcopy(binning["SR_1b"])

binning["WmunuCR_2b"]=copy.deepcopy(binning["SR_1b"])

binning["WmunuCR_1b"]=copy.deepcopy(binning["SR_1b"])



