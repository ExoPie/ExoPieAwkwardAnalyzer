import copy 
#regions=["SR_1b", "SR_2b", "ZeeCR_2b", "ZeeCR_1b", "ZmumuCR_2b", "ZmumuCR_1b", "TopenuCR_2b", "TopenuCR_1b", "TopmunuCR_2b", "TopmunuCR_1b", "WenuCR_2b", "WenuCR_1b", "WmunuCR_2b", "WmunuCR_1b"]
regions=["SR_1b", "SR_2b", "ZeeCR_2b", "ZmumuCR_2b", "TopenuCR_2b", "TopmunuCR_2b",  "WenuCR_1b", "WmunuCR_1b","ZeeCR_1b","ZmumuCR_1b"]

vardict={"metpt":"MET",
         "metphi":"METPhi",
         "jetpt0":"Jet1Pt",
         "jetpt1":"Jet2Pt",
         "jeteta0":"Jet1Eta",
         "jeteta1":"Jet2Eta",
         "csv0":"Jet1deepCSV",
         "csv1":"Jet2deepCSV",
         "recoil_Wmunu0":"Recoil",
         "recoil_Wenu0":"Recoil",
         "recoil_Zmumu0":"Recoil",
         "recoil_Zee0":"Recoil",
         "recoil_WmunuPhi0":"RecoilPhi",
         "recoil_WenuPhi0":"RecoilPhi",
         "Zee_recoilPhi":"RecoilPhi",
         "Zmumu_recoilPhi":"RecoilPhi",
         "nTrueInt":"nPV",
         "nJetLoose":"nJets",
         "nEleLoose":"NEle",
         
         "min_dphi_jet_met":"min_dPhi",
         
         "mupt0":"lep1_pT",
         "mupt1":"lep2_pT",
         "mueta0":"lep1_eta",
         "mueta1":"lep2_eta",
         
         "elept0":"lep1_pT",
         "elept1":"lep2_pT",
         "eleeta0":"lep1_eta",
         "eleeta1":"lep2_eta",
         
         "mt_Wmunu0":"Wmass",
         "mt_Wenu0":"Wmass",
         "pt_Wmunu0":"WpT",
         "pt_Wenu0":"WpT",
         "Zee_mass":"Zmass",
         "Zmumu_mass":"Zmass",
         "Zmumu_pt":"ZpT",
         "Zee_pt":"ZpT"
}


variables_common={"SR_1b":["metpt", "metphi", "jetpt0", "jeteta0", "csv0", "nTrueInt", "nJetLoose", "nEleLoose", "min_dphi_jet_met"]}

for ireg in regions:
    variables_common[ireg] = copy.deepcopy(variables_common["SR_1b"])


#variables_common["SR_2b"].append("csv1")
#variables_common["SR_2b"].append("jetpt1")
#variables_common["SR_2b"].append("jeteta1")

