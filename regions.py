import awkward1 as ak

def get_mask_wmunu1b(cms_events):
    mask_wmunu1b = ( (cms_events.mettrig ) &
                     (cms_events.nEleLoose==0 ) &
                     #(cms_events.npho==0) &
                     (cms_events.metpt>100.) &
                     (cms_events.ncleanpho==0) &
                     (cms_events.ntau==0) &
                     (cms_events.nMuLoose==1) &
                     (cms_events.nMuTight==1) &
                     (cms_events.delta_met_topmu < 0.5 ) &
                     (cms_events.recoil_Wmunu0>250.) &
                     (cms_events.min_dphi_jet_met > 0.5) &
                     (cms_events.nJetLoose==1 ) &
                     (cms_events.nJetTight==1 ) &
                     (cms_events.nJetb ==1 ) &
                     (cms_events.mt_Wmunu0 >=0 ) & (cms_events.mt_Wmunu0 < 160 )
                 )
    
    return mask_wmunu1b
    
def get_mask_wenu1b (cms_events):
    mask_wenu1b = ( (cms_events.eletrig) &
                    (cms_events.nEleTight==1 ) &
                    (cms_events.nEleLoose==1 ) &
                    (cms_events.metpt>100.) &
                    (cms_events.ncleanpho==0) &
                    (cms_events.ntau==0) &
                    (cms_events.nMuLoose==0) &
                    (cms_events.recoil_Wenu0>250.) &
                    (cms_events.delta_met_tope < 0.5 ) &
                    (cms_events.min_dphi_jet_met > 0.5) &
                    (cms_events.nJetLoose==1 ) &
                    (cms_events.nJetTight==1 ) &
                    (cms_events.nJetb ==1 ) &
                    (cms_events.mt_Wenu0 >=0 ) & (cms_events.mt_Wenu0 < 160 )
                )

    return mask_wenu1b

def get_mask_topmunu2b (cms_events):
    mask_topmunu2b = ( (cms_events.mettrig) &
                       (cms_events.nEleLoose==0 ) &
                       (cms_events.ncleanpho==0) &
                       #(cms_events.npho==0) &
                       (cms_events.ntau==0) &
                       (cms_events.nMuLoose==1) &
                       (cms_events.nMuTight==1) &
                       (cms_events.delta_met_topmu < 0.5 ) &
                       (cms_events.recoil_Wmunu0>250.) &
                       (cms_events.mt_Wmunu0 >=0 ) & (cms_events.mt_Wmunu0 < 160 )&
                       (cms_events.min_dphi_jet_met > 0.5) &
                       ((cms_events.nJetLoose>=3 ) & (cms_events.nJetLoose<=6 )) &
                       (cms_events.nJetTight>=1 ) &
                       ((cms_events.nJetb ==2 ) & (cms_events.nJetLoose<=600 )) &
                       (cms_events.jet_sel_tight0==True) &
                       (cms_events.jet_sel_b_0 == True) &
                       (cms_events.jet_sel_b_1 == True) 
    )

    return mask_topmunu2b

def get_mask_topenu2b (cms_events):
    mask_topenu2b = ( (cms_events.eletrig) &
                      (cms_events.eletrig)& 
                      (cms_events.nMuLoose==0) &
                      (cms_events.ntau==0) &
                      (cms_events.nEleTight==1 ) &
                      (cms_events.nEleLoose==1 ) &
                      (cms_events.ncleanpho==0) &
                      #(cms_events.npho==0) &
                      (cms_events.min_dphi_jet_met > 0.5) &
                      (cms_events.delta_met_tope < 0.5 ) & 
                      (cms_events.mt_Wenu0 >=0 ) & (cms_events.mt_Wenu0 < 160 ) &
                      (cms_events.recoil_Wenu0>250.) &
                      ((cms_events.nJetLoose>=3 ) & (cms_events.nJetLoose<=100 )) &
                      (cms_events.nJetTight>=1 ) &
                      (cms_events.nJetb ==2 ) &
                      (cms_events.jet_sel_tight0==True) &
                      (cms_events.jet_sel_b_0 == True) &
                      (cms_events.jet_sel_b_1 == True) 
    )

    return mask_topenu2b

def get_mask_Zmumu1b (cms_events):
    mask_Zmumu1b = ( (cms_events.mettrig) &
                     (cms_events.nEleLoose==0 ) &
                     (cms_events.ncleanpho==0) &
                     #(cms_events.npho==0) &
                     (cms_events.ntau==0) &
                     (cms_events.nMuLoose==2) &
                     (cms_events.nMuTight>=1) &
                     (cms_events.mu_sel_tight0==True) &
                     (cms_events.delta_met_zmumu < 0.5) &
                     (cms_events.Zmumu_recoil>250.) &
                     (cms_events.min_dphi_jet_met > 0.5) &
                     ((cms_events.nJetLoose>=1 ) &  (cms_events.nJetLoose<=2 )) &
                     (cms_events.nJetTight>=1 ) &
                     #(cms_events.nJetb ==1 ) &
                     (cms_events.jet_sel_tight0 == True) &
                     #(cms_events.jet_sel_b_0 == True) &
                     ( (cms_events.mu_q0 * cms_events.mu_q1 ) < 0) &
                     (cms_events.Zmumu_mass >=70 ) & (cms_events.Zmumu_mass <= 110 )
                 )

    return mask_Zmumu1b

def get_mask_Zmumu2b (cms_events):
    mask_Zmumu2b = ( (cms_events.mettrig) &
                     (cms_events.nEleLoose==0 ) &
                     (cms_events.ncleanpho==0) &
                     #(cms_events.npho==0) &
                     (cms_events.ntau==0) &
                     (cms_events.nMuLoose==2) &
                     (cms_events.nMuTight>=1) &
                     (cms_events.mu_sel_tight0==True) &
                     (cms_events.Zmumu_recoil>250.) &
                     (cms_events.delta_met_zmumu < 0.5) &
                     (cms_events.min_dphi_jet_met > 0.5) &
                     ((cms_events.nJetLoose>=2 ) &  (cms_events.nJetLoose<=3 )) &
                     (cms_events.nJetTight>=1 ) &
                     #(cms_events.nJetb ==2 ) &
                     (cms_events.jet_sel_tight0==True) &
                     # (cms_events.jet_sel_b_0 == True) &
                     # (cms_events.jet_sel_b_1 == True) &
                     ( (cms_events.mu_q0 * cms_events.mu_q1 ) < 0) &
                     (cms_events.Zmumu_mass >=70 ) & (cms_events.Zmumu_mass <= 110 )
                 )

    return mask_Zmumu2b

def get_mask_Zee1b (cms_events):
    mask_Zee1b = ( (cms_events.eletrig) &
                   (cms_events.nEleLoose==2 ) &
                   (cms_events.nEleTight>=1 ) &
                   (cms_events.ele_sel_tight0==True)&
                   (cms_events.ncleanpho==0) &
                   #(cms_events.npho==0) &
                   (cms_events.ntau==0) &
                   (cms_events.nMuLoose==0) &
                   (cms_events.Zee_recoil>250.) &
                   (cms_events.min_dphi_jet_met > 0.5) &
                   (cms_events.delta_met_zee < 0.5) &
                   ((cms_events.nJetLoose>=1 ) &  (cms_events.nJetLoose<=2 )) &
                   (cms_events.nJetTight>=1 ) &
                   #(cms_events.nJetb ==1 ) &
                   (cms_events.jet_sel_tight0 == True) &
                   #(cms_events.jet_sel_b_0 == True) &
                   ( (cms_events.ele_q0 * cms_events.ele_q1 ) < 0) &
                   (cms_events.Zee_mass >=70 ) & (cms_events.Zee_mass <= 110 )
                 )

    return mask_Zee1b

def get_mask_Zee2b (cms_events):
    mask_Zee2b = ( (cms_events.eletrig) &
                   (cms_events.nEleLoose==2 ) &
                   (cms_events.nEleTight>=1 ) &
                   (cms_events.ele_sel_tight0==True) &
                   (cms_events.ncleanpho==0) &
                   #(cms_events.npho==0) &
                   (cms_events.ntau==0) &
                   (cms_events.nMuLoose==0) &
                   (cms_events.Zee_recoil>250.) &
                   (cms_events.min_dphi_jet_met > 0.5) &
                   (cms_events.delta_met_zee < 0.5) &
                   ((cms_events.nJetLoose>=2 ) &  (cms_events.nJetLoose<=3 )) &
                   (cms_events.nJetTight>=1 ) &
                   #(cms_events.nJetb ==2 ) &
                   (cms_events.jet_sel_tight0==True) &
                   #(cms_events.jet_sel_b_0 == True) &
                   #(cms_events.jet_sel_b_1 == True) &
                   ( (cms_events.ele_q0 * cms_events.ele_q1 ) < 0) &
                   (cms_events.Zee_mass >=70 ) & (cms_events.Zee_mass <= 110 )
               )

    return mask_Zee2b

def get_mask_SR1b (cms_events):
    mask_SR1b =  ( (cms_events.mettrig) &
                   (cms_events.metpt > 250) &
                   (cms_events.nEleLoose==0 ) &
                   (cms_events.ncleanpho==0) &
                   (cms_events.ntau==0) &
                   (cms_events.nMuLoose==0) &
                   ((cms_events.min_dphi_jet_met) > 0.5) &
                   (cms_events.delta_met_sr < 0.5 ) &
                   ((cms_events.nJetLoose>=1 ) &  (cms_events.nJetLoose<=2 )) &
                   (cms_events.nJetTight>=1 ) &
                   (cms_events.nJetb ==1 ) &
                   (cms_events.jet_sel_tight0 == True) &
                   (cms_events.jet_sel_b_0 == True)
               )

    return mask_SR1b

def get_mask_SR2b (cms_events):
    mask_SR2b =  ( (cms_events.mettrig) & 
                   (cms_events.metpt > 250) &
                   (cms_events.nEleLoose==0 ) &
                   #(cms_events.npho==0) &
                   (cms_events.ncleanpho==0) &
                   (cms_events.ntau==0) &
                   (cms_events.nMuLoose==0) &
                   (cms_events.delta_met_sr < 0.5 ) &
                   (cms_events.min_dphi_jet_met > 0.5) &
                   ((cms_events.nJetLoose>=2 ) &  (cms_events.nJetLoose<=3 )) &
                   (cms_events.nJetTight>=1 ) &
                   (cms_events.nJetb ==2 ) &
                   (cms_events.jet_sel_tight0 == True) &
                   (cms_events.jet_sel_b_0 == True) &
                   (cms_events.jet_sel_b_1 == True)

               )

    return mask_SR2b 





def get_mask_wmunu2b (cms_events):
    mask_wmunu2b = ( (cms_events.nEleLoose==0 ) &
                     (cms_events.npho==0) &
                     (cms_events.ntau==0) &
                     (cms_events.nMuLoose==1) &
                     (cms_events.nMuTight==1) &
                     (cms_events.recoil_Wmunu0>250.) &
                     (cms_events.min_dphi_jet_met > 0.5) &
                     (cms_events.nJetLoose==2 ) &
                     (cms_events.nJetTight>=1 ) &
                     (cms_events.nJetb ==2 ) &
                     (cms_events.mt_Wmunu0 >=0 ) & (cms_events.mt_Wmunu0 < 160 )
                 )
    
    return mask_wmunu2b


def get_mask_wenu2b (cms_events):
    mask_wenu2b = ( (cms_events.nEleTight==1 ) &
                    (cms_events.nEleLoose==1 ) &
                    (cms_events.npho==0) &
                    (cms_events.ntau==0) &
                    (cms_events.nMuLoose==0) &
                    (cms_events.recoil_Wenu0>250.) &
                    (cms_events.min_dphi_jet_met > 0.5) &
                    (cms_events.nJetLoose==2 ) &
                    (cms_events.nJetTight>=1 ) &
                    (cms_events.nJetb ==2 ) &
                    (cms_events.mt_Wenu0 >=0 ) & (cms_events.mt_Wenu0 < 160 )
                )

    return mask_wenu2b

def get_mask_topmunu1b (cms_events):
    mask_topmunu1b = ( (cms_events.nEleLoose==0 ) &
                       (cms_events.npho==0) &
                       (cms_events.ntau==0) &
                       (cms_events.nMuLoose==1) &
                       (cms_events.nMuTight==1) &
                       (cms_events.recoil_Wmunu0>250.) &
                       (cms_events.min_dphi_jet_met > 0.5) &
                       ((cms_events.nJetLoose>=2 ) & (cms_events.nJetLoose<=600 )) &
                       (cms_events.nJetTight>=1 ) &
                       (cms_events.nJetb ==1 ) &
                       (cms_events.jet_sel_tight0 == True) &
                       (cms_events.jet_sel_b_0 == True) &
                       (cms_events.mt_Wmunu0 >0 ) & (cms_events.mt_Wmunu0 < 160 )
                   )
    
    return mask_topmunu1b

def get_mask_topenu1b (cms_events):
    mask_topenu1b = ( (cms_events.nEleTight==1 ) &
                      (cms_events.nEleLoose==1 ) &
                      (cms_events.npho==0) &
                      (cms_events.ntau==0) &
                      (cms_events.nMuLoose==0) &
                      (cms_events.recoil_Wenu0>250.) &
                      (cms_events.min_dphi_jet_met > 0.5) &
                      ((cms_events.nJetLoose>=2 ) & (cms_events.nJetLoose<=60 )) &
                      (cms_events.nJetTight>=1 ) &
                      (cms_events.nJetb ==1 ) &
                      (cms_events.jet_sel_tight0 == True) &
                      (cms_events.jet_sel_b_0 == True) &
                      (cms_events.mt_Wenu0 >0 ) & (cms_events.mt_Wenu0 < 160 )
                  )

    return mask_topenu1b

