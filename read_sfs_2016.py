from coffea.btag_tools import BTagScaleFactor
btag_sf = BTagScaleFactor("data/DeepCSV_2016LegacySF_V1.csv.gz", "medium")

## ele tight ID weights                                                                                                                                                                                 
from coffea.lookup_tools import extractor
ext = extractor()

ext.add_weight_sets(["EGamma_SF2D_T EGamma_SF2D data/2016LegacyReReco_ElectronTight_Fall17V2.root",
                     "EGamma_SF2D_T_err EGamma_SF2D_error data/2016LegacyReReco_ElectronTight_Fall17V2.root"])

ext.add_weight_sets(["EGamma_SF2D_L EGamma_SF2D data/2016LegacyReReco_ElectronLoose_Fall17V2.root",
                     "EGamma_SF2D_L_err EGamma_SF2D_error data/2016LegacyReReco_ElectronLoose_Fall17V2.root"])

ext.add_weight_sets(["EGamma_SF2D_Trig hEffEtaPt data/electron_Trigger_eleTrig.root",
                     "EGamma_SF2D_Trig_err hEffEtaPt_error data/electron_Trigger_eleTrig.root" ])
ext.add_weight_sets(["EGamma_SF2D_Reco EGamma_SF2D data/EGM2D_BtoH_GT20GeV_RecoSF_Legacy2016.root",
                     "EGamma_SF2D_Reco_err EGamma_SF2D_error data/EGM2D_BtoH_GT20GeV_RecoSF_Legacy2016.root"])

ext.add_weight_sets(["EGamma_SF2D_Reco_lowpt EGamma_SF2D data/EGM2D_BtoH_low_RecoSF_Legacy2016.root",
                     "EGamma_SF2D_Reco_lowpt_err EGamma_SF2D_error data/EGM2D_BtoH_low_RecoSF_Legacy2016.root"])


ext.add_weight_sets(["met_trig hden_monojet_recoil_clone_passed data/metTriggerEfficiency_zmm_recoil_monojet_TH1F.root"])



ext.add_weight_sets(["muon_lowpt_BCDEF_LooseID NUM_LooseID_DEN_genTracks_pt_abseta data/Muon_low-pT_RunBCDEF_SF_ID.root",
                     "muon_lowpt_BCDEF_TightID NUM_TightID_DEN_genTracks_pt_abseta data/Muon_low-pT_RunBCDEF_SF_ID.root",
                     "muon_lowpt_GH_LooseID NUM_LooseID_DEN_genTracks_pt_abseta data/Muon_low-pT_RunGH_SF_ID.root",
                     "muon_lowpt_GH_TightID NUM_TightID_DEN_genTracks_pt_abseta data/Muon_low-pT_RunGH_SF_ID.root",

                     "muon_highpt_BCDEF_LooseISO NUM_LooseRelIso_DEN_LooseID_eta_pt data/Muon_RunBCDEF_SF_ISO.root",
                     "muon_highpt_BCDEF_TightISO NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt data/Muon_RunBCDEF_SF_ISO.root",
                     "muon_highpt_GH_LooseISO NUM_LooseRelIso_DEN_LooseID_eta_pt data/Muon_RunGH_SF_ISO.root",
                     "muon_highpt_GH_TightISO NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt data/Muon_RunGH_SF_ISO.root",

                     "muon_highpt_BCDEF_LooseID NUM_LooseID_DEN_genTracks_eta_pt data/Muon_RunBCDEF_SF_ID.root",
                     "muon_highpt_BCDEF_TightID NUM_TightID_DEN_genTracks_eta_pt data/Muon_RunBCDEF_SF_ID.root",
                     "muon_highpt_GH_LooseID NUM_LooseID_DEN_genTracks_eta_pt data/Muon_RunGH_SF_ID.root",
                     "muon_highpt_GH_TightID NUM_TightID_DEN_genTracks_eta_pt data/Muon_RunGH_SF_ID.root"])


ext.add_weight_sets(["pu_weight puweight data/PU_Reweight_2016.root",
                     "pu_weight_Up puweight_Up data/PU_Reweight_2016.root",
                     "pu_weight_Down puweight_Down data/PU_Reweight_2016.root"])


ext.add_weight_sets(["btag_eff_lwp efficiency_btag_lwp data/bTagEffs_2016.root",
                     "btag_eff_mwp efficiency_btag_mwp data/bTagEffs_2016.root",
                     "ctag_eff_lwp efficiency_ctag_lwp data/bTagEffs_2016.root",
                     "ctag_eff_mwp efficiency_ctag_mwp data/bTagEffs_2016.root",
                     "ltag_eff_lwp efficiency_lighttag_lwp data/bTagEffs_2016.root",
                     "ltag_eff_mwp efficiency_lighttag_mwp data/bTagEffs_2016.root" ])

## following are same for 17 and 18                                                                                                                                                                     
ext.add_weight_sets(["k_dy_qcd_nlo kfac_dy_filter data/kfac_dy_filter.root",
                     "k_znn_qcd_nlo kfac_znn_filter data/kfac_znn_filter.root",
                     "k_wj_qcd_nlo wjet_dress_monojet data/2017_gen_v_pt_qcd_sf.root",
                     "k_dy_zj_ewk_nlo kfactor_monojet_ewk data/merged_kfactors_zjets.root",
                     "k_gj_ewk_nlo kfactor_monojet_ewk data/merged_kfactors_gjets.root",
                     "k_wj_ewk_nlo kfactor_monojet_ewk data/merged_kfactors_wjets.root"
                 ])
ext.finalize()
evaluator = ext.make_evaluator()
