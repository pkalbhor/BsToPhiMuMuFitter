#!/usr/bin/env python
# -*- coding: utf-8 -*-
# vim: set sw=4 ts=4 fdm=indent foldnestmax=3 ft=python et:

import re
import types
import functools
import itertools
from array import array
from copy import copy
import math, time

import BsToPhiMuMuFitter.cpp

from v2Fitter.Fitter.DataReader import DataReader
from v2Fitter.Fitter.ObjProvider import ObjProvider
from BsToPhiMuMuFitter.varCollection import dataArgs, Bmass, CosThetaL, CosThetaK, Phi, Phimass, dataArgsGEN
from BsToPhiMuMuFitter.anaSetup import q2bins, bMassRegions, cuts, cuts_noResVeto,  modulePath, dataFilePath, sigMC, UnfilteredMC

import ROOT
from ROOT import TChain
from ROOT import TEfficiency, TH2D, TH3D, TCanvas
from ROOT import RooArgList
from ROOT import RooDataHist

from BsToPhiMuMuFitter.StdProcess import p

CFG = DataReader.templateConfig()
CFG.update({
    'argset': dataArgs,
    'lumi': -1,  # Keep a record, useful for mixing simulations samples
    #'ifriendIndex': ["Bmass", "Mumumass"],
})

# dataReader
def customizeOne(self, targetBMassRegion=None, extraCuts=None):
    """Define datasets with arguments."""
    if targetBMassRegion is None:
        targetBMassRegion = []
    if not self.process.cfg['binKey'] in q2bins.keys():
        self.logger.logERROR("Bin {0} is not defined.\n".format(self.process.cfg['binKey']))
        raise ValueError

    # With shallow copied CFG, have to bind cfg['dataset'] to a new object.
    self.cfg['dataset'] = []
    for key, val in bMassRegions.items():
        if any([re.match(pat, key) for pat in targetBMassRegion]):
            self.cfg['dataset'].append(
                (
                    "{0}.{1}".format(self.cfg['name'], key),
                    "({0}) && ({1}) && ({2}) && ({3})".format(
                        val['cutString'],
                        q2bins[self.process.cfg['binKey']]['cutString'],
                        cuts[-1] if self.process.cfg['binKey'] not in ['jpsi', 'psi2s'] else cuts_noResVeto,
                        "1" if not extraCuts else extraCuts,
                    )
                )
            )

    # Customize preload TFile
    if self.cfg['preloadFile']:
        self.cfg['preloadFile'] = self.cfg['preloadFile'].format(binLabel=q2bins[self.process.cfg['binKey']]['label'])

dataReaderCfg = copy(CFG)
dataReaderCfg.update({
    'name': "dataReader",
    'ifile':  [dataFilePath], #["/afs/cern.ch/work/p/pkalbhor/BFitter/BPhysicsData/data/Modified_BDT_Modified_sel_BsToPhiMuMu_2016_combine_data_cut0_s0.root"],
    #'ifriend': ["/afs/cern.ch/work/p/pchen/public/BuToKstarMuMu/v2Fitter/BsToPhiMuMuFitter/script/plotMatchCandPreSelector.root"],
    'preloadFile': modulePath + "/data/preload_dataReader_{binLabel}.root",
    'lumi': 19.98,
})
dataReader = DataReader(dataReaderCfg)
customizeData = functools.partial(customizeOne, targetBMassRegion=['^Fit$', '^SR$', '^.{0,1}SB$'])
dataReader.customize = types.MethodType(customizeData, dataReader)

# sigMCReader
sigMCReaderCfg = copy(CFG)
sigMCReaderCfg.update({
    'name': "sigMCReader",
    'ifile': [sigMC], #["/afs/cern.ch/work/p/pkalbhor/BFitter/BPhysicsData/data/Modified_BDT_Modified_sel_BsToPhiMuMu_combine_MC_2016_mc.lite_cut0.root"],
    'preloadFile': modulePath + "/data/preload_sigMCReader_{binLabel}.root",
    'lumi': 16281.440 + 21097.189,
})
sigMCReader = DataReader(sigMCReaderCfg)
customizeSigMC = functools.partial(customizeOne, targetBMassRegion=['^Fit$'])  # Assuming cut_kshortWindow makes no impact
sigMCReader.customize = types.MethodType(customizeSigMC, sigMCReader)

# sigMCGENReader
def customizeGEN(self):
    """Define datasets with arguments."""
    if not self.process.cfg['binKey'] in q2bins.keys():
        print("ERROR\t: Bin {0} is not defined.\n".format(self.process.cfg['binKey']))
        raise AttributeError

    # With shallow copied CFG, have to bind cfg['dataset'] to a new object.
    self.cfg['dataset'] = []
    self.cfg['dataset'].append(
        (
            "{0}.Fit".format(self.cfg['name']),
            re.sub("Mumumass", "sqrt(genQ2)", q2bins[self.process.cfg['binKey']]['cutString'])
        )
    )

    # Customize preload TFile
    if self.cfg['preloadFile']:
        self.cfg['preloadFile'] = self.cfg['preloadFile'].format(binLabel=q2bins[self.process.cfg['binKey']]['label'])

sigMCGENReaderCfg = copy(CFG)
sigMCGENReaderCfg.update({
    'name': "sigMCGENReader",
    #'ifile': ["/eos/cms/store/user/pchen/BToKstarMuMu/dat/sel/v3p5/unfilteredSIG_genonly/sel_*.root"],
    'ifile': [UnfilteredMC], #["/afs/cern.ch/work/p/pkalbhor/BFitter/BPhysicsData/data/Modified_sel_BsToPhiMuMu_NofilterMC_signal_2016_mc.lite_cut0_s.root"],
    'preloadFile': modulePath + "/data/preload_sigMCGENReader_{binLabel}.root",
    'argset': dataArgsGEN,
})
sigMCGENReader = DataReader(sigMCGENReaderCfg)
sigMCGENReader.customize = types.MethodType(customizeGEN, sigMCGENReader)

# effiHistReader
accXEffThetaLBins = array('d', [-1, -0.7, -0.3, 0., 0.3, 0.7, 1.])
accXEffThetaKBins = array('d', [-1, -0.7, 0., 0.4, 0.8, 1.])
accXEffPhiBins = array('d', [-3.15, -2.0, -1.0, 0.0, 1.0, 2.0, 3.15])

def buildAccXRecEffiHist(self):
    """Build efficiency histogram for later fitting/plotting"""
    print("Now I am Here in buildAccXRecEffiHist")
    fin = self.process.filemanager.open("buildAccXRecEffiHist", modulePath + "/data/accXrecEffHists_Run2016.root", "UPDATE")

    # Build acceptance, reco efficiency, and accXrec
    forceRebuild = False
    for binKey in q2bins.keys():
        print("Q2Bins: ", q2bins.keys())#Pritam
        if binKey in ['jpsi', 'psi2s', 'peaks']:
            continue
        h3_accXrec = fin.Get("h3_accXrec_{0}".format(binKey))
        if h3_accXrec == None or forceRebuild:
            h3_acc = fin.Get("h3_acc_{0}".format(binKey))
            h3_rec = fin.Get("h3_rec_{0}".format(binKey))

            # Fill histograms
            setupEfficiencyBuildProcedure = {}
            setupEfficiencyBuildProcedure['acc'] = {
                'ifiles': [UnfilteredMC], # ["/afs/cern.ch/work/p/pkalbhor/BFitter/BPhysicsData/data/Modified_sel_BsToPhiMuMu_NofilterMC_signal_2016_mc.lite_cut0_s.root"],# sigMCGENReaderCfg.cfg['ifile'],
                'baseString': re.sub("Mumumass", "sqrt(genQ2)", q2bins[binKey]['cutString']),
                'cutString': "fabs(genMupEta)<2.5 && fabs(genMumEta)<2.5 && genMupPt>2.5 && genMumPt>2.5",
                'fillXYZ': "genPhi:genCosThetaK:genCosThetaL"  # Z:Y:X
            }
            setupEfficiencyBuildProcedure['rec'] = {
                'ifiles': sigMCReader.cfg['ifile'],
                'baseString': "{0}".format(setupEfficiencyBuildProcedure['acc']['baseString']),
                'cutString': "Bmass > 0.5 && ({0})".format(cuts[-1]),
                'fillXYZ': "genPhi:genCosThetaK:genCosThetaL"  # Z:Y:X
            }
            # print(setupEfficiencyBuildProcedure['rec'])
            for h3, label in (h3_acc, 'acc'), (h3_rec, 'rec'):
                if h3 == None or forceRebuild:
                    print("h3: ", h3, label)
                    treein = TChain("events")
                    for f in setupEfficiencyBuildProcedure[label]['ifiles']:
                        treein.Add(f)
                        print(setupEfficiencyBuildProcedure[label]['ifiles'])

                    treein.Draw(">>totEvtList", setupEfficiencyBuildProcedure[label]['baseString'])
                    totEvtList = ROOT.gDirectory.Get("totEvtList")
                    totEvtList.Print() 
                    treein.SetEventList(totEvtList)
                    
                    treein.Draw(">>accEvtList", setupEfficiencyBuildProcedure[label]['cutString'])
                    accEvtList = ROOT.gDirectory.Get("accEvtList")
                    accEvtList.Print()

                    h3_total = TH3D("h3_{0}_{1}_total".format(label, binKey), "", len(accXEffThetaLBins) - 1, accXEffThetaLBins, len(accXEffThetaKBins) - 1, accXEffThetaKBins, len(accXEffPhiBins)-1, accXEffPhiBins)
                    h3_passed = h3_total.Clone("h3_{0}_{1}_passed".format(label, binKey))

                    h3_fine_total = TH3D("h3_{0}_fine_{1}_total".format(label, binKey), "", 20, -1, 1, 20, -1, 1, 20, -3.15, 3.15)
                    h3_fine_passed = h3_fine_total.Clone("h3_{0}_fine_{1}_passed".format(label, binKey))

                    treein.SetEventList(totEvtList)
                    for hist in h3_total, h3_fine_total:
                        treein.Draw("{0}>>{1}".format(setupEfficiencyBuildProcedure[label]['fillXYZ'], hist.GetName()), "", "goff")

                    treein.SetEventList(accEvtList)
                    for hist in h3_passed, h3_fine_passed:
                        treein.Draw("{0}>>{1}".format(setupEfficiencyBuildProcedure[label]['fillXYZ'], hist.GetName()), "", "goff")

                    h3_eff = TEfficiency(h3_passed, h3_total)
                    h3_eff_fine = TEfficiency(h3_fine_passed, h3_fine_total)

                    fin.cd()
                    for proj, var in [("ProjectionX", CosThetaL), ("ProjectionY", CosThetaK), ("ProjectionZ", Phi)]:
                        proj_fine_total = getattr(h3_fine_total, proj)("{0}_{1}".format(h3_fine_total.GetName(), proj), 0, -1, 0, -1, "e")
                        proj_fine_passed = getattr(h3_fine_passed, proj)("{0}_{1}".format(h3_fine_passed.GetName(), proj), 0, -1, 0, -1, "e")
                        h_eff = TEfficiency(proj_fine_passed, proj_fine_total)
                        h_eff.Write("h_{0}_fine_{1}_{2}".format(label, binKey, proj), ROOT.TObject.kOverwrite)

                    h3_eff.Write("h3_{0}_{1}".format(label, binKey), ROOT.TObject.kOverwrite)
                    h3_eff_fine.Write("h3_{0}_fine_{1}".format(label, binKey), ROOT.TObject.kOverwrite)
            
            # Merge acc and rec to accXrec
            fin.cd()
            for proj in ["ProjectionX", "ProjectionY", "ProjectionZ"]:
                h_acc_fine = fin.Get("h_acc_fine_{0}_{1}".format(binKey, proj))
                h_rec_fine = fin.Get("h_rec_fine_{0}_{1}".format(binKey, proj))
                h_accXrec_fine = h_acc_fine.GetPassedHistogram().Clone("h_accXrec_fine_{0}_{1}".format(binKey, proj))
                h_accXrec_fine.Reset("ICESM")
                for b in range(1, h_accXrec_fine.GetNbinsX() + 1):
                    h_accXrec_fine.SetBinContent(b, h_acc_fine.GetEfficiency(b) * h_rec_fine.GetEfficiency(b))
                    # print("news: ", h_accXrec_fine.GetBinContent(b))
                    h_accXrec_fine.SetBinError(b, h_accXrec_fine.GetBinContent(b) * math.sqrt(1 / h_acc_fine.GetTotalHistogram().GetBinContent(b) + 1 / h_acc_fine.GetPassedHistogram().GetBinContent(b) + 1 / h_rec_fine.GetTotalHistogram().GetBinContent(b) + 1 / h_rec_fine.GetPassedHistogram().GetBinContent(b)))
                h_accXrec_fine.Write("h_accXrec_{0}_{1}".format(binKey, proj), ROOT.TObject.kOverwrite)

            h3_acc = fin.Get("h3_acc_{0}".format(binKey)) #3D Effificiency
            h3_rec = fin.Get("h3_rec_{0}".format(binKey))
            h3_accXrec = h3_acc.GetPassedHistogram().Clone("h3_accXrec_{0}".format(binKey))
            h3_accXrec.Reset("ICESM")
            for iL, iK, iP in itertools.product(range(1, len(accXEffThetaLBins)), range(1, len(accXEffThetaKBins)), range(1, len(accXEffPhiBins))):
                if h3_acc.GetTotalHistogram().GetBinContent(iL, iK, iP) == 0 or h3_acc.GetPassedHistogram().GetBinContent(iL, iK, iP) == 0 or h3_rec.GetTotalHistogram().GetBinContent(iL, iK, iP) == 0 or h3_rec.GetPassedHistogram().GetBinContent(iL, iK, iP) == 0:
                    h3_accXrec.SetBinContent(iL, iK, iP, 0)
                    h3_accXrec.SetBinError(iL, iK, iP, 1)
                else:
                    iLKP = h3_acc.GetGlobalBin(iL, iK, iP)
                    h3_accXrec.SetBinContent(iL, iK, iP, h3_acc.GetEfficiency(iLKP) * h3_rec.GetEfficiency(iLKP))
                    #print("Errors: ", h3_acc.GetTotalHistogram().GetBinContent(iLKP), h3_acc.GetPassedHistogram().GetBinContent(iLKP), h3_rec.GetTotalHistogram().GetBinContent(iLKP), h3_rec.GetPassedHistogram().GetBinContent(iLKP) )
                    h3_accXrec.SetBinError(iL, iK, iP, h3_accXrec.GetBinContent(iL, iK, iP) * math.sqrt(1 / h3_acc.GetTotalHistogram().GetBinContent(iLKP) + 1 / h3_acc.GetPassedHistogram().GetBinContent(iLKP) + 1 / h3_rec.GetTotalHistogram().GetBinContent(iLKP) + 1 / h3_rec.GetPassedHistogram().GetBinContent(iLKP)))
                    # print("live: ", h3_accXrec.GetBinContent(iL, iK, iP))
            h3_accXrec.SetXTitle("cos#theta_{l}")
            h3_accXrec.SetYTitle("cos#theta_{K}")
            h3_accXrec.SetZTitle("#phi")

            h3_accXrec.Write("h3_accXrec_{0}".format(binKey), ROOT.TObject.kOverwrite)
            self.logger.logINFO("Overall efficiency is built.")

    # Register the chosen one to sourcemanager
    h3_accXrec = fin.Get("h3_accXrec_{0}".format(self.process.cfg['binKey']))
    self.cfg['source']['effiHistReader.h3_accXrec'] = h3_accXrec
    self.cfg['source']['effiHistReader.accXrec'] = RooDataHist("accXrec", "", RooArgList(CosThetaL, CosThetaK, Phi), ROOT.RooFit.Import(h3_accXrec))
    self.cfg['source']['effiHistReader.h_accXrec_fine_ProjectionX'] = fin.Get("h_accXrec_{0}_ProjectionX".format(self.process.cfg['binKey']))
    self.cfg['source']['effiHistReader.h_accXrec_fine_ProjectionY'] = fin.Get("h_accXrec_{0}_ProjectionY".format(self.process.cfg['binKey']))

effiHistReader = ObjProvider({
    'name': "effiHistReader",
    'obj': {
        'effiHistReader.h3_accXrec': [buildAccXRecEffiHist, ],
    }
})

if __name__ == '__main__':
    #  p.setSequence([dataReader])
    #  p.setSequence([sigMCReader])
    p.setSequence([effiHistReader])
    p.beginSeq()
    p.runSeq()
    p.endSeq()
