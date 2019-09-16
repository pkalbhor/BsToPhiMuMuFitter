#!/usr/bin/env python
# -*- coding: utf-8 -*-
# vim: set sw=4 ts=4 fdm=indent fdl=2 ft=python et:

from v2Fitter.Fitter.FitterCore import FitterCore

import BsToPhiMuMuFitter.cpp
from BsToPhiMuMuFitter.anaSetup import q2bins
from BsToPhiMuMuFitter.StdProcess import setStyle
from BsToPhiMuMuFitter.varCollection import CosThetaL, CosThetaK, Phi
from BsToPhiMuMuFitter.FitDBPlayer import FitDBPlayer

import re
import itertools

import ROOT

class EfficiencyFitter(FitterCore):
    """Implementation to standard efficiency fitting procdeure to BuToKstarMuMu angular analysis"""

    @classmethod
    def templateConfig(cls):
        cfg = FitterCore.templateConfig()
        cfg.update({
            'name': "EfficiencyFitter",
            'data': "effiHistReader.accXrec",
            'dataX': "effiHistReader.h_accXrec_fine_ProjectionX",
            'dataY': "effiHistReader.h_accXrec_fine_ProjectionY",
            'dataZ': "effiHistReader.h_accXrec_fine_ProjectionZ",
            'pdf': "effi_sigA",
            'pdfX': "effi_cosl",
            'pdfY': "effi_cosK",
            'pdfZ': "effi_Phi",
            'updateArgs': True,
        })
        del cfg['createNLLOpt']
        return cfg

    def _bookMinimizer(self):
        print("""Pass complicate fitting control.""")
        pass

    def _preFitSteps(self):
        print("""Prefit uncorrelated term""")
        args = self.pdf.getParameters(self.data)
        FitDBPlayer.initFromDB(self.process.dbplayer.odbfile, args)
        self.ToggleConstVar(args, isConst=True)

        # Disable xTerm correction and fit to 1-D
        args.find('hasXTerm').setVal(0)
        print("args in _preFitSteps: ", args.Print())

        h_accXrec_fine_ProjectionX = self.process.sourcemanager.get(self.cfg['dataX'])
        h_accXrec_fine_ProjectionY = self.process.sourcemanager.get(self.cfg['dataY'])
        h_accXrec_fine_ProjectionZ = self.process.sourcemanager.get(self.cfg['dataZ'])
        effi_cosl = self.process.sourcemanager.get(self.cfg['pdfX'])
        effi_cosK = self.process.sourcemanager.get(self.cfg['pdfY'])
        effi_Phi = self.process.sourcemanager.get(self.cfg['pdfZ'])
        for proj, pdf, var, argPats in [(h_accXrec_fine_ProjectionX, effi_cosl, CosThetaL, [r"^l\d+$"]), (h_accXrec_fine_ProjectionY, effi_cosK, CosThetaK, [r"^k\d+$"]), (h_accXrec_fine_ProjectionZ, effi_Phi, Phi, [r"^p\d+$"])]:
            hdata = ROOT.RooDataHist("hdata", "", ROOT.RooArgList(var), ROOT.RooFit.Import(proj))
            self.ToggleConstVar(args, isConst=False, targetArgs=argPats)
            pdf.chi2FitTo(hdata, ROOT.RooLinkedList())
            self.ToggleConstVar(args, isConst=True, targetArgs=argPats)

        args.find('effi_norm').setConstant(False)
        self.pdf.chi2FitTo(self.data, ROOT.RooFit.Minos(True))
        args.find('effi_norm').setVal(args.find('effi_norm').getVal() / 4.)
        args.find('effi_norm').setConstant(True)

        # Fix uncorrelated term and for later update with xTerms in main fit step
        args.find('hasXTerm').setVal(0)
        self.ToggleConstVar(args, isConst=False, targetArgs=[r"^x\d+$"])

    def _postFitSteps(self):
        print("""Post-processing""")
        args = self.pdf.getParameters(self.data)
        self.ToggleConstVar(args, True)
        FitDBPlayer.UpdateToDB(self.process.dbplayer.odbfile, args)

    def _runFitSteps(self):
        h3_accXrec = self.process.sourcemanager.get("effiHistReader.h3_accXrec")

        effi_sigA_formula = self.pdf.formula().GetExpFormula().Data()
        args = self.pdf.getParameters(self.data)
        args_it = args.createIterator()
        arg = args_it.Next()
        nPar = 0
        print("args in _runFitSteps: ", args.Print())
        print("Effi Formula1: ", effi_sigA_formula)
        while arg:
            if any(re.match(pat, arg.GetName()) for pat in ["effi_norm", "hasXTerm", r"^l\d+$", r"^k\d+$", r"^p\d+$"]):
                effi_sigA_formula = re.sub(arg.GetName(), "({0})".format(arg.getVal()), effi_sigA_formula)
            elif re.match(r"^x\d+$", arg.GetName()):
                nPar = nPar + 1
            arg = args_it.Next()
        effi_sigA_formula = re.sub(r"x(\d{1,2})", r"[\1]", effi_sigA_formula)
        effi_sigA_formula = re.sub(r"CosThetaL", r"x", effi_sigA_formula)
        effi_sigA_formula = re.sub(r"CosThetaK", r"y", effi_sigA_formula)
        effi_sigA_formula = re.sub(r"Phi", r"z", effi_sigA_formula)
        f3_effi_sigA = ROOT.TF3("f3_effi_sigA", effi_sigA_formula, -1, 1, -1, 1, -3.15, 3.15)

        print("Effi Formula: ", effi_sigA_formula)
        fitter = ROOT.EfficiencyFitter()
        minuit = fitter.Init(nPar, h3_accXrec, f3_effi_sigA)
        for xIdx in range(nPar):
            minuit.DefineParameter(xIdx, "x{0}".format(xIdx), 0., 1E-4, -1E+1, 1E+1)
        minuit.Command("MINI")
        minuit.Command("MINI")
        minuit.Command("MINOS")

        parVal = ROOT.Double(0)
        parErr = ROOT.Double(0)
        for xIdx in range(nPar):
            minuit.GetParameter(xIdx, parVal, parErr)
            arg = args.find("x{0}".format(xIdx))
            arg.setVal(parVal)
            arg.setError(parErr)

        # Check if efficiency is positive definite
        f3_max_x, f3_max_y, f3_max_z = ROOT.Double(0), ROOT.Double(0), ROOT.Double(0)
        f3_min_x, f3_min_y, f3_min_z = ROOT.Double(0), ROOT.Double(0), ROOT.Double(0)
        f3_effi_sigA.GetMaximumXYZ(f3_max_x, f3_max_y, f3_max_z)
        f3_effi_sigA.GetMinimumXYZ(f3_min_x, f3_min_y, f3_min_z)
        self.logger.logINFO("Sanitary check: Efficiency ranges from {0:.2e} to {1:.2e}".format(f3_effi_sigA.Eval(f3_min_x, f3_min_y, f3_min_z), f3_effi_sigA.Eval(f3_max_x, f3_max_y, f3_max_z)))

        # Plot comparison between fitting result to data
        setStyle()
        canvas = ROOT.TCanvas()
        latex = ROOT.TLatex()
        h3_effi_3D_comp = h3_accXrec.Clone("h3_effi_3D_comp")
        h3_effi_3D_comp.Reset("ICESM")
        for lBin, KBin, PBin in itertools.product(list(range(1, h3_effi_3D_comp.GetNbinsX() + 1)), list(range(1, h3_effi_3D_comp.GetNbinsY() + 1)), list(range(1, h3_effi_3D_comp.GetNbinsZ() + 1))):
            if (h3_accXrec.GetBinContent(lBin, KBin, PBin) != 0.0):
                h3_effi_3D_comp.SetBinContent(lBin, KBin, PBin, f3_effi_sigA.Eval(h3_accXrec.GetXaxis().GetBinCenter(lBin), h3_accXrec.GetYaxis().GetBinCenter(KBin), h3_accXrec.GetZaxis().GetBinCenter(PBin)) / h3_accXrec.GetBinContent(lBin, KBin, PBin))
            else:
                h3_effi_3D_comp.SetBinContent(lBin, KBin, PBin, 0.0)
            print(f3_effi_sigA.Eval(h3_accXrec.GetXaxis().GetBinCenter(lBin), h3_accXrec.GetYaxis().GetBinCenter(KBin), h3_accXrec.GetZaxis().GetBinCenter(PBin)), h3_accXrec.GetBinContent(lBin, KBin, PBin))
        h3_effi_3D_comp.SetMinimum(0)
        h3_effi_3D_comp.SetMaximum(1.5)
        h3_effi_3D_comp.SetTitleOffset(1.6, "X")
        h3_effi_3D_comp.SetTitleOffset(1.8, "Y")
        h3_effi_3D_comp.SetTitleOffset(1.5, "Z")
        h3_effi_3D_comp.SetZTitle("#phi") #("#varepsilon_{fit}/#varepsilon_{measured}")
        h3_effi_3D_comp.GetXaxis().CenterTitle()
        h3_effi_3D_comp.GetYaxis().CenterTitle()
        h3_effi_3D_comp.GetZaxis().CenterTitle()
        canvas.SetRightMargin(1.0)
        h3_effi_3D_comp.Draw("LEGO2 COLZ")
        latex.DrawLatexNDC(.08, .93, "#font[61]{CMS} #font[52]{#scale[0.8]{Simulation}}")
        latex.DrawLatexNDC(.08, .89, "#chi^{{2}}={0:.2f}".format(fitter.GetChi2()))
        canvas.Print("effi_3D_comp_{0}.pdf".format(q2bins[self.process.cfg['binKey']]['label']))

    @staticmethod
    def isPosiDef(formula2D):
        f2_min_x, f2_min_y = ROOT.Double(0), ROOT.Double(0)
        formula2D.GetMinimumXY(f2_min_x, f2_min_y)
        f2_min = formula2D.Eval(f2_min_x, f2_min_y)
        if f2_min > 0:
            return True
        else:
            print("WARNING\t: Sanitary check failed: Minimum of efficiency map is {0:.2e} at {1}, {2}".format(f2_min, f2_min_x, f2_min_y))
        return False
