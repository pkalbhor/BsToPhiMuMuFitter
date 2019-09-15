#!/usr/bin/env python
# -*- coding: utf-8 -*-
# vim: set sw=4 ts=4 fdm=indent fdl=2 ft=python et:

# Description     : Shared object definition.

from ROOT import RooRealVar
from ROOT import RooArgSet, RooArgList

Bmass = RooRealVar("Bmass","m_{K^{*}#mu#mu} [GeV]", 4.9, 5.9)
CosThetaK = RooRealVar("CosThetaK", "cos#theta_{K}", -1., 1.)
CosThetaL = RooRealVar("CosThetaL", "cos#theta_{l}", -1., 1.)
Phi = RooRealVar("Phi","\phi", -3.15, 3.15)
Mumumass = RooRealVar("Mumumass", "m_{#mu#mu} [GeV]", 0., 10.)
Mumumasserr = RooRealVar("Mumumasserr", "Error of m_{#mu#mu} [GeV]", 0., 10.)
Phimass = RooRealVar("Phimass", "m_{K^{+} K^{-}} [GeV]", 1.01, 1.03)
# Kshortmass = RooRealVar("Kshortmass", "m_{K_{S}} [GeV]", 0.427, 0.577)
Q2 = RooRealVar("Q2", "q^{2} [GeV^{2}]", 0.5, 20.)
Triggers = RooRealVar("Triggers", "", 0, 100)
Bdt = RooRealVar("Bdt", "", 0, 1)
arglist = RooArgList(
    Bmass,
    CosThetaK,
    CosThetaL,
    Phi,
    Mumumass,
    Mumumasserr,
    Phimass,
    Q2,
    Triggers
   )

dataArgs = RooArgSet(
        arglist,
        Bdt
        )
genCosThetaK = RooRealVar("genCosThetaK", "cos#theta_{K}", -1., 1.)
genCosThetaL = RooRealVar("genCosThetaL", "cos#theta_{l}", -1., 1.)
genPhi = RooRealVar("genPhi", "\Phi^{0}", -3.14, 3.14)
genQ2 = RooRealVar("genQ2", "q^{2} [GeV^{2}]", 0.5, 20.)
dataArgsGEN = RooArgSet(
    genQ2,
    genCosThetaK,
    genCosThetaL,
    genPhi
    )
