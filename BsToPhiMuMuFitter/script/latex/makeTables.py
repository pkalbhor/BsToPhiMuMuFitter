#!/usr/bin/env python
# vim: set sw=4 sts=4 fdm=indent fdl=0 fdn=3 nopaste et:

from __future__ import print_function

import os
import math
import shelve
from collections import OrderedDict

import ROOT
import BsToPhiMuMuFitter.cpp

from BsToPhiMuMuFitter.anaSetup import q2bins, modulePath
from BsToPhiMuMuFitter.StdFitter import unboundFlToFl, unboundAfbToAfb
from BsToPhiMuMuFitter.StdProcess import p

# For developers:
#   * Input db is forced to be StdProcess.dbplayer.absInputDir
#   * function name for labelled table is table_label1[_label2]


db_dir = p.dbplayer.absInputDir

indent = "  "

def table_sysFL_sysAFB():
    baseIndentLevel = 2
    for var in ["fl", "afb"]:
        dbKeyToLine = OrderedDict()
        dbKeyToLine['syst_randEffi'] = [r"Limited MC size"]
        dbKeyToLine['syst_altEffi'] = [r"Eff.\ mapping"]
        dbKeyToLine['syst_simMismodel'] = [r"Simu.\ mismodel"]
        dbKeyToLine['syst_altSP'] = [r"$S$ - $P$ wave interf.\ "]
        dbKeyToLine['syst_altBkgCombA'] = [r"Comb.\ Bkg.\ shape"]
        dbKeyToLine['syst_vetoJpsiX'] = [r"$J/\psi + X$ contrib.\ "]
        dbKeyToLine['syst_altFitRange'] = [r"$B$ mass range"]
        totalErrorLine = ["Total"]
        for binKey in ['belowJpsi', 'betweenPeaks', 'abovePsi2s', 'summary']:
            db = shelve.open("{0}/fitResults_{1}.db".format(db_dir, q2bins[binKey]['label']))
            totalSystErr = 0.
            for systKey, latexLine in dbKeyToLine.items():
                err = db["{0}_{1}".format(systKey, var)]['getError']
                latexLine.append("{0:.03f}".format(err))
                totalSystErr += pow(err, 2)
            db.close()
            totalErrorLine.append("{0:.03f}".format(math.sqrt(totalSystErr)))

        print("[table_sysFL_sysAFB] Printing table of syst. unc. for {0}".format(var))
        print("")
        print(indent * (baseIndentLevel + 0) + r"\begin{tabular}{|l|c|c|c|c|}")
        print(indent * (baseIndentLevel + 1) + r"\hline")
        print(indent * (baseIndentLevel + 1) + r"Syst.\ err.\ $\backslash$ $q^2$ bin & 1 & 3 & 5 & 0 \\")
        print(indent * (baseIndentLevel + 1) + r"\hline")
        print(indent * (baseIndentLevel + 1) + r"\hline")
        print(indent * (baseIndentLevel + 1) + r"\multicolumn{5}{|c|}{Uncorrelated systematic uncertainties} \\")
        print(indent * (baseIndentLevel + 1) + r"\hline")
        for systKey, latexLine in dbKeyToLine.items():
            print(indent * (baseIndentLevel + 1) + " & ".join(latexLine) + r" \\")
        print(indent * (baseIndentLevel + 1) + r"\hline")
        print(indent * (baseIndentLevel + 1) + " & ".join(totalErrorLine) + r" \\")
        print(indent * (baseIndentLevel + 1) + r"\hline")
        print(indent * (baseIndentLevel + 0) + r"\end{tabular}")
        print("")

def table_yields():
    baseIndentLevel = 2
    print("[table_yields] Printing table of yields")
    print("")
    print(indent * (baseIndentLevel + 0) + r"\begin{tabular}{|c|c|c|}")
    print(indent * (baseIndentLevel + 1) + r"\hline")
    print(indent * (baseIndentLevel + 1) + r"$q^2$ bin & $Y_S$ & $Y^C_B$ \\")
    print(indent * (baseIndentLevel + 1) + r"\hline")
    print(indent * (baseIndentLevel + 1) + r"\hline")

    binKeyToLine = OrderedDict()
    binKeyToLine['belowJpsi'] = ["1"]
    binKeyToLine['betweenPeaks'] = ["3"]
    binKeyToLine['abovePsi2s'] = ["5"]
    binKeyToLine['summary'] = ["0"]
    for binKey, latexLine in binKeyToLine.items():
        db = shelve.open("{0}/fitResults_{1}.db".format(db_dir, q2bins[binKey]['label']))
        latexLine.append("${0:.01f} \pm {1:.01f}$".format(db['nSig']['getVal'], db['nSig']['getError']))
        latexLine.append("${0:.01f} \pm {1:.01f}$".format(db['nBkgComb']['getVal'], db['nBkgComb']['getError']))
        db.close()
        print(indent * (baseIndentLevel + 1) + " & ".join(latexLine) + r" \\")

    print(indent * (baseIndentLevel + 1) + r"\hline")
    print(indent * (baseIndentLevel + 0) + r"\end{tabular}")
    print("")

def table_coverageAFBFL():
    baseIndentLevel = 2
    print("[table_coverageAFBFL] Printing table of stat error coverage")
    print("")
    print(indent * (baseIndentLevel + 0) + r"\begin{tabular}{|c|c|c|}")
    print(indent * (baseIndentLevel + 1) + r"\hline")
    print(indent * (baseIndentLevel + 1) + r"$q^2$ bin & $A_{FB}$ & $F_{L}$ \\")
    print(indent * (baseIndentLevel + 1) + r"\hline")
    print(indent * (baseIndentLevel + 1) + r"\hline")

    binKeyToLine = OrderedDict()
    binKeyToLine['belowJpsi'] = ["1"]
    binKeyToLine['betweenPeaks'] = ["3"]
    binKeyToLine['abovePsi2s'] = ["5"]
    binKeyToLine['summary'] = ["0"]
    for binKey, latexLine in binKeyToLine.items():
        db = shelve.open("{0}/fitResults_{1}.db".format(db_dir, q2bins[binKey]['label']))
        latexLine.append("{0:.1f}\%".format(db['stat_FC_afb']['coverage'] * 100.))
        latexLine.append("{0:.1f}\%".format(db['stat_FC_fl']['coverage'] * 100.))
        db.close()
        print(indent * (baseIndentLevel + 1) + " & ".join(latexLine) + r" \\")

    print(indent * (baseIndentLevel + 1) + r"\hline")
    print(indent * (baseIndentLevel + 0) + r"\end{tabular}")
    print("")

def table_dataresAFBFL():
    baseIndentLevel = 2

    print("[table_dataresAFBFL] Printing table of final result")
    print("")
    print(indent * (baseIndentLevel + 0) + r"\begin{tabular}{|c|c|c|c|c|}")
    print(indent * (baseIndentLevel + 1) + r"\hline")
    print(indent * (baseIndentLevel + 1) + r"$q^2$ bin index & $q^2$ range (in $\GeV^2$) & Signal Yield & $A_{FB}$ & $F_{L}$ \\")
    print(indent * (baseIndentLevel + 1) + r"\hline")
    print(indent * (baseIndentLevel + 1) + r"\hline")

    binKeyToLine = OrderedDict()
    binKeyToLine['belowJpsi'] = ["1", r"1.00 -- 8.68"]
    binKeyToLine['jpsi'] = ["2", r"8.68 -- 10.09", r"\multicolumn{3}{|c|} {$\JPsi$ resonance region}"]
    binKeyToLine['betweenPeaks'] = ["3", r"10.09 -- 12.86"]
    binKeyToLine['psi2s'] = ["4", r"12.86 -- 14.18", r"\multicolumn{3}{|c|} {$\psi'$ resonance region}"]
    binKeyToLine['abovePsi2s'] = ["5", r"14.18 -- 19.00"]
    binKeyToLine['summary'] = ["0", r"bin\#1 $+$ bin\#3 $+$ bin\#5"]

    syst_sources = [
        'syst_randEffi',
        'syst_altEffi',
        'syst_simMismodel',
        'syst_altSP',
        'syst_altBkgCombA',
        'syst_vetoJpsiX',
        #  'syst_altFitRange',
    ]
    for binKey, latexLine in binKeyToLine.items():
        if binKey not in ['jpsi', 'psi2s']:
            db = shelve.open(r"{0}/fitResults_{1}.db".format(db_dir, q2bins[binKey]['label']))
            latexLine.append(r"${0:.01f} \pm {1:.01f}$".format(db['nSig']['getVal'], db['nSig']['getError']))
            fl = unboundFlToFl(db['unboundFl']['getVal'])
            latexLine.append("${0:.2f}^{{{1:+.2f}}}_{{{2:+.2f}}} \pm {3:.2f}$".format(
                unboundAfbToAfb(db['unboundAfb']['getVal'], fl),
                db['stat_FC_afb']['getErrorHi'],
                db['stat_FC_afb']['getErrorLo'],
                math.sqrt(sum([pow(db[v + '_afb']['getError'], 2) for v in syst_sources]))))
            latexLine.append("${0:.2f}^{{{1:+.2f}}}_{{{2:+.2f}}} \pm {3:.2f}$".format(
                fl,
                db['stat_FC_fl']['getErrorHi'],
                db['stat_FC_fl']['getErrorLo'],
                math.sqrt(sum([pow(db[v + '_fl']['getError'], 2) for v in syst_sources]))))
            db.close()
        print(indent * (baseIndentLevel + 1) + " & ".join(latexLine) + r" \\")
        print(indent * (baseIndentLevel + 1) + r"\hline")

    print(indent * (baseIndentLevel + 0) + r"\end{tabular}")
    print("")

def table_FinalDataresAFBFL():
    table_dataresAFBFL()

if __name__ == '__main__':
    #  table_sysFL_sysAFB()
    #  table_yields()
    #  table_coverageAFBFL()
    #  table_dataresAFBFL()
    #  table_FinalDataresAFBFL()
