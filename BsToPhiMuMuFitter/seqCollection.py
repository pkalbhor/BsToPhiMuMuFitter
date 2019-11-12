#!/usr/bin/env python
# -*- coding: utf-8 -*-
# vim: set sw=4 ts=4 fdm=indent fdl=1 fdn=3 ft=python et:

import sys
import BsToPhiMuMuFitter.cpp
import BsToPhiMuMuFitter.dataCollection as dataCollection
import BsToPhiMuMuFitter.toyCollection as toyCollection
import BsToPhiMuMuFitter.pdfCollection as pdfCollection
import BsToPhiMuMuFitter.fitCollection as fitCollection

from BsToPhiMuMuFitter.StdProcess import p

# Standard fitting procedures
predefined_sequence = {}
predefined_sequence['loadData'] = [dataCollection.dataReader]
predefined_sequence['buildAllPdfs'] = [dataCollection.dataReader, pdfCollection.stdWspaceReader, pdfCollection.stdPDFBuilder]
predefined_sequence['buildEfficiecyHist'] = [dataCollection.effiHistReader]

predefined_sequence['fitEfficiency'] = [dataCollection.effiHistReader, pdfCollection.stdWspaceReader, fitCollection.effiFitter]
predefined_sequence['fitSigM'] = [dataCollection.sigMCReader, pdfCollection.stdWspaceReader, fitCollection.sigMFitter]
predefined_sequence['fitBkgCombA'] = [dataCollection.dataReader, pdfCollection.stdWspaceReader, fitCollection.bkgCombAFitter]
predefined_sequence['fitFinal3D'] = [dataCollection.dataReader, pdfCollection.stdWspaceReader, fitCollection.finalFitter]

predefined_sequence['stdFit'] = [dataCollection.effiHistReader, dataCollection.sigMCReader, dataCollection.dataReader, pdfCollection.stdWspaceReader, fitCollection.effiFitter, fitCollection.sigMFitter, fitCollection.bkgCombAFitter, fitCollection.sig2DFitter, dataCollection.sigMCGENReader, fitCollection.sigAFitter, fitCollection.finalFitter]

# For fitter validation and syst
predefined_sequence['fitSig2D'] = [dataCollection.sigMCReader, pdfCollection.stdWspaceReader, fitCollection.sig2DFitter]
predefined_sequence['fitSigMCGEN'] = [dataCollection.sigMCGENReader, pdfCollection.stdWspaceReader, fitCollection.sigAFitter]

if __name__ == '__main__':
    binKey = ['belowJpsiA']#, 'belowJpsiB', 'belowJpsiC', 'betweenPeaks', 'abovePsi2sA', 'abovePsi2sB', 'summary', 'summaryLowQ2']
    binKey = [sys.argv[1]]
    for b in binKey:
        p.cfg['binKey'] = b
        #  p.setSequence(predefined_sequence['buildEfficiecyHist'])
        # p.setSequence(predefined_sequence['fitEfficiency'])
        # p.setSequence(predefined_sequence['fitBkgCombA'])
        # p.setSequence(predefined_sequence['fitFinal3D'])
        # p.setSequence(predefined_sequence['fitSigM'])
        # p.setSequence(predefined_sequence['fitSig2D']) #Fitting MC distributions
        # p.setSequence(predefined_sequence['fitSigMCGEN'])
        p.setSequence(predefined_sequence['stdFit'])
        try:
            p.beginSeq()
            p.runSeq()
        finally:
            p.endSeq()
