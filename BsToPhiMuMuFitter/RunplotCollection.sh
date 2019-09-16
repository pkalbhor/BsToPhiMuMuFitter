#!/bin/bash
for arg in 'belowJpsiA' 'belowJpsiB' 'belowJpsiC' 'betweenPeaks' 'abovePsi2sA' 'abovePsi2sB' 'summary' 'summaryLowQ2' 
do
    python plotCollection.py $arg
done
