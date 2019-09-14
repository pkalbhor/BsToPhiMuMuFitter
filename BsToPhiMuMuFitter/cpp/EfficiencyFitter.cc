#include "TH2.h"
#include "TH3.h"
#include "TF2.h"
#include "TF3.h"
#include "TMinuit.h"

#ifndef EFFICIENCYFITTER_H
#define EFFICIENCYFITTER_H

TH3 *h3_fcn = 0;
TF3 *f3_fcn = 0;
int chi2Val = 0;

void fcn_binnedChi2_3D(int &npar, double *gin, double &f, double *par, int iflag)
{//{{{
    f=0;
    for (int i = 1; i <= h3_fcn->GetNbinsX(); i++) {
        for (int j = 1; j <= h3_fcn->GetNbinsY(); j++) {
            for (int m = 1; m <= h3_fcn->GetNbinsZ(); m++) {
                int gBin = h3_fcn->GetBin(i,j, m);
                double measure  = h3_fcn->GetBinContent(gBin);
                double error    = h3_fcn->GetBinError(gBin);
                for (int k = 0; k < f3_fcn->GetNpar(); k++){
                    f3_fcn->SetParameter(k,par[k]);
                }
                double xi = h3_fcn->GetXaxis()->GetBinLowEdge(i);
                double xf = h3_fcn->GetXaxis()->GetBinUpEdge(i);
                double yi = h3_fcn->GetYaxis()->GetBinLowEdge(j);
                double yf = h3_fcn->GetYaxis()->GetBinUpEdge(j);
                double zi = h3_fcn->GetZaxis()->GetBinLowEdge(m);
                double zf = h3_fcn->GetZaxis()->GetBinUpEdge(m);
                f += pow( (f3_fcn->Integral(xi,xf,yi,yf,zi,zf)/(xf-xi)/(yf-yi)/(zf-zi)-measure)/error,2);
            }
        }
    }

    chi2Val = f;

    // Prevent from negative function
    double f3_minX, f3_minY, f3_minZ;
    f3_fcn->GetMinimumXYZ(f3_minX, f3_minY, f3_minZ);
    if (f3_fcn->Eval(f3_minX, f3_minY, f3_minZ) < 0){
        f += 100*h3_fcn->GetNbinsX()*h3_fcn->GetNbinsY()*h3_fcn->GetNbinsZ();
    }

}//}}}

class EfficiencyFitter{
public:
    EfficiencyFitter();
    virtual ~EfficiencyFitter();
    TH3* GetH3(){return h3_fcn;}
    TF3* GetF3(){return f3_fcn;}
    int  GetChi2(){return chi2Val;}
    TMinuit* Init(int, TH3*, TF3*);
private:
    TMinuit *minuit = 0;
};

EfficiencyFitter::EfficiencyFitter(){}
EfficiencyFitter::~EfficiencyFitter(){
    delete minuit;
    minuit = 0;
}
TMinuit* EfficiencyFitter::Init(int nPar, TH3 *h3, TF3 *f3){
    h3_fcn = h3;
    f3_fcn = f3;
    minuit = new TMinuit(nPar);
    minuit->SetFCN(fcn_binnedChi2_3D);
    return minuit;
}
#endif
