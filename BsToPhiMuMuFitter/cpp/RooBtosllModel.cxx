/***************************************************************************** 
 * Project: RooFit                                                           * 
 *                                                                           * 
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/ 

// Your description goes here... 

#include "Riostream.h" 

#include "RooBtosllModel.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include <math.h> 
#include "TMath.h" 

ClassImp(RooBtosllModel); 

RooBtosllModel::RooBtosllModel(const char *name, const char *title, 
                       RooAbsReal& _CosThetaL,
                       RooAbsReal& _CosThetaK,
                       RooAbsReal& _unboundAfb,
                       RooAbsReal& _unboundFl/*,
                       RooAbsReal& _fs,
                       RooAbsReal& _transAs*/) :
  RooAbsPdf(name,title), 
  CosThetaL("CosThetaL","CosThetaL",this,_CosThetaL),
  CosThetaK("CosThetaK","CosThetaK",this,_CosThetaK),
  unboundAfb("unboundAfb","unboundAfb",this,_unboundAfb),
  unboundFl("unboundFl","unboundFl",this,_unboundFl)/*,
  fs("fs","fs",this,_fs),
  transAs("transAs","transAs",this,_transAs)*/
{ 
} 


RooBtosllModel::RooBtosllModel(const RooBtosllModel& other, const char* name) :  
  RooAbsPdf(other,name), 
  CosThetaL("CosThetaL",this,other.CosThetaL),
  CosThetaK("CosThetaK",this,other.CosThetaK),
  unboundAfb("unboundAfb",this,other.unboundAfb),
  unboundFl("unboundFl",this,other.unboundFl)/*,
  fs("fs",this,other.fs),
  transAs("transAs",this,other.transAs)*/
{ 
} 


Double_t RooBtosllModel::evaluate() const 
{ 
  // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE 
  // Remark: The tranAs formula in older AN is wrong, unboundFl should be replaced with fl
  Double_t fl = 0.5+TMath::ATan(unboundFl)/TMath::Pi() ; // [0.,1.]
  Double_t afb = 2.0*(1.-fl)*TMath::ATan(unboundAfb)/TMath::Pi() ; // [-0.75, 0.75]
  //Double_t as = 1.78*TMath::Sqrt(3.*fs*(1.-fs)*fl)*transAs ; // [-(fs+3*fl(1-fs)), fs+3*fs*(1-fs)]
  
  //Double_t afb = (0.5-TMath::ATan(unboundFl)/TMath::Pi())*2*TMath::ATan(unboundAfb)/TMath::Pi(); 
  Double_t result = (9.0/16.0)*((0.5*(1.0-fl)*(1.0-CosThetaK*CosThetaK)*(1.0+CosThetaL*CosThetaL)) + (2.0*fl*CosThetaK*CosThetaK*(1.0-CosThetaL*CosThetaL)) + (afb*(1.0-CosThetaK*CosThetaK)*CosThetaL)) ;
  
  //Double_t result = 0.5625*((0.666667*fs+1.333333*as*CosThetaK)*(1.-pow(CosThetaL,2))+(1.-fs)*(2.*fl*pow(CosThetaK,2)*(1.-pow(CosThetaL,2))+0.5*(1.-fl)*(1.-pow(CosThetaK,2))*(1.+pow(CosThetaL,2))+1.333333*afb*(1.-pow(CosThetaK,2))*CosThetaL)) ;
  // Double_t result = 0.5625*((0.666667*fs+2.666667*transAs*TMath::Sqrt(3.*fs*(1.-fs)*(0.5+TMath::ATan(unboundFl)/TMath::Pi()))*CosThetaK)*(1-pow(CosThetaL,2))+(1-fs)*(2*(0.5+TMath::ATan(unboundFl)/TMath::Pi())*pow(CosThetaK,2)*(1.-pow(CosThetaL,2))+0.5*(0.5-TMath::ATan(unboundFl)/TMath::Pi())*(1.-pow(CosThetaK,2))*(1.+pow(CosThetaL,2))+(2.*(0.5-TMath::ATan(unboundFl)/TMath::Pi())*TMath::ATan(unboundAfb)/TMath::Pi())*(1.-pow(CosThetaK,2))*CosThetaL)) ;

  return result;
  // return result > 0 ? result : 0;
} 
