/*****************************************************************************
 * Project: CMS detector at the CERN
 *
 * Package: PhysicsTools/TagAndProbe/RooCMSShape
 *
 *
 * Authors:
 *   Nadia Adam, Princeton - neadam@princeton.edu
 *   Adam Hunt, Princeton  - ahunt@princeton.edu
 *   Kalanand Mishra, Fermilab - kalanand@fnal.gov
 *
 * Description:
 *   Defines a probability density function which has exponential decay 
 *   distribution at high mass beyond the pole position (say, Z peak)  
 *   but turns over (i.e., error function) at low mass due to threshold 
 *   effect. We use this to model the background shape in Z->ll invariant 
 *   mass.
 * History:
 *   
 *
 *****************************************************************************/

#include "../include/RooBkgShape.h"

ClassImp(RooBkgShape) 

 RooBkgShape::RooBkgShape(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _alpha,
                        RooAbsReal& _beta,
                        RooAbsReal& _gamma,
                        RooAbsReal& _peak,
                        RooAbsReal& _gamma2,
                        RooAbsReal& _peak2) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   alpha("alpha","alpha",this,_alpha),
   beta("beta","beta",this,_beta),
   gamma("gamma","gamma",this,_gamma),
   peak("peak","peak",this,_peak),
   gamma2("gamma2","gamma2",this,_gamma2),
   peak2("peak2","peak2",this,_peak2)
 { } 


 RooBkgShape::RooBkgShape(const RooBkgShape& other, const char* name):
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   alpha("alpha",this,other.alpha),
   beta("beta",this,other.beta),
   gamma("gamma",this,other.gamma),
   peak("peak",this,other.peak),
   gamma2("gamma2",this,other.gamma2),
   peak2("peak2",this,other.peak2)
 { } 



 Double_t RooBkgShape::evaluate() const 
 { 
  // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE 

  //Double_t erf = TMath::Erfc((alpha - x) * beta);
  Double_t erf = RooMath::erfc((alpha - x) * beta);
  Double_t u = (x - peak)*gamma;
  Double_t u2 = (x - peak2)*gamma2;
  
  if(u < -70) {u = 1e20;u2= 1e20;}
  else if( u>70 ) {u = 0;u2=0;}
  else {u = exp(-u); u2=exp(-u2);}   //exponential decay
  return erf*(u+u2);
 } 
