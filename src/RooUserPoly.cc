 /*****************************************************************************
  * Project: RooFit                                                           *
  * Package: RooFitModels                                                     *
  * @(#)root/roofit:$Id$
  * Authors:                                                                  *
  *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
  *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
  *                                                                           *
  * Copyright (c) 2000-2005, Regents of the University of California          *
  *                          and Stanford University. All rights reserved.    *
  *                                                                           *
  * Redistribution and use in source and binary forms,                        *
  * with or without modification, are permitted according to the terms        *
  * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
  *****************************************************************************/
 
 /**
 \file RooUserPoly.cxx
 \class RooUserPoly
 \ingroup Roofit
 
 P.d.f implementing the Crystall Ball line shape
 **/
 
 #include "RooFit.h"
 
 #include "Riostream.h"
 #include "Riostream.h"
 #include <math.h>
 
 #include "../include/RooUserPoly.h"
 #include "RooAbsReal.h"
 #include "RooRealVar.h"
 #include "RooMath.h"
 #include "TMath.h"
 
 #include "TError.h"
 
 using namespace std;
 
 ClassImp(RooUserPoly)
 ;
 
 
 ////////////////////////////////////////////////////////////////////////////////
 
 RooUserPoly::RooUserPoly(const char *name, const char *title, RooAbsReal& _m, RooAbsReal& _slope) :
   RooAbsPdf(name, title),
   m("m", "Dependent", this, _m),
   slope("slope", "Slope", this, _slope)
 {
 }
 
 
 ////////////////////////////////////////////////////////////////////////////////
 
 RooUserPoly::RooUserPoly(const RooUserPoly& other, const char* name) :
   RooAbsPdf(other, name), m("m", this, other.m), slope("slope", this, other.slope)
 {
 }
 
 
 ////////////////////////////////////////////////////////////////////////////////
 
 Double_t RooUserPoly::evaluate() const {
   Double_t constant(1);
   Double_t minvalue = slope*40.0;
   Double_t maxvalue = slope*140.0;
   if(minvalue <= 0 || maxvalue <= 0)constant = minvalue < maxvalue? fabs(minvalue)+1: fabs(maxvalue)+1;
   
   return slope*m + constant;

 }
 
 
