#ifndef ROONLLVARDEBUG_H
#define ROONLLVARDEBUG_H

#include "RooNLLVar.h"
#include "RooAbsReal.h"
#include "RooAbsData.h"
#include "RooAbsDataStore.h"
#include "RooRealSumPdf.h"
#include "RooRealVar.h"

class RooNLLVarDebug : public RooNLLVar {

public:

RooNLLVarDebug(const RooNLLVar& other, const char* name=0) : RooNLLVar(other,name) {}
RooNLLVarDebug(const RooAbsReal& other, const char* name=0) : RooNLLVar(*static_cast<const RooNLLVar*>(&other),name) {}

void setPoi(RooRealVar* poi) { poi_ = poi; }

double getVal() const;

protected:
Double_t evaluatePartition(Int_t firstEvent, Int_t lastEvent, Int_t stepSize) const override;

RooRealVar * poi_ = nullptr;

};

#endif
