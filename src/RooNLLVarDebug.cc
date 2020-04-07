#include "HiggsAnalysis/CombinedLimit/interface/RooNLLVarDebug.h"

#include "TMath.h"

#include <iostream>
#include <sstream>

Double_t RooNLLVarDebug::getVal() const {

    // Evaluate as straight FUNC
    Int_t nFirst(0), nLast(_nEvents), nStep(1) ;
    
    switch (_mpinterl) {
    case RooFit::BulkPartition:
      nFirst = _nEvents * _setNum / _numSets ;
      nLast  = _nEvents * (_setNum+1) / _numSets ;
      nStep  = 1 ;
      break;
      
    case RooFit::Interleave:
      nFirst = _setNum ;
      nLast  = _nEvents ;
      nStep  = _numSets ;
      break ;
      
    case RooFit::SimComponents:
      nFirst = 0 ;
      nLast  = _nEvents ;
      nStep  = 1 ;
      break ;
      
    case RooFit::Hybrid:
      throw(std::string("this should never happen")) ;
      break ;
    }

    Double_t ret = evaluatePartition(nFirst,nLast,nStep);

    if (numSets()==1) {
      const Double_t norm = globalNormalization();
      ret /= norm;
      _evalCarry /= norm;
    }
    
    return ret ;
}

Double_t RooNLLVarDebug::evaluatePartition(Int_t firstEvent, Int_t lastEvent, Int_t stepSize) const
{
  std::cout << "DEBUGNLLSTART" << std::endl;
  // Throughout the calculation, we use Kahan's algorithm for summing to
  // prevent loss of precision - this is a factor four more expensive than
  // straight addition, but since evaluating the PDF is usually much more
  // expensive than that, we tolerate the additional cost...
  Int_t i ;
  Double_t result(0), carry(0);

  RooAbsPdf* pdfClone = (RooAbsPdf*) _funcClone ;

  // cout << "RooNLLVar::evaluatePartition(" << GetName() << ") projDeps = " << (_projDeps?*_projDeps:RooArgSet()) << endl ;

  _dataClone->store()->recalculateCache( _projDeps, firstEvent, lastEvent, stepSize,(_binnedPdf?kFALSE:kTRUE) ) ;

  Double_t sumWeight(0), sumWeightCarry(0);

  // If pdf is marked as binned - do a binned likelihood calculation here (sum of log-Poisson for each bin)
  Double_t resultTerm(0);
  if (_binnedPdf) {

    std::stringstream debug_msg;

    for (i=firstEvent ; i<lastEvent ; i+=stepSize) {

      _dataClone->get(i) ;

      if (!_dataClone->valid()) continue;

      Double_t eventWeight = _dataClone->weight();

      // debug data format: r bin data r*s+b nll
      // Calculate log(Poisson(N|mu) for this bin
      Double_t N = eventWeight ;
      Double_t mu = _binnedPdf->getVal()*_binw[i] ;
      //cout << "RooNLLVar::binnedL(" << GetName() << ") N=" << N << " mu = " << mu << endl ;

      Double_t term = 0.;

      if (mu<=0 && N>0) {

	// Catch error condition: data present where zero events are predicted
	logEvalError(Form("Observed %f events in bin %d with zero event yield",N,i)) ;

      } else if (fabs(mu)<1e-10 && fabs(N)<1e-10) {

	// Special handling of this case since log(Poisson(0,0)=0 but can't be calculated with usual log-formula
	// since log(mu)=0. No update of result is required since term=0.

      } else {

	term = -1*(-mu + N*log(mu) - TMath::LnGamma(N+1)) ;

	// Kahan summation of sumWeight
	Double_t y = eventWeight - sumWeightCarry;
	Double_t t = sumWeight + y;
	sumWeightCarry = (t - sumWeight) - y;
	sumWeight = t;

	// Kahan summation of result
	y = term - carry;
	t = result + y;
	carry = (t - result) - y;
	result = t;
      }
	  
	  debug_msg << "DEBUGNLLHIST: " << (poi_ ? poi_->getVal() : TMath::QuietNaN()) << " " << i << " " << N << " " << mu << " " << term << "\n";
      resultTerm += term;
    }
    std::cout << debug_msg.str() << std::endl;


  } else {

    for (i=firstEvent ; i<lastEvent ; i+=stepSize) {

      _dataClone->get(i) ;

      if (!_dataClone->valid()) continue;

      Double_t eventWeight = _dataClone->weight();
      if (0. == eventWeight * eventWeight) continue ;
      if (_weightSq) eventWeight = _dataClone->weightSquared() ;

      Double_t term = -eventWeight * pdfClone->getLogVal(_normSet);


      Double_t y = eventWeight - sumWeightCarry;
      Double_t t = sumWeight + y;
      sumWeightCarry = (t - sumWeight) - y;
      sumWeight = t;

      y = term - carry;
      t = result + y;
      carry = (t - result) - y;
      result = t;
    }

    // include the extended maximum likelihood term, if requested
    if(_extended && _setNum==_extSet) {
      if (_weightSq) {

	// Calculate sum of weights-squared here for extended term
	Double_t sumW2(0), sumW2carry(0);
	for (i=0 ; i<_dataClone->numEntries() ; i++) {
	  _dataClone->get(i);
	  Double_t y = _dataClone->weightSquared() - sumW2carry;
	  Double_t t = sumW2 + y;
	  sumW2carry = (t - sumW2) - y;
	  sumW2 = t;
	}

	Double_t expected= pdfClone->expectedEvents(_dataClone->get());

	// Adjust calculation of extended term with W^2 weighting: adjust poisson such that
	// estimate of Nexpected stays at the same value, but has a different variance, rescale
        // both the observed and expected count of the Poisson with a factor sum[w] / sum[w^2] which is
        // the effective weight of the Poisson term.
	// i.e. change Poisson(Nobs = sum[w]| Nexp ) --> Poisson( sum[w] * sum[w] / sum[w^2] | Nexp * sum[w] / sum[w^2] )
        // weighted by the effective weight  sum[w^2]/ sum[w] in the likelihood.
        // Since here we compute the likelihood with the weight square we need to multiply by the
        // square of the effective weight
        // expectedW = expected * sum[w] / sum[w^2]   : effective expected entries
        // observedW =  sum[w]  * sum[w] / sum[w^2]   : effective observed entries
        // The extended term for the likelihood weighted by the square of the weight will be then:
        //  (sum[w^2]/ sum[w] )^2 * expectedW -  (sum[w^2]/ sum[w] )^2 * observedW * log (expectedW)  and this is
        //  using the previous expressions for expectedW and observedW
        //  sum[w^2] / sum[w] * expected - sum[w^2] * log (expectedW)
        //  and since the weights are constants in the likelihood we can use log(expected) instead of log(expectedW)

	Double_t expectedW2 = expected * sumW2 / _dataClone->sumEntries() ;
	Double_t extra= expectedW2 - sumW2*log(expected );

	// Double_t y = pdfClone->extendedTerm(sumW2, _dataClone->get()) - carry;

	Double_t y = extra - carry ;

	Double_t t = result + y;
	carry = (t - result) - y;
	result = t;
      } else {
	Double_t y = pdfClone->extendedTerm(_dataClone->sumEntries(), _dataClone->get()) - carry;
	Double_t t = result + y;
	carry = (t - result) - y;
	result = t;
      }
    }
  }


  // If part of simultaneous PDF normalize probability over
  // number of simultaneous PDFs: -sum(log(p/n)) = -sum(log(p)) + N*log(n)
  if (_simCount>1) {
    Double_t y = sumWeight*log(1.0*_simCount) - carry;
    Double_t t = result + y;
    carry = (t - result) - y;
    result = t;
  }

  //timer.Stop() ;
  //cout << "RooNLLVar::evalPart(" << GetName() << ") SET=" << _setNum << " first=" << firstEvent << ", last=" << lastEvent << ", step=" << stepSize << ") result = " << result << " CPU = " << timer.CpuTime() << endl ;

  // At the end of the first full calculation, wire the caches
  if (_first) {
    _first = kFALSE ;
    _funcClone->wireAllCaches() ;
  }


  // Check if value offset flag is set.
  if (_doOffset) {

    // If no offset is stored enable this feature now
    if (_offset==0 && result !=0 ) {
      coutI(Minimization) << "RooNLLVar::evaluatePartition(" << GetName() << ") first = "<< firstEvent << " last = " << lastEvent << " Likelihood offset now set to " << result << std::endl ;
      _offset = result ;
      _offsetCarry = carry;
    }

    // Substract offset
    Double_t y = -_offset - (carry + _offsetCarry);
    Double_t t = result + y;
    carry = (t - result) - y;
    result = t;
  }


  _evalCarry = carry;
  
  std::cout << "DEBUGNLLSUM: " << resultTerm << " " << result << std::endl;
  return result ;
}
