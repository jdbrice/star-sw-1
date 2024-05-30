/* Copyright 2013
 *   Authors: Sergey Yashchenko and Tadeas Bilka
 *
 *   This file is part of GENFIT.
 *
 *   GENFIT is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU Lesser General Public License as published
 *   by the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   GENFIT is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU Lesser General Public License for more details.
 *
 *   You should have received a copy of the GNU Lesser General Public License
 *   along with GENFIT.  If not, see <http://www.gnu.org/licenses/>.
 */
/** @addtogroup genfit
 * @{
 */

/* July 8th, 2022
 * Gavin Wilks
 * Modify GFGbl.h/cc class to work with Forward Tracking System at STAR.
 *
 * This is used for the alignment of the forward system, which includes
 * the Forward Silicon Tracker (FST) and Forward Thickgem Tracker (FTT).
 *
 * This class is used to create the mille.dat file, which is required
 * for alignment with Millepede II.
 */

#ifndef STFWDGBL_H
#define STFWDGBL_H

#include "GenFit/GblTrajectory.h"
#include "GenFit/AbsFitter.h"

#include <map>
#include <iostream>

#include <TMatrixD.h>
#include <assert.h>
#include <sstream>

#include <TMath.h>
#include <TVector3.h>

#include <TH1D.h>
#include <TFile.h>


//namespace genfit {
  
  
/** @brief Generic GBL implementation
 * 
 * The interface class to GBL track fit
 *
 */
class StFwdGbl : public genfit::AbsFitter {
  
private:
  StFwdGbl(const StFwdGbl&);
  StFwdGbl& operator=(StFwdGbl const&);
  
  std::string m_milleFileName;
  std::string m_gblInternalIterations;
  double m_pValueCut;
  int m_minNdf;
  double m_chi2Cut;
  bool m_enableScatterers;
  bool m_enableIntermediateScatterer;

  bool m_includeVertex = false;
  bool m_successfulFit = false;  

  bool m_allFST = false;  

  // constants to set length scale for derivatives
  double lu_iSi = 16.5;    // length in u of inner silicon sensor
  double lv_iSi = 4.37062; // width  in v of inner silicon sensor
  double lu_oSi = 28;      // length in u of outer silicon sensor
  double lv_oSi = 7.67030; // width  in v of outer silicon sensor
  double lu_pent = 60.2361; // length and width of FTT pentagon modules
  double lv_pent = 60.2361;

  // Millepede Binary File for output of GBL trajectories for alignment
  gbl::MilleBinary* milleFile;
  // Minimum scattering sigma (will be squared and inverted...)
  const double scatEpsilon = 1.e-8;

  std::string rootFileName = "gbl.root";

  TFile* diag;
  TH1F* chi2OndfHistoGBL;
  TH1F* pValueHistoGBL;
  TH1F* chi2OndfHisto;
  TH1F* pValueHisto;

public:
  
  /**
   * Constructor
   */
  StFwdGbl();
  
  /**
   * Destructor
   */
  virtual ~StFwdGbl() {;}
  
  /**
   * Creates the mille binary file for output of
   * data for Millepede II alignment, can be set by setMP2Options
   */
  void beginRun();
  
  /**
   * Required to write and close ROOT file
   * with debug output. Destructor cannot be used.
   * To be called from endRun function of a module
   */
  void endRun();
  
  
  /**
   * @brief Sets internal GBL down-weighting
   * @param internalIterations GBL internal down-weighting settings, see GBL doc
   * @return void
   */
  void setGBLOptions(std::string internalIterations = "THC", bool enableScatterers = true, bool enableIntermediateScatterer = true) {
    m_gblInternalIterations = internalIterations;
    if (!enableScatterers)
      enableIntermediateScatterer = false;
    m_enableScatterers = enableScatterers;
    m_enableIntermediateScatterer = enableIntermediateScatterer;
  }
  
  /**
   * @brief Sets GBL & Millepede settings
   * @param pValueCut minimum track p-value for MP2 output
   * @param minNdf minimum track NDF for MP2 output
   * @param mille_file_name name of MP2 binary file for output
   * @return void
   */
  void setMP2Options(double pValueCut = 0., unsigned int minNdf = 1, std::string mille_file_name = "millefile.dat", double chi2Cut = 0.) {
    m_pValueCut = pValueCut;
    m_minNdf = minNdf;
    m_milleFileName = mille_file_name;
    m_chi2Cut = chi2Cut;
  }
  
  /**
   * Performs fit on a Track.
   * Hit resorting currently NOT supported.
   */
  void processTrackWithRep(genfit::Track* trk, const genfit::AbsTrackRep* rep, bool resortHits = false) ;//override;
  

  void setIncludeVertexFlag(bool flag) {m_includeVertex = flag;}
  
  void setSuccessfulFitFlag(bool flag) {m_successfulFit = flag;} 
  bool getSuccessfulFitFlag() {return m_successfulFit;} 

  void setAllFSTFlag(bool flag) {m_allFST = flag;} 
  
public:
  
  ClassDef(StFwdGbl, 1)
  
};

//}  /* End of namespace genfit */
/** @} */

#endif // STFWDGBL_H

