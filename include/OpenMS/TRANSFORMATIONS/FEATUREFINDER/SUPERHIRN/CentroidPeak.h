/*
 *  CentroidPeak.h
 *  PeakDetection
 *
 *  Created by Markus Mueller on 10/19/06.
 *  Copyright 2006 ETHZ, IMSB, Switzerland. All rights reserved.
 *  
 *  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
 *  December 2010
 *
 */

#ifndef _CENTROIDPEAK_H_
#define _CENTROIDPEAK_H_
//using namespace std;
#include <ostream>
#include <cmath>
#include <vector>
// Class for centroid peaks
class CentroidPeak{
public:
  
  static int	sfCentroidWindowWidth; // Centroid window width
  
  
  CentroidPeak();
  CentroidPeak(double,double);
  CentroidPeak(double,double, double);
  CentroidPeak(const CentroidPeak&);
  // Copy constructor
  CentroidPeak(const CentroidPeak*);
  
  CentroidPeak& operator=(const CentroidPeak&);
  
  bool operator<(const CentroidPeak&);
  
  virtual ~CentroidPeak();
  
  inline double getMass(){return fMass;}
  inline double getIntensity(){return fIntensity;}
  inline int	  getIsotopIdx(){return fIsotopIdx;}
  inline double getSignalToNoise(){return fSignalToNoise;}	
  inline double getFittedIntensity(){return fFittedIntensity;}
  inline double getOrgIntensity(){return fOrgIntensity;}
  inline std::string getExtraPeakInfo( ){ return extraPeakInfo;};
  inline double getRetentionTime( ){ return fTr;};
  
  inline void setMass(double pMass){fMass = pMass;}
  inline void setIntensity(double pIntensity){fIntensity = pIntensity;}
  inline void setIsotopIdx(double pIsotopIdx){fIsotopIdx = (int)pIsotopIdx;}
  inline void setSignalToNoise(double in){fSignalToNoise = in;}	
  inline void setFittedIntensity(double pFittedIntensity){fFittedIntensity = pFittedIntensity;}
  inline void setOrgIntensity(double pOrgIntensity){fOrgIntensity = pOrgIntensity;}
  inline void setExtraPeakInfo( std::string in){ extraPeakInfo = in;};
  inline void setRetentionTime(double in ){ fTr = in;};

  
  // shows the info of the peak:
  void show_info();
  void subtractIntensity(double);
  
protected:
  
  int		fIsotopIdx;
  double	fMass;
  double	fIntensity;
  double	fFittedIntensity;
  double	fOrgIntensity;
  double fTr;
  
  double fSignalToNoise;
  
  std::string extraPeakInfo;
};

std::ostream& operator<<(std::ostream&, CentroidPeak&);

// Class for deconvoluted isotopic patterns 
class DeconvPeak : public CentroidPeak{
public:
  
  DeconvPeak();
  DeconvPeak(double,double,int,int,double,double);
  DeconvPeak(const DeconvPeak&);
  DeconvPeak(const DeconvPeak*);
  
  DeconvPeak& operator=(const DeconvPeak&);
  
  virtual ~DeconvPeak();
  
  inline int getCharge(){return fCharge;}
  inline int getNrIsotopes(){return fNrIsotopes;}
  inline double getC13MassError(){return fC13MassError;}
  inline double getScore(){return fScore;}
  // shows the info of the peak:
  void show_info();
  
  inline std::vector<CentroidPeak> getIsotopicPeaks() {return fIsotopicPeaks;}
  
  inline void setCharge(int pCharge){fCharge = pCharge;}	
  inline void setC13MassError(double pC13MassError){fC13MassError = pC13MassError;}
  inline void setNrIsotopes(int pNrIsotopes){fNrIsotopes = pNrIsotopes;}
  inline void setScore(double pScore){fScore = pScore;}
  inline void setIsotopicPeaks(std::vector<CentroidPeak> pIsotopicPeaks) {
    fIsotopicPeaks = pIsotopicPeaks;
  }
  
protected:
  
  int						fCharge;
  int						fNrIsotopes;
  double					fC13MassError;
  double					fScore;
  std::vector<CentroidPeak>	fIsotopicPeaks;
};

std::ostream& operator<<(std::ostream&, DeconvPeak&);


#endif
