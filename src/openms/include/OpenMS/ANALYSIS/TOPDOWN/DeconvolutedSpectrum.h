//
// Created by Kyowon Jeong on 4/22/20.
//
#pragma once

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <iomanip>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvHelperStructs.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/ANALYSIS/TOPDOWN/QScore.h>

namespace OpenMS
{
  class PeakGroup;
  class OPENMS_DLLAPI DeconvolutedSpectrum
  {

  public: // TODO protect and public devide..
    typedef FLASHDeconvHelperStructs::Parameter Parameter;
    typedef FLASHDeconvHelperStructs::LogMzPeak LogMzPeak;
    //typedef FLASHDeconvHelperStructs::hash_LogMzPeak hash_LogMzPeak;

    DeconvolutedSpectrum();

    explicit DeconvolutedSpectrum(MSSpectrum &s, int n);

    ~DeconvolutedSpectrum();

    void writeDeconvolutedMasses(std::fstream &fs, Parameter &param);
    void writeAttCsv(std::fstream &fs, int msLevel, double qScoreThreshold, int numMaxMS2);
    void writeMassList(std::fstream &fs, double retDelta, double qScoreThreshold, int numMaxMS2);

    void writeTopFD(std::fstream &fs, int id);
    MSSpectrum toSpectrum();

    static void writeDeconvolutedMassesHeader(std::fstream &fs, int &n, bool detail);
    static void writeAttCsvHeader(std::fstream &fs);
    static void writeThermoInclusionHeader(std::fstream &fs);

    bool empty() const;
    void clearPeakGroupsChargeInfo();

    MSSpectrum *spec;
    std::vector<PeakGroup> peakGroups;
    //std::vector<LogMzPeak> peaks;
    //std::vector<double> mzs; // sorted peaks from the original spectrum

    PeakGroup *precursorPeakGroup = nullptr;
    Precursor precursorPeak;
    //double originalPrecursorIntensity = 0;
    std::string activationMethod;

    bool registerPrecursor(DeconvolutedSpectrum &precursorSpectrum);

    int scanNumber, precursorScanNumber;
  };
}

