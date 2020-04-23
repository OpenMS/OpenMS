//
// Created by Kyowon Jeong on 4/22/20.
//
#pragma once

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvHelperStructs.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>

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

    explicit DeconvolutedSpectrum(MSSpectrum &s);

    ~DeconvolutedSpectrum();

    void writeDeconvolutedMasses(std::fstream &fs, Parameter &param);
    void writeAttCsv(std::fstream &fs, int msLevel);
    void writeTopFD(std::fstream &fs, int id);

    static void writeDeconvolutedMassesHeader(std::fstream &fs, int &n, bool detail);
    static void writeAttCsvHeader(std::fstream &fs);

    bool empty() const;
    MSSpectrum *spec;
    std::vector<PeakGroup> peakGroups;
    std::vector<LogMzPeak> peaks;
    //std::vector<double> mzs; // sorted peaks from the original spectrum

    PeakGroup *precursorPeakGroup = nullptr;
    LogMzPeak *precursorPeak = nullptr;
    std::string activationMethod;

    bool registerPrecursor(DeconvolutedSpectrum &precursorSpectrum);

    int specIndex, massCntr, scanNumber;


    //fs << "MS_ONE_ID=" << pg.precursorSpecIndex << "\n"
    //       << "MS_ONE_SCAN=" << pg.precursorScanNumber << "\n"
    //       << "PRECURSOR_MZ="
    //       << std::to_string(pg.precursorMz) << "\n"
    //       << "PRECURSOR_CHARGE=" << pg.precursorCharge << "\n"
    //       << "PRECURSOR_MASS=" << std::to_string(pg.precursorMonoMass) << "\n"
    //       << "PRECURSOR_INTENSITY=" << pg.precursorIntensity << "\n";*/

  };
}

