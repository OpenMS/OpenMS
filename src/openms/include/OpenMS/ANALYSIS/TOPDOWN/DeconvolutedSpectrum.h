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

  public:
    typedef FLASHDeconvHelperStructs::Parameter Parameter;
    //typedef FLASHDeconvHelperStructs::PeakGroup PeakGroup;

    explicit DeconvolutedSpectrum(MSSpectrum &s);

    ~DeconvolutedSpectrum();

    void write(std::fstream &fs, Parameter &param);

    void writeTopFD(std::fstream &fs, int id);

    static void writeHeader(std::fstream &fs, int &n, bool detail);

    bool empty() const;
    MSSpectrum *spec;
    DeconvolutedSpectrum *precursorSpectrum;


    std::vector<PeakGroup> peakGroups;
    int specIndex, massCntr, scanNumber;
  };
}

