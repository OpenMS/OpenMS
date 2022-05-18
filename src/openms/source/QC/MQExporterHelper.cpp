// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow$
// $Authors: Virginia Rossow, Lenny Kovac, Hendrik Beschorner$
// --------------------------------------------------------------------------

#include <OpenMS/QC/MQExporterHelper.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/KERNEL/MSExperiment.h>


#include <QtCore/QDir>
#include <cmath> // isnan
#include <fstream>

using namespace OpenMS;



Size MQExporterHelper::proteinGroupID_(std::map<OpenMS::String, OpenMS::Size>& protein_id,
                                       const String& protein_accession)
{
  auto it = protein_id.find(protein_accession);
  if (it == protein_id.end())
  {
    protein_id.emplace(protein_accession, protein_id.size() + 1);
    return protein_id.size();
  }
  else
  {
    return it->second;
  }
}

std::map<Size, Size> MQExporterHelper::makeFeatureUIDtoConsensusMapIndex_(const ConsensusMap& cmap)
{
  std::map<Size, Size> f_to_ci;
  for (Size i = 0; i < cmap.size(); ++i)
  {
    for (const auto& fh : cmap[i].getFeatures())
    {
      auto[it, was_created_newly] = f_to_ci.emplace(fh.getUniqueId(), i);
      if (!was_created_newly)
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Adding [" + String(it->first) + "," + String(it->second) +  "] failed. FeatureHandle exists twice in ConsensusMap!");
      }
      f_to_ci[fh.getUniqueId()] = i;
    }
  }
  return f_to_ci;
}


bool MQExporterHelper::isValid(std::string filename)
{
  return File::writable(filename);
}