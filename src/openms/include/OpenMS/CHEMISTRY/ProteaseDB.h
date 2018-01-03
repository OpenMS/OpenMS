// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Xiao Liang $
// $Authors: Xiao Liang, Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_CHEMISTRY_PROTEASEDB_H
#define OPENMS_CHEMISTRY_PROTEASEDB_H

#include <OpenMS/CHEMISTRY/DigestionEnzymeDB.h>
#include <OpenMS/CHEMISTRY/DigestionEnzymeProtein.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <vector>

namespace OpenMS
{
  /**
    @ingroup Chemistry

    @brief Database for enzymes that digest proteins (proteases)

    The enzymes stored in this DB are defined in an XML file under share/CHEMISTRY/Enzymes.xml.
  */
  class OPENMS_DLLAPI ProteaseDB: public DigestionEnzymeDB<DigestionEnzymeProtein, ProteaseDB>
  {
    // allow access to constructor in DigestionEnzymeDB::getInstance():
    friend class DigestionEnzymeDB<DigestionEnzymeProtein, ProteaseDB>;

  protected:
    /// constructor
    ProteaseDB();

  public:
    /// returns all the enzyme names available for XTandem
    void getAllXTandemNames(std::vector<String>& all_names) const;

    /// returns all the enzyme names available for Comet
    void getAllCometNames(std::vector<String>& all_names) const;

    /// returns all the enzyme names available for Crux
    void getAllCruxNames(std::vector<String>& all_names) const;

    /// returns all the enzyme names available for OMSSA
    void getAllOMSSANames(std::vector<String>& all_names) const;

    /// returns all the enzyme names available for MSGFPlus
    void getAllMSGFNames(std::vector<String>& all_names) const;
  };
}

#endif // OPENMS_CHEMISTRY_PROTEASEDB_H

