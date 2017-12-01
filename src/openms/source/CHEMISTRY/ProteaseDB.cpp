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
// $Maintainer: Xiao Liang  $
// $Authors: Xiao Liang, Chris Bielow $
// --------------------------------------------------------------------------
//

#include <OpenMS/CHEMISTRY/ProteaseDB.h>

using namespace std;

namespace OpenMS
{
  ProteaseDB::ProteaseDB():
    DigestionEnzymeDB<DigestionEnzymeProtein, ProteaseDB>("CHEMISTRY/Enzymes.xml")
  {
  }

  void ProteaseDB::getAllXTandemNames(vector<String>& all_names) const
  {
    all_names.clear();
    for (ConstEnzymeIterator it = const_enzymes_.begin(); it != const_enzymes_.end(); ++it)
    {
      if ((*it)->getXTandemID() != "")
      {
        all_names.push_back((*it)->getName());
      }
    }
  }

  void ProteaseDB::getAllCometNames(vector<String>& all_names) const
  {
    all_names.clear();
    all_names.push_back("unspecific cleavage");
    for (ConstEnzymeIterator it = const_enzymes_.begin(); it != const_enzymes_.end(); ++it)
    {
      if ((*it)->getCometID() != 0)
      {
        all_names.push_back((*it)->getName());
      }
    }
  }

  void ProteaseDB::getAllOMSSANames(vector<String>& all_names) const
  {
    all_names.clear();
    all_names.push_back("Trypsin");
    for (ConstEnzymeIterator it = const_enzymes_.begin(); it != const_enzymes_.end(); ++it)
    {
      if ((*it)->getOMSSAID() != 0)
      {
        all_names.push_back((*it)->getName());
      }
    }
  }

  void ProteaseDB::getAllMSGFNames(vector<String>& all_names) const
  {
    all_names.clear();
    for (ConstEnzymeIterator it = const_enzymes_.begin(); it != const_enzymes_.end(); ++it)
    {
      if ((*it)->getMSGFID() != -1) // MS-GF+ starts enzyme numbering at 0
      {
        all_names.push_back((*it)->getName());
      }
    }
  }

} // namespace OpenMS
