// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Hendrik Weisser $
// $Authors: Julianus Pfeuffer, Oliver Alka, Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/ID/IdentificationData.h>

#include <QtSql/QSqlDatabase>
#include <QtSql/QSqlError>

namespace OpenMS
{
  class IdentificationData;
  /**
      @brief This class supports reading and writing of OMS files.

      OMS files are SQLite databases consisting of several tables.
  */
  class OPENMS_DLLAPI OMSFile
  {

  public:

    /** @brief Write out an IdentificationData object to SQL-based OMS file
     *
     * @param filename The output file
     * @param id_data The IdentificationData object
     *
    */
    static void store(const String& filename, const IdentificationData& id_data);

    /** @brief Read in a OMS file and construct an IdentificationData
     *
     * @param filename The input file
     * @param id_data The IdentificationData object
     *
    */
    static void load(const String& filename, IdentificationData& id_data);

  protected:

    using Key = qint64;

    // convenience functions:

    static bool tableExists_(QSqlDatabase& db, const String& name);

    static void raiseDBError_(const QSqlError& error, QSqlDatabase& db,
                              int line, const char* function,
                              const String& context);

    // store helper functions:

    static void createTable_(const String& name, const String& definition,
                             QSqlDatabase& db, bool may_exist = false);

    static void storeVersionAndDate_(QSqlDatabase& db);

    static void createTableDataValue_(QSqlDatabase& db);

    static Key storeDataValue_(const DataValue& value, QSqlDatabase& db);

    static void createTableCVTerm_(QSqlDatabase& db);

    static Key storeCVTerm_(const CVTerm& cv_term, QSqlDatabase& db);

    static void storeScoreTypes_(const IdentificationData& id_data,
                                 QSqlDatabase& db);

    static void storeDataProcessingSoftwares_(const IdentificationData& id_data,
                                              QSqlDatabase& db);

    static void storeParentMolecules_(const IdentificationData& id_data,
                                      QSqlDatabase& db);

    // load helper functions:

    // static CVTerm loadCVTerm_(int id, QSqlDatabase& db);

    static void loadScoreTypes_(IdentificationData& id_data, QSqlDatabase& db);

    static void loadParentMolecules_(IdentificationData& id_data,
                                     QSqlDatabase& db);
  };
} // namespace OpenMS


