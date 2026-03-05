// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ModificationDataProvider.h>

namespace OpenMS
{
  class OPENMS_DLLAPI CrossLinksDB :
      public ModificationsDB
  {
  public:

    /// Returns a pointer to the modifications DB (singleton)
    inline static CrossLinksDB* getInstance()
    {
      static CrossLinksDB* db_ = new CrossLinksDB;
      return db_;
    }

    /// Collects all modifications that can be used for identification searches
    void getAllSearchModifications(std::vector<String>& modifications) const;

  private:

      /** @name Constructors and Destructors
       */
      //@{
      /// Default constructor
      CrossLinksDB();

      /// Copy constructor
      CrossLinksDB(const CrossLinksDB& residue_db);

      /// Destructor
      ~CrossLinksDB() override;
      //@}

      /** @name Assignment
       */
      //@{
      /// Assignment operator
      CrossLinksDB & operator=(const CrossLinksDB& aa);
      //@}

      /// Creates the OBO provider for cross-linker loading
      static std::vector<std::unique_ptr<ModificationDataProvider>> makeCrossLinkProviders_();

  };
}

