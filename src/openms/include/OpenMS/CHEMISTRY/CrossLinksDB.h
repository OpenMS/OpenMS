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
  /**
    @brief Process-wide singleton database of cross-linking modifications.

    CrossLinksDB is a specialization of @ref ModificationsDB that loads
    cross-linker definitions from the PSI XLMOD ontology
    (@c CHEMISTRY/XLMOD.obo, shipped with OpenMS) instead of the
    general-purpose UniMod/PSI-MOD modifications. It exposes the same
    lookup, iteration, and matching interface as the base class but its
    underlying pool contains @em only entries flagged as cross-links by
    the OBO provider.

    @section CrossLinksDB_singleton Singleton model

    Like @ref ModificationsDB, CrossLinksDB is a singleton: there is at
    most one instance per process. Obtain it via #getInstance — the
    constructors, destructor, and assignment operator are intentionally
    @c private so callers cannot create additional copies. Loading of the
    XLMOD ontology happens once on the first call to #getInstance.

    @section CrossLinksDB_typical_use Typical use

    @code
      auto* db = CrossLinksDB::getInstance();
      // Use any ModificationsDB API to look up cross-linkers by name,
      // origin, mass, or accession:
      const ResidueModification* mod = db->getModification("DSS");

      // Or enumerate the searchable subset for ID-search UIs:
      std::vector<String> ids;
      db->getAllSearchModifications(ids);
    @endcode

    Currently in OpenMS this database is used by @c MzIdentMLHandler
    when serialising cross-link identifications, by the pyOpenMS chemistry
    bindings, and by tests; new cross-link search algorithms should pull
    their candidate modification list from here rather than from the
    general ModificationsDB.

    @ingroup Chemistry
    @see ModificationsDB
  */
  class OPENMS_DLLAPI CrossLinksDB :
      public ModificationsDB
  {
  public:

    /**
      @brief Returns the process-wide singleton instance.

      First call lazily constructs the instance and populates it from
      @c CHEMISTRY/XLMOD.obo. Subsequent calls return the same pointer.
      The returned pointer is owned by the database and must @em not be
      deleted by the caller.
    */
    inline static CrossLinksDB* getInstance()
    {
      static CrossLinksDB* db_ = new CrossLinksDB;
      return db_;
    }

    /**
      @brief Returns the IDs of all cross-linker modifications that are
      usable for identification searches (i.e. carry a PSI-MOD accession).

      The list is sorted ascending and contains @ref
      ResidueModification::getFullId values, suitable for offering to the
      user in a search-engine GUI dropdown or for passing to a
      cross-link-aware search engine.

      @param[out] modifications Sorted list of search-eligible modification
                                full-IDs. Any previous contents are cleared.
    */
    void getAllSearchModifications(std::vector<String>& modifications) const;

  private:

      /** @name Constructors and Destructors
       */
      //@{
      /// @brief Private default constructor: parses the XLMOD ontology
      /// and forwards the resulting providers to @ref ModificationsDB.
      /// Called once via #getInstance.
      CrossLinksDB();

      /// @brief Private copy constructor (singleton — copying is not
      /// permitted).
      CrossLinksDB(const CrossLinksDB& residue_db);

      /// @brief Destructor. Never invoked in practice because the
      /// singleton is allocated with @c new and outlives the process.
      ~CrossLinksDB() override;
      //@}

      /** @name Assignment
       */
      //@{
      /// @brief Private assignment operator (singleton — assignment is
      /// not permitted).
      CrossLinksDB & operator=(const CrossLinksDB& aa);
      //@}

      /// @brief Builds the OBO data provider list used by the base-class
      /// constructor. Loads @c CHEMISTRY/XLMOD.obo with the
      /// @c cross_links_only flag set.
      static std::vector<std::unique_ptr<ModificationDataProvider>> makeCrossLinkProviders_();

  };
}

