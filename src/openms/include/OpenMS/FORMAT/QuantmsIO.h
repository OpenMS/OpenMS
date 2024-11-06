// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------
#include "OpenMS/KERNEL/ConsensusMap.h"
//#include "OpenMS/KERNEL/FeatureMap.h"
#include "OpenMS/KERNEL/MSExperiment.h"
#include "arrow/io/api.h"
#include "arrow/table.h"
#include "arrow/builder.h"
#include "arrow/ipc/writer.h"
#include "parquet/api/reader.h"
#include "parquet/arrow/reader.h"
#include "parquet/arrow/writer.h"
#include "arrow/array.h"
#include "arrow/scalar.h"

#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

namespace OpenMS
{
    class QuantmsIOWriter
    {
    public:
        QuantmsIOWriter();
        ~QuantmsIOWriter();
        // Can write PSM table
        void write(const ProteinIdentification& prot_ids,
                   const std::vector<PeptideIdentification>& pep_ids,
                   const std::string& output_root,
                   bool partitioned = false);

        // Can write PSM table with spectra
        void write(const ProteinIdentification& prot_ids,
                   const std::vector<PeptideIdentification>& pep_ids,
                   const MSExperiment& exp,
                   const std::string& output_root,
                   bool partitioned = false);
                   
        // Can write psm, feature, and protein group tables
        void write(ConsensusMap& consensus_map, const ProteinIdentification& proteins, const std::string &output_root);

        // Can write psm, feature, and protein group tables for a single file
        //void write(FeatureMap& consensus_map, const ProteinIdentification& proteins, const std::string &output_root);
    };
}