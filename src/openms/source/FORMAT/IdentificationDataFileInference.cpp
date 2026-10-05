// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Timo Sachsenberg $
#include "IdentificationDataFileSupport.h"

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/METADATA/DataArrays.h>
#include <bit>
#include <limits>
#include <set>

namespace OpenMS::Internal::IdentificationDataIO
{
namespace
{
  using Field = std::shared_ptr<arrow::Field>;
  Field required(const std::string& name, const std::shared_ptr<arrow::DataType>& type)
  { return arrow::field(name, type, false); }
  // Integer physical columns preserve distinct NaN payloads and signed zero even
  // when Parquet dictionary encoding is enabled for these legacy numeric values.
  void appendReal(arrow::ArrayBuilder& builder, double item)
  { append<arrow::UInt64Builder>(builder, std::bit_cast<UInt64>(item)); }
  double readReal(const arrow::Array& array, int64_t row)
  { return std::bit_cast<double>(number<arrow::UInt64Array>(array, row)); }
  std::shared_ptr<arrow::DataType> listOf(const std::shared_ptr<arrow::DataType>& type, bool nullable = false)
  { return arrow::list(arrow::field("item", type, nullable)); }
  const arrow::StructArray& structure(const arrow::Array& array, int64_t row)
  {
    if (array.IsNull(row)) invalid("Unexpected null structure");
    return static_cast<const arrow::StructArray&>(array);
  }
  struct ListView
  {
    const arrow::Array& values;
    int64_t begin;
    int64_t end;
    ListView(const arrow::Array& array, int64_t row):
        values(*static_cast<const arrow::ListArray&>(array).values()),
        begin(static_cast<const arrow::ListArray&>(array).value_offset(row)),
        end(begin + static_cast<const arrow::ListArray&>(array).value_length(row))
    {
      if (array.IsNull(row)) invalid("Unexpected null list");
    }
  };
  arrow::ArrayBuilder& beginList(arrow::ArrayBuilder& builder)
  {
    auto& list = static_cast<arrow::ListBuilder&>(builder);
    check(list.Append());
    return *list.value_builder();
  }
  arrow::StructBuilder& beginStruct(arrow::ArrayBuilder& builder)
  {
    auto& item = static_cast<arrow::StructBuilder&>(builder);
    check(item.Append());
    return item;
  }
  void requireOrdinal(UInt64 actual, UInt64 expected)
  {
    if (actual != expected) invalid("Inference child rows are missing, duplicated or out of order");
  }
  void requirePayload(Size bytes, const Options& options)
  {
    if (bytes > options.max_record_bytes) invalid("Record exceeds max_record_bytes");
  }
  std::string dateText(const DateTime& date)
  { return date.toString("yyyy-MM-ddThh:mm:ss.zzz"); }
  DateTime parseDate(const std::string& value)
  {
    DateTime date;
    if (! value.empty()) date.set(value);
    return date;
  }
  Json formulaJson(const EmpiricalFormula& formula)
  { return {{"composition", formula.toString()}, {"charge", formula.getCharge()}}; }
  EmpiricalFormula readFormulaJson(const Json& item)
  {
    EmpiricalFormula result(item.at("composition").get<std::string>());
    result.setCharge(integer<Int>(item.at("charge")));
    return result;
  }
  void validateJsonStrings(const Json& item)
  {
    if (item.is_string()) validateText(item.get_ref<const std::string&>());
    else if (item.is_array() || item.is_object())
    {
      for (const auto& child : item)
        validateJsonStrings(child);
    }
  }
  std::shared_ptr<arrow::DataType> identityType()
  { return arrow::struct_({required("database", arrow::utf8()), required("accession", arrow::utf8())}); }
  void appendIdentity(arrow::ArrayBuilder& builder, const ID::QualifiedAccession& identity)
  {
    if (identity.accession.empty()) invalid("Empty qualified parent accession");
    auto& fields = beginStruct(builder);
    appendText(*fields.field_builder(0), identity.database);
    appendText(*fields.field_builder(1), identity.accession);
  }
  ID::QualifiedAccession readIdentity(const arrow::Array& array, int64_t row)
  {
    const auto& fields = structure(array, row);
    ID::QualifiedAccession identity {text(*fields.field(0), row), text(*fields.field(1), row)};
    if (identity.accession.empty()) invalid("Empty qualified parent accession");
    return identity;
  }
  void appendAlias(arrow::ArrayBuilder& builder, const std::string& alias, const ID::InferenceResult& result, std::set<std::string>& emitted_aliases)
  {
    const auto found = result.qualified_accessions.find(alias);
    if (found == result.qualified_accessions.end()) check(builder.AppendNull());
    else
    {
      appendIdentity(builder, found->second);
      emitted_aliases.insert(alias);
    }
  }
  void readAlias(const arrow::Array& array, int64_t row, const std::string& alias, ID::InferenceResult& result)
  {
    if (array.IsNull(row)) return;
    const auto identity = readIdentity(array, row);
    const auto [found, inserted] = result.qualified_accessions.emplace(alias, identity);
    if (! inserted && found->second != identity) invalid("Conflicting qualified identities for protein alias");
  }
  std::shared_ptr<arrow::DataType> formulaType()
  { return arrow::struct_({required("composition", arrow::utf8()), required("charge", arrow::int32())}); }
  void appendFormula(arrow::ArrayBuilder& builder, const EmpiricalFormula& formula)
  {
    auto& fields = beginStruct(builder);
    appendText(*fields.field_builder(0), formula.toString());
    append<arrow::Int32Builder>(*fields.field_builder(1), formula.getCharge());
  }
  EmpiricalFormula readFormula(const arrow::Array& array, int64_t row)
  {
    const auto& fields = structure(array, row);
    EmpiricalFormula result(text(*fields.field(0), row));
    result.setCharge(number<arrow::Int32Array>(*fields.field(1), row));
    return result;
  }
  std::shared_ptr<arrow::DataType> modificationType()
  {
    return listOf(arrow::struct_({required("position", arrow::uint64()),
                                  required("id", arrow::utf8()),
                                  required("full_id", arrow::utf8()),
                                  required("psi_mod_accession", arrow::utf8()),
                                  required("unimod_record_id", arrow::int32()),
                                  required("full_name", arrow::utf8()),
                                  required("name", arrow::utf8()),
                                  required("term_specificity", arrow::uint8()),
                                  required("origin", arrow::uint8()),
                                  required("classification", arrow::uint8()),
                                  required("provenance", arrow::uint8()),
                                  required("average_mass_bits", arrow::uint64()),
                                  required("mono_mass_bits", arrow::uint64()),
                                  required("diff_average_mass_bits", arrow::uint64()),
                                  required("diff_mono_mass_bits", arrow::uint64()),
                                  required("formula", arrow::utf8()),
                                  required("diff_formula", formulaType()),
                                  required("synonyms", listOf(arrow::utf8())),
                                  required("neutral_loss_formulas", listOf(formulaType())),
                                  required("neutral_loss_mono_masses_bits", listOf(arrow::uint64())),
                                  required("neutral_loss_average_masses_bits", listOf(arrow::uint64()))}));
  }
  Size appendModifications(arrow::ArrayBuilder& builder, const ProteinHit& protein)
  {
    auto& items = beginList(builder);
    Size bytes = 0;
    for (const auto& [position, mod] : protein.getModifications())
    {
      if (static_cast<unsigned>(mod.getTermSpecificity()) >= ResidueModification::NUMBER_OF_TERM_SPECIFICITY
          || static_cast<unsigned>(mod.getSourceClassification()) >= ResidueModification::NUMBER_OF_SOURCE_CLASSIFICATIONS
          || static_cast<unsigned>(mod.getProvenance()) >= ResidueModification::NUMBER_OF_PROVENANCE)
        invalid("Unknown protein modification enum value");
      auto& fields = beginStruct(items);
      append<arrow::UInt64Builder>(*fields.field_builder(0), position);
      appendText(*fields.field_builder(1), mod.getId());
      appendText(*fields.field_builder(2), mod.getFullId());
      appendText(*fields.field_builder(3), mod.getPSIMODAccession());
      append<arrow::Int32Builder>(*fields.field_builder(4), mod.getUniModRecordId());
      appendText(*fields.field_builder(5), mod.getFullName());
      appendText(*fields.field_builder(6), mod.getName());
      append<arrow::UInt8Builder>(*fields.field_builder(7), mod.getTermSpecificity());
      append<arrow::UInt8Builder>(*fields.field_builder(8), static_cast<unsigned char>(mod.getOrigin()));
      append<arrow::UInt8Builder>(*fields.field_builder(9), mod.getSourceClassification());
      append<arrow::UInt8Builder>(*fields.field_builder(10), mod.getProvenance());
      appendReal(*fields.field_builder(11), mod.getAverageMass());
      appendReal(*fields.field_builder(12), mod.getMonoMass());
      appendReal(*fields.field_builder(13), mod.getDiffAverageMass());
      appendReal(*fields.field_builder(14), mod.getDiffMonoMass());
      appendText(*fields.field_builder(15), mod.getFormula());
      appendFormula(*fields.field_builder(16), mod.getDiffFormula());
      auto& synonyms = beginList(*fields.field_builder(17));
      for (const auto& synonym : mod.getSynonyms())
      {
        appendText(synonyms, synonym);
        bytes += synonym.size();
      }
      auto& formulas = beginList(*fields.field_builder(18));
      for (const auto& formula : mod.getNeutralLossDiffFormulas())
      {
        appendFormula(formulas, formula);
        bytes += formula.toString().size() + 8;
      }
      auto& mono = beginList(*fields.field_builder(19));
      for (double mass : mod.getNeutralLossMonoMasses())
        appendReal(mono, mass);
      auto& average = beginList(*fields.field_builder(20));
      for (double mass : mod.getNeutralLossAverageMasses())
        appendReal(average, mass);
      bytes += 128 + mod.getId().size() + mod.getFullId().size() + mod.getPSIMODAccession().size() + mod.getFullName().size() + mod.getName().size()
               + mod.getFormula().size() + mod.getDiffFormula().toString().size()
               + 8 * (mod.getNeutralLossMonoMasses().size() + mod.getNeutralLossAverageMasses().size());
    }
    return bytes;
  }
  void readModifications(const arrow::Array& array, int64_t row, ProteinHit& protein)
  {
    ListView list(array, row);
    std::set<std::pair<Size, ResidueModification>> modifications;
    for (int64_t i = list.begin; i < list.end; ++i)
    {
      const auto& fields = structure(list.values, i);
      auto field = [&](int n) -> const arrow::Array& { return *fields.field(n); };
      ResidueModification mod;
      mod.setId(text(field(1), i));
      mod.setPSIMODAccession(text(field(3), i));
      mod.setUniModRecordId(number<arrow::Int32Array>(field(4), i));
      mod.setFullName(text(field(5), i));
      mod.setName(text(field(6), i));
      const auto specificity = number<arrow::UInt8Array>(field(7), i);
      const auto classification = number<arrow::UInt8Array>(field(9), i);
      const auto provenance = number<arrow::UInt8Array>(field(10), i);
      if (specificity >= ResidueModification::NUMBER_OF_TERM_SPECIFICITY || classification >= ResidueModification::NUMBER_OF_SOURCE_CLASSIFICATIONS
          || provenance >= ResidueModification::NUMBER_OF_PROVENANCE)
        invalid("Unknown protein modification enum value");
      mod.setTermSpecificity(static_cast<ResidueModification::TermSpecificity>(specificity));
      mod.setOrigin(static_cast<char>(number<arrow::UInt8Array>(field(8), i)));
      mod.setSourceClassification(static_cast<ResidueModification::SourceClassification>(classification));
      mod.setProvenance(static_cast<ResidueModification::Provenance>(provenance));
      mod.setAverageMass(readReal(field(11), i));
      mod.setMonoMass(readReal(field(12), i));
      mod.setDiffAverageMass(readReal(field(13), i));
      mod.setDiffMonoMass(readReal(field(14), i));
      mod.setFormula(text(field(15), i));
      mod.setDiffFormula(readFormula(field(16), i));
      ListView names(field(17), i);
      std::set<std::string> synonyms;
      for (int64_t j = names.begin; j < names.end; ++j)
        if (! synonyms.insert(text(names.values, j)).second) invalid("Duplicate modification synonym");
      mod.setSynonyms(synonyms);
      ListView formulas(field(18), i);
      std::vector<EmpiricalFormula> losses;
      for (int64_t j = formulas.begin; j < formulas.end; ++j)
        losses.push_back(readFormula(formulas.values, j));
      mod.setNeutralLossDiffFormulas(losses);
      ListView mono(field(19), i);
      std::vector<double> masses;
      for (int64_t j = mono.begin; j < mono.end; ++j)
        masses.push_back(readReal(mono.values, j));
      mod.setNeutralLossMonoMasses(masses);
      ListView average(field(20), i);
      masses.clear();
      for (int64_t j = average.begin; j < average.end; ++j)
        masses.push_back(readReal(average.values, j));
      mod.setNeutralLossAverageMasses(masses);
      const auto full_id = text(field(2), i);
      // An empty stored full ID is distinct from asking ResidueModification to derive one.
      if (! full_id.empty()) mod.setFullId(full_id);
      if (! modifications.emplace(number<arrow::UInt64Array>(field(0), i), std::move(mod)).second) invalid("Duplicate protein modification");
    }
    protein.setModifications(modifications);
  }
  std::shared_ptr<arrow::DataType> dataProcessingType()
  {
    return listOf(arrow::struct_({required("software_name", arrow::utf8()), required("software_version", arrow::utf8()),
                                  required("software_metadata", metadataType()), required("actions", listOf(arrow::uint32())),
                                  required("completion_time", arrow::utf8()), required("metadata", metadataType())}),
                  true);
  }
  std::shared_ptr<arrow::DataType> dataArrayType(const std::shared_ptr<arrow::DataType>& type)
  {
    return listOf(arrow::struct_({required("name", arrow::utf8()), required("metadata", metadataType()), required("processing", dataProcessingType()),
                                  required(type->id() == arrow::Type::UINT32 ? "value_bits" : "values", listOf(type))}));
  }
  Size modificationBytes(const ProteinHit& protein)
  {
    Size bytes = 0;
    for (const auto& [position, mod] : protein.getModifications())
    {
      bytes += 128 + mod.getId().size() + mod.getFullId().size() + mod.getPSIMODAccession().size() + mod.getFullName().size() + mod.getName().size()
               + mod.getFormula().size() + mod.getDiffFormula().toString().size()
               + 8 * (mod.getNeutralLossMonoMasses().size() + mod.getNeutralLossAverageMasses().size());
      for (const auto& synonym : mod.getSynonyms())
        bytes += synonym.size();
      for (const auto& formula : mod.getNeutralLossDiffFormulas())
        bytes += formula.toString().size() + 8;
    }
    return bytes;
  }
  template<class Arrays>
  Size arrayBytes(const Arrays& arrays)
  {
    Size bytes = 0;
    for (const auto& array : arrays)
    {
      bytes += 32 + array.getName().size() + metadataBytes(array);
      for (const auto& processing : array.getDataProcessing())
      {
        if (! processing) continue;
        const auto& software = processing->getSoftware();
        bytes += 64 + software.getName().size() + software.getVersion().size() + metadataBytes(software) + metadataBytes(*processing)
                 + 4 * processing->getProcessingActions().size();
      }
      for (const auto& item : array)
      {
        if constexpr (std::is_same_v<typename Arrays::value_type::value_type, std::string>) bytes += item.size() + 4;
        else
          bytes += sizeof(item);
      }
    }
    return bytes;
  }
  template<class Arrays>
  void collectArrayMetadata(const Arrays& arrays, Dictionary& dictionary)
  {
    for (const auto& array : arrays)
    {
      dictionary.collect(array);
      for (const auto& processing : array.getDataProcessing())
      {
        if (! processing) continue;
        dictionary.collect(*processing);
        dictionary.collect(processing->getSoftware());
      }
    }
  }
  Size appendProcessing(arrow::ArrayBuilder& builder, const MetaInfoDescription& description, const Dictionary& dictionary)
  {
    auto& items = beginList(builder);
    Size bytes = 0;
    for (const auto& processing : description.getDataProcessing())
    {
      if (! processing)
      {
        check(items.AppendNull());
        continue;
      }
      const auto& software = processing->getSoftware();
      auto& fields = beginStruct(items);
      appendText(*fields.field_builder(0), software.getName());
      appendText(*fields.field_builder(1), software.getVersion());
      appendMetadata(*fields.field_builder(2), software, dictionary);
      auto& actions = beginList(*fields.field_builder(3));
      for (const auto action : processing->getProcessingActions())
      {
        if (static_cast<unsigned>(action) >= DataProcessing::SIZE_OF_PROCESSINGACTION) invalid("Invalid processing action");
        append<arrow::UInt32Builder>(actions, action);
      }
      appendText(*fields.field_builder(4), dateText(processing->getCompletionTime()));
      appendMetadata(*fields.field_builder(5), *processing, dictionary);
      bytes += 64 + software.getName().size() + software.getVersion().size() + metadataBytes(software) + metadataBytes(*processing)
               + 4 * processing->getProcessingActions().size();
    }
    return bytes;
  }
  void readProcessing(const arrow::Array& array, int64_t row, MetaInfoDescription& description, const Dictionary& dictionary)
  {
    ListView list(array, row);
    std::vector<DataProcessingPtr> values;
    for (int64_t i = list.begin; i < list.end; ++i)
    {
      if (list.values.IsNull(i))
      {
        values.push_back(nullptr);
        continue;
      }
      const auto& fields = structure(list.values, i);
      auto processing = std::make_shared<DataProcessing>();
      auto& software = processing->getSoftware();
      software.setName(text(*fields.field(0), i));
      software.setVersion(text(*fields.field(1), i));
      readMetadata(*fields.field(2), i, software, dictionary);
      ListView actions(*fields.field(3), i);
      std::set<DataProcessing::ProcessingAction> kinds;
      for (int64_t j = actions.begin; j < actions.end; ++j)
      {
        const auto action = number<arrow::UInt32Array>(actions.values, j);
        if (action >= DataProcessing::SIZE_OF_PROCESSINGACTION || ! kinds.insert(static_cast<DataProcessing::ProcessingAction>(action)).second)
          invalid("Invalid or duplicate processing action");
      }
      processing->setProcessingActions(kinds);
      processing->setCompletionTime(parseDate(text(*fields.field(4), i)));
      readMetadata(*fields.field(5), i, *processing, dictionary);
      values.push_back(std::move(processing));
    }
    description.setDataProcessing(values);
  }
  template<class Builder, class Arrays>
  Size appendArrays(arrow::ArrayBuilder& builder, const Arrays& arrays, const Dictionary& dictionary)
  {
    auto& items = beginList(builder);
    Size bytes = 0;
    for (const auto& array : arrays)
    {
      auto& fields = beginStruct(items);
      appendText(*fields.field_builder(0), array.getName());
      appendMetadata(*fields.field_builder(1), array, dictionary);
      bytes += appendProcessing(*fields.field_builder(2), array, dictionary);
      auto& values = beginList(*fields.field_builder(3));
      for (const auto& item : array)
      {
        if constexpr (std::is_same_v<Builder, arrow::StringBuilder>)
        {
          appendText(values, item);
          bytes += item.size() + 4;
        }
        else
        {
          if constexpr (std::is_same_v<Builder, arrow::FloatBuilder>) append<arrow::UInt32Builder>(values, std::bit_cast<UInt32>(item));
          else
            append<Builder>(values, item);
          bytes += sizeof(item);
        }
      }
      bytes += 32 + array.getName().size() + metadataBytes(array);
    }
    return bytes;
  }
  template<class Array, class Arrays>
  void readArrays(const arrow::Array& array, int64_t row, Arrays& arrays, const Dictionary& dictionary)
  {
    ListView list(array, row);
    for (int64_t i = list.begin; i < list.end; ++i)
    {
      const auto& fields = structure(list.values, i);
      typename Arrays::value_type output;
      output.setName(text(*fields.field(0), i));
      readMetadata(*fields.field(1), i, output, dictionary);
      readProcessing(*fields.field(2), i, output, dictionary);
      ListView values(*fields.field(3), i);
      for (int64_t j = values.begin; j < values.end; ++j)
      {
        if constexpr (std::is_same_v<Array, arrow::StringArray>) output.push_back(text(values.values, j));
        else if constexpr (std::is_same_v<Array, arrow::FloatArray>)
          output.push_back(std::bit_cast<float>(number<arrow::UInt32Array>(values.values, j)));
        else
          output.push_back(number<Array>(values.values, j));
      }
      arrays.push_back(std::move(output));
    }
  }
  std::shared_ptr<arrow::Schema> parentsSchema()
  {
    return arrow::schema({required("ordinal", arrow::uint64()), required("identity", identityType()), required("target_decoy", arrow::uint8()),
                          required("sequence", arrow::utf8()), required("description", arrow::utf8()), required("metadata", metadataType())});
  }
  std::shared_ptr<arrow::Schema> inputsSchema()
  {
    return arrow::schema({required("input_id", arrow::uint64()), required("run_uuid", arrow::utf8()), required("run_identifier", arrow::utf8()),
                          arrow::field("score_definition", arrow::uint32()), required("selection", arrow::utf8())});
  }
  std::shared_ptr<arrow::Schema> proteinsSchema()
  {
    return arrow::schema({required("protein_id", arrow::uint64()), required("alias", arrow::utf8()), arrow::field("identity", identityType()),
                          required("score_bits", arrow::uint64()), required("rank", arrow::uint32()), required("sequence", arrow::utf8()),
                          required("coverage_bits", arrow::uint64()), required("modifications", modificationType()),
                          required("metadata", metadataType())});
  }
  std::shared_ptr<arrow::Schema> groupsSchema()
  {
    const auto members = listOf(arrow::struct_({required("alias", arrow::utf8()), arrow::field("identity", identityType())}));
    return arrow::schema({required("group_id", arrow::uint64()), required("kind", arrow::uint8()), required("score_bits", arrow::uint64()),
                          required("members", members), required("float_arrays", dataArrayType(arrow::uint32())),
                          required("integer_arrays", dataArrayType(arrow::int32())), required("string_arrays", dataArrayType(arrow::utf8()))});
  }
  // Membership order and duplicate aliases are meaningful and are retained by the list.
  void readGroupMembers(const arrow::Array& array, int64_t row, ProteinIdentification::ProteinGroup& group,
                        ID::InferenceResult* result, Size bytes, const Options& options)
  {
    ListView members(array, row);
    requirePayload(bytes, options);
    for (int64_t i = members.begin; i < members.end; ++i)
    {
      const auto& fields = structure(members.values, i);
      const auto alias = text(*fields.field(0), i);
      bytes += 16 + alias.size();
      if (!fields.field(1)->IsNull(i))
      {
        const auto identity = readIdentity(*fields.field(1), i);
        bytes += identity.database.size() + identity.accession.size();
      }
      requirePayload(bytes, options);
      if (result)
      {
        readAlias(*fields.field(1), i, alias, *result);
        group.accessions.push_back(alias);
      }
    }
  }
} // namespace

Json processingJson(const ProteinIdentification& processing)
{
  const auto& sp = processing.getSearchParameters();
  const auto specificity = static_cast<unsigned>(sp.enzyme_term_specificity);
  if (static_cast<unsigned>(sp.mass_type) >= static_cast<unsigned>(ProteinIdentification::PeakMassType::SIZE_OF_PEAKMASSTYPE)
      || (specificity > 3 && specificity != 8 && specificity != 9))
    invalid("Invalid search parameter enum");
  const auto& enzyme = sp.digestion_enzyme;
  Json result = {{"identifier", processing.getIdentifier()},
                 {"search_engine", processing.getSearchEngine()},
                 {"search_engine_version", processing.getSearchEngineVersion()},
                 {"date", dateText(processing.getDateTime())},
                 {"score_type", processing.getScoreType()},
                 {"higher_better", processing.isHigherScoreBetter()},
                 {"significance_threshold_bits", std::bit_cast<UInt64>(processing.getSignificanceThreshold())},
                 {"metadata", metadataJson(processing)},
                 {"search_parameters",
                  {{"db", sp.db},
                   {"db_version", sp.db_version},
                   {"taxonomy", sp.taxonomy},
                   {"charges", sp.charges},
                   {"mass_type", sp.mass_type},
                   {"fixed_modifications", sp.fixed_modifications},
                   {"variable_modifications", sp.variable_modifications},
                   {"missed_cleavages", sp.missed_cleavages},
                   {"fragment_tolerance_bits", std::bit_cast<UInt64>(sp.fragment_mass_tolerance)},
                   {"fragment_tolerance_ppm", sp.fragment_mass_tolerance_ppm},
                   {"precursor_tolerance_bits", std::bit_cast<UInt64>(sp.precursor_mass_tolerance)},
                   {"precursor_tolerance_ppm", sp.precursor_mass_tolerance_ppm},
                   {"specificity", sp.enzyme_term_specificity},
                   {"metadata", metadataJson(sp)},
                   {"enzyme",
                    {{"name", enzyme.getName()},
                     {"regex", enzyme.getRegEx()},
                     {"synonyms", enzyme.getSynonyms()},
                     {"description", enzyme.getRegExDescription()},
                     {"n_term_gain", formulaJson(enzyme.getNTermGain())},
                     {"c_term_gain", formulaJson(enzyme.getCTermGain())},
                     {"psi_id", enzyme.getPSIID()},
                     {"xtandem_id", enzyme.getXTandemID()},
                     {"comet_id", enzyme.getCometID()},
                     {"msgf_id", enzyme.getMSGFID()},
                     {"omssa_id", enzyme.getOMSSAID()}}}}}};
  validateJsonStrings(result);
  return result;
}

ProteinIdentification readProcessingJson(const Json& json)
{
  validateJsonStrings(json);
  ProteinIdentification result;
  result.setIdentifier(json.at("identifier").get<std::string>());
  result.setSearchEngine(json.at("search_engine").get<std::string>());
  result.setSearchEngineVersion(json.at("search_engine_version").get<std::string>());
  result.setDateTime(parseDate(json.at("date").get<std::string>()));
  result.setScoreType(json.at("score_type").get<std::string>());
  result.setHigherScoreBetter(json.at("higher_better").get<bool>());
  result.setSignificanceThreshold(std::bit_cast<double>(integer<UInt64>(json.at("significance_threshold_bits"))));
  readMetadataJson(json.at("metadata"), result);
  const auto& parameters = json.at("search_parameters");
  auto& sp = result.getSearchParameters();
  sp.db = parameters.at("db").get<std::string>();
  sp.db_version = parameters.at("db_version").get<std::string>();
  sp.taxonomy = parameters.at("taxonomy").get<std::string>();
  sp.charges = parameters.at("charges").get<std::string>();
  const auto mass_type = integer<unsigned>(parameters.at("mass_type"));
  const auto specificity = integer<unsigned>(parameters.at("specificity"));
  if (mass_type >= static_cast<unsigned>(ProteinIdentification::PeakMassType::SIZE_OF_PEAKMASSTYPE)
      || (specificity > 3 && specificity != 8 && specificity != 9))
    invalid("Invalid search parameter enum");
  sp.mass_type = static_cast<ProteinIdentification::PeakMassType>(mass_type);
  sp.enzyme_term_specificity = static_cast<EnzymaticDigestion::Specificity>(specificity);
  sp.fixed_modifications = parameters.at("fixed_modifications").get<std::vector<std::string>>();
  sp.variable_modifications = parameters.at("variable_modifications").get<std::vector<std::string>>();
  sp.missed_cleavages = integer<UInt>(parameters.at("missed_cleavages"));
  sp.fragment_mass_tolerance = std::bit_cast<double>(integer<UInt64>(parameters.at("fragment_tolerance_bits")));
  sp.fragment_mass_tolerance_ppm = parameters.at("fragment_tolerance_ppm").get<bool>();
  sp.precursor_mass_tolerance = std::bit_cast<double>(integer<UInt64>(parameters.at("precursor_tolerance_bits")));
  sp.precursor_mass_tolerance_ppm = parameters.at("precursor_tolerance_ppm").get<bool>();
  readMetadataJson(parameters.at("metadata"), sp);
  const auto& enzyme = parameters.at("enzyme");
  sp.digestion_enzyme
    = Protease(enzyme.at("name").get<std::string>(), enzyme.at("regex").get<std::string>(), enzyme.at("synonyms").get<std::set<std::string>>(),
               enzyme.at("description").get<std::string>(), readFormulaJson(enzyme.at("n_term_gain")), readFormulaJson(enzyme.at("c_term_gain")),
               enzyme.at("psi_id").get<std::string>(), enzyme.at("xtandem_id").get<std::string>(), integer<Int>(enzyme.at("comet_id")),
               integer<Int>(enzyme.at("msgf_id")), integer<Int>(enzyme.at("omssa_id")));
  return result;
}

Json writeParents(const std::filesystem::path& path, const std::vector<ID::ParentRecord>& parents, Dictionary& dictionary, const Options& options)
{
  for (const auto& parent : parents)
    dictionary.collect(parent);
  TableWriter writer(path, parentsSchema(), options);
  UInt64 ordinal = 0;
  for (const auto& parent : parents)
  {
    if (static_cast<unsigned>(parent.target_decoy) > static_cast<unsigned>(ID::TargetDecoy::BOTH)) invalid("Invalid parent target/decoy value");
    append<arrow::UInt64Builder>(writer.column(0), ordinal++);
    appendIdentity(writer.column(1), parent.identity);
    append<arrow::UInt8Builder>(writer.column(2), static_cast<unsigned>(parent.target_decoy));
    appendText(writer.column(3), parent.sequence);
    appendText(writer.column(4), parent.description);
    appendMetadata(writer.column(5), parent, dictionary);
    writer.finishRow(64 + parent.identity.database.size() + parent.identity.accession.size() + parent.sequence.size() + parent.description.size()
                     + metadataBytes(parent));
  }
  writer.close();
  return writer.reference();
}

std::vector<ID::ParentRecord> readParents(const std::filesystem::path& root, const Json& reference, const Dictionary& dictionary, const Options& options)
{
  TableReader reader(root, reference, parentsSchema(), options);
  std::vector<ID::ParentRecord> parents;
  while (reader.next())
  {
    const auto row = reader.row();
    requireOrdinal(number<arrow::UInt64Array>(reader.column(0), row), parents.size());
    ID::ParentRecord parent;
    parent.identity = readIdentity(reader.column(1), row);
    const auto state = number<arrow::UInt8Array>(reader.column(2), row);
    if (state > static_cast<unsigned>(ID::TargetDecoy::BOTH)) invalid("Invalid parent target/decoy value");
    parent.target_decoy = static_cast<ID::TargetDecoy>(state);
    parent.sequence = text(reader.column(3), row);
    parent.description = text(reader.column(4), row);
    readMetadata(reader.column(5), row, parent, dictionary);
    requirePayload(64 + parent.identity.database.size() + parent.identity.accession.size() + parent.sequence.size() + parent.description.size()
                     + metadataBytes(parent),
                   options);
    parents.push_back(std::move(parent));
  }
  return parents;
}

Json writeInference(const std::filesystem::path& directory, const ID::InferenceResult& result, const Options& options)
{
  Dictionary dictionary;
  for (const auto& protein : result.proteins.getHits())
    dictionary.collect(protein);
  for (const auto* groups : {&result.proteins.getProteinGroups(), &result.proteins.getIndistinguishableProteins()})
    for (const auto& group : *groups)
    {
      collectArrayMetadata(group.getFloatDataArrays(), dictionary);
      collectArrayMetadata(group.getIntegerDataArrays(), dictionary);
      collectArrayMetadata(group.getStringDataArrays(), dictionary);
    }
  Json descriptor = {{"identifier", result.identifier},
                     {"processing", processingJson(result.proteins)},
                     {"metadata_fields", dictionary.toJson()},
                     {"parent_score", result.parent_score ? scoreJson(*result.parent_score) : Json(nullptr)},
                     {"group_score", result.group_score ? scoreJson(*result.group_score) : Json(nullptr)},
                     {"input_scores", Json::array()},
                     {"tables", Json::object()},
                     {"counts", Json::object()}};
  std::set<std::string> emitted_aliases;
  auto recordTable = [&](const std::string& name, const TableWriter& writer) {
    descriptor["tables"][name] = writer.reference();
    descriptor["counts"][name] = writer.rows();
  };
  {
    TableWriter inputs(directory / "inputs.parquet", inputsSchema(), options);
    UInt64 input_id = 0;
    for (const auto& input : result.inputs)
    {
      append<arrow::UInt64Builder>(inputs.column(0), input_id);
      appendText(inputs.column(1), input.run_uuid);
      appendText(inputs.column(2), input.run_identifier);
      if (input.score)
      {
        if (descriptor["input_scores"].size() > std::numeric_limits<UInt32>::max()) invalid("Too many inference input score definitions");
        append<arrow::UInt32Builder>(inputs.column(3), descriptor["input_scores"].size());
        descriptor["input_scores"].push_back(scoreJson(*input.score));
      }
      else
        check(inputs.column(3).AppendNull());
      appendText(inputs.column(4), input.selection);
      inputs.finishRow(64 + input.run_uuid.size() + input.run_identifier.size() + input.selection.size());
      ++input_id;
    }
    inputs.close();
    recordTable("inputs", inputs);
  }
  {
    TableWriter proteins(directory / "proteins.parquet", proteinsSchema(), options);
    UInt64 ordinal = 0;
    for (const auto& protein : result.proteins.getHits())
    {
      append<arrow::UInt64Builder>(proteins.column(0), ordinal++);
      appendText(proteins.column(1), protein.getAccession());
      appendAlias(proteins.column(2), protein.getAccession(), result, emitted_aliases);
      appendReal(proteins.column(3), protein.getScore());
      append<arrow::UInt32Builder>(proteins.column(4), protein.getRank());
      appendText(proteins.column(5), protein.getSequence());
      appendReal(proteins.column(6), protein.getCoverage());
      const Size modifications = appendModifications(proteins.column(7), protein);
      appendMetadata(proteins.column(8), protein, dictionary);
      Size bytes = 64 + protein.getAccession().size() + protein.getSequence().size() + modifications + metadataBytes(protein);
      const auto alias = result.qualified_accessions.find(protein.getAccession());
      if (alias != result.qualified_accessions.end()) bytes += alias->second.database.size() + alias->second.accession.size();
      proteins.finishRow(bytes);
    }
    proteins.close();
    recordTable("proteins", proteins);
  }
  {
    TableWriter groups(directory / "groups.parquet", groupsSchema(), options);
    UInt64 group_id = 0;
    unsigned kind = 0;
    for (const auto* collection : {&result.proteins.getProteinGroups(), &result.proteins.getIndistinguishableProteins()})
    {
      for (const auto& group : *collection)
      {
        append<arrow::UInt64Builder>(groups.column(0), group_id++);
        append<arrow::UInt8Builder>(groups.column(1), kind);
        appendReal(groups.column(2), group.probability);
        Size bytes = 32;
        bytes += appendArrays<arrow::FloatBuilder>(groups.column(4), group.getFloatDataArrays(), dictionary);
        bytes += appendArrays<arrow::Int32Builder>(groups.column(5), group.getIntegerDataArrays(), dictionary);
        bytes += appendArrays<arrow::StringBuilder>(groups.column(6), group.getStringDataArrays(), dictionary);
        auto& members = beginList(groups.column(3));
        for (const auto& alias : group.accessions)
        {
          bytes += 16 + alias.size();
          const auto identity = result.qualified_accessions.find(alias);
          if (identity != result.qualified_accessions.end()) bytes += identity->second.database.size() + identity->second.accession.size();
          requirePayload(bytes, options);
          auto& fields = beginStruct(members);
          appendText(*fields.field_builder(0), alias);
          appendAlias(*fields.field_builder(1), alias, result, emitted_aliases);
        }
        groups.finishRow(bytes);
      }
      ++kind;
    }
    groups.close();
    recordTable("groups", groups);
  }
  if (emitted_aliases.size() != result.qualified_accessions.size())
    invalid("Qualified protein alias is unused by both hits and groups; this format cannot represent detached aliases");
  validateJsonStrings(descriptor);
  return descriptor;
}

void validateParents(const std::filesystem::path& root, const Json& reference, const Dictionary& dictionary, const Options& options)
{
  TableReader reader(root, reference, parentsSchema(), options);
  UInt64 ordinal = 0;
  while (reader.next())
  {
    const auto row = reader.row();
    requireOrdinal(number<arrow::UInt64Array>(reader.column(0), row), ordinal++);
    const auto identity = readIdentity(reader.column(1), row);
    if (number<arrow::UInt8Array>(reader.column(2), row) > static_cast<unsigned>(ID::TargetDecoy::BOTH)) invalid("Invalid parent target/decoy value");
    const auto sequence = text(reader.column(3), row);
    const auto description = text(reader.column(4), row);
    MetaInfoInterface metadata;
    readMetadata(reader.column(5), row, metadata, dictionary);
    requirePayload(64 + identity.database.size() + identity.accession.size() + sequence.size() + description.size() + metadataBytes(metadata),
                   options);
  }
}

void validateInferenceTables(const std::filesystem::path& directory, const Json& descriptor, const Options& options)
{
  validateJsonStrings(descriptor);
  readProcessingJson(descriptor.at("processing"));
  if (! descriptor.at("parent_score").is_null()) readScoreJson(descriptor.at("parent_score"));
  if (! descriptor.at("group_score").is_null()) readScoreJson(descriptor.at("group_score"));
  for (const auto& score : descriptor.at("input_scores"))
    readScoreJson(score);
  Dictionary dictionary;
  dictionary.load(descriptor.at("metadata_fields"));
  auto open = [&](const std::string& name, const std::shared_ptr<arrow::Schema>& schema) {
    TableReader table(directory, descriptor.at("tables").at(name), schema, options);
    if (table.rows() != integer<UInt64>(descriptor.at("counts").at(name))) invalid("Inference row count does not match its manifest");
    return table;
  };
  auto validateUuid = [](const std::string& uuid) {
    if (uuid.size() != 36) invalid("Invalid inference run UUID");
    for (Size i = 0; i < uuid.size(); ++i)
    {
      if (i == 8 || i == 13 || i == 18 || i == 23)
      {
        if (uuid[i] != '-') invalid("Invalid inference run UUID");
      }
      else if (! ((uuid[i] >= '0' && uuid[i] <= '9') || (uuid[i] >= 'a' && uuid[i] <= 'f')))
        invalid("Invalid inference run UUID");
    }
  };
  {
    auto inputs = open("inputs", inputsSchema());
    UInt64 ordinal = 0;
    while (inputs.next())
    {
      const auto row = inputs.row();
      const auto id = number<arrow::UInt64Array>(inputs.column(0), row);
      requireOrdinal(id, ordinal++);
      const auto uuid = text(inputs.column(1), row);
      validateUuid(uuid);
      const auto identifier = text(inputs.column(2), row);
      if (const auto score = optionalNumber<arrow::UInt32Array>(inputs.column(3), row))
        if (*score >= descriptor.at("input_scores").size()) invalid("Unknown inference input score definition");
      const auto selection = text(inputs.column(4), row);
      requirePayload(64 + uuid.size() + identifier.size() + selection.size(), options);
    }
  }
  {
    auto proteins = open("proteins", proteinsSchema());
    UInt64 ordinal = 0;
    while (proteins.next())
    {
      const auto row = proteins.row();
      requireOrdinal(number<arrow::UInt64Array>(proteins.column(0), row), ordinal++);
      const auto alias = text(proteins.column(1), row);
      Size bytes = 64 + alias.size();
      if (! proteins.column(2).IsNull(row))
      {
        const auto identity = readIdentity(proteins.column(2), row);
        bytes += identity.database.size() + identity.accession.size();
      }
      ProteinHit protein;
      protein.setAccession(alias);
      if (protein.getAccession() != alias) invalid("Protein alias cannot be represented without trimming");
      protein.setScore(readReal(proteins.column(3), row));
      protein.setRank(number<arrow::UInt32Array>(proteins.column(4), row));
      const auto sequence = text(proteins.column(5), row);
      protein.setSequence(sequence);
      if (protein.getSequence() != sequence) invalid("Protein sequence cannot be represented without trimming");
      protein.setCoverage(readReal(proteins.column(6), row));
      readModifications(proteins.column(7), row, protein);
      readMetadata(proteins.column(8), row, protein, dictionary);
      requirePayload(bytes + sequence.size() + modificationBytes(protein) + metadataBytes(protein), options);
    }
  }
  {
    auto groups = open("groups", groupsSchema());
    UInt64 next_group = 0;
    unsigned previous_kind = 0;
    while (groups.next())
    {
      const auto row = groups.row();
      const auto id = number<arrow::UInt64Array>(groups.column(0), row);
      requireOrdinal(id, next_group++);
      const auto kind = number<arrow::UInt8Array>(groups.column(1), row);
      if (kind > 1 || kind < previous_kind) invalid("Invalid inference group kind or order");
      previous_kind = kind;
      ProteinIdentification::ProteinGroup group;
      group.probability = readReal(groups.column(2), row);
      readArrays<arrow::FloatArray>(groups.column(4), row, group.getFloatDataArrays(), dictionary);
      readArrays<arrow::Int32Array>(groups.column(5), row, group.getIntegerDataArrays(), dictionary);
      readArrays<arrow::StringArray>(groups.column(6), row, group.getStringDataArrays(), dictionary);
      readGroupMembers(groups.column(3), row, group, nullptr,
                       32 + arrayBytes(group.getFloatDataArrays()) + arrayBytes(group.getIntegerDataArrays()) + arrayBytes(group.getStringDataArrays()), options);
    }
  }
}

ID::InferenceResult readInference(const std::filesystem::path& directory, const Json& descriptor, const Options& options)
{
  validateJsonStrings(descriptor);
  ID::InferenceResult result;
  result.identifier = descriptor.at("identifier").get<std::string>();
  result.proteins = readProcessingJson(descriptor.at("processing"));
  if (! descriptor.at("parent_score").is_null()) result.parent_score = readScoreJson(descriptor.at("parent_score"));
  if (! descriptor.at("group_score").is_null()) result.group_score = readScoreJson(descriptor.at("group_score"));
  Dictionary dictionary;
  dictionary.load(descriptor.at("metadata_fields"));
  std::vector<ID::ScoreDefinition> scores;
  for (const auto& item : descriptor.at("input_scores"))
    scores.push_back(readScoreJson(item));
  auto reference = [&](const std::string& table) { return descriptor.at("tables").at(table); };
  auto count = [&](const std::string& table, UInt64 actual) {
    if (integer<UInt64>(descriptor.at("counts").at(table)) != actual) invalid("Inference row count does not match its manifest");
  };
  {
    TableReader inputs(directory, reference("inputs"), inputsSchema(), options);
    while (inputs.next())
    {
      const auto row = inputs.row();
      const auto input_id = number<arrow::UInt64Array>(inputs.column(0), row);
      requireOrdinal(input_id, result.inputs.size());
      ID::InferenceInput input;
      input.run_uuid = text(inputs.column(1), row);
      input.run_identifier = text(inputs.column(2), row);
      if (const auto score = optionalNumber<arrow::UInt32Array>(inputs.column(3), row))
      {
        if (*score >= scores.size()) invalid("Unknown inference input score definition");
        input.score = scores[*score];
      }
      input.selection = text(inputs.column(4), row);
      requirePayload(64 + input.run_uuid.size() + input.run_identifier.size() + input.selection.size(), options);
      result.inputs.push_back(std::move(input));
    }
    count("inputs", result.inputs.size());
  }
  {
    TableReader proteins(directory, reference("proteins"), proteinsSchema(), options);
    while (proteins.next())
    {
      const auto row = proteins.row();
      requireOrdinal(number<arrow::UInt64Array>(proteins.column(0), row), result.proteins.getHits().size());
      ProteinHit protein;
      const auto alias = text(proteins.column(1), row);
      protein.setAccession(alias);
      if (protein.getAccession() != alias) invalid("Protein alias cannot be represented without trimming");
      readAlias(proteins.column(2), row, protein.getAccession(), result);
      protein.setScore(readReal(proteins.column(3), row));
      protein.setRank(number<arrow::UInt32Array>(proteins.column(4), row));
      const auto sequence = text(proteins.column(5), row);
      protein.setSequence(sequence);
      if (protein.getSequence() != sequence) invalid("Protein sequence cannot be represented without trimming");
      protein.setCoverage(readReal(proteins.column(6), row));
      readModifications(proteins.column(7), row, protein);
      readMetadata(proteins.column(8), row, protein, dictionary);
      Size bytes = 64 + protein.getAccession().size() + protein.getSequence().size() + modificationBytes(protein) + metadataBytes(protein);
      const auto identity = result.qualified_accessions.find(protein.getAccession());
      if (identity != result.qualified_accessions.end()) bytes += identity->second.database.size() + identity->second.accession.size();
      requirePayload(bytes, options);
      result.proteins.insertHit(std::move(protein));
    }
    count("proteins", result.proteins.getHits().size());
  }
  {
    TableReader groups(directory, reference("groups"), groupsSchema(), options);
    UInt64 group_count = 0;
    unsigned previous_kind = 0;
    while (groups.next())
    {
      const auto row = groups.row();
      const auto group_id = number<arrow::UInt64Array>(groups.column(0), row);
      requireOrdinal(group_id, group_count++);
      const auto kind = number<arrow::UInt8Array>(groups.column(1), row);
      if (kind > 1 || kind < previous_kind) invalid("Invalid inference group kind or order");
      previous_kind = kind;
      ProteinIdentification::ProteinGroup group;
      group.probability = readReal(groups.column(2), row);
      readArrays<arrow::FloatArray>(groups.column(4), row, group.getFloatDataArrays(), dictionary);
      readArrays<arrow::Int32Array>(groups.column(5), row, group.getIntegerDataArrays(), dictionary);
      readArrays<arrow::StringArray>(groups.column(6), row, group.getStringDataArrays(), dictionary);
      readGroupMembers(groups.column(3), row, group, &result,
                       32 + arrayBytes(group.getFloatDataArrays()) + arrayBytes(group.getIntegerDataArrays()) + arrayBytes(group.getStringDataArrays()), options);
      if (kind == 0) result.proteins.insertProteinGroup(group);
      else
        result.proteins.insertIndistinguishableProteins(group);
    }
    count("groups", group_count);
  }
  return result;
}
} // namespace OpenMS::Internal::IdentificationDataIO
