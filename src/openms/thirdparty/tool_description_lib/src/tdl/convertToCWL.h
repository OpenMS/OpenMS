// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "ToolInfo.h"
#include "cwl_v1_2.h"

#include <cassert>
#include <filesystem>
#include <functional>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>


#if (defined(_WIN32) || defined(WIN32)) && !defined(__GNUG__)
#   // Maybe visual studio will get there one day to support the c++ standard...
#   // Until then we have to live with this:
#   define and &&
#   define or ||
#   define not !
#endif


namespace tdl {

//!\brief overload structure allowing fancy 'std::visit` syntax
template<class... Ts> struct overloaded : Ts... { using Ts::operator()...; };

//!\brief required deduction guide for c++17 (not required for c++20)
template<class... Ts> overloaded(Ts...) -> overloaded<Ts...>;

namespace cwl = https___w3id_org_cwl_cwl;

/**!\brief a global callback function to adjust the exporting for cwl
 *
 * This callback allows to adjust the exported yaml file to add/change/remove
 * cwl entries, which currently aren't controllable via tdl itself.
 */
inline std::function<void(YAML::Node&)> post_process_cwl;



namespace detail {

inline auto simplifyType(YAML::Node type) -> YAML::Node {
    // 1. Collapsing optional scalar types into one option
    if (type.IsSequence() and type.size() == 2) {
        if (type[0].IsScalar()
            and type[0].as<std::string>() == "null"
            and type[1].IsScalar()) {
            type = type[1].as<std::string>() + "?";
        }
    }

    // 2. Collapsing array types into one option
    if (type.IsMap()
        and type["type"].as<std::string>("") == "array"
        and type["items"].IsScalar()) {
        type = type["items"].as<std::string>() + "[]";
    }

    // 3. Collapsing optional array types into one option
    if (type.IsSequence() and type.size() == 2) {
        if (type[0].IsScalar()
            and type[0].as<std::string>() == "null"
            and type[1].IsMap()
            and type[1]["type"].as<std::string>() == "array"
            and type[1]["items"].IsScalar()) {
            type = type[1]["items"].as<std::string>() + "[]?";
        }
    }
    return type;
}

inline auto findCLIMapping(std::string const& referenceName, ToolInfo const& doc) -> CLIMapping const* {
    for (auto const& mapping : doc.cliMapping) {
        if (mapping.referenceName == referenceName) {
            return &mapping;
        }
    }
    return nullptr;
};

template <typename InputType>
void setIdOrName(InputType& input, std::string name) {
    if constexpr (std::is_same_v<InputType, cwl::CommandInputRecordField>
                    || std::is_same_v<InputType, cwl::CommandOutputRecordField>) {
        input.name         = name;
    } else {
        input.id           = name;
    }
}
template<typename InputType = cwl::CommandInputParameter, typename OutputType = cwl::CommandOutputParameter, size_t deep = 5, typename InputCB, typename OutputCB, typename BaseCommandCB>
inline void f(Node::Children const& children, ToolInfo const& doc, InputCB const& inputCB, OutputCB const& outputCB, BaseCommandCB const& baseCommandCB) {
    if constexpr (deep > 0) {
    for (auto child : children) {
        // find CLIMapping
        auto cliMapping = findCLIMapping(child.name, doc);

        auto addInput = [&](auto type) {
            auto input         = InputType{};
            setIdOrName(input, child.name);
            if (child.tags.count("required")) {
                input.type         = type;
            } else {
                using namespace https___w3id_org_cwl_cwl;
                input.type         = std::vector<std::variant<CWLType, CommandInputRecordSchema, CommandInputEnumSchema, CommandInputArraySchema, std::string>>{cwl::CWLType::null, type};
            }
            input.doc          = child.description;
            if (cliMapping) {
                auto binding       = cwl::CommandLineBinding{};
                binding.prefix     = cliMapping->optionIdentifier;
                input.inputBinding = binding;
            }
            inputCB(std::move(input));
        };
        auto addInputArray = [&](auto type) {
            auto input     = InputType{};
            setIdOrName(input, child.name);
            auto arrayType = cwl::CommandInputArraySchema{};
            arrayType.items = type;
            if (child.tags.count("required")) {
                input.type = arrayType;
            } else {
                using namespace https___w3id_org_cwl_cwl;
                input.type         = std::vector<std::variant<CWLType, CommandInputRecordSchema, CommandInputEnumSchema, CommandInputArraySchema, std::string>>{cwl::CWLType::null, arrayType};
            }
            input.doc          = child.description;

            if (cliMapping) {
                auto binding       = cwl::CommandLineBinding{};
                binding.prefix     = cliMapping->optionIdentifier;
                input.inputBinding = binding;
            }
            inputCB(std::move(input));
        };
        auto addOutput = [&](auto type) {
            auto input         = InputType{};
            setIdOrName(input, child.name);
            if (child.tags.count("required")) {
                input.type         = cwl::CWLType::string;
            } else {
                using namespace https___w3id_org_cwl_cwl;
                input.type         = std::vector<std::variant<CWLType, CommandInputRecordSchema, CommandInputEnumSchema, CommandInputArraySchema, std::string>>{cwl::CWLType::null, cwl::CWLType::string};
            }
            input.doc          = child.description;

            if (cliMapping) {
                auto binding       = cwl::CommandLineBinding{};
                binding.prefix     = cliMapping->optionIdentifier;
                input.inputBinding = binding;
            }
            inputCB(std::move(input));

            auto output         = OutputType{};
            setIdOrName(output, child.name);
            if (child.tags.count("required")) {
                output.type         = type;
            } else {
                using namespace https___w3id_org_cwl_cwl;
                output.type         = std::vector<std::variant<CWLType, CommandOutputRecordSchema, CommandOutputEnumSchema, CommandOutputArraySchema, std::string>>{cwl::CWLType::null, type};
            }
            auto binding         = cwl::CommandOutputBinding{};
            binding.glob         = "$(inputs." + child.name + ")";
            output.outputBinding = binding;

            outputCB(std::move(output));
        };
        auto addOutputPrefixed = [&](auto type, bool multipleFiles) {
            auto input         = InputType{};
            setIdOrName(input, child.name);
            if (child.tags.count("required")) {
                input.type         = cwl::CWLType::string;
            } else {
                using namespace https___w3id_org_cwl_cwl;
                input.type         = std::vector<std::variant<CWLType, CommandInputRecordSchema, CommandInputEnumSchema, CommandInputArraySchema, std::string>>{cwl::CWLType::null, cwl::CWLType::string};
            }
            input.doc          = child.description;

            if (cliMapping) {
                auto binding       = cwl::CommandLineBinding{};
                binding.prefix     = cliMapping->optionIdentifier;
                input.inputBinding = binding;
            }
            inputCB(std::move(input));

            auto output         = OutputType{};
            setIdOrName(output, child.name);

            // Setting the correct value type
            output.type         = type;
            using namespace https___w3id_org_cwl_cwl;
            // Add an array, if a list of files
            if (multipleFiles) {
                auto arrayType  = cwl::CommandOutputArraySchema{};
                arrayType.items = type;
                output.type     = arrayType;
            }

            // Add a null, if not required and an array, if a list of files
            if (!child.tags.count("required")) {
                output.type         = std::vector<std::variant<CWLType, CommandOutputRecordSchema, CommandOutputEnumSchema, CommandOutputArraySchema, std::string>>{cwl::CWLType::null, type };

                if (multipleFiles) {
                    auto arrayType  = cwl::CommandOutputArraySchema{};
                    arrayType.items = type;
                    output.type     = std::vector<std::variant<CWLType, CommandOutputRecordSchema, CommandOutputEnumSchema, CommandOutputArraySchema, std::string>>{cwl::CWLType::null, arrayType };
                }
            }
            auto binding         = cwl::CommandOutputBinding{};
            binding.glob         = "$(inputs." + child.name + ")*";
            output.outputBinding = binding;
            outputCB(std::move(output));
        };

        std::visit(overloaded{
            [&](BoolValue const&) {
                addInput(cwl::CWLType::boolean);
            },
            [&](IntValue const&) {
                addInput(cwl::CWLType::long_);
            },
            [&](DoubleValue const&) {
                addInput(cwl::CWLType::double_);
            },
            [&](StringValue const&) {
                if (child.tags.count("output")) {
                    if (child.tags.count("file")) {
                        addOutput(cwl::CWLType::File);
                    } else if (child.tags.count("directory")) {
                        addOutput(cwl::CWLType::Directory);
                    } else if (child.tags.count("prefixed")) {
                        addOutputPrefixed(cwl::CWLType::File, /*.mutliplieFiles = */ false);
                    }
                } else if (child.tags.count("file")) {
                    addInput(cwl::CWLType::File);
                } else if (child.tags.count("directory")) {
                    addInput(cwl::CWLType::Directory);
                } else {
                    addInput(cwl::CWLType::string);
                }
            },
            [&](IntValueList const&) {
                addInputArray(cwl::CWLType::long_);
            },
            [&](DoubleValueList const&) {
                addInputArray(cwl::CWLType::double_);
            },
            [&](StringValueList const&) {
                if (child.tags.count("output")) {
                    if (child.tags.count("prefixed")) {
                        addOutputPrefixed(cwl::CWLType::File, /*.multipleFiles =*/ true);
                    } else {
                        //!TODO not implemented
                    }
                } else if (child.tags.count("file")) {
                    addInputArray(cwl::CWLType::File);
                } else if (child.tags.count("directory")) {
                    addInputArray(cwl::CWLType::Directory);
                } else {
                    addInputArray(cwl::CWLType::string);
                }
            },
            [&](Node::Children const& v) {
                if (child.tags.count("basecommand")) {
                    baseCommandCB(child.name);

                    f<InputType, OutputType, deep-1>(v, doc, inputCB, outputCB, baseCommandCB);
                    return;
                }

                auto inputs  = std::vector<cwl::CommandInputRecordField>{};
                auto outputs = std::vector<cwl::CommandOutputRecordField>{};
                f<cwl::CommandInputRecordField, cwl::CommandOutputRecordField, (deep-1)>(v, doc, [&](auto input) {
                    inputs.push_back(input);
                } , [&](auto output) {
                    outputs.push_back(output);
                }, baseCommandCB);

                auto inputType  = cwl::CommandInputRecordSchema{};
                auto outputType = cwl::CommandOutputRecordSchema{};

                inputType.fields  = inputs;
                outputType.fields = outputs;
//                child.tags.insert("required"); //!TODO, is this required?
                addInput(inputType);
            },
        }, child.value);
    }
    }
}
}

/*!\brief converts a ToolInfo into a string that
 * holds the CWL representation of the given tool
 */
inline auto convertToCWL(ToolInfo const& doc) -> std::string {
    auto& tool_info = doc.metaInfo;
    auto const schema_location = std::string{"/SCHEMAS/Param_1_7_0.xsd"};
    auto const schema_version  = std::string{"1.7.0"};


    auto tool = cwl::CommandLineTool{};
    tool.cwlVersion  = cwl::CWLVersion::v1_2;
    tool.label       = tool_info.name;
    tool.doc         = tool_info.description;
    // = tool_info.category; //!TODO
    // = tool_info.docurl; //!TODO
    // = tool_info.version; //!TODO


    //!TODO Add citation information
    // for (auto& [doi, url] : tool_info.citations) {
    //     citationNode.children.push_back({/*.tag = */"citation",
    //                                      /*.attr = */{{"doi", doi}, {"url", url}}});
    // }
    // toolNode.children.push_back(citationNode);
    auto baseCommand = std::vector<std::string>{};
    baseCommand.push_back(std::filesystem::path{tool_info.executableName}.filename().string());

    detail::f(doc.params, doc, [&](auto input) {
        tool.inputs->push_back(std::move(input));
    }, [&](auto output) {
        tool.outputs->push_back(std::move(output));
    }, [&](auto command) {
        baseCommand.push_back(command);
    });

    tool.baseCommand = baseCommand;

    auto y = toYaml(tool);

    // function to traverse yaml tree and executes 'simplifyType' on all nodes with name 'type'
    auto traverseTree = std::function<void(YAML::Node)>{};
    traverseTree = [&traverseTree](YAML::Node node) {
        if (node.IsMap()) {
            for (auto n : node) {
                if (n.first.as<std::string>("") == "type") {
                    n.second = detail::simplifyType(n.second);
                    traverseTree(n.second);
                } else {
                    traverseTree(n.second);
                }
            }
        } else if (node.IsSequence()) {
            for (auto n : node) {
                traverseTree(n);
            }
        }
    };
    // Post procssing inputs and outputs of the yaml object
    for (auto param : {"inputs", "outputs"}) {
        traverseTree(y[param]);
    }

    // post process generated cwl yaml file
    if (post_process_cwl) {
        post_process_cwl(y);
    }

    YAML::Emitter out;
    out << y;
    return out.c_str();
}

}

#if (defined(_WIN32) || defined(WIN32)) && !defined(__GNUG__)
#undef and
#undef or
#undef not
#endif
