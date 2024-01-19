// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include <numeric>
#include <optional>
#include <set>
#include <string>
#include <tuple>
#include <variant>
#include <vector>

namespace tdl {

namespace detail {

//!\brief Helper struct to detect if a datastructure is a vector
template <typename T>
struct is_vector : std::false_type {};

template <typename T>
struct is_vector<std::vector<T>> : std::true_type {};

template <typename T>
inline static constexpr bool is_vector_v = is_vector<T>::value;

/*!\brief Stores a value with optional limits
 */
template <typename T, typename ListType=T>
struct TValue {
    ListType value{};
    std::optional<T> minLimit{};
    std::optional<T> maxLimit{};

    enum class State {Ok, LimitsInvalid, ValueToLow, ValueToHigh};
    auto state() const {
        if (minLimit and maxLimit and *minLimit > *maxLimit) {
            return State::LimitsInvalid;
        }
        if constexpr (is_vector_v<T>) {
            for (auto e : value) {
                if (minLimit and *minLimit > e) {
                    return State::ValueToLow;
                }
                if (maxLimit and *maxLimit < e) {
                    return State::ValueToHigh;
                }
            }
        } else {
            if (minLimit and *minLimit > value) {
                return State::ValueToLow;
            }
            if (maxLimit and *maxLimit < value) {
                return State::ValueToHigh;
            }
        }
        return State::Ok;
    }
};

/*!\brief Stores a string with optional list of valid strings.
 */
template <typename T, typename ListType=T>
struct TStringValue {
    ListType value{};
    std::optional<std::vector<T>> validValues{};

    enum class State {Valid, Invalid};
    auto state() const {
        if (validValues) {
            auto checkSingleValue = [this](T const& value) {
                return std::accumulate(begin(*validValues), end(*validValues), false, [&](auto acc, auto const& pattern) {
                    if (value == pattern) {
                        return true;
                    }
                    return acc;
                });
            };

            if constexpr (is_vector_v<T>) {
                for (auto const& e : value) {
                    auto isValid = checkSingleValue(e);
                    if (not isValid) {
                        return State::Invalid;;
                    }
                }
            } else {
                auto isValid = checkSingleValue(value);
                if (not isValid) {
                    return State::Invalid;;
                }
            }
        }
        return State::Valid;
    }
};

}

// Value types that are valid entries in the Node
using BoolValue       = bool;
using IntValue        = detail::TValue<int>;
using DoubleValue     = detail::TValue<double>;
using StringValue     = detail::TStringValue<std::string>;
using IntValueList    = detail::TValue<int, std::vector<int>>;
using DoubleValueList = detail::TValue<double, std::vector<double>>;
using StringValueList = detail::TStringValue<std::string, std::vector<std::string>>;

/*!\brief represents a parameter tree or a subtree of the parameter tree.
 *
 * This represents values that are structured in a tree and strongly typed.
 */
struct Node {
    using Children = std::vector<Node>;
    using Value = std::variant<BoolValue,          // just a single bool value
                               IntValue,           // single int, double or string value
                               DoubleValue,
                               StringValue,
                               IntValueList,       // list of int, double or string values
                               DoubleValueList,
                               StringValueList,
                               Children>;          // not a value, but a node with children

    std::string name{};           //!< Name of the entry.
    std::string description{};    //!< Entry description.
    std::set<std::string> tags{}; //!< List of tags, e.g.: advanced parameter tag.
    Value value{Children{}};      //!< Current value of this entry
};

//! A pair of mapping from tree parameter names to cli names
struct CLIMapping {
    std::string optionIdentifier; //!< full name on the command line (including '-' or '--')
    std::string referenceName;    //!< name of the option inside the parameter tree
};

//\brief Citation information of the app
struct Citation {
    std::string doi; //!\brief the doi (document object identifier)
    std::string url; //!\brief an url for direct access.
};

//!\brief Meta data of the tool
struct MetaInfo {
    std::string version{};              //!\brief version as a string
    std::string name{};                 //!\brief name of the app
    std::string docurl{};               //!\brief url to the documentation of the app
    std::string category{};             //!\brief category of the app
    std::string description{};          //!\brief a brief description of the app
    std::string executableName{};       //!\brief the actual call of this app
    std::vector<Citation> citations{};  //!\brief list publication integrated into this app
};

//! A full parameter tree document with cli mappings
struct ToolInfo {
    MetaInfo                metaInfo{};
    Node::Children          params{};
    std::vector<CLIMapping> cliMapping{};
};
}
