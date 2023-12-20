<!-- SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin -->
<!-- SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik -->
<!-- SPDX-License-Identifier: CC0-1.0 -->
```
#include <tdl/tdl.h>

int main(int argc, char** argv) {

  auto toolInfo = tdl::ToolInfo {
    .meta = {
      .version = "7.6.5",                             // version of your tool
      .name    = "testApp",                           // name of your tool
      .docurl  = "example.com",                       // url to the documentation of your tool
      .category = "sequence analysis",                //!TODO what are good categories?
      .description = "awesome tool! does everything", // a description of your tool
      .executableName = argv[0],                      // executable name of your tool, we advise to fill it with what ever is given in argv[0]

       // A list of citations you would like people to look at when they use your tool for research
      .citations      = {{"doi:123", "https://en.wikipedia.org/wiki/Meaning_of_life"},
                         {"doi:456", "https://en.wikipedia.org/wiki/Turing_completeness"}},
    }
    .params = {
      // allow a call like `./myTool build --input <file1> <file2>`
      tdl::Node {
        .name        = "build",
        .description = "builds some index for search",
        .tags        = {}, // no tags
        .value       = tdl::Node::Children { tdl::Node {
            .name        = "input", // This must be the same as the referenceName of CLIMapping
            .description = "input file",
            .value       = tdl::StringValueList{}, // indicates that we are accepting a list of strings as values
        }}
      },

      // allow a call like `./myTool search --index myindex.db --queries <file1> <file2>`
      tdl::Node {
        .name        = "search",
        .description = "using index to search",
        .value       = tdl::Node::Children {
          tdl::Node {
            .name        = "queries", // This must be the same as the referenceName of CLIMapping
            .description = "files with search queries",
            .value       = tdl::StringValueList{}, // indicates that we are accepting a list of strings as values
          }, tdl::Node {
            .name        = "index", // This must be the same as the referenceName of CLIMapping
            .description = "path to an index file",
            .value       = tdl::StringValue{}, // indicates that we are accepting a single string value
          }
        }
     }
   },
   .cliMapping = {
       tdl::CLIMapping {.optionIdentifier = "--input",   .referenceName = "input"},
       tdl::CLIMapping {.optionIdentifier = "--index",   .referenceName = "index"},
       tdl::CLIMapping {.optionIdentifier = "--queries", .referenceName = "queries"},
  });
};
```

