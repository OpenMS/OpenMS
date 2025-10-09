# Boolean Flag Migration Guide

This guide explains how to migrate TOPP tools with boolean string parameters (true/false restrictions) to use proper command-line flags.

## Background

String parameters with `true/false` restrictions behave inconsistently:
- Parameters with default=`false` and order `true,false` work as flags (no argument needed)
- Parameters with default=`true` or order `false,true` require arguments
- This confuses users and violates the principle of least surprise

## Migration Pattern

### For Parameters with Default = "false"

Convert directly to a flag:

**Before:**
```cpp
registerStringOption_("my_param", "<choice>", "false", "Enable feature X", false);
setValidStrings_("my_param", {"true", "false"});
```

**After:**
```cpp
registerFlag_("my_param", "Enable feature X", false);
```

**Usage:**
```cpp
// Before
bool enable_x = getStringOption_("my_param") == "true";

// After
bool enable_x = getFlag_("my_param");
```

### For Parameters with Default = "true"

Create a negated flag and keep old parameter for backward compatibility:

**Before:**
```cpp
registerStringOption_("enable_feature", "<choice>", "true", "Enable feature X", false);
setValidStrings_("enable_feature", {"true", "false"});
```

**After:**
```cpp
registerFlag_("disable_feature", "Disable feature X (enabled by default)", false);
registerStringOption_("enable_feature", "<choice>", "true", 
                      "(DEPRECATED: Use -disable_feature flag instead) Enable feature X", 
                      false, true);  // Note: true = advanced parameter
setValidStrings_("enable_feature", {"true", "false"});
```

**Usage:**
```cpp
// Before
bool enable_x = getStringOption_("enable_feature") == "true";

// After
bool enable_x = true;  // default
if (getFlag_("disable_feature"))
{
  enable_x = false;
}
else if (!getParam_("enable_feature").isEmpty())
{
  // Old parameter was used
  enable_x = getStringOption_("enable_feature") == "true";
  writeLogWarn_("Warning: Parameter 'enable_feature' is deprecated. Use flag '-disable_feature' instead.");
}
```

## Complete Example

See `src/topp/TICCalculator.cpp` for a complete example:

**Registration:**
```cpp
registerFlag_("skip_data", "Do not load and decode binary data", false);
registerStringOption_("loadData", "<method>", "true", 
                      "(DEPRECATED: Use -skip_data flag instead) Whether to load binary data", 
                      false, true);
setValidStrings_("loadData", {"true", "false"});
```

**Usage:**
```cpp
bool load_data = true;  // default
if (getFlag_("skip_data"))
{
  load_data = false;
}
else if (!getParam_("loadData").isEmpty())
{
  load_data = getStringOption_("loadData") == "true";
  writeLogWarn_("Warning: Parameter 'loadData' is deprecated. Use flag '-skip_data' instead.");
}
```

## Naming Conventions

For negated flags, use clear prefixes/patterns:
- `disable_X` for `enable_X`
- `skip_X` for `load_X` or `do_X`
- `no_X` for `X`
- `keep_X` for `remove_X`

Examples:
- `enable_ms1 → disable_ms1`
- `loadData → skip_data`
- `indexed_file → no_index`
- `remove_proteins → keep_proteins`

## Tools Already Updated

- TICCalculator
- FileFilter
- OpenSwathWorkflow
- IDMapper
- CometAdapter
- MSGFPlusAdapter
- FalseDiscoveryRate

## Tools Still Needing Updates

Run this command to find remaining boolean string parameters:

```bash
grep -rn "setValidStrings_.*true.*false" src/topp/ | grep -v "DEPRECATED"
```

Known tools with remaining work:
- ProteinInference (protein_fdr, conservative_fdr, picked_fdr)
- Epifany (protein_fdr, conservative_fdr, picked_fdr)
- ProteomicsLFQ (multiple parameters)
- FileConverter (write_scan_index)
- OpenSwathMzMLFileCacher (lossy_compression, full_meta)
- MapAlignerTreeGuided (copy_data)
- OpenSwathDecoyGenerator (switchKR)
- IDMerger (annotate_file_origin)
- ProteinQuantifier (greedy_group_resolution, file_and_channel_level_output)
- OpenNuXL (NuXL:decoys)

## Testing

For each updated tool, test:

1. **New syntax works:**
   ```bash
   Tool -new_flag
   ```

2. **Old syntax shows deprecation warning:**
   ```bash
   Tool -old_param false
   # Should show: "Warning: Parameter 'old_param' is deprecated..."
   ```

3. **INI file compatibility:**
   - Generate INI with `-write_ini`
   - Modify old parameter in INI
   - Run with `-ini` - should work with warning

4. **CTD output remains compatible:**
   ```bash
   Tool -write_ctd /tmp
   # Check that CTD still has string with true/false restrictions
   ```

## Notes

- The CTD output intentionally keeps string parameters with true/false for ParamXML 1.6.2 compatibility
- The deprecation warning is shown by TOPPBase automatically when old syntax is used
- Old parameters are marked as "advanced" to hide them from default help but keep them accessible
- This migration maintains full backward compatibility while guiding users to the new syntax
