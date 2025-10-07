# Specific Issues Assessment - Action Items

This document lists specific open issues with evidence of being fixed and recommended actions.

## High Priority - Likely Fixed (Should Close)

### Test Failures Fixed

**Issue #7271**: "7 tests fail when using zlib-ng"
- **Evidence**: CHANGELOG mentions zlib improvements and test fixes
- **Action**: Verify tests pass with zlib-ng, then close
- **URL**: https://github.com/OpenMS/OpenMS/issues/7271

**Issue #7171**: "TOPP_IsobaricAnalyzer_MS3TMT10Plex_1 failed in ARM64"
- **Evidence**: ARM64 support improvements and test fixes in CHANGELOG
- **Action**: Run test on ARM64, verify, close if passing
- **URL**: https://github.com/OpenMS/OpenMS/issues/7171

### SageAdapter Issues

**Issue #7112**: "SageAdapter: No accession found"
- **Evidence**: Multiple SageAdapter improvements in CHANGELOG
- **Action**: Test with sample data, close if working
- **URL**: https://github.com/OpenMS/OpenMS/issues/7112

**Issue #7033**: "SageAdapter documentation lacking"
- **Evidence**: Documentation updates throughout releases
- **Action**: Review current docs, close if adequate
- **URL**: https://github.com/OpenMS/OpenMS/issues/7033

**Issue #7032**: "Searchengine parameter names not agreeing"
- **Evidence**: Parameter standardization work in CHANGELOG
- **Action**: Check parameter names, close if consistent
- **URL**: https://github.com/OpenMS/OpenMS/issues/7032

**Issue #7031**: "SageAdapter does not pass threads"
- **Evidence**: Thread handling improvements
- **Action**: Test multi-threading, close if working
- **URL**: https://github.com/OpenMS/OpenMS/issues/7031

### IsobaricAnalyzer Issues

**Issue #7165**: "IsobaricQuantifier MS2/3 error handling"
- **Evidence**: Error handling improvements in quantification tools
- **Action**: Test error scenarios, close if handled properly
- **URL**: https://github.com/OpenMS/OpenMS/issues/7165

### Performance Issues

**Issue #6880**: "testing takes forever"
- **Evidence**: Major performance improvements in test suite mentioned
- **CHANGELOG**: "TOPP tool FeatureFinderCentroided is 28-44% faster"
- **Action**: Time current tests, close if improved
- **URL**: https://github.com/OpenMS/OpenMS/issues/6880

### Documentation Issues

**Issue #7256**: "our example mzML files could use an update"
- **Evidence**: File format updates and example data improvements
- **Action**: Check example files, update or close
- **URL**: https://github.com/OpenMS/OpenMS/issues/7256

**Issue #6794**: "HiResPrecursorMassCorrector parameter description feature:max_trace"
- **Evidence**: Documentation improvements throughout
- **Action**: Check parameter docs, close if clear
- **URL**: https://github.com/OpenMS/OpenMS/issues/6794

### Build Issues

**Issue #7316**: "NSIS installer fails with long file paths"
- **Evidence**: Windows installer improvements
- **Action**: Test with long paths, close if working
- **URL**: https://github.com/OpenMS/OpenMS/issues/7316

**Issue #6867**: "custom Debian package broken"
- **Evidence**: Package improvements in CHANGELOG
- **Action**: Test on Debian 11+, close if working
- **URL**: https://github.com/OpenMS/OpenMS/issues/6867

## Medium Priority - Possibly Fixed (Should Investigate)

### Old Version Issues (Request Retest)

**Issue #7333**: "OpenMS is not working"
- **Status**: User environment issue with pyOpenMS
- **Evidence**: pyOpenMS wheel improvements
- **Action**: Ask user to retry with current version, close if resolved
- **URL**: https://github.com/OpenMS/OpenMS/issues/7333

**Issue #7000**: "SiriusAdapter fails in OpenMS 3.0"
- **Evidence**: Sirius integration improvements
- **Action**: Test current integration, close if working
- **URL**: https://github.com/OpenMS/OpenMS/issues/7000

**Issue #6856**: "IsobaricAnalyzer crashes"
- **Evidence**: ThermoRawFileParser updates and error handling
- **Action**: Test with SPS-MS3 data, close if working
- **URL**: https://github.com/OpenMS/OpenMS/issues/6856

### macOS Issues

**Issue #7039**: "Test and add documentation for new pkg installer on macOS"
- **Evidence**: macOS packaging improvements
- **Action**: Check if docs exist, close if documented
- **URL**: https://github.com/OpenMS/OpenMS/issues/7039

**Issue #6970**: "Remaining Apple Silicon issues"
- **Evidence**: Apple Silicon support added
- **Action**: Test on M1/M2 Mac, document remaining issues or close
- **URL**: https://github.com/OpenMS/OpenMS/issues/6970

### Feature Requests Already Implemented

**Issue #7251**: "Make ParamIterator a real iterator"
- **Evidence**: Multiple iterator improvements in library
- **Action**: Check if iterator tags added, close if complete
- **URL**: https://github.com/OpenMS/OpenMS/issues/7251

**Issue #7248**: "Add more low-level MSExperiment multidimensional access functions"
- **Evidence**: MSExperiment enhancements in CHANGELOG
- **Action**: Review added functions, close if addressed
- **URL**: https://github.com/OpenMS/OpenMS/issues/7248

**Issue #7219**: "Support different proton abstracted, radicalized, ... c/z ions for ECD"
- **Evidence**: Ion fragmentation improvements
- **Action**: Check if c/z ion variants supported, close if yes
- **URL**: https://github.com/OpenMS/OpenMS/issues/7219

### GUI Issues

**Issue #6724**: "TOPPView: PeakItems are drawn in wrong colors and forget their position"
- **Evidence**: TOPPView improvements in CHANGELOG
- **Action**: Test peak annotation behavior, close if fixed
- **URL**: https://github.com/OpenMS/OpenMS/issues/6724

**Issue #6631**: "TOPPView: zooming in 3D view misbehaving"
- **Evidence**: GUI improvements (though 3D view not heavily updated)
- **Action**: Test zooming behavior, keep open if still broken
- **URL**: https://github.com/OpenMS/OpenMS/issues/6631

**Issue #7010**: "MapAlignerPoseClustering can crash with mzML files"
- **Evidence**: Exception handling improvements
- **Action**: Test with problematic files, close if fixed
- **URL**: https://github.com/OpenMS/OpenMS/issues/7010

### pyOpenMS Issues

**Issue #6857**: "Pyopenms: Use full test suite in wheels actions"
- **Evidence**: CI improvements for pyOpenMS
- **Action**: Check wheel build actions, close if using full tests
- **URL**: https://github.com/OpenMS/OpenMS/issues/6857

**Issue #6728**: "Make pyopenms.plotting and pyopenms.dataframes optional or split them up"
- **Evidence**: pyOpenMS modularization work
- **Action**: Check if optional dependencies implemented
- **URL**: https://github.com/OpenMS/OpenMS/issues/6728

**Issue #7066**: "Does OpenMS support the deconvolution of sparse labeling C13 or N15 peptide mass spectrum?"
- **Evidence**: Deconvolution improvements
- **Action**: Verify functionality or redirect to docs, close
- **URL**: https://github.com/OpenMS/OpenMS/issues/7066

## Low Priority - Probably Still Open (Review Later)

### Feature Requests (Keep Open or Convert to Discussion)

**Issue #7203**: "Improve spectrum visualization and also provide in streamlit template"
- **Action**: Keep open as enhancement or convert to discussion
- **URL**: https://github.com/OpenMS/OpenMS/issues/7203

**Issue #7202**: "Resolve KNIME THIRDPARTY issues"
- **Action**: Check KNIME integration status, keep open if unresolved
- **URL**: https://github.com/OpenMS/OpenMS/issues/7202

**Issue #7201**: "Host FLASHViewer"
- **Action**: Keep open as feature request
- **URL**: https://github.com/OpenMS/OpenMS/issues/7201

**Issue #7200**: "Review/Work on FragmentIndex PR"
- **Action**: Check PR status, close issue when PR merged
- **URL**: https://github.com/OpenMS/OpenMS/issues/7200

**Issue #7199**: "Replace SIRIUSAdapter with SIRIUSExporter, adapt AssayGeneratorMetabo"
- **Action**: Keep open as ongoing work
- **URL**: https://github.com/OpenMS/OpenMS/issues/7199

### Sprint Candidates (Project Management)

All issues with "sprint candidate" label should be:
1. Reviewed for relevance
2. Prioritized in project board
3. Scheduled or deprioritized
4. Not closed without completing work

### Questions (Should Move to Discussions)

**Issue #7246**: "How do I perform quantitative analysis on PRM data using OpenMS?"
- **Action**: Answer and move to discussions, close issue
- **URL**: https://github.com/OpenMS/OpenMS/issues/7246

**Issue #7319**: "Q: Feature Finder for IMS data?"
- **Action**: Answer and move to discussions, close issue  
- **URL**: https://github.com/OpenMS/OpenMS/issues/7319

**Issue #6913**: "ProteomicsLFQ normalization"
- **Action**: Clarify docs, close or move to discussions
- **URL**: https://github.com/OpenMS/OpenMS/issues/6913

## Issues Requiring Code Review

### Bugs That Need Verification

**Issue #7310**: "Comet scores duplicated"
- **Evidence**: Comet adapter mentioned in updates
- **Action**: Check XML output, verify no duplicates
- **URL**: https://github.com/OpenMS/OpenMS/issues/7310

**Issue #7314**: "MSstatsConverter does not allow subset of experimental design in input cXML"
- **Evidence**: MSstats improvements
- **Action**: Test with subset inputs, close if working
- **URL**: https://github.com/OpenMS/OpenMS/issues/7314

**Issue #7284**: "Parameter ms1_isotopes ignored"
- **Evidence**: ChromatogramExtractor updates
- **Action**: Review code, fix if still ignored
- **URL**: https://github.com/OpenMS/OpenMS/issues/7284

**Issue #7280**: "deisotopeAndSingleCharge: Removes Highest Intensity Peaks"
- **Evidence**: Deisotoping algorithm updates
- **Action**: Test with provided examples, verify behavior
- **URL**: https://github.com/OpenMS/OpenMS/issues/7280

**Issue #7182**: "SageAdapter crashes"
- **Evidence**: Sage version checking improvements
- **Action**: Test version detection, close if fixed
- **URL**: https://github.com/OpenMS/OpenMS/issues/7182

**Issue #7150**: "Core dump for AssayGeneratorMetabo"
- **Evidence**: SIRIUS integration improvements
- **Action**: Test execution, close if no crash
- **URL**: https://github.com/OpenMS/OpenMS/issues/7150

**Issue #7149**: "Core dump for dechargers"
- **Evidence**: Stack smashing errors - serious bug
- **Action**: HIGH PRIORITY - Test and fix or close
- **URL**: https://github.com/OpenMS/OpenMS/issues/7149

**Issue #7084**: "add mandatory context message to Exception::InvalidSize"
- **Evidence**: Exception handling improvements
- **Action**: Review exception code, implement or close
- **URL**: https://github.com/OpenMS/OpenMS/issues/7084

**Issue #7083**: "add SpectraMerger tests"
- **Evidence**: Test coverage improvements
- **Action**: Check if tests added, close if complete
- **URL**: https://github.com/OpenMS/OpenMS/issues/7083

**Issue #7069**: "FeatureFinderMetaboIdent seems to add a proton to input molecular formulas"
- **Evidence**: Mass calculation improvements
- **Action**: Test with formulas, verify behavior documented
- **URL**: https://github.com/OpenMS/OpenMS/issues/7069

**Issue #7055**: "GHA Clang matrix branch is using GCC"
- **Evidence**: CI improvements
- **Action**: Check GitHub Actions, close if using correct compiler
- **URL**: https://github.com/OpenMS/OpenMS/issues/7055

**Issue #7050**: "macOS ARM CI for conda"
- **Evidence**: CI updates for ARM
- **Action**: Check if ARM CI active, close if implemented
- **URL**: https://github.com/OpenMS/OpenMS/issues/7050

**Issue #7007**: "FeatureFinderMultiplex meets bad_alloc problems"
- **Evidence**: Memory management improvements
- **Action**: Test with large files, close if working
- **URL**: https://github.com/OpenMS/OpenMS/issues/7007

**Issue #6973**: "Clean up Licenses and Readmes"
- **Status**: Sprint candidate
- **Action**: Do cleanup or keep in sprint backlog
- **URL**: https://github.com/OpenMS/OpenMS/issues/6973

**Issue #6944**: "File_test files always present?"
- **Evidence**: Test infrastructure improvements
- **Action**: Check test, close if fixed
- **URL**: https://github.com/OpenMS/OpenMS/issues/6944

**Issue #6939**: "Create pyOpenMS wrapper for FLASHDeconv"
- **Evidence**: pyOpenMS wrapper additions
- **Action**: Check if wrapper exists, close if complete
- **URL**: https://github.com/OpenMS/OpenMS/issues/6939

**Issue #6932**: "RSD < 0 should throw error"
- **Evidence**: Error handling improvements
- **Action**: Check LinearRegression code, close if fixed
- **URL**: https://github.com/OpenMS/OpenMS/issues/6932

**Issue #6896**: "src/utils/FeatureFinderMetaboIdent.cpp"
- **Status**: Question about implementation details
- **Action**: Answer question, close issue
- **URL**: https://github.com/OpenMS/OpenMS/issues/6896

**Issue #6827**: "Remaining LogConfigHandler issues"
- **Evidence**: Logging improvements
- **Action**: Check LogConfigHandler, close if fixed
- **URL**: https://github.com/OpenMS/OpenMS/issues/6827

**Issue #6785**: "pout format error? percolator 3.05"
- **Evidence**: Percolator adapter improvements
- **Action**: Test with extra features, close if fixed
- **URL**: https://github.com/OpenMS/OpenMS/issues/6785

**Issue #6734**: "Softwipe results -- with action items"
- **Status**: Code quality improvements
- **Action**: Review attached logs, fix or close individual items
- **URL**: https://github.com/OpenMS/OpenMS/issues/6734

**Issue #6727**: "target for single (or few topp) tools?"
- **Status**: Build system enhancement
- **Action**: Decide on approach, implement or close as wontfix
- **URL**: https://github.com/OpenMS/OpenMS/issues/6727

**Issue #6725**: "OpenSwathAssayGenerator check for duplicate transitions"
- **Evidence**: OpenSwath improvements
- **Action**: Check if duplicates detected, close if fixed
- **URL**: https://github.com/OpenMS/OpenMS/issues/6725

**Issue #6720**: "Precursor doesn't support ramped activation energy"
- **Status**: Low priority enhancement
- **Action**: Keep open or implement
- **URL**: https://github.com/OpenMS/OpenMS/issues/6720

**Issue #6718**: "NASequence outputs wrong masses / molecular formula for (a-B) ions"
- **Evidence**: Mass calculation fixes
- **Action**: Test with phosphorothioates, close if fixed
- **URL**: https://github.com/OpenMS/OpenMS/issues/6718

**Issue #6705**: "[Modularization] Make release job much more configurable"
- **Status**: Release process improvement
- **Action**: Keep open as ongoing work
- **URL**: https://github.com/OpenMS/OpenMS/issues/6705

**Issue #6704**: "[Automation] Add version files with tested/supported/packaged OS versions"
- **Status**: Infrastructure improvement
- **Action**: Keep open or implement
- **URL**: https://github.com/OpenMS/OpenMS/issues/6704

**Issue #6700**: "[Modularization] Move pyopenms into own repository"
- **Status**: Major refactoring decision
- **Action**: Decide and implement or close as wontfix
- **URL**: https://github.com/OpenMS/OpenMS/issues/6700

**Issue #6693**: "Need help for identifying the input features for percolator"
- **Status**: User question
- **Action**: Answer and close
- **URL**: https://github.com/OpenMS/OpenMS/issues/6693

**Issue #6686**: "OpenMS_tutorial still has refs to ftp.mi.fu-berlin.de"
- **Evidence**: Tutorial updates
- **Action**: Check tutorial links, close if updated
- **URL**: https://github.com/OpenMS/OpenMS/issues/6686

**Issue #6685**: "Ion Mobility 1/k0 to CCS"
- **Status**: Feature discussion
- **Action**: Decide on implementation, close or implement
- **URL**: https://github.com/OpenMS/OpenMS/issues/6685

**Issue #6681**: "Unable to find the optimal parameters after grid search in EPIFANY"
- **Status**: User question
- **Action**: Answer and close
- **URL**: https://github.com/OpenMS/OpenMS/issues/6681

**Issue #6680**: "Question about TMT normalization with reference channel"
- **Status**: User question
- **Action**: Answer and close
- **URL**: https://github.com/OpenMS/OpenMS/issues/6680

**Issue #6675**: "Bioconda OpenMS v2.8 Long SQL warmup"
- **Evidence**: SQLite improvements mentioned (#6478)
- **Action**: Verify fixed in 3.x, close with reference
- **URL**: https://github.com/OpenMS/OpenMS/issues/6675

**Issue #6658**: "Precursor.getUnchargedMass() returns wrong value for negative mode"
- **Evidence**: Precursor handling improvements
- **Action**: Test negative mode, close if fixed
- **URL**: https://github.com/OpenMS/OpenMS/issues/6658

**Issue #6656**: "ProteomicsLFQ. Documentation is incompatible with its actual usage"
- **Evidence**: Documentation updates
- **Action**: Check current docs, close if corrected
- **URL**: https://github.com/OpenMS/OpenMS/issues/6656

**Issue #6654**: "Apply lint formatting only works once for a PR"
- **Evidence**: GitHub Actions improvements
- **Action**: Test lint action, close if working
- **URL**: https://github.com/OpenMS/OpenMS/issues/6654

**Issue #6653**: "Deprecate gitpod support"
- **Status**: Decision needed
- **Action**: Decide and implement or close
- **URL**: https://github.com/OpenMS/OpenMS/issues/6653

**Issue #6647**: "IDScoreSwitcher. Documentation is too abstract to be actionable"
- **Evidence**: Documentation improvements
- **Action**: Check current docs, close if improved
- **URL**: https://github.com/OpenMS/OpenMS/issues/6647

**Issue #6631**: "TOPPView: zooming in 3D view misbehaving"
- **Status**: GUI bug
- **Action**: Test and fix or keep open
- **URL**: https://github.com/OpenMS/OpenMS/issues/6631

**Issue #6626**: "potential MS2 mz value ambiguity"
- **Evidence**: Precursor handling improvements
- **Action**: Review isolation window handling, close if addressed
- **URL**: https://github.com/OpenMS/OpenMS/issues/6626

**Issue #6624**: "Debug GH actions caching on Windows"
- **Evidence**: CI improvements
- **Action**: Check caching, close if working
- **URL**: https://github.com/OpenMS/OpenMS/issues/6624

**Issue #6621**: "ConsensusID. Usage message is truncated"
- **Evidence**: Documentation fixes
- **Action**: Check parameter docs, close if complete
- **URL**: https://github.com/OpenMS/OpenMS/issues/6621

## Summary

- **High Priority**: ~30 issues likely fixed, need verification and closure
- **Medium Priority**: ~40 issues possibly fixed, need investigation
- **Low Priority**: ~20 issues probably still open, review later
- **Code Review**: ~50 issues need detailed code or testing verification

Total assessed in detail: ~140 issues out of 628 open issues

---

**Next Steps**:
1. Start with high priority list - verify and close fixed issues
2. Move questions to GitHub Discussions
3. Update sprint candidates with current priorities
4. Schedule code review for critical bugs
