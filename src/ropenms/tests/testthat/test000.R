library(purrr)
context("Running main test for pyopenms R bindings")

# perform MetaInfoInterface tests for subclass ----------------------------

check_MetaInfoInterface = function(x) {

  x <- RichPeak2D$new()
  x$setMetaValue("key", 42)
  x$setMetaValue("key2", 42)

  # both character() & list() work
  keys = character()
  x$getKeys(keys)
  expect_true(is.character(keys),length(keys) > 0)
  expect_equivalent(x$getMetaValue(keys[1]),42L)
  expect_true(x$metaValueExists("key"))
  x$removeMetaValue("key")

  keys = list()
  x$getKeys(keys)
  expect_equivalent(x$getMetaValue(keys[1]),42L)

  x$clearMetaInfo()
  keys = character()
  x$getKeys(keys)
  expect_length(keys,0)

}

# on_failure(check_MetaInfoInterface) <- function(call,env){
#   paste0(deparse(class(call$x)[1])," does not satisfy MetaInfoInterface tests")
# }

# Test MetaInfoInterface --------------------------------------------------

test_that("test MetaInfoInterface",{

  md <- MetaInfoInterface$new()
  md$setMetaValue("key", 42)
  md$setMetaValue("key2", 42)

  # both character() or list() would work
  keys = character()
  md$getKeys(keys)
  expect_true(is.character(keys),length(keys) > 0)
  expect_equivalent(md$getMetaValue(keys[1]),42L)
  expect_true(md$metaValueExists("key"))
  md$removeMetaValue("key")

  keys = list()
  md$getKeys(keys)
  expect_equivalent(md$getMetaValue(keys[1]),42L)

  md$clearMetaInfo()
  keys = character()
  md$getKeys(keys)
  expect_length(keys,0)

})

# perform UniqueIdInterface tests for subclass--------------------------------------------------
check_UniqueIdInterface <- function(x){

  expect_true(as.logical(x$hasInvalidUniqueId()))
  expect_false(as.logical(x$hasValidUniqueId()))
  expect_true(as.logical(x$ensureUniqueId()))
  expect_true(is_scalar_integer(x$getUniqueId()))

  # tentative as UniqueId is 64 bit signed/unsigned integer
  # if we attempt to convert 64 bit integer from python -> R, -1 is returned
  # so, casting all 64 bit int result to character.

  expect_true(as.character(py_call(x$.__enclos_env__$private$py_obj$getUniqueId)) > 0)
  expect_false(as.logical(x$hasInvalidUniqueId()))
  expect_true(as.logical(x$hasValidUniqueId()))

  x$clearUniqueId()
  expect_equal(x$getUniqueId(),0)
  expect_true(as.logical(x$hasInvalidUniqueId()))
  expect_false(as.logical(x$hasValidUniqueId()))

  expect_true(as.logical(x$ensureUniqueId()))
  expect_true(is_scalar_integer(x$getUniqueId()))
  # See above
  expect_true(as.character(py_call(x$.__enclos_env__$private$py_obj$getUniqueId)) > 0)
  expect_false(as.logical(x$hasInvalidUniqueId()))
  expect_true(as.logical(x$hasValidUniqueId()))

  x$setUniqueId(1234)
  expect_equivalent(x$getUniqueId(),1234L)

}

# on_failure(check_UniqueIdInterface) <- function(call,env){
#   paste0(deparse(class(call$x)[1])," does not satisfy UniqueIdInterface tests")
# }


# DataProcessing check function -------------------------------------------

check_DataProcessing <- function(dp){

  check_MetaInfoInterface(dp)

  expect_true(dp==dp)
  expect_false(dp!=dp)

  dp$clearMetaInfo()
  k <- list()
  dp$getKeys(k)
  expect_equal(k,list())

  dp$getMetaValue
  ac = dp$getProcessingActions()
  expect_false(T %in% duplicated(ac))

  ac = c(ProcessingAction()$PEAK_PICKING, ProcessingAction()$BASELINE_REDUCTION)
  dp$setProcessingActions(ac)

  expect_length(dp$getProcessingActions(), 2)
  expect_true(is_scalar_character(dp$getSoftware()$getName()))

  expect_true(is_scalar_character(dp$getSoftware()$getVersion()))

  dp$isMetaEmpty()
  dp$metaValueExists
  dp$removeMetaValue
  s = dp$getSoftware()
  s$setName("pyopenms")
  dp$setSoftware(s)

  expect_equal(dp$getSoftware()$getName(), "pyopenms")

}

# on_failure(check_DataProcessing) <- function(call, env){
#   paste0(deparse(class(call$x)[1])," does not satisfy DataProcessing tests")
# }

# Test ProgressLogger -----------------------------------------------------

test_that("test ProgressLogger",{

  ff <- MzXMLFile$new()
  # Change here : Enum
  ff$setLogType(LogType()$NONE)
  expect_equal(ff$getLogType(),LogType()$NONE)
  ff$startProgress(0,3, "label")
  ff$setProgress(0)
  ff$setProgress(1)
  ff$setProgress(2)
  ff$setProgress(3)
  ff$endProgress()

})

check_ProgressLogger <- function(ff){

  # Change here: Enum
  ff$setLogType(LogType()$NONE)
  expect_equal(ff$getLogType(),LogType()$NONE)
  ff$startProgress(0,3, "label")
  ff$setProgress(0)
  ff$setProgress(1)
  ff$setProgress(2)
  ff$setProgress(3)
  ff$endProgress()

}

# on_failure(check_ProgressLogger) <- function(call,env){
#   paste0(deparse(class(call$x)[1])," does not satisfy ProgressLogger tests")
# }

# Test SpectrumAlignment --------------------------------------------------

test_that("test SpectrumAlignment",{

  spec <- MSSpectrum$new()
  p <- Peak1D$new()
  p$setMZ(1000.0)
  p$setIntensity(200.0)
  spec$push_back(p)
  p$setMZ(2000.0)
  p$setIntensity(200.0)
  spec$push_back(p)

  rich_spec <- MSSpectrum$new()
  p = Peak1D$new()
  p$setMZ(1000.001)
  p$setIntensity(200.0)
  rich_spec$push_back(p)
  p$setMZ(2000.001)
  p$setIntensity(200.0)
  rich_spec$push_back(p)
  p$setMZ(3000.001)
  p$setIntensity(200.0)
  rich_spec$push_back(p)

  aligner <- SpectrumAlignment$new()
  result <- list()

  aligner$getSpectrumAlignment(result, spec, spec)
  expect_equivalent(result,list(c(0L,0L),c(1L,1L)))
  aligner$getSpectrumAlignment(result, rich_spec, spec)
  expect_equivalent(result,list(c(0L,0L),c(1L,1L)))
  aligner$getSpectrumAlignment(result,spec,rich_spec)
  expect_equivalent(result,list(c(0L,0L),c(1L,1L)))
  aligner$getSpectrumAlignment(result, rich_spec, rich_spec)
  expect_equivalent(result,list(c(0L,0L),c(1L,1L),c(2L,2L)))

  aligner <- SpectrumAlignmentScore$new()

  # N.B. To provide alternative for "aligner(spec)"
  # & "aligner(spec,spec)" two functions score_1,
  # score_2 are added as addons.

  expect_true(is_scalar_double(aligner$score_1(spec)))
  expect_true(is_scalar_double(aligner$score_1(rich_spec)))

  expect_true(is_scalar_double(aligner$score_2(spec, rich_spec)))
  expect_true(is_scalar_double(aligner$score_2(rich_spec,spec)))

  expect_true(is_scalar_double(aligner$score_2(spec,spec)))
  expect_true(is_scalar_double(aligner$score_2(rich_spec,rich_spec)))
})


# Test AASequence ---------------------------------------------------------

test_that("test AASequence",{
  aas <- AASequence$new()

  # += operator is not supported
  aas = aas + aas

  aas = AASequence$fromString("DFPIANGER")
  expect_equal(aas$getCTerminalModificationName(),"")
  expect_equal(aas$getNTerminalModificationName(),"")
  aas$setCTerminalModification("")
  aas$setNTerminalModification("")
  expect_equal(aas$toString(),"DFPIANGER")
  expect_equal(aas$toUnmodifiedString(),"DFPIANGER")
  aas = AASequence$fromStringPermissive("DFPIANGER", T)
  expect_equal(aas$toUnmodifiedString(),"DFPIANGER")

  seq = AASequence$fromString("PEPTIDESEKUEM(Oxidation)CER")
  expect_equal(seq$toString(),"PEPTIDESEKUEM(Oxidation)CER")
  expect_equal(seq$toUnmodifiedString(),"PEPTIDESEKUEMCER")
  expect_equal(seq$toBracketString(),"PEPTIDESEKUEM[147]CER")
  expect_equal(seq$toBracketString(TRUE),"PEPTIDESEKUEM[147]CER")

  expect_true(seq$toBracketString(FALSE) == "PEPTIDESEKUEM[147.03540001709996]CER" || seq$toBracketString(FALSE) == "PEPTIDESEKUEM[147.035400017100017]CER")
  expect_true(seq$toBracketString(FALSE) == "PEPTIDESEKUEM[147.03540001709996]CER" || seq$toBracketString(F) == "PEPTIDESEKUEM[147.035400017100017]CER")

  # Change here: bool
  expect_equal(seq$toUniModString(),"PEPTIDESEKUEM(UniMod:35)CER")
  expect_true(seq$isModified())
  expect_false(seq$hasCTerminalModification())
  expect_false(seq$hasNTerminalModification())
  expect_false(seq$empty())

  # has selenocysteine
  expect_true(!is.null(seq$getResidue(1)))
  expect_equal(seq$size(),16)
  expect_equal(seq$getFormula(Residue$ResidueType()$Full,0), EmpiricalFormula$new("C75H122N20O32S2Se1"))

  expect_lt((seq$getMonoWeight(Residue$ResidueType()$Full, 0) - 1952.7200317517998),1e-5)

})


# Test testElement --------------------------------------------------------

test_that("",{
  ins = Element$new()
  ins$setAtomicNumber(6)
  ins$getAtomicNumber()
  ins$setAverageWeight(12.011)
  ins$getAverageWeight()
  ins$setMonoWeight(12)
  ins$getMonoWeight()

  iso = IsotopeDistribution$new()
  ins$setIsotopeDistribution(iso)
  ins$getIsotopeDistribution()
  ins$setName("Carbon")
  ins$getName()
  ins$setSymbol("C")
  ins$getSymbol()

  e = Element$new()
  e$setSymbol("blah")
  e$setSymbol(c("blah"))
  oms_string <- String$new("blu")
  e$setSymbol(oms_string)
  expect_true(!is.null(oms_string))
  expect_equal(oms_string$toString(),"blu")

  evil = "blü"
  e$setSymbol(evil)
  expect_equal(e$getSymbol(),"blü")
  String$new(e$getSymbol()) == String$new("blü")
  expect_equal(String$new(e$getSymbol())$toString(),"blü")

})


# Test Residue ------------------------------------------------------------

test_that("test Residue",{
  ins <- Residue$new()
  # Change here: Enum
  rt <- Residue$ResidueType()
  rt$Full
  rt$Internal
  rt$NTerminal
  rt$CTerminal
  rt$AIon
  rt$BIon
  rt$CIon
  rt$XIon
  rt$YIon
  rt$ZIon
  rt$SizeOfResidueType
  expect_true(T)
})


# Test IsotopeDistribution ------------------------------------------------

test_that("test IsotopeDistribution",{

  # N.B. iteration is not supported currently.
  ins <- IsotopeDistribution$new()
  ins$getMax()
  ins$getMin()
  ins$size()
  ins$clear()
  ins$renormalize()
  ins$trimLeft(6.0)
  ins$trimRight(8.0)

  ins$clear()
  ins$insert(1, 2)
  ins$insert(6, 5)

  expect_equal(ins$size(), 2)

  for(i in ins$getContainer()){
    print(paste("MZ:",i$getMZ(),"Intensity:",i$getIntensity()))
  }

})


# Test FineIsotopePatternGenerator -----------------------------------------
test_that("test FineIsotopePatternGenerator",{

  iso = FineIsotopePatternGenerator$new()
  iso$setThreshold(1e-5)
  iso$setAbsolute(TRUE)
  iso$getAbsolute()

  methanol = EmpiricalFormula$new("CH3OH")
  water = EmpiricalFormula$new("H2O")
  mw = methanol + water
  iso_dist = mw$getIsotopeDistribution(FineIsotopePatternGenerator$new(1e-20, F, F))
  expect_equal(length(iso_dist$getContainer()), 56)
  iso_dist = mw$getIsotopeDistribution(FineIsotopePatternGenerator$new(1e-200, F, F))
  expect_equal(length(iso_dist$getContainer()), 84)

  c100 = EmpiricalFormula$new("C100")
  iso_dist = c100$getIsotopeDistribution(FineIsotopePatternGenerator$new(1e-200, F, F))
  expect_equal(length(iso_dist$getContainer()), 101)

  expect_equal(c100$getIsotopeDistribution(FineIsotopePatternGenerator$new(1e-2, F, F))$size(), 6)
  expect_equal(c100$getIsotopeDistribution(FineIsotopePatternGenerator$new(1e-2, F, T))$size(), 5)
  expect_equal(c100$getIsotopeDistribution(FineIsotopePatternGenerator$new(1e-2, T, F))$size(), 5)
  expect_equal(c100$getIsotopeDistribution(FineIsotopePatternGenerator$new(1e-2, T, T))$size(), 5)

  expect_equal(c100$getIsotopeDistribution(FineIsotopePatternGenerator$new(1e-10, F, F))$size(), 14)
  expect_equal(c100$getIsotopeDistribution(FineIsotopePatternGenerator$new(1e-10, F, T))$size(), 13)
  expect_equal(c100$getIsotopeDistribution(FineIsotopePatternGenerator$new(1e-10, T, F))$size(), 10)
  expect_equal(c100$getIsotopeDistribution(FineIsotopePatternGenerator$new(1e-10, T, T))$size(), 10)

  iso <- FineIsotopePatternGenerator$new(1e-5, F, F)
  isod <- iso$run(methanol)
  expect_equal(length(isod$getContainer()), 6)
  expect_lt( abs(isod$getContainer()[[1]]$getMZ() - 32.0262151276), 1e-5)
  expect_lt( abs(isod$getContainer()[[1]]$getIntensity() - 0.986442089081), 1e-5)
})


#  Test CoarseIsotopePatternGenerator -------------------------------------

test_that("test CoarseIsotopePatternGenerator",{

  iso <- CoarseIsotopePatternGenerator$new()
  iso$setMaxIsotope(5)
  expect_equivalent(iso$getMaxIsotope(),5L)
  res <- iso$estimateFromPeptideWeight(500)

  methanol <- EmpiricalFormula$new("CH3OH")
  water <- EmpiricalFormula$new("H2O")
  mw <- methanol + water
  iso_dist <- mw$getIsotopeDistribution(CoarseIsotopePatternGenerator$new(3))
  expect_length(iso_dist$getContainer(),3)
  iso_dist <- mw$getIsotopeDistribution(CoarseIsotopePatternGenerator$new(0))
  expect_length(iso_dist$getContainer(),18)

  iso <- CoarseIsotopePatternGenerator$new(10)
  isod <- iso$run(methanol)
  expect_length(isod$getContainer(),10)
  expect_true(abs(isod$getContainer()[[1]]$getMZ() - 32.0262151276) < 1e-5)
  expect_true(abs(isod$getContainer()[[1]]$getIntensity() - 0.986442089081) < 1e-5)
})


# test EmpiricalFormula --------------------------------------------------------

test_that("tests EmpiricalFormula",{

  ins <- EmpiricalFormula$new()
  ins$getMonoWeight()
  ins$getAverageWeight()
  ins$getIsotopeDistribution(CoarseIsotopePatternGenerator$new(0))
  ins$getNumberOfAtoms()
  ins$setCharge(2)
  ins$getCharge()
  ins$toString()
  ins$isEmpty()
  ins$isCharged()
  ins$hasElement(Element$new())

  ef <- EmpiricalFormula$new("C2H5")
  s <- ef$toString()
  expect_equal(s,"C2H5")
  m <- ef$getElementalComposition()
  # dict is returned as collections::dict()
  expect_equal(m$get("C"),2)
  expect_equal(m$get("H"),5)
  expect_equal(ef$getNumberOfAtoms(),7)
})


# Test testIdentificationHit ----------------------------------------------

test_that("test testIdentificationHit",{
  f = IdentificationHit$new()
  check_MetaInfoInterface(f)

  expect_true(!is.null(IdentificationHit$new()$setId))
  expect_true(!is.null(IdentificationHit$new()$getId))
  expect_true(!is.null(IdentificationHit$new()$setCharge))
  expect_true(!is.null(IdentificationHit$new()$getCharge))
  expect_true(!is.null(IdentificationHit$new()$setCalculatedMassToCharge))
  expect_true(!is.null(IdentificationHit$new()$getCalculatedMassToCharge))
  expect_true(!is.null(IdentificationHit$new()$setExperimentalMassToCharge))
  expect_true(!is.null(IdentificationHit$new()$getExperimentalMassToCharge))
  expect_true(!is.null(IdentificationHit$new()$setName))
  expect_true(!is.null(IdentificationHit$new()$getName))
  expect_true(!is.null(IdentificationHit$new()$setPassThreshold))
  expect_true(!is.null(IdentificationHit$new()$getPassThreshold))
  expect_true(!is.null(IdentificationHit$new()$setRank))
  expect_true(!is.null(IdentificationHit$new()$getRank))


  f$setId("test_id")
  expect_equal(f$getId(),"test_id")

  f$setId("test_id")
  expect_equal(f$getId(),"test_id")

  f$setCharge(5)
  expect_equal(f$getCharge(),5)

  f$setCalculatedMassToCharge(5.0)
  expect_equal(f$getCalculatedMassToCharge(),5.0)

  f$setExperimentalMassToCharge(5.0)
  expect_equal(f$getExperimentalMassToCharge(),5.0)

  f$setName("test")
  expect_equal(f$getName(),"test")

  f$setPassThreshold(T)
  expect_true(f$getPassThreshold())

  f$setRank(42)
  expect_equal(f$getRank(),42)

})


# Test SpectrumIdentification ---------------------------------------------

test_that("test SpectrumIdentification",{
  f = SpectrumIdentification$new()
  check_MetaInfoInterface(f)

  expect_true(!is.null(SpectrumIdentification$new()$setHits))
  expect_true(!is.null(SpectrumIdentification$new()$addHit))
  expect_true(!is.null(SpectrumIdentification$new()$getHits))

  hit <- IdentificationHit$new()
  hit$setName("test1")
  f$addHit(hit)
  hit <- IdentificationHit$new()
  hit$setName("test2")
  f$addHit(hit)
  all_hits <- f$getHits()
  expect_length(all_hits,2)
  expect_true(any(sapply(all_hits,function(h) h$getName()=="test1")))
  expect_true(any(sapply(all_hits,function(h) h$getName()=="test2")))
})


# Test Identification -----------------------------------------------------

test_that("test Identification",{
  f <- Identification$new()
  check_MetaInfoInterface(f)

  expect_true(!is.null(Identification$new()$setCreationDate))
  expect_true(!is.null(Identification$new()$getCreationDate))
  expect_true(!is.null(Identification$new()$setSpectrumIdentifications))
  expect_true(!is.null(Identification$new()$addSpectrumIdentification))
  expect_true(!is.null(Identification$new()$getSpectrumIdentifications))

  id1 <- SpectrumIdentification$new()
  f$addSpectrumIdentification(id1)
  expect_length(f$getSpectrumIdentifications(),1)
  id2 = SpectrumIdentification$new()
  f$addSpectrumIdentification(id2)
  expect_length(f$getSpectrumIdentifications(),2)
})


# Test ModificationDefinitionsSet -----------------------------------------

test_that("test ModificationDefinitionsSet",{
  empty <- ModificationDefinitionsSet$new()
  fixed <- list("Carbamidomethyl")
  variable <- list("Oxidation")
  full <- ModificationDefinitionsSet$new(fixed,variable)
  expect_true(T)
})



# Test AcquisitionInfo ----------------------------------------------------

test_that("test AcquisitionInfo",{
  ai = AcquisitionInfo$new()
  ai == ai
  !(ai != ai)
  ai$setMethodOfCombination("ABC")
  expect_equal(ai$getMethodOfCombination(),"ABC")
})



# Test BaseFeature --------------------------------------------------------

test_that("test BaseFeature",{
  bf <- BaseFeature$new()
  check_MetaInfoInterface(bf)
  check_UniqueIdInterface(bf)
  bf$clearUniqueId()
  expect_true(as.logical(bf$ensureUniqueId()))
  expect_equal(bf$getCharge(),0)
  expect_true(is_scalar_double(bf$getQuality()))
  expect_true(is_scalar_integer(bf$getUniqueId()))
  expect_true(is_scalar_double(bf$getWidth()))

  expect_true(!(bf$hasInvalidUniqueId()))
  expect_true(as.logical(bf$hasValidUniqueId()))

  check_MetaInfoInterface(bf)
  bf$setCharge(1)
  bf$setQuality(0.0)
  bf$setWidth(1.0)
})


# Test AnnotationState ----------------------------------------------------

test_that("test AnnotationState",{
  state <- AnnotationState()
  expect_true(!is.null(state$FEATURE_ID_NONE))
  expect_true(!is.null(state$FEATURE_ID_SINGLE))
  expect_true(!is.null(state$FEATURE_ID_MULTIPLE_SAME))
  expect_true(!is.null(state$FEATURE_ID_MULTIPLE_DIVERGENT))
  expect_true(!is.null(state$SIZE_OF_ANNOTATIONSTATE))
})



# Test ChecksumType -------------------------------------------------------

test_that("test ChecksumType",{
  expect_true(is_scalar_integer(ChecksumType()$MD5))
  expect_true(is_scalar_integer(ChecksumType()$SHA1))
  expect_true(is_scalar_integer(ChecksumType()$SIZE_OF_CHECKSUMTYPE))
  expect_true(is_scalar_integer(ChecksumType()$UNKNOWN_CHECKSUM))
})


# Test ChromatogramPeak ---------------------------------------------------

test_that("test ChromatogramPeak",{
  p <- ChromatogramPeak$new()
  expect_true(p == p)
  expect_true(!(p!=p))
  p$setIntensity(12.0)
  p$setRT(34.0)
  expect_equal(p$getIntensity(), 12.0)
  expect_equal(p$getRT(),34.0)
})


# Test testChromatogramToosl ----------------------------------------------

test_that("test testChromatogramToosl",{
  ct <- ChromatogramTools$new()
  expect_visible(ct$convertChromatogramsToSpectra)
  expect_visible(ct$convertSpectraToChromatograms)
})


# Test ConsensusFeature ---------------------------------------------------

test_that("test ConsensusFeature",{
  f <- ConsensusFeature$new()
  check_UniqueIdInterface(f)
  check_MetaInfoInterface(f)

  f$setCharge(1)
  f$setQuality(2.0)
  f$setWidth(4.0)
  expect_equal(f$getCharge(),1)
  expect_equal(f$getQuality(),2.0)
  expect_equal(f$getWidth(),4.0)

  f$insert(0, Peak2D$new(), 1)
  f$insert(1, BaseFeature$new())
  f$insert(2, ConsensusFeature$new())

  f$computeConsensus()
  f$computeDechargeConsensus
  f$computeMonoisotopicConsensus()

  expect_gt(f$size(),0)
  p = f$getPeptideIdentifications()
  f$setPeptideIdentifications(p)

})


# Test ConsensusMap -------------------------------------------------------

test_that("test ConsensusMap",{
  m = ConsensusMap$new()
  m$clear()
  m$clearUniqueId()
  m$ensureUniqueId()
  m$getDataProcessing()
  m$getColumnHeaders()$as_list()
  m$getProteinIdentifications()
  m$getUnassignedPeptideIdentifications()
  # Returns -1 as Unique Id is a 64 bit signed integer.
  m$getUniqueId()
  m$hasInvalidUniqueId()
  m$hasValidUniqueId()
  m$setDataProcessing
  m$setColumnHeaders
  m$setProteinIdentifications
  m$setUnassignedPeptideIdentifications
  m$setUniqueId
  m$setUniqueIds
  m$size()
  m$sortByIntensity()
  m$sortByMZ()
  m$sortByMaps()
  m$sortByPosition()
  m$sortByQuality()
  m$sortByRT()
  m$sortBySize()
  m$updateRanges()

  expect_true(is_scalar_double(m$getMin()[[1]]))
  expect_true(is_scalar_double(m$getMin()[[1]]))
  expect_true(is_scalar_double(m$getMax()[[2]]))
  expect_true(is_scalar_double(m$getMax()[[2]]))
  expect_true(is_scalar_double(m$getMinInt()))
  expect_true(is_scalar_double(m$getMaxInt()))

  m$getIdentifier()
  m$getLoadedFileType()
  m$getLoadedFilePath()
  expect_true(m == m)
  expect_true(!(m != m))
})


# Test ConsensusXMLFile ---------------------------------------------------

test_that("test ConsensusXMLFile",{
  f = ConsensusXMLFile$new()
  f$getOptions()
  expect_true(!is.null(f$load))
  expect_true(!is.null(f$store))
})


# Test XTandemXMLFile -----------------------------------------------------

test_that("test XTandemXMLFile",{
  f = XTandemXMLFile$new()
  expect_true(!is.null(f$load))
})


# Test XTandemInfile ------------------------------------------------------

test_that("test XTandemInfile",{
  f <- XTandemInfile$new()
  expect_true(!is.null(f$setFragmentMassTolerance))
  expect_true(!is.null(f$getFragmentMassTolerance))
  expect_true(!is.null(f$setPrecursorMassTolerancePlus))
  expect_true(!is.null(f$getPrecursorMassTolerancePlus))
  expect_true(!is.null(f$setPrecursorMassToleranceMinus))
  expect_true(!is.null(f$getPrecursorMassToleranceMinus))
  expect_true(!is.null(f$setPrecursorErrorType))
  expect_true(!is.null(f$getPrecursorErrorType))
  expect_true(!is.null(f$setFragmentMassErrorUnit))
  expect_true(!is.null(f$getFragmentMassErrorUnit))
  expect_true(!is.null(f$setPrecursorMassErrorUnit))
  expect_true(!is.null(f$getPrecursorMassErrorUnit))
  expect_true(!is.null(f$setNumberOfThreads))
  expect_true(!is.null(f$getNumberOfThreads))
  expect_true(!is.null(f$setModifications))
  expect_true(!is.null(f$getModifications))
  expect_true(!is.null(f$setOutputFilename))
  expect_true(!is.null(f$getOutputFilename))
  expect_true(!is.null(f$setInputFilename))
  expect_true(!is.null(f$getInputFilename))
  expect_true(!is.null(f$setTaxonomyFilename))
  expect_true(!is.null(f$getTaxonomyFilename))
  expect_true(!is.null(f$setDefaultParametersFilename))
  expect_true(!is.null(f$getDefaultParametersFilename))

  f$setTaxon("testTaxon")
  expect_equal(f$getTaxon(),"testTaxon")

  expect_true(!is.null(f$setMaxPrecursorCharge))
  expect_true(!is.null(f$getMaxPrecursorCharge))
  expect_true(!is.null(f$setNumberOfMissedCleavages))
  expect_true(!is.null(f$getNumberOfMissedCleavages))
  expect_true(!is.null(f$setMaxValidEValue))
  expect_true(!is.null(f$getMaxValidEValue))
  expect_true(!is.null(f$setSemiCleavage))
  expect_true(!is.null(f$setAllowIsotopeError))
  expect_true(!is.null(f$write))
  expect_true(!is.null(f$setCleavageSite))
  expect_true(!is.null(f$getCleavageSite))
})


# Test SignalToNoiseEstimatorMedian ---------------------------------------

test_that("test SignalToNoidseEstimatorMedian",{
  f <- SignalToNoiseEstimatorMedian$new()
  expect_true(!is.null(f$init))
  expect_true(!is.null(f$getSignalToNoise))
})


# Test SignalToNoiseEstimatorMedianChrom ---------------------------------------

test_that("test SignalToNoiseEstimatorMedianChrom",{
  f <- SignalToNoiseEstimatorMedianChrom$new()
  expect_true(!is.null(f$init))
  expect_true(!is.null(f$getSignalToNoise))
})


# Test ConvexHull2D -------------------------------------------------------

test_that("test ConvexHull2D",{
  ch <- ConvexHull2D$new()
  ch$clear()
  expect_true(ch==ch)
})


# Test DataProcessing -----------------------------------------------------
test_that("test DataProcessing",{

    dp = DataProcessing$new()
    check_MetaInfoInterface(dp)

    expect_true(dp==dp)
    expect_false(dp!=dp)

    dp$clearMetaInfo()
    k <- list()
    dp$getKeys(k)
    expect_equal(k,list())

    dp$getMetaValue
    ac = dp$getProcessingActions()
    expect_false(T %in% duplicated(ac))

    ac = c(ProcessingAction()$PEAK_PICKING, ProcessingAction()$BASELINE_REDUCTION)
    dp$setProcessingActions(ac)

    expect_length(dp$getProcessingActions(), 2)
    expect_true(is_scalar_character(dp$getSoftware()$getName()))

    expect_true(is_scalar_character(dp$getSoftware()$getVersion()))

    dp$isMetaEmpty()
    dp$metaValueExists
    dp$removeMetaValue
    s = dp$getSoftware()
    s$setName("pyopenms")
    dp$setSoftware(s)

    expect_equal(dp$getSoftware()$getName(), "pyopenms")
})


# Test DataType -----------------------------------------------------------

test_that("test DataType",{
  expect_true(is_scalar_integer(DataType()$DOUBLE_LIST))
  expect_true(is_scalar_integer(DataType()$DOUBLE_VALUE))
  expect_true(is_scalar_integer(DataType()$EMPTY_VALUE))
  expect_true(is_scalar_integer(DataType()$INT_LIST))
  expect_true(is_scalar_integer(DataType()$INT_VALUE))
  expect_true(is_scalar_integer(DataType()$STRING_LIST))
  expect_true(is_scalar_integer(DataType()$STRING_VALUE))
})


# Test DataValue ----------------------------------------------------------

test_that("test DataValue",{
  a = DataValue$new()
  expect_false(!a$isEmpty())

  a = DataValue$new(1)
  expect_true(!a$isEmpty())
  expect_equivalent(a$toInt(), 1L)
  expect_equal(a$valueType(), DataType()$INT_VALUE)

  # Note if you provide 1.0 while creating the object,
  # the constructor for Integer is called.
  a = DataValue$new(1.01)
  expect_true(!a$isEmpty())
  expect_equivalent(a$toDouble(), 1.01)
  expect_equal(a$valueType(), DataType()$DOUBLE_VALUE)

  a = DataValue$new("1")
  expect_true(!a$isEmpty())
  expect_equal(a$toString(), "1")
  expect_equal(a$valueType(), DataType()$STRING_VALUE)

  a = DataValue$new(list(1L))
  expect_true(!a$isEmpty())
  expect_equivalent(a$toIntList(), c(1L))
  expect_equal(a$valueType(), DataType()$INT_LIST)

  # 1L is integer while 1 is double in R. It is important to add L to specify integer,
  # because reticulate converts 1L -> 1 (int) while 1 -> 1.0 (float) in python.
  a = DataValue$new(list(1.1))
  expect_true(!a$isEmpty())
  expect_equivalent(a$toDoubleList(), c(1.1))
  expect_equal(a$valueType(), DataType()$DOUBLE_LIST)

  a = DataValue$new(list("1.0"))
  expect_true(!a$isEmpty())
  expect_equivalent(a$toStringList(), list("1.0"))
  expect_equal(a$valueType(), DataType()$STRING_LIST)

  expect_true(is.null(MSSpectrum$new()$getMetaValue("nonexisingkey")))
})


# Test Adduct -------------------------------------------------------------

test_that("test Adduct",{
  a = Adduct$new()
  expect_true(T)
})

# Test GaussFitter -------------------------------------------------------------

test_that("test GaussFitter",{
  a = GaussFitter$new()
  expect_true(T)
})

# Test GaussFitResult -------------------------------------------------------------

test_that("test GaussFitResult",{
  ins = GaussFitResult$new(0,0,0)
  ins$A <- 5
  ins$x0 <- 5
  ins$sigma <- 5
  expect_true(T)
})

# Test ChargePair -------------------------------------------------------------

test_that("test ChargePair",{
  a = ChargePair$new()
  expect_true(T)
})


# Test Compomer -----------------------------------------------------------
test_that("test Compomer",{
  a = Compomer$new()
  expect_true(T)
})

# Test CVMappings -----------------------------------------------------------
test_that("test CVMappings",{
  a = CVMappings$new()
  expect_true(T)
})

# Test CVMappingFile -----------------------------------------------------------
test_that("test CVMappingFile",{
  a = CVMappingFile$new()
  expect_true(!is.null(CVMappingFile$new()$load))
})

# test ControlledVocabulary -----------------------------------------------
test_that("test ControlledVocabulary",{
  val = ControlledVocabulary$new()
  expect_true(!is.null(ControlledVocabulary$new()$loadFromOBO))
})


# test SemanticValidator --------------------------------------------------
test_that("test SemanticValidator",{
  m = CVMappings$new()
  cv = ControlledVocabulary$new()

  val = SemanticValidator$new(m,cv)
  expect_true(!is.null(val$validate))
  expect_true(!is.null(val$setCheckTermValueTypes))
  expect_true(!is.null(val$setCheckUnits))
})


# test Feature ------------------------------------------------------------

test_that("Feature",{
  f = Feature$new()
  check_MetaInfoInterface(f)
  check_UniqueIdInterface(f)

  f$setConvexHulls(f$getConvexHulls())
  f$setSubordinates(f$getSubordinates())
  f$setUniqueId(12345)

  expect_true(f == f)
  expect_false(f != f)

  f$setCharge(-1)
  expect_equal(f$getCharge(), -1)
  f$setIntensity(10.0)
  expect_equivalent(f$getIntensity(), 10)
  f$setQuality(0, 20.0)
  expect_equal(f$getQuality(0), 20)
  f$setRT(30.0)
  expect_equal(f$getRT(), 30)
  f$setWidth(40.0)
  expect_equal(f$getWidth(), 40)

  p = f$getPeptideIdentifications()
  f$setPeptideIdentifications(p)
})

# Test Param --------------------------------------------------------------
check_Param <- function(p){

  expect_true(p == p)

  dd = p$asDict()
  expect_equal(dd$size(), p$size())
  expect_true(is.environment(dd) && identical(parent.env(dd), asNamespace("collections")))


  for(k in p$keys()){
    value <- p[k]
    p[k] <- value
    p$update(p)
    p$update(p$asDict())
    expect_equal(p[k], value)
    desc = p$getDescription(k)
    tags = p$getTags(k)
    # print("********************")
    # print(value)
    # print(desc)
    # print(tags)
    # print("********************")
    p$setValue(k, value, desc, tags)
    p$setValue(k, value, desc)
    expect_true(p$exists(k))

    if (length(strsplit(k,":")[[1]]) < 2) {next}
    f = strsplit(k,":")[[1]][1]
    p$setSectionDescription(f, k)
    expect_equal(p$getSectionDescription(f), k)
    expect_true(!is.null(p$get(k)))
  }

  expect_equal(length(p$values()), length(lapply(p$keys(),function(k) p[k])))
  expect_equal(p$items(), lapply(p$keys(), function(k) list(k, p[k])))

  expect_true(!p$exists("asdflkj01231321321v"))
  p$addTag(k, "a")
  p$addTags(k, c("","c"))
  expect_equal(p$getTags(k), c("","a","c"))
  p$clearTags(k)
  expect_equal(p$getTags(k), list())

  pn = Param$new()
  pn$insert("master:", p)
  expect_true(pn$exists(paste0("master:",k)))

  p1 = pn$copy("master:", T)
  expect_equal(p1, p)
  p1$update(p)
  p1$update(p,0)
  p1$update(p,1)
  p1$update(dd)

  p$setValidStrings
  p$setMinFloat
  p$setMaxFloat
  p$setMinInt
  p$setMaxInt

  ph = ParamXMLFile$new()
  ph$store("test.ini", p)
  p1 = Param$new()
  ph$load("test.ini", p1)
  expect_equal(p, p1)

  e1 = p1$getEntry(k)

  for (f in c("name","description","value","tags","valid_strings","min_float","max_float","min_int","max_int")) {
    expect_true(!is.null(reticulate::py_get_attr(e1$.__enclos_env__$private$py_obj,f)))
  }

  expect_true(e1 == e1)
  expect_equivalent(p1$get("abcde", 7), 7L)
}

# on_failure(check_Param) <- function(call,env){
#   paste0(deparse(class(call$x)[1]),"does not satisfy Param tests")
# }

# Test FeatureFinder ------------------------------------------------------
test_that("test FeatureFinder",{
  ff = FeatureFinder$new()
  name = FeatureFinderAlgorithmPicked$getProductName()
  ff$run(name, MSExperiment$new(), FeatureMap$new(), Param$new(), FeatureMap$new())

  check_ProgressLogger(ff)

  p = ff$getParameters(name)
  check_Param(p)
})


# Test FeatureFileOptions -------------------------------------------------

test_that("test FeatureFileOptions",{
  fo <- FeatureFileOptions$new()
  fo$getLoadConvexHull()
  fo$getLoadSubordinates()
  fo$getSizeOnly()
  expect_true(!is.null(fo$setLoadConvexHull))
  expect_true(!is.null(fo$setLoadSubordinates))
  expect_true(!is.null(fo$setMetadataOnly))
  expect_true(!is.null(fo$setSizeOnly))
})



# Test FeatureFinderAlgorithmPicked ---------------------------------------

test_that("test FeatureFinderAlgorithmPicked",{
  ff = FeatureFinderAlgorithmPicked$new()
  p = ff$getDefaults()
  check_Param(p)

  check_Param(ff$getParameters())
  expect_equal(ff$getName(), "FeatureFinderAlgorithm")
  expect_equal(FeatureFinderAlgorithmPicked$getProductName(), "centroided")
  ff$setName("test")
  expect_equal(ff$getName(), "test")
})


# test FeatureFinderAlgorithmSH -------------------------------------------

test_that("test FeatureFinderAlgorithmSH",{

  ff = FeatureFinderAlgorithmSH$new()
  p = ff$getDefaults()
  check_Param(p)

  expect_equal(ff$getName(),"FeatureFinderAlgorithm")
  expect_equal(FeatureFinderAlgorithmSH$getProductName(),"superhirn")

  ff$setParameters(Param$new())

  ff$setName("test")
  expect_equal(ff$getName(),"test")
})

# test FeatureFinderAlgorithmIsotopeWavelet -------------------------------------------

test_that("test FeatureFinderAlgorithmIsotopeWavelet",{

  ff = FeatureFinderAlgorithmIsotopeWavelet$new()
  p = ff$getDefaults()
  check_Param(p)

  expect_equal(ff$getName(),"FeatureFinderAlgorithm")
  expect_equal(FeatureFinderAlgorithmIsotopeWavelet$getProductName(),"isotope_wavelet")

  ff$setParameters(Param$new())

  ff$setName("test")
  expect_equal(ff$getName(),"test")
})


# test CompNovoIdentification ---------------------------------------------

test_that("test CompNovoIdentification",{
  ff = CompNovoIdentification$new()
  p = ff$getDefaults()
  check_Param(p)

  expect_true(!is.null(CompNovoIdentification$new()$getIdentification))
  expect_true(!is.null(CompNovoIdentification$new()$getIdentifications))
})

# test CompNovoIdentificationCID ---------------------------------------------

test_that("test CompNovoIdentificationCID",{
  ff = CompNovoIdentificationCID$new()
  p = ff$getDefaults()
  check_Param(p)

  expect_true(!is.null(CompNovoIdentificationCID$new()$getIdentification))
  expect_true(!is.null(CompNovoIdentificationCID$new()$getIdentifications))
})


# test ExperimentalSettings -----------------------------------------------

test_that("test ExperimentalSettings",{
  ff = ExperimentalSettings$new()
  expect_true(T)
})


# Test FeatureDeconvolution -----------------------------------------------

test_that("test FeatureDeconvolution",{
  ff = FeatureDeconvolution$new()
  p = ff$getDefaults()
  check_Param(p)

  expect_true(!is.null(FeatureDeconvolution$new()$compute))
})

# Test InternalCalibration -----------------------------------------------

test_that("test InternalCalibration",{
  ff = InternalCalibration$new()
  # getDefaults() also not present in pyopenms.
  # p = ff$getDefaults()
  # check_Param(p)

  # expect_true(!is.null(InternalCalibration$compute))
  expect_true(T)
})


# Test ItraqConstants -----------------------------------------------------
test_that("test ItraqConstants",{
  constants = ItraqConstants$new()

  expect_true(!is.null(ITRAQ_TYPES()$FOURPLEX))
  expect_true(!is.null(ITRAQ_TYPES()$EIGHTPLEX))
  expect_true(!is.null(ITRAQ_TYPES()$TMT_SIXPLEX))

  expect_true(!is.null(constants$getIsotopeMatrixAsStringList))
  expect_true(!is.null(constants$updateIsotopeMatrixFromStringList))
  expect_true(!is.null(constants$translateIsotopeMatrix))

})


# test LinearResampler ----------------------------------------------------

test_that("test LinearResampler",{
  ff = LinearResampler$new()
  p = ff$getDefaults()
  check_Param(p)

  expect_true(!is.null(LinearResampler$new()$raster))
  expect_true(!is.null(LinearResampler$new()$rasterExperiment))
})


# test PeptideAndProteinQuant ---------------------------------------------

test_that("test PeptideAndProteinQuant",{
  ff = PeptideAndProteinQuant$new()
  p = ff$getDefaults()
  check_Param(p)

  expect_true(!is.null(PeptideAndProteinQuant$new()$quantifyPeptides))
  expect_true(!is.null(PeptideAndProteinQuant$new()$quantifyProteins))

})

# test SeedListGenerator ---------------------------------------------

test_that("test SeedListGenerator",{
  ff = SeedListGenerator$new()
  # p = ff$getDefaults()
  # check_Param(p)
  expect_true(T)
})


# test TOFCalibration -----------------------------------------------------

test_that("test TOFCalibration",{
  ff = TOFCalibration$new()
  p = ff$getDefaults()

  expect_true(!is.null(TOFCalibration$new()$calibrate))
  expect_true(!is.null(TOFCalibration$new()$pickAndCalibrate))
})

# test FalseDiscoveryRate -----------------------------------------------------

test_that("test FalseDiscoveryRate",{
  ff = FalseDiscoveryRate$new()
  p = ff$getDefaults()
  check_Param(p)

  expect_true(!is.null(FalseDiscoveryRate$new()$apply))
})

# test IDFilter -----------------------------------------------------

test_that("test IDFilter",{
  ff = IDFilter$new()
  expect_true(T)
})

# test ProteinResolver -----------------------------------------------------

test_that("test ProteinResolver",{
  ff = ProteinResolver$new()

  expect_true(!is.null(ProteinResolver$new()$resolveConsensus))
  expect_true(!is.null(ProteinResolver$new()$resolveID))
  expect_true(!is.null(ProteinResolver$new()$setProteinData))
  expect_true(!is.null(ProteinResolver$new()$getResults))
})

# test SvmTheoreticalSpectrumGeneratorTrainer -----------------------------------------------------

test_that("test SvmTheoreticalSpectrumGeneratorTrainer",{
  ff = SvmTheoreticalSpectrumGeneratorTrainer$new()

  expect_true(!is.null(SvmTheoreticalSpectrumGeneratorTrainer$new()$trainModel))
  expect_true(!is.null(SvmTheoreticalSpectrumGeneratorTrainer$new()$normalizeIntensity))
})


# test PosteriorErrorProbabilityModel -------------------------------------

test_that("test PosteriorErrorProbabilityModel",{

  model = PosteriorErrorProbabilityModel$new()
  p = model$getDefaults()
  check_Param(p)

  expect_true(!is.null(PosteriorErrorProbabilityModel$new()$fit))
  expect_true(!is.null(PosteriorErrorProbabilityModel$new()$computeProbability))

  scores = as.double(0:9)
  model$fit(scores, "none")
  model$fit(scores, scores, "none")

  model$fillLogDensities(scores, scores, scores)

  expect_true(!is.null(model$computeLogLikelihood))
  expect_true(!is.null(model$pos_neg_mean_weighted_posteriors))

  Gfr = model$getCorrectlyAssignedFitResult()
  Gfr = model$getIncorrectlyAssignedFitResult()
  model$getNegativePrior()
  model$computeProbability(5.0)

  # model.InitPlots

  target = as.double(0:9)
  model$getGumbelGnuplotFormula(Gfr)
  model$getGaussGnuplotFormula(Gfr)
  model$getBothGnuplotFormula(Gfr, Gfr)
  model$plotTargetDecoyEstimation(target, target)
  model$getSmallestScore()

})


# Test SeedListGenerator --------------------------------------------------

test_that("test SeedListGenerator",{
  ff = SeedListGenerator$new()
  expect_true(T)
})


# test ConsensusMapNormalizerAlgorithmMedian ------------------------------

test_that("test ConsensusMapNormalizerAlgorithmMedian",{
  ff = ConsensusMapNormalizerAlgorithmMedian$new()
  expect_true(!is.null(ConsensusMapNormalizerAlgorithmMedian$new()$normalizeMaps))
})


# Test ConsensusMapNormalizerAlgorithmQuantile ----------------------------

test_that("test ConsensusMapNormalizerAlgorithmQuantile",{
  ff = ConsensusMapNormalizerAlgorithmQuantile$new()
  expect_true(!is.null(ConsensusMapNormalizerAlgorithmQuantile$new()$normalizeMaps))
})

# Test ConsensusMapNormalizerAlgorithmThreshold ----------------------------

test_that("test ConsensusMapNormalizerAlgorithmThreshold",{
  ff = ConsensusMapNormalizerAlgorithmThreshold$new()
  expect_true(!is.null(ConsensusMapNormalizerAlgorithmThreshold$new()$normalizeMaps))
})


# test FeatureFinderAlgorithmSH -------------------------------------------

test_that("test FeatureFinderAlgorithmSH",{
  ff = FeatureFinderAlgorithmSH$new()
  expect_true(!is.null(FeatureFinderAlgorithmSH$new()$setData))
  expect_true(!is.null(FeatureFinderAlgorithmSH$new()$run))
})


# test FeatureFinderAlgorithmIsotopeWavelet -------------------------------

test_that("test FeatureFinderAlgorithmIsotopeWavelet",{
  ff = FeatureFinderAlgorithmIsotopeWavelet$new()

  expect_true(!is.null(FeatureFinderAlgorithmIsotopeWavelet$new()$setData))
  expect_true(!is.null(FeatureFinderAlgorithmIsotopeWavelet$new()$run))
})


# test AScore -------------------------------------------------------------

test_that("test AScore",{
  ff = AScore$new()
  hit = PeptideHit$new()
  spectrum = MSSpectrum$new()
  ff$compute(hit, spectrum)
  expect_true(T)
})


# test IDRipper -----------------------------------------------------------

test_that("test IDRipper",{
  ff = IDRipper$new()
  expect_true(!is.null(IDRipper$new()$rip))
})


# test FASTAFile ----------------------------------------------------------

test_that("test FASTAFile",{
  ff = FASTAFile$new()
  expect_true(!is.null(FASTAFile$new()$load))
  expect_true(!is.null(FASTAFile$new()$store))
})


# test FASTAEntry ---------------------------------------------------------

test_that("test FASTAEntry",{
  ff = FASTAEntry$new()
  expect_true(T)
})


# test InternalCalibration ------------------------------------------------

test_that("test InternalCalibration",{
  ff = InternalCalibration$new()
  expect_true(!is.null(InternalCalibration$new()$fillCalibrants))
  expect_true(!is.null(InternalCalibration$new()$getCalibrationPoints))
  expect_true(!is.null(InternalCalibration$new()$calibrate))
})


# test TransitionTSVFile --------------------------------------------------

test_that("test TransitionTSVFile",{
  ff = TransitionTSVFile$new()
  expect_true(!is.null(TransitionTSVFile$new()$convertTargetedExperimentToTSV))
  expect_true(!is.null(TransitionTSVFile$new()$convertTSVToTargetedExperiment))
  expect_true(!is.null(TransitionTSVFile$new()$validateTargetedExperiment))

})


# test ProteaseDigestion --------------------------------------------------

test_that("test ProteaseDigestion",{
  ff = ProteaseDigestion$new()

  expect_true(!is.null(ProteaseDigestion$new()$getMissedCleavages))
  expect_true(!is.null(ProteaseDigestion$new()$setMissedCleavages))
  expect_true(!is.null(ProteaseDigestion$new()$digest))
  expect_true(!is.null(ProteaseDigestion$new()$peptideCount))

  ff$setMissedCleavages(5)
  expect_equal(ff$getMissedCleavages(), 5L)
})


# test EnzymaticDigestionLogModel -----------------------------------------

test_that("test EnzymaticDigestionLogModel",{
  ff = EnzymaticDigestionLogModel$new()

  expect_true(!is.null(EnzymaticDigestionLogModel$new()$getLogThreshold))
  expect_true(!is.null(EnzymaticDigestionLogModel$new()$setLogThreshold))
  expect_true(!is.null(EnzymaticDigestionLogModel$new()$digest))
  expect_true(!is.null(EnzymaticDigestionLogModel$new()$peptideCount))

  ff$setLogThreshold(0.25)
  expect_equivalent(ff$getLogThreshold(),0.25)
})


# test IDDecoyProbability -------------------------------------------------

test_that("test IDDecoyProbability",{
  ff = IDDecoyProbability$new()
  expect_true(!is.null(IDDecoyProbability$new()$apply))
})


# test FeatureGrouping ----------------------------------------------------

test_that("test FeatureGroupingAlgorithm",{
  expect_true(!is.null(FeatureGroupingAlgorithm$new()$getDefaults))
  expect_true(!is.null(FeatureGroupingAlgorithm$new()$getName))
  expect_true(!is.null(FeatureGroupingAlgorithm$new()$getParameters))
  expect_true(!is.null(FeatureGroupingAlgorithm$new()$setName))
  expect_true(!is.null(FeatureGroupingAlgorithm$new()$setParameters))
  expect_true(!is.null(FeatureGroupingAlgorithm$new()$transferSubelements))

  qt = FeatureGroupingAlgorithmQT$new()
  qt$getDefaults()
  qt$getParameters()
  qt$getName()
  expect_true(!is.null(qt$group))
  expect_true(!is.null(qt$setName))
  expect_true(!is.null(qt$setParameters))
  expect_true(!is.null(qt$transferSubelements))
})


# test FeatureMap ---------------------------------------------------------

test_that("test FeatureMap",{
  fm = FeatureMap$new()

  # clone and deepclone not supported for now
  fm_ = FeatureMap$new(fm)
  expect_true(fm_ == fm)

  check_UniqueIdInterface(fm)
  fm$clear()
  fm$clearUniqueId()
  fm$getIdentifier()
  fm$getLoadedFileType()
  fm$getLoadedFilePath()

  f = Feature$new()
  fm$push_back(f)

  expect_length(list(fm),1)
  expect_equal(fm$size(), 1)
  expect_true(fm[1] == f)

  fm$sortByIntensity()
  expect_equal(fm$size(), 1)
  expect_true(fm[1] == f)

  fm$sortByIntensity(FALSE)
  expect_equal(fm$size(), 1)
  expect_true(fm[1] == f)

  fm$sortByPosition()
  expect_equal(fm$size(), 1)
  expect_true(fm[1] == f)

  fm$sortByRT()
  expect_equal(fm$size(), 1)
  expect_true(fm[1] == f)

  fm$sortByMZ()
  expect_equal(fm$size(), 1)
  expect_true(fm[1] == f)

  fm$sortByOverallQuality()
  expect_equal(fm$size(), 1)
  expect_true(fm[1] == f)

  fm2 = FeatureMap$new()

  fm$swap(fm2)
  expect_equal(fm2$size(), 1)
  expect_true(fm2[1] == f)

  expect_equal(fm$size(), 0)

  fm2$updateRanges()

  expect_true(is_scalar_double(fm2$getMin()[1]))
  expect_true(is_scalar_double(fm2$getMin()[2]))
  expect_true(is_scalar_double(fm2$getMax()[1]))
  expect_true(is_scalar_double(fm2$getMinInt()))
  expect_true(is_scalar_double(fm2$getMaxInt()))

  expect_equal(fm2$getProteinIdentifications(),list())
  fm2$setUnassignedPeptideIdentifications(list())

  fm2$clear()
  expect_equal(fm2$size(), 0)

  dp = DataProcessing$new()
  fm2$setDataProcessing(c(dp))
  expect_equal(fm2$getDataProcessing(), c(dp))
  check_DataProcessing(dp)

  fm2$setUniqueIds()
  fm = (fm + fm)
  expect_false((fm+fm2) == fm)
})


# test FeatureXMLFile -----------------------------------------------------

test_that("test FeatureXMLFile",{
  fm = FeatureMap$new()
  fm$setUniqueIds()
  fh = FeatureXMLFile$new()
  fh$store("test.featureXML",fm)
  fh$load("test.featureXML",fm)
  fh = FileHandler$new()
  expect_true(fh$loadFeatures("test.featureXML", fm))
})


# test ColumnHeader ----------------------------------------------------

test_that("test ColumnHeader",{
  fd = ColumnHeader$new()
  expect_true(is_scalar_character(fd$filename))
  expect_true(is_scalar_character(fd$label))
  expect_true(is_scalar_integer(fd$size))
})


# test FileHandler --------------------------------------------------------

test_that("test FileHandler",{
  mse = MSExperiment$new()
  fh = FileHandler$new()

  expect_invisible(fh$storeExperiment("test1.mzML", mse))
  expect_true(fh$loadExperiment("test1.mzML", mse))
  expect_invisible(fh$storeExperiment("test1.mzML", mse))
  expect_true(fh$loadExperiment("test1.mzML", mse))
  expect_invisible(fh$storeExperiment("test1.mzData", mse))
  expect_true(fh$loadExperiment("test1.mzData", mse))
})


# test CachedMzML ---------------------------------------------------------

test_that("tests CachedMzML",{

  mse = MSExperiment$new()
  s = MSSpectrum$new()
  mse$addSpectrum(s)

  # First load data and cache to disk
  CachedmzML$store("myCache.mzML", mse)

  # Now load data
  cfile = CachedmzML$new()
  CachedmzML$load("myCache.mzML", cfile)

  meta_data = cfile$getMetaData()
  expect_equal(cfile$getNrChromatograms(), 0)
  expect_equal(cfile$getNrSpectra(), 1)
})


# test IndexedMzMLFile ----------------------------------------------------

test_that("test IndexedMzMLFile",{

  mse = MSExperiment$new()
  s = MSSpectrum$new()
  mse$addSpectrum(s)

  # First load data and cache to disk
  MzMLFile$new()$store("tfile_idx.mzML", mse)

  # Now load data
  ih = IndexedMzMLHandler$new("tfile_idx.mzML")

  expect_equal(ih$getNrChromatograms(), 0)
  expect_equal(ih$getNrSpectra(), 1)

  s = ih$getMSSpectrumById(0)
  s2 = ih$getSpectrumById(0)
})


# test IDMapper -----------------------------------------------------------

test_that("test IDMapper",{

  idm = IDMapper$new()
  expect_true(!is.null(idm$annotate))
  idm$getDefaults()
  idm$setName("x")
  expect_equal(idm$getName(), "x")
  idm$setParameters(idm$getParameters())
})


# test IdXMLFile ----------------------------------------------------------

test_that("test IdXMLFile",{
  expect_true(!is.null(IdXMLFile$new()$load))
  expect_true(!is.null(IdXMLFile$new()$store))
})


# test PepXMLFile ---------------------------------------------------------

test_that("test PepXMLFile",{
  f = PepXMLFile$new()

  expect_true(!is.null(PepXMLFile$new()$load))
  expect_true(!is.null(PepXMLFile$new()$store))
})


# test ProtXMLFile --------------------------------------------------------
test_that("test ProtXMLFile",{
  f = ProtXMLFile$new()

  expect_true(!is.null(ProtXMLFile$new()$load))
  expect_true(!is.null(ProtXMLFile$new()$store))

})


# test DTA2DFile ----------------------------------------------------------
test_that("test DTA2DFile",{
  f = DTA2DFile$new()

  expect_true(!is.null(DTA2DFile$new()$load))
  expect_true(!is.null(DTA2DFile$new()$store))
})

# test DTAFile ----------------------------------------------------------
test_that("test DTAFile",{
  f = DTAFile$new()

  expect_true(!is.null(DTAFile$new()$load))
  expect_true(!is.null(DTAFile$new()$store))
})


# test EDTAFile ----------------------------------------------------------
test_that("test EDTAFile",{
  f = EDTAFile$new()

  expect_true(!is.null(EDTAFile$new()$load))
  expect_true(!is.null(EDTAFile$new()$store))
})


# test KroenikFile ----------------------------------------------------------
test_that("test KroenikFile",{
  f = KroenikFile$new()

  expect_true(!is.null(KroenikFile$new()$load))
  expect_true(!is.null(KroenikFile$new()$store))
})


# test MSPFile ------------------------------------------------------------
test_that("test MSPFile",{
  f = MSPFile$new()

  expect_true(!is.null(MSPFile$new()$load))
  expect_true(!is.null(MSPFile$new()$store))
})


# test MzIdentMLFile ------------------------------------------------------

test_that("test MzIdentMLFile",{
  f = MzIdentMLFile$new()

  expect_true(!is.null(MzIdentMLFile$new()$load))
  expect_true(!is.null(MzIdentMLFile$new()$store))
  expect_true(!is.null(MzIdentMLFile$new()$isSemanticallyValid))
})


# test MzTabFile ----------------------------------------------------------

test_that("test MzTabFile",{
  f = MzTabFile$new()

  expect_true(!is.null(MzTabFile$new()$store))
})


# test MzTab ----------------------------------------------------------
test_that("test MzTab",{
  f = MzTab$new()
  expect_true(!is.null(f))
})


# test InstrumentSettings -------------------------------------------------
test_that("test InstrumentSettings",{

  ins = InstrumentSettings$new()
  check_MetaInfoInterface(ins)
  ins$setPolarity(IonSource$Polarity()$NEGATIVE)
  expect_equal(ins$getPolarity(), IonSource$Polarity()$NEGATIVE)
  expect_true(ins == ins)
  expect_false(ins != ins)
})


# test ContactPerson ------------------------------------------------------
test_that("test ContactPerson",{

  ins = ContactPerson$new()

  expect_equal(ins$getFirstName(), "")
  ins$setFirstName("test")
  ins$getLastName()
  ins$setLastName("test")
  ins$setName("Testy Test")
  ins$getInstitution()
  ins$setInstitution("test")
  ins$getEmail()
  ins$setEmail("test")
  ins$getURL()
  ins$setURL("test")
  ins$getAddress()
  ins$setAddress("test")
  ins$getContactInfo()
  ins$setContactInfo("test")

})


# test DocumentIdentifier -------------------------------------------------

test_that("test DocumentIdentifier",{
  ins = DocumentIdentifier$new()

  ins$setIdentifier("test")
  expect_equal(ins$getIdentifier(),"test")
  ins$getLoadedFilePath()
  ins$getLoadedFileType()
})


# test Gradient -----------------------------------------------------------
test_that("test Gradient",{
  ins = Gradient$new()

  ins$addEluent("test")
  ins$clearEluents()
  expect_length(ins$getEluents(), 0)
  ins$addEluent("test")
  expect_length(ins$getEluents(), 1)

  ins$clearTimepoints()
  ins$addTimepoint(5)
  expect_equal(ins$getTimepoints(), 5)

  ins$setPercentage("test", 5, 20)
  expect_equal(ins$getPercentage("test", 5), 20)
  ins$clearPercentages()
  ins$isValid()
})



# test HPLC ---------------------------------------------------------------

test_that("test HPLC",{
  ins = HPLC$new()

  ins$setInstrument("test")
  expect_equal(ins$getInstrument(), "test")
  ins$setColumn("test")
  expect_equal(ins$getColumn(), "test")
  ins$setTemperature(6)
  expect_equal(ins$getTemperature(), 6)
  ins$setPressure(6)
  expect_equal(ins$getPressure(),6)
  ins$setFlux(8)
  expect_equal(ins$getFlux(), 8)
  ins$setComment("test")
  expect_equal(ins$getComment(), "test")

  g = Gradient$new()
  ins$setGradient(g)
  ins$getGradient()

})


# test Instrument ---------------------------------------------------------

test_that("test Instrument",{
  ins = Instrument$new()

  ins$setName("test")
  expect_equal(ins$getName(), "test")
  ins$setVendor("test")
  expect_equal(ins$getVendor(), "test")
  ins$setModel("test")
  expect_equal(ins$getModel(),"test")
  ins$setCustomizations("test")
  expect_equal(ins$getCustomizations(), "test")

  ion_sources = lapply(1:5, function(x) IonSource$new())
  ins$setIonSources(ion_sources)
  src <- ins$getIonSources()
  expect_true(all(sapply(src, function(s) is.R6(s) && class(s)[1]=="IonSource")))
  expect_length(src,5)

  mass_analyzers = lapply(1:5, function(x) MassAnalyzer$new())
  ins$setMassAnalyzers(mass_analyzers)
  ma <- ins$getMassAnalyzers()
  expect_true(all(sapply(ma, function(s) is.R6(s) && class(s)[1]=="MassAnalyzer")))

  ion_detectors = lapply(1:5, function(x) IonDetector$new())
  ins$setIonDetectors(ion_detectors)
  expect_length(ins$getIonDetectors(), 5)

  s = Software$new()
  ins$setSoftware(s)
  ins$getSoftware()

})


# test IonDetector --------------------------------------------------------

test_that("test IonDetector",{
  ins = IonDetector$new()

  m = IonDetector$AcquisitionMode()$ACQMODENULL
  ins$setAcquisitionMode(m)
  expect_equal(ins$getAcquisitionMode(), m)

  ins$setResolution(8.0)
  expect_equal(ins$getResolution(), 8)

  ins$setADCSamplingFrequency(8.0)
  expect_equal(ins$getADCSamplingFrequency(), 8)

  ins$setOrder(8)
  expect_equal(ins$getOrder(), 8)
})



# test IonSource ----------------------------------------------------------

test_that("test IonSource",{
  ins = IonSource$new()

  p = IonSource$Polarity()$POSITIVE
  ins$setPolarity(p)
  expect_equal(ins$getPolarity(), p)

  i = IonSource$InletType()$INLETNULL
  ins$setInletType(i)
  expect_equal(ins$getInletType(), i)

  i = IonSource$IonizationMethod()$ESI
  ins$setIonizationMethod(i)
  expect_equal(ins$getIonizationMethod(), i)

  ins$setOrder(5)
  expect_equal(ins$getOrder(), 5)

})


# test MassAnalyzer -------------------------------------------------------

test_that("test MassAnalyzer",{

  ins = MassAnalyzer$new()

  ma = MassAnalyzer$AnalyzerType()$QUADRUPOLE
  ins$setType(ma)
  expect_equal(ins$getType(), ma)

  res = MassAnalyzer$ResolutionMethod()$FWHM
  ins$setResolutionMethod(res)
  expect_equal(ins$getResolutionMethod(), res)

  res = MassAnalyzer$ResolutionType()$CONSTANT
  ins$setResolutionType(res)
  expect_equal(ins$getResolutionType(), res)

  res = MassAnalyzer$ScanDirection()$UP
  ins$setScanDirection(res)
  expect_equal(ins$getScanDirection(), res)

  res = MassAnalyzer$ScanLaw()$LINEAR
  ins$setScanLaw(res)
  expect_equal(ins$getScanLaw(), res)

  res = MassAnalyzer$ReflectronState()$ON
  ins$setReflectronState(res)
  expect_equal(ins$getReflectronState(), res)

  ins$setResolution(5.0)
  ins$getResolution()
  ins$setAccuracy(5.0)
  ins$getAccuracy()
  ins$setScanRate(5.0)
  ins$getScanRate()
  ins$setScanTime(5.0)
  ins$getScanTime()
  ins$setTOFTotalPathLength(5.0)
  ins$getTOFTotalPathLength()
  ins$setIsolationWidth(5.0)
  ins$getIsolationWidth()
  ins$setFinalMSExponent(5)
  ins$getFinalMSExponent()
  ins$setMagneticFieldStrength(5.0)
  ins$getMagneticFieldStrength()
  ins$setOrder(5)
  ins$getOrder()

})


# test Sample -------------------------------------------------------------

test_that("test Sample",{

  ins = Sample$new()

  ins$setName("test")
  ins$getName()
  ins$setOrganism("test")
  ins$getOrganism()
  ins$setNumber("test")
  ins$getNumber()
  ins$setComment("test")
  ins$getComment()

  state = Sample$SampleState()$LIQUID
  ins$setState(state)
  expect_equal(ins$getState(), state)
  ins$setMass(42.0)
  expect_equal(ins$getMass(), 42)
  ins$setVolume(42.0)
  expect_equal(ins$getVolume(), 42)
  ins$setConcentration(42.0)
  expect_equal(ins$getConcentration(), 42)

  a = ins$getSubsamples()
  ins$setSubsamples(a)

  has_exception = F

  tryCatch({
    ins$removeTreatment(0)
  }, error = function(e){
    has_exception <<- T
  })

  expect_true(has_exception)
  expect_equal(ins$countTreatments(), 0)

})


# test LogType ------------------------------------------------------------

test_that("test LogType",{
  expect_true(is_scalar_integer(LogType()$CMD))
  expect_true(is_scalar_integer(LogType()$GUI))
  expect_true(is_scalar_integer(LogType()$NONE))
})


# test MSExperiment -------------------------------------------------------
test_that("test MSExperiment",{

  mse = MSExperiment$new()
  mse_ = MSExperiment$new(mse)
  expect_true(mse_ == mse)

  check_MetaInfoInterface(mse)
  mse$updateRanges()
  mse$sortSpectra(T)
  expect_true(is_scalar_double(mse$getMaxRT()))
  expect_true(is_scalar_double(mse$getMinRT()))
  expect_true(is_scalar_double(mse$getMaxMZ()))
  expect_true(is_scalar_double(mse$getMinMZ()))
  expect_true(is_scalar_character(mse$getLoadedFilePath()))
  expect_true(is_scalar_double(mse$getMinInt()))
  expect_true(is_scalar_double(mse$getMaxInt()))

  expect_true(is_scalar_double(mse$getMin()[1]))
  expect_true(is_scalar_double(mse$getMin()[2]))
  expect_true(is_scalar_double(mse$getMax()[1]))
  expect_true(is_scalar_double(mse$getMax()[2]))
  mse$setLoadedFilePath("")
  expect_equal(mse$size(), 0)

  mse$getIdentifier()
  mse$getLoadedFileType()
  mse$getLoadedFilePath()

  mse$addSpectrum(MSSpectrum$new())
  expect_equal(mse$size(), 1)

  expect_true(!is.null(mse[1]))

  expect_true(typeof(list(mse)) == "list")

  expect_true(mse == mse)
  expect_false(mse != mse)

  expect_true(mse$getSize() >= 0)
  mse$isSorted()
})


# test MSQuantifications --------------------------------------------------

test_that("test MSQuantifications",{

  msq = MSQuantifications$new()
  expect_true(msq == msq)

  expect_false(msq != msq)
  msq$setConsensusMaps(msq$getConsensusMaps())
  summary = msq$getAnalysisSummary()
  msq$setDataProcessingList(msq$getDataProcessingList())
  expect_true(typeof(msq$getAssays()) == "list")
  expect_true(typeof(msq$getFeatureMaps()) == "list")
  msq$setAnalysisSummaryQuantType(MSQuantifications$QUANT_TYPES()$LABELFREE)

  msq$addConsensusMap(ConsensusMap$new())
  msq$assignUIDs()

})


# tests for SpectrumSettings ----------------------------------------------
check_SpectrumSettings <- function(s){

  s <- SpectrumSettings$new()
  expect_true(s$getType() %in% c(SpectrumSettings$SpectrumType()$UNKNOWN))
  expect_equal(class(s$getInstrumentSettings())[1], "InstrumentSettings")
  expect_equal(class(s$getSourceFile())[1], "SourceFile")
  expect_equal(typeof(s$getPeptideIdentifications()), "list")
  expect_equal(typeof(s$getDataProcessing()), "list")

  s$setAcquisitionInfo(s$getAcquisitionInfo())
  s$setInstrumentSettings(s$getInstrumentSettings())
  s$setSourceFile(s$getSourceFile())
  s$setPeptideIdentifications(s$getPeptideIdentifications())
  s$setDataProcessing(s$getDataProcessing())
  s$setComment(s$getComment())
  s$setPrecursors(s$getPrecursors())
  s$setProducts(s$getProducts())
  s$setType(s$getType())
  s$setNativeID(s$getNativeID())
  s$setType(s$getType())
  if (py_to_r(py_builtin$isinstance(r_to_py(s), Pymod$SpectrumSettings))){
    s$unify(s)
  }
}



# test MSSpectrum ---------------------------------------------------------

test_that("test MSSpectrum",{

  spec = MSSpectrum$new()
  spec_ = MSSpectrum$new(spec)
  expect_true(spec == spec_)

  check_MetaInfoInterface(spec)

  check_SpectrumSettings(spec)

  spec$setRT(3.0)
  expect_true(spec$getRT() == 3.0)
  spec$setMSLevel(2)
  expect_true(spec$getMSLevel() == 2)
  spec$setName("spec")
  expect_true(spec$getName() == "spec")

  p = Peak1D$new()
  p$setMZ(1000.0)
  p$setIntensity(200.0)
  spec$push_back(p)

  expect_equal(spec$size(), 1)
  expect_true(spec[1] == p)

  spec$updateRanges()
  expect_true(is_scalar_integer(spec$findNearest(0.0)))

  expect_true(is_scalar_double(spec$getMin()[1]))
  expect_true(is_scalar_double(spec$getMax()[1]))
  expect_true(is_scalar_double(spec$getMinInt()))
  expect_true(is_scalar_double(spec$getMaxInt()))

  expect_true(spec == spec)
  expect_false(spec != spec)

  peaks = spec$get_peaks()
  mz = peaks[[1]]; ii = peaks[[2]]
  expect_true(length(mz) == length(ii))
  expect_true(length(mz) == 1)

  spec$set_peaks(list(mz, ii))
  peaks = spec$get_peaks()
  mz0 = peaks[[1]]; ii0 = peaks[[2]]
  expect_true(mz0 == mz)
  expect_true(ii0 == ii)

  expect_true(as.integer(spec$isSorted()) %in% c(0,1))

  spec$clear(F)
  p = Peak1D$new()
  p$setMZ(1000.0)
  p$setIntensity(200.0)
  spec$push_back(p)
  p = Peak1D$new()
  p$setMZ(2000.0)
  p$setIntensity(400.0)
  spec$push_back(p)

  peaks = spec$get_peaks()
  mz = peaks[[1]]; ii = peaks[[2]]
  expect_true(spec[1]$getMZ() == 1000.0)
  expect_true(spec[2]$getMZ() == 2000.0)
  expect_true(spec[1]$getIntensity() == 200.0)
  expect_true(spec[2]$getIntensity() == 400.0)
  expect_true(mz[1] == 1000.0)
  expect_true(mz[2] == 2000.0)
  expect_true(ii[1] == 200.0)
  expect_true(ii[2] == 400.0)

  spec$clear(F)
  data_mz = c(5.0, 8.0)
  data_i = c(50.0, 80.0)
  spec$set_peaks(list(data_mz,data_i))

  peaks = spec$get_peaks()
  mz = peaks[[1]]; ii = peaks[[2]]
  expect_true(spec[1]$getMZ() == 5.0)
  expect_true(spec[2]$getMZ() == 8.0)
  expect_true(spec[1]$getIntensity() == 50.0)
  expect_true(spec[2]$getIntensity() == 80.0)
  expect_true(mz[1] == 5.0)
  expect_true(mz[2] == 8.0)
  expect_true(ii[1] == 50.0)
  expect_true(ii[2] == 80.0)

  # Fast
  spec$clear(F)
  data_mz = c(5.0, 8.0)
  data_i = c(50.0, 80.0)
  spec$set_peaks( list(data_mz,data_i) )

  peaks = spec$get_peaks()
  mz = peaks[[1]]; ii = peaks[[2]]
  expect_true(spec[1]$getMZ() == 5.0)
  expect_true(spec[2]$getMZ() == 8.0)
  expect_true(spec[1]$getIntensity() == 50.0)
  expect_true(spec[2]$getIntensity() == 80.0)
  expect_true(mz[1] == 5.0)
  expect_true(mz[2] == 8.0)
  expect_true(ii[1] == 50.0)
  expect_true(ii[2] == 80.0)

  ###################################
  # get data arrays
  ###################################
  expect_true(length(spec$getStringDataArrays()) == 0)
  string_da = c(StringDataArray$new())
  string_da[[1]]$push_back("hello")
  string_da[[1]]$push_back("world")

  string_da <- c(string_da,StringDataArray$new())

  string_da[[2]]$push_back("other")
  spec$setStringDataArrays( string_da )
  expect_true(length(spec$getStringDataArrays()) == 2)
  expect_true(spec$getStringDataArrays()[[1]][1] == "hello")
  expect_true(spec$getStringDataArrays()[[2]][1] == "other")


  spec = MSSpectrum$new()
  expect_true(length(spec$getIntegerDataArrays()) == 0)
  int_da = list(IntegerDataArray$new())
  int_da[[1]]$push_back(5)
  int_da[[1]]$push_back(6)
  int_da <- c(int_da, IntegerDataArray$new() )
  int_da[[2]]$push_back(8)
  spec$setIntegerDataArrays( int_da )
  expect_true(length(spec$getIntegerDataArrays()) == 2)
  expect_true(spec$getIntegerDataArrays()[[1]][1] == 5)
  expect_true(spec$getIntegerDataArrays()[[2]][1] == 8)

  spec = MSSpectrum$new()
  # c(5,8,42) would also work fine.
  data = c(5L, 8L, 42L)
  int_da = list(IntegerDataArray$new())
  int_da[[1]]$set_data(data)
  spec$setIntegerDataArrays( int_da )
  expect_true(length(spec$getIntegerDataArrays()) == 1)
  expect_equivalent(spec$getIntegerDataArrays()[[1]][1], 5L)
  expect_equivalent(spec$getIntegerDataArrays()[[1]][3], 42L)
  expect_true(length(int_da[[1]]$get_data() ) == 3)

  spec = MSSpectrum$new()
  expect_true(length(spec$getFloatDataArrays()) == 0)
  f_da = list(FloatDataArray$new())
  f_da[[1]]$push_back(5.0)
  f_da[[1]]$push_back(6.0)
  f_da <- c(f_da,FloatDataArray$new())
  f_da[[2]]$push_back(8.0)
  spec$setFloatDataArrays( f_da )
  expect_true( length(spec$getFloatDataArrays()) == 2)
  expect_true( spec$getFloatDataArrays()[[1]][1] == 5.0 )
  expect_true( spec$getFloatDataArrays()[[2]][1] == 8.0 )

  spec = MSSpectrum$new()
  data = c(5, 8, 42)
  f_da = c(FloatDataArray$new())
  f_da[[1]]$set_data(data)
  spec$setFloatDataArrays( f_da )
  expect_equal(length(spec$getFloatDataArrays()), 1)
  expect_equivalent(spec$getFloatDataArrays()[[1]][1], 5.0)
  expect_equivalent(spec$getFloatDataArrays()[[1]][3], 42.0)
  expect_true(length(f_da[[1]]$get_data() ) == 3 )
})


# test StringDataArray ----------------------------------------------------

test_that("test StringDataArray",{
  da = StringDataArray$new()
  expect_equal(da$size(), 0)

  da$push_back("hello")
  da$push_back("world")
  expect_equal(da$size(), 2)
  expect_equal(da[1], "hello")
  expect_equal(da[2], "world")
  da[2] = "hello world"
  expect_equal(da[2], "hello world")
  da$clear()

  expect_true(da$size() == 0)
  da$push_back("hello")
  expect_true(da$size() == 1)
  da$resize(3)
  da[1] = "hello"
  da[2] = ""
  da[3] = "world"
  expect_true(da$size() == 3)

})

# test IntegerDataArray ---------------------------------------------------

test_that("test IntegerDataArray",{
  da = IntegerDataArray$new()
  expect_true(da$size() == 0)
  da$push_back(1)
  da$push_back(4)
  expect_true(da$size() == 2)

  expect_true(da[1] == 1)
  expect_true(da[2] == 4)
  da[2] = 7
  expect_equivalent(da[2],7L)

  da$clear()

  expect_true(da$size() == 0)
  da$push_back(1)
  expect_true(da$size() == 1)
  da$resize(3)

  da[1] = 1
  da[2] = 2
  da[3] = 3
  expect_true(da$size() == 3)


  q = da$get_data()
  q = c(q, 4L)
  da$set_data(q)
  expect_equal(da$size(), 4)
})

# test FloatDataArray -----------------------------------------------------
test_that("",{

  da = FloatDataArray$new()
  expect_true(da$size() == 0)
  da$push_back(1.0)
  da$push_back(4.0)

  expect_true(da$size() == 2)
  expect_true(da[1] == 1.0)
  expect_true(da[2] == 4.0)
  da[1] = 7.0
  expect_equivalent(da[1],7.0)

  da$clear()
  expect_true(da$size() == 0)
  da$push_back(1.0)
  expect_true(da$size() == 1)

  da$resize(3)
  da[1] = 1.0
  da[2] = 2.0
  da[3] = 3.0
  expect_equal(da$size(),3)

  q = da$get_data()
  q <- c(q,4.0)
  da$set_data(q)
  expect_true(da$size() == 4)
})

# test MSChromatogram ---------------------------------------------------------

test_that("test MSChromatogram",{

  chrom = MSChromatogram$new()
  chrom_ = MSChromatogram$new(chrom)
  expect_true(chrom == chrom_)

  check_MetaInfoInterface(chrom)

  chrom$setName("chrom")
  expect_equal(chrom$getName(), "chrom")

  p = ChromatogramPeak$new()

  p$setRT(1000.0)
  p$setIntensity(200.0)

  chrom$push_back(p)
  expect_equal(chrom$size(), 1)
  expect_true(chrom[1] == p)

  chrom$updateRanges()
  expect_true(is_scalar_integer(chrom$findNearest(0.0)))

  expect_true(is_scalar_double(chrom$getMin()[1]))
  expect_true(is_scalar_double(chrom$getMax()[1]))
  expect_true(is_scalar_double(chrom$getMinInt()))
  expect_true(is_scalar_double(chrom$getMaxInt()))

  expect_true(chrom == chrom)
  expect_false(chrom != chrom)

  peaks = chrom$get_peaks()
  mz = peaks[[1]]; ii = peaks[[2]]
  expect_true( length(mz) == length(ii) )
  expect_true(length(mz) == 1)

  chrom$set_peaks(list(mz, ii))
  peaks = chrom$get_peaks()
  mz0 = peaks[[1]]; ii0 = peaks[[2]]
  expect_true(mz0 == mz)
  expect_true(ii0 == ii)

  expect_true(as.integer(chrom$isSorted()) %in%  c(0,1))

  chrom$clear(as.integer(F))
  p = ChromatogramPeak$new()
  p$setRT(1000.0)
  p$setIntensity(200.0)
  chrom$push_back(p)
  p = ChromatogramPeak$new()
  p$setRT(2000.0)
  p$setIntensity(400.0)
  chrom$push_back(p)

  peaks = chrom$get_peaks()
  mz = peaks[[1]]; ii = peaks[[2]]
  expect_true(chrom[1]$getRT() == 1000.0)
  expect_true(chrom[2]$getRT() == 2000.0)
  expect_true(chrom[1]$getIntensity() == 200.0)
  expect_true(chrom[2]$getIntensity() == 400.0)
  expect_true(mz[1] == 1000.0)
  expect_true(mz[2] == 2000.0)
  expect_true(ii[1] == 200.0)
  expect_true(ii[2] == 400.0)

  chrom$clear(as.integer(F))
  data_mz =  c(5.0, 8.0)
  data_i = c(50.0, 80.0)
  chrom$set_peaks( list(data_mz,data_i) )

  peaks = chrom$get_peaks()
  mz = peaks[[1]]; ii = peaks[[2]]
  expect_true(chrom[1]$getRT() == 5.0)
  expect_true(chrom[2]$getRT() == 8.0)
  expect_true(chrom[1]$getIntensity() == 50.0)
  expect_true(chrom[2]$getIntensity() == 80.0)
  expect_true(mz[1] == 5.0)
  expect_true(mz[2] == 8.0)
  expect_true(ii[1] == 50.0)
  expect_true(ii[2] == 80.0)

})






# test MRMFeature ---------------------------------------------------------

test_that("test MRMFeature",{

  mrmfeature = MRMFeature$new()
  f = Feature$new()

  fs = mrmfeature$getFeatures()
  expect_true(length(fs) == 0)

  mrmfeature$addFeature(f, "myFeature")
  fs = mrmfeature$getFeatures()
  expect_length(fs, 1)
  expect_true(!is.null(mrmfeature$getFeature("myFeature")))
  slist = list()
  mrmfeature$getFeatureIDs(slist)
  expect_length(slist, 1)

  mrmfeature$addPrecursorFeature(f, "myFeature_Pr0")
  slist = list()
  mrmfeature$getPrecursorFeatureIDs(slist)
  expect_length(slist, 1)
  expect_true(!is.null(mrmfeature$getPrecursorFeature("myFeature_Pr0")))

  s = mrmfeature$getScores()
  expect_lt(abs(s$yseries_score - 0.0), 1e-4)
  s$yseries_score = 4.0
  mrmfeature$setScores(s)
  s2 = mrmfeature$getScores()
  expect_lt(abs(s2$yseries_score - 4.0), 1e-4)

})


# test ConfidenceScoring --------------------------------------------------

test_that("test ConfidenceScoring",{
  scoring = ConfidenceScoring$new()
  expect_true(!is.null(scoring))
})


# test MRMDecoy -----------------------------------------------------------

test_that("test MRMDecoy",{
  mrmdecoy = MRMDecoy$new()
  expect_true(!is.null(mrmdecoy))
  expect_true(!is.null(MRMDecoy$new()$generateDecoys))
})


# test MRMTransitionGroup -------------------------------------------------

test_that("test MRMTransitionGroup",{
  mrmgroup = MRMTransitionGroupCP$new()
  expect_true(!is.null(mrmgroup))

  mrmgroup$setTransitionGroupID("this_id")
  expect_equal(mrmgroup$getTransitionGroupID(), "this_id")

  expect_length(mrmgroup$getTransitions(), 0)
  mrmgroup$addTransition(ReactionMonitoringTransition$new(), "tr1")
  expect_length(mrmgroup$getTransitions(), 1)
})


# test ReactionMonitoringTransition ---------------------------------------

test_that("test ReactionMonitoringTransition",{
  tr = ReactionMonitoringTransition$new()
  expect_true(!is.null(tr))
})



# test TargetedExperiment -------------------------------------------------

test_that("test TargetedExperiment",{

  m = TargetedExperiment$new()
  m_ = TargetedExperiment$new(m)
  expect_true(m_ == m)

  m$clear(T)
  m$setCVs(m$getCVs())

  targeted = m

  targeted$setCVs(targeted$getCVs())
  targeted$setTargetCVTerms(targeted$getTargetCVTerms())
  targeted$setPeptides(targeted$getPeptides())
  targeted$setProteins(targeted$getProteins())
  targeted$setTransitions(targeted$getTransitions())

  expect_true(m == m)
  expect_false(m != m)
})


# test TargetedExperimentHelper -------------------------------------------

test_that("test TargetedExperimentHelper",{

  rtu = RetentionTime$RTUnit()
  rtu = RetentionTime$RTUnit()$SECOND
  rtu = RetentionTime$RTUnit()$MINUTE
  rtt = RetentionTime$RTType()
  rtt = RetentionTime$RTType()$LOCAL
  rtt = RetentionTime$RTType()$NORMALIZED
  rtt = RetentionTime$RTType()$IRT

  rt = RetentionTime$new()
  expect_true(!is.null(rt$software_ref))
  expect_false(rt$isRTset())
  rt$setRT(5.0)
  rt$retention_time_unit = RetentionTime$RTUnit()$SECOND
  rt$retention_time_type = RetentionTime$RTType()$NORMALIZED
  expect_true(rt$isRTset())
  expect_equal(rt$getRT(), 5.0)

  p = Peptide$new()
  expect_false(is.null(p$rts))
  expect_false(is.null(p$id))
  expect_false(is.null(p$protein_refs))
  expect_false(is.null(p$evidence))
  expect_false(is.null(p$sequence))
  expect_false(is.null(p$mods))

  expect_false(p$hasCharge())
  p$setChargeState(5)
  expect_true(!is.null(p$hasCharge()))
  expect_equal(p$getChargeState(), 5)

  expect_false(p$hasRetentionTime())
  p$rts = list(rt)
  expect_true(p$hasRetentionTime())
  expect_equal(p$getRetentionTime(), 5.0)
  expect_equal(p$getRetentionTimeUnit(), RetentionTime$RTUnit()$SECOND)
  expect_equal(p$getRetentionTimeType(), RetentionTime$RTType()$NORMALIZED)

  c = Compound$new()
  expect_false(is.null(c$rts))
  expect_false(is.null(c$id))
  expect_false(is.null(c$molecular_formula))
  expect_false(is.null(c$smiles_string))
  expect_false(is.null(c$theoretical_mass))

  expect_false(c$hasCharge())
  c$setChargeState(5)
  expect_true(c$hasCharge())
  expect_equal(c$getChargeState(), 5)

  expect_false(c$hasRetentionTime())
  c$rts = list(rt)
  expect_true(c$hasRetentionTime())
  expect_equal(c$getRetentionTime(), 5.0)
  expect_true(c$getRetentionTimeUnit() == RetentionTime$RTUnit()$SECOND)
  expect_true(c$getRetentionTimeType() == RetentionTime$RTType()$NORMALIZED)

})



# test MapAlignment -------------------------------------------------------

test_that("test MapAlignment",{

  ma = MapAlignmentAlgorithmPoseClustering$new()
  expect_equal(ma$getDefaults(), Param$new())
  expect_equal(ma$getParameters(), Param$new())
  expect_true(is_scalar_character(ma$getName()))

  ma$setName(ma$getName())

  ma$getDefaults()
  ma$getParameters()

  ma$setParameters(ma$getDefaults())

  ma$setReference
  ma$align

  MapAlignmentTransformer$transformRetentionTimes
})


# test MatrixDouble -------------------------------------------------------

test_that("test MatrixDouble",{

  m = MatrixDouble$new()
  n = 90
  m$resize(n-1,n+2, 5.0)
  expect_equal(m$rows(), 89)
  expect_equal(m$cols(), 92)

  rows = n-1
  cols = n+2
  test = c()

  for(i in 0:(rows-1)){
    for(j in 0:(cols-1)){
      test <- c(test, m$getValue(i, j))
    }
  }

  # testm = np.asarray(test)
  # testm = testm.reshape(rows, cols)
  testm = array(test,c(rows,cols))

  expect_equal(sum(testm), 40940.0)
  expect_equal(sum(testm), (n-1)*(n+2)*5)

  matrix = m$get_matrix()
  expect_equal(sum(matrix), 40940.0)
  expect_equal(sum(matrix), (n-1)*(n+2)*5)

  matrix_view = m$get_matrix_as_view()
  expect_equal(sum(matrix_view), 40940.0)
  expect_equal(sum(matrix_view), (n-1)*(n+2)*5)


  # Column = 3 / Row = 5
  ## Now change a value:

  expect_equal(m$getValue(3, 5), 5.0)
  m$setValue(3, 5, 8.0)
  expect_equal(m$getValue(3, 5), 8.0)

  mat = m$get_matrix_as_view()

  # P.S. mat is an R matrix; indexing starts from (1,1)
  expect_equal(mat[4, 6], 8.0)

  mat = m$get_matrix()
  expect_true(m$getValue(3, 5) == 8.0)
  expect_true(mat[4, 6] == 8.0)

  # Whatever we change here gets changed in the raw data as well
  m$setValue(1, 6, 11.0)
  expect_equal(m$getValue(1, 6), 11.0)
  matrix_view = m$get_matrix_as_view()
  expect_equal(matrix_view[2, 7], 11.0)

  m$clear()
  expect_equal(m$rows(), 0)
  expect_true(m$cols() == 0)

  mat[4, 7] = 9.0
  m$set_matrix(mat)
  expect_true(m$getValue(3, 5) == 8.0)
  expect_true(m$getValue(3, 6) == 9.0)
})



# test MapAlignmentAlgorithmIdentification --------------------------------

test_that("test MapAlignmentAlgorithmIdentification",{
  ma = MapAlignmentAlgorithmIdentification$new()
  expect_false(is.null(MapAlignmentAlgorithmIdentification$new()$align))
  expect_false(is.null(MapAlignmentAlgorithmIdentification$new()$setReference))
})

# test MapAlignmentTransformer --------------------------------

test_that("test MapAlignmentTransformer",{
  ma = MapAlignmentTransformer$new()
  expect_false(is.null(MapAlignmentTransformer$new()$transformRetentionTimes))
})


# test testMxxxFile -------------------------------------------------------

test_that("test testMxxxFile",{

  mse = MSExperiment$new()
  s = MSSpectrum$new()
  mse$addSpectrum(s)

  fh = MzDataFile$new()

  check_ProgressLogger(fh)
  fh$store("test.mzData", mse)
  fh$load("test.mzData", mse)

  fh$setOptions(fh$getOptions())

  fh = MzMLFile$new()
  check_ProgressLogger(fh)

  fh$store("test.mzML", mse)
  fh$load("test.mzML", mse)
  fh$setOptions(fh$getOptions())

  myStr = String$new()
  fh$storeBuffer(myStr, mse)
  expect_equal(nchar(myStr$toString()), 5269)
  mse2 = MSExperiment$new()
  fh$loadBuffer(myStr, mse2)
  expect_true(mse2 == mse)
  expect_equal(mse2$size(), 1)

  fh = MzXMLFile$new()
  check_ProgressLogger(fh)
  fh$store("test.mzXML", mse)
  fh$load("test.mzXML", mse)
  fh$setOptions(fh$getOptions())

  fh = MzQuantMLFile$new()
  fh$isSemanticallyValid
  fh$load
  fh$store
})


# test ParamXMLFile -------------------------------------------------------

test_that("test ParamXMLFile",{
  fh = ParamXMLFile$new()
  p = Param$new()
  expect_invisible(fh$store("test.ini", p))
  expect_invisible(fh$load("test.ini", p))
})


# test Peak ---------------------------------------------------------------

test_that("test Peak1D/2D",{

  p1 = Peak1D$new()
  p1$setIntensity(12.0)
  expect_true(p1$getIntensity() == 12.0)
  p1$setMZ(13.0)
  p1$getMZ() == 13.0

  expect_true(p1 == p1)
  expect_false(p1 != p1)

  p2 = Peak2D$new()
  expect_true(p2 == p2)
  expect_false(p2 != p2)
  p2$setIntensity(22.0)
  expect_equal(p2$getIntensity(), 22.0)
  p2$setMZ(23.0)
  expect_equal(p2$getMZ(), 23.0)
  p2$setRT(45.0)
  expect_true(p2$getRT() == 45.0)
})


# test MSNumpressCoder ----------------------------------------------------

test_that("test NumpressCoder",{

  np = MSNumpressCoder$new()

  nc = NumpressConfig$new()
  nc$np_compression = MSNumpressCoder$NumpressCompression()$LINEAR
  nc$estimate_fixed_point = T
  tmp = String$new()
  out = vector()
  inp =  c(1.0, 2.0, 3.0)
  np$encodeNP(inp, tmp, T, nc)

  res = tmp$toString()
  expect_true(length(res) != 0)
  expect_true(res != "")
  np$decodeNP(res, out, T, nc)
  expect_length(out, 3)
  expect_equal(out, inp)


  res = ""
  has_error = F

  tryCatch({
    np$decodeNP(inp, res, T, nc)
  }, error = function(e){
    has_error <<- T
  })

  expect_true(has_error)
})


# test NumpressConfig -----------------------------------------------------
test_that("test MSNumpressCoder",{

  n = MSNumpressCoder$new()
  np = NumpressConfig$new()
  np$np_compression = MSNumpressCoder$NumpressCompression()$LINEAR
  expect_equal(np$np_compression, MSNumpressCoder$NumpressCompression()$LINEAR)
  np$numpressFixedPoint = 4.2
  np$numpressErrorTolerance = 4.2
  np$estimate_fixed_point = T
  np$linear_fp_mass_acc = 4.2
  np$setCompression("linear")

})


# test Base64 -------------------------------------------------------------

test_that("test Base64",{

  b = Base64$new()
  out = String$new()
  inp =  c(1.0, 2.0, 3.0)
  # Here python expects an integer(0/1) for argument zlib_compression.
  b$encode64(inp, Base64$ByteOrder()$BYTEORDER_LITTLEENDIAN, out, F)
  res = out$toString()
  expect_true(length(res) != 0)
  expect_true(res != "")

  convBack = vector()
  b$decode64(res, Base64$ByteOrder()$BYTEORDER_LITTLEENDIAN, convBack, F)
  expect_true(all.equal(convBack,inp))

  # For 32 bit
  out = String$new()
  b$encode32(inp, Base64$ByteOrder()$BYTEORDER_LITTLEENDIAN, out, F)
  res = out$toString()
  expect_true(length(res) != 0)
  expect_true(res != "")

  convBack = list()
  b$decode32(res, Base64$ByteOrder()$BYTEORDER_LITTLEENDIAN, convBack, F)
  expect_true(all.equal(convBack, inp))

})


# test PeakFileOptions ----------------------------------------------------

test_that("test PeakFileOptions",{
  pfo = PeakFileOptions$new()
  pfo$addMSLevel
  pfo$clearMSLevels()
  expect_false(pfo$containsMSLevel(1))
  pfo$getCompression()
  pfo$getMSLevels()
  pfo$getMetadataOnly()
  pfo$getWriteSupplementalData()
  expect_false(pfo$hasMSLevels())
  pfo$setCompression
  pfo$setMSLevels
  pfo$setMetadataOnly
  pfo$setWriteSupplementalData

})


# test MRMMapping ---------------------------------------------------------
test_that("test MRMMapping",{

  p = MRMMapping$new()
  expect_true(!is.null(p$mapExperiment))
  e = MSExperiment$new()
  c = MSChromatogram$new()
  e$addChromatogram(c)
  expect_equal(e$getNrChromatograms(), 1)

  o = MSExperiment$new()
  t = TargetedExperiment$new()
  p$mapExperiment(e, t, o)
  expect_true(o$getNrChromatograms() == 0) # not so easy to test
})


# test PeakPickerMRM ------------------------------------------------------

test_that("test PeakPickerMRM",{
  p = PeakPickerMRM$new()
  p$pickChromatogram
  expect_true(T)
})


# test PeakPickerHiRes ----------------------------------------------------

test_that("PeakPickerHiRes",{
  p = PeakPickerHiRes$new()
  expect_true(!is.null(p$pick))
  expect_true(!is.null(p$pickExperiment))
})


# test PeakTypeEstimator --------------------------------------------------

test_that("test PeakTypeEstimator",{
  expect_equal(PeakTypeEstimator$new()$estimateType(MSSpectrum$new()), 0)
})


# test PeptideHit ---------------------------------------------------------
test_that("test PeptideHit",{

  ph = PeptideHit$new()
  expect_true(ph == ph)
  expect_false(ph != ph)

  ph = PeptideHit$new(1.0, 1, 0, AASequence$fromString("A"))
  check_MetaInfoInterface(ph)

  expect_true(length(ph$getPeptideEvidences()) == 0)
  expect_equal(ph$getPeptideEvidences(),list())

  pe = PeptideEvidence$new()
  pe$setProteinAccession('B_id')

  ph$addPeptideEvidence(pe)
  expect_true(length(ph$getPeptideEvidences()) == 1)
  expect_true(ph$getPeptideEvidences()[[1]]$getProteinAccession() == 'B_id')

  ph$setPeptideEvidences(c(pe,pe))
  expect_true(length(ph$getPeptideEvidences()) == 2)
  expect_true(ph$getPeptideEvidences()[[1]]$getProteinAccession() == 'B_id')

  expect_true(ph$getScore() == 1.0)
  expect_true(ph$getRank() == 1)
  expect_true(ph$getSequence()$toString() == "A")

  ph$setScore(2.0)
  expect_true(ph$getScore() == 2.0)
  ph$setRank(30)
  expect_true(ph$getRank() == 30)
  ph$setSequence(AASequence$fromString("AAA"))
  expect_true(ph$getSequence()$toString() == "AAA")

  expect_true(ph == ph)
  expect_false(ph != ph)
})


# test PeptideEvidence ----------------------------------------------------

test_that("test PeptideEvidence",{

  pe = PeptideEvidence$new()
  expect_true(pe == pe)
  expect_false(pe != pe)

  pe$setProteinAccession('B_id')
  expect_true(pe$getProteinAccession() == "B_id")

  pe$setAABefore('A')
  expect_true(pe$getAABefore() == 'A')
  pe$setAAAfter('C')
  expect_true(pe$getAAAfter() == 'C')

  pe$setStart(5)
  expect_true(pe$getStart() == 5)
  pe$setEnd(9)
  expect_true(pe$getEnd() == 9)

  expect_true(pe == pe)
  expect_false( pe != pe)
})


# test PeptideIdentification ----------------------------------------------

test_that("test PeptideIdentification",{

  pi = PeptideIdentification$new()
  check_MetaInfoInterface(pi)
  expect_true(pi == pi)
  expect_false(pi != pi)

  pe = PeptideEvidence$new()
  pe$setProteinAccession('B_id')

  ph = PeptideHit$new(1.0, 1, 0, AASequence$fromString("A"))
  ph$addPeptideEvidence(pe)
  pi$insertHit(ph)
  phx = pi$getHits()
  expect_equal(phx[[1]], ph)

  pi$setHits(c(ph))
  phx = pi$getHits()
  expect_equal(phx[[1]], ph)

  rv = list()
  peptide_hits = pi$getReferencingHits(pi$getHits(), rv)
  expect_equal(rv,list())
  # assert len(peptide_hits) == 1

  expect_true(is_scalar_double(pi$getSignificanceThreshold()))
  is_scalar_character(pi$getScoreType())
  pi$setScoreType("A")
  is.logical(pi$isHigherScoreBetter())
  is_scalar_character(pi$getIdentifier())
  pi$setIdentifier("id")
  pi$assignRanks()
  pi$sort()
  expect_false(pi$empty())

  pi$setSignificanceThreshold(6.0)
})


# test Polarity -----------------------------------------------------------

test_that("testPolarity",{
  expect_true(is_scalar_integer(IonSource$Polarity()$NEGATIVE))
  expect_true(is_scalar_integer(IonSource$Polarity()$POLNULL))
  expect_true(is_scalar_integer(IonSource$Polarity()$POSITIVE))
})


# Precursor ------------------------------------------------------------
test_that("Precursor",{

  pc = Precursor$new()
  pc$setMZ(123.0)
  pc$setIntensity(12.0)
  expect_true(pc$getMZ() == 123.0)
  expect_true(pc$getIntensity() == 12.0)

  pc$setActivationMethods(pc$getActivationMethods())
  pc$setActivationEnergy(6.0)
  pc$getActivationEnergy()

  pc$setIsolationWindowUpperOffset(500.0)
  pc$getIsolationWindowUpperOffset()
  pc$setIsolationWindowLowerOffset(600.0)
  pc$getIsolationWindowLowerOffset()

  pc$setCharge(2)
  pc$getCharge()

  pc$setPossibleChargeStates(pc$getPossibleChargeStates())

  pc$getUnchargedMass()
})


# ProcessingAction --------------------------------------------------------
test_that("test ProcessingAction",{
  pa = ProcessingAction()
  expect_true(is_scalar_integer(pa$ALIGNMENT))
  expect_true(is_scalar_integer(pa$BASELINE_REDUCTION))
})



# test Product ------------------------------------------------------------
test_that("test Product",{
  p = Product$new()
  p$setMZ(12.0)
  p$setIsolationWindowLowerOffset(10.0)
  p$setIsolationWindowUpperOffset(15.0)
  expect_true(p$getMZ() == 12.0)
  expect_true(p$getIsolationWindowLowerOffset() == 10.0)
  expect_true(p$getIsolationWindowUpperOffset() == 15.0)

  expect_true(p == p)
  expect_false(p != p)
})



# test ProteinHit ---------------------------------------------------------
test_that("test ProteinHit",{

  ph = ProteinHit$new()
  ph == ph
  ph != ph
  check_MetaInfoInterface(ph)
  ph$setAccession("A")
  ph$setCoverage(0.5)
  ph$setRank(2)
  ph$setScore(1.5)
  ph$setSequence("ABA")
  expect_true(ph$getAccession() == "A")
  expect_true(ph$getCoverage() == 0.5)
  ph$getRank() == 2
  ph$getScore() == 1.5
  ph$getSequence() == "ABA"
})


# test ProteinIdentification ----------------------------------------------
test_that("ProteinIdentification",{

  pi = ProteinIdentification$new()
  check_MetaInfoInterface(pi)
  expect_true(pi == pi)
  expect_false(pi != pi)

  expect_equal(pi$getHits(),list())
  ph = ProteinHit$new()
  pi$insertHit(ph)
  ph2 = pi$getHits()
  expect_true(ph2[[1]] == ph)

  pi$setHits(c(ph))
  ph2 = pi$getHits()
  expect_equal(ph2[[1]], ph)
})



# test RichPeak -----------------------------------------------------------

test_that("test RichPeak",{

  p2 = RichPeak2D$new()
  check_MetaInfoInterface(p2)
  check_UniqueIdInterface(p2)
  p2 == p2
  p2 != p2
  p2$setMZ(22.0)
  p2$setIntensity(23.0)
  p2$setRT(43.0)
  p2$getMZ() == 22.0
  expect_equal(p2$getIntensity(),23.0)
  p2$getRT() == (43.0)

})


# test Software -----------------------------------------------------------
test_that("test Software",{

  sw = Software$new()
  sw$setName("name")
  sw$setVersion("1.0.0")
  expect_true(sw$getName() == "name")
  expect_true(sw$getVersion() == "1.0.0")

})


# test SourceFile ---------------------------------------------------------

test_that("SourceFile",{

  sf = SourceFile$new()
  sf$setNameOfFile("file.txt")
  expect_true(sf$getNameOfFile() == "file.txt")
  sf$setPathToFile("file.txt")
  expect_true(sf$getPathToFile() == "file.txt")
  sf$setFileType(".txt")
  sf$getFileType() == ".txt"
  sf$setChecksum("abcde000", ChecksumType()$UNKNOWN_CHECKSUM)
  expect_true(sf$getChecksum() == "abcde000")

  expect_true(sf$getChecksumType() %in% c(ChecksumType()$UNKNOWN_CHECKSUM,ChecksumType()$SHA1,ChecksumType()$MD5))
})



# TEST TransformationDescription ------------------------------------------

test_that("test TransformationDescription",{

  td = TransformationDescription$new()
  expect_equal(td$getDataPoints(), list())
  expect_true(is_scalar_double(td$apply(0.0)))

  td$fitModel
  p = td$getModelParameters()
  expect_equal(td$getModelType(), "none")
  td$invert
})


# test TransformationModelInterpolated ------------------------------------

test_that("test TransformationModelInterpolated",{

  for(clz in c(TransformationModelLinear,TransformationModelInterpolated)){

    p = Param$new()
    data = c(TM_DataPoint$new(9.0, 8.9),
             TM_DataPoint$new(5.0, 6.0),
             TM_DataPoint$new(8.0, 8.0))
    mod = clz$new(data, p)
    mod$evaluate(7.0)
    mod$getDefaultParameters(p)
  }
  expect_true(T)
})


# test TransformationXMLFile ----------------------------------------------

test_that("test TransformationXMLFile",{

  fh = TransformationXMLFile$new()
  td = TransformationDescription$new()
  fh$store("test.transformationXML", td)
  fh$load("test.transformationXML", td, T)
  expect_equal(td$getDataPoints(), list())
})



# test IBSpectraFile ------------------------------------------------------
test_that("test IBSpectraFile",{

  fh = IBSpectraFile$new()
  cmap = ConsensusMap$new()
  correctError = F

  tryCatch({
    fh$store( String$new("test.ibspectra.file"), cmap)
  },error = function(e){
    correctError <<- T
  })

  expect_true(correctError)
})


# test Type ---------------------------------------------------------------
test_that("test Type",{
  for(ti in  c(
    FileType()$CONSENSUSXML
    ,FileType()$DTA
    ,FileType()$DTA2D
    ,FileType()$EDTA
    ,FileType()$FASTA
    ,FileType()$FEATUREXML
    ,FileType()$GELML
    ,FileType()$HARDKLOER
    ,FileType()$IDXML
    ,FileType()$INI
    ,FileType()$KROENIK
    ,FileType()$MASCOTXML
    ,FileType()$MGF
    ,FileType()$MS2
    ,FileType()$MSP
    ,FileType()$MZDATA
    ,FileType()$MZIDENTML
    ,FileType()$MZML
    ,FileType()$MZXML
    ,FileType()$OMSSAXML
    ,FileType()$PEPLIST
    ,FileType()$PEPXML
    ,FileType()$PNG
    ,FileType()$PROTXML
    ,FileType()$SIZE_OF_TYPE
    ,FileType()$TOPPAS
    ,FileType()$TRAML
    ,FileType()$TRANSFORMATIONXML
    ,FileType()$TSV
    ,FileType()$UNKNOWN
    ,FileType()$XMASS)){
    expect_true(is_scalar_integer(ti))
    }
})



# test VersionDetails -----------------------------------------------------

test_that("test VersionDetails",{

  is_scalar_character(VersionInfo$getVersion())
  is_scalar_character(VersionInfo$getRevision())
  is_scalar_character(VersionInfo$getTime())

  vd = VersionDetails$create("19.2.1")
  expect_true(vd$version_major == 19)
  expect_true(vd$version_minor == 2)
  expect_true(vd$version_patch == 1)

  vd = VersionDetails$create("19.2.1-alpha")
  expect_true(vd$version_major == 19)
  expect_true(vd$version_minor == 2)
  expect_true(vd$version_patch == 1)
  vd$pre_release_identifier == "alpha"

  vd == vd
  expect_false(vd < vd)
  expect_false(vd > vd)

})


# test InspectInfile ------------------------------------------------------
test_that("test InspectInfile",{

  inst = InspectInfile$new()

  expect_true(!is.null(inst$getModifications))
  mods = inst$getModifications()
  # P.S. mods is collections::dict(). It is an environment.
  mods$print()
  expect_equal(mods$size(), 0)
})


# test IsotopeMarker ------------------------------------------------------

test_that("test IsotopeMarker",{

  inst = IsotopeMarker$new()
  ptr = inst$create()

  expect_false(is.null(ptr$apply))

  # This should be a collections::dict
  res = collections::dict()
  spec = MSSpectrum$new()
  ptr$apply(res, spec)

})


# test Attachment ---------------------------------------------------------
test_that("test Attachment",{
  inst = Attachment$new()

  inst$name
  inst$value
  inst$cvRef
  inst$cvAcc
  inst$unitRef
  inst$unitAcc
  inst$binary
  inst$qualityRef
  inst$colTypes
  inst$tableRows

  inst$toXMLString
  inst$toCSVString

  inst$name = "test"
  inst$value = "test"
  inst$cvRef = "test"
  inst$cvAcc = "test"
  inst$unitRef = "test"
  inst$unitAcc = "test"
  inst$binary = "test"
  inst$qualityRef = "test"
  inst$colTypes = c("test","test2")
  inst$tableRows = list(c("test", "test2"), c("otherTest") )
  expect_equal(inst$tableRows, list(c("test","test2"), c("otherTest")))

  expect_true(inst$tableRows[[2]][[1]] == "otherTest")
})


# test OptimizePeakDeconvolution ------------------------------------------

test_that("test OptimizePeakDeconvolution",{

  inst = OptimizePeakDeconvolution$new()
  inst$getParameters

  inst$getPenalties
  inst$setPenalties
  inst$getCharge
  inst$setCharge
  inst$optimize


  inst = PenaltyFactorsIntensity$new()
  expect_true(!is.null(inst$height))

  inst = OptimizePeakDeconvolution_Data$new()
  inst$peaks
  inst$peaks
  inst$signal
  inst$penalties
  expect_equal(inst$charge, 0)

})



# test testKernelMassTrace ------------------------------------------------
test_that("test KernelMassTrace",{

  trace = Kernel_MassTrace$new()

  expect_true(!is.null(trace$getSize))
  expect_true(!is.null(trace$getLabel))
  expect_true(!is.null(trace$setLabel))

  trace$getCentroidMZ
  trace$getCentroidRT
  trace$getCentroidSD
 trace$getFWHM
  trace$getTraceLength
  trace$getFWHMborders
  trace$getSmoothedIntensities
  trace$getAverageMS1CycleTime

  trace$computeSmoothedPeakArea
  trace$computePeakArea
  trace$findMaxByIntPeak
  trace$estimateFWHM
  trace$computeFwhmArea
  trace$computeFwhmAreaSmooth
  # assert trace.computeFwhmAreaRobust
  # assert trace.computeFwhmAreaSmoothRobust
  trace$getIntensity
  trace$getMaxIntensity

  trace$getConvexhull

  trace$setCentroidSD
  trace$setSmoothedIntensities
  trace$updateSmoothedMaxRT
  trace$updateWeightedMeanRT
  trace$updateSmoothedWeightedMeanRT
  trace$updateMedianRT
  trace$updateMedianMZ
  trace$updateMeanMZ
  trace$updateWeightedMeanMZ
  trace$updateWeightedMZsd

  s = trace$getSize()
  expect_equal(s, 0)

})



# test testElutionPeakDetection -------------------------------------------

test_that("test ElutionPeakDetection",{

  detection = ElutionPeakDetection$new()

  expect_true(!is.null(detection$detectPeaks))
  expect_true(!is.null(detection$filterByPeakWidth))
  expect_true(!is.null(detection$computeMassTraceNoise))
  expect_true(!is.null(detection$computeMassTraceSNR))
  expect_true(!is.null(detection$computeApexSNR))
  expect_true(!is.null(detection$findLocalExtrema))
  expect_true(!is.null(detection$smoothData))

  trace = Kernel_MassTrace$new()
  expect_invisible(detection$smoothData(trace, 4))

})


# test IndexedMzMLDecoder -------------------------------------------------

test_that("test IndexedMzMLDecoder",{

  decoder = IndexedMzMLDecoder$new()

  tryCatch({
    pos = decoder$findIndexListOffset("abcde", 100)
  }, error = function(e){invisible()
  })

  expect_true(T)
})


test_that("test_streampos",{

  p = streampos$new()
  expect_true(is_scalar_integer(as.integer(p)))
})


test_that("test_MapConversion",{

  feature = Feature$new()
  feature$setRT(99)

  cmap = ConsensusMap$new()
  fmap = FeatureMap$new()
  fmap$push_back(feature)
  MapConversion$new()$convert(0, fmap, cmap, 1)

  expect_true(cmap$size() == 1)
  expect_true(cmap[1]$getRT() == 99.0)

  fmap = FeatureMap$new()
  MapConversion$new()$convert(cmap, T, fmap)

  expect_true(fmap$size() == 1)
  expect_true(fmap[1]$getRT() == 99.0)

  exp = MSExperiment$new()
  sp = MSSpectrum$new()
  peak = Peak1D$new()
  peak$setIntensity(10)
  peak$setMZ(20)
  sp$push_back(peak)
  exp$addSpectrum(sp)
  exp$addSpectrum(sp)

  cmap = ConsensusMap$new()
  MapConversion$new()$convert(0, exp, cmap, 2)

  expect_true(cmap$size() == 2)
  expect_true(cmap[1]$getIntensity() == 10.0)
  expect_true(cmap[1]$getMZ() == 20.0)

})

test_that("test_BSpline2d",{

  x = c(1.0, 6.0, 8.0, 10.0, 15.0)
  y = c(2.0, 5.0, 6.0, 12.0, 13.0)
  spline = BSpline2d$new(x,y,0, BoundaryCondition()$BC_ZERO_ENDPOINTS, 0)

  expect_true(spline$ok())
  expect_equal(abs(spline$eval(6.0) - 5.0 < 0.01), 1)
  expect_equal(abs(spline$derivative(6.0) - 5.0 < 0.01), 1)

  y_new = c(4.0, 5.0, 6.0, 12.0, 13.0)
  spline$solve(y_new)

  spline$ok()
  abs(spline$eval(6.0) - 5.0 < 0.01)

})

test_that("testConsensusIDAlgorithmAverage",{
  algo = ConsensusIDAlgorithmAverage$new()
  expect_true(!is.null(algo$apply))
})

test_that("testConsensusIDAlgorithmBest",{
  algo = ConsensusIDAlgorithmBest$new()
  expect_true(!is.null(algo$apply))
})

test_that("testConsensusIDAlgorithmIdentity",{
  algo = ConsensusIDAlgorithmIdentity$new()
  expect_true(!is.null(algo$apply))
})

test_that("testConsensusIDAlgorithmPEPIons",{
  algo = ConsensusIDAlgorithmPEPIons$new()
  expect_true(!is.null(algo$apply))
})

test_that("testConsensusIDAlgorithmPEPMatrix",{
  algo = ConsensusIDAlgorithmPEPMatrix$new()
  expect_true(!is.null(algo$apply))
})

test_that("testConsensusIDAlgorithmRanks",{
  algo = ConsensusIDAlgorithmRanks$new()
  expect_true(!is.null(algo$apply))
})

test_that("testConsensusIDAlgorithmSimilarity",{
  algo = ConsensusIDAlgorithmSimilarity$new()
  expect_true(!is.null(algo$apply))
})

test_that("testConsensusIDAlgorithmWorst",{
  algo = ConsensusIDAlgorithmWorst$new()
  expect_true(!is.null(algo$apply))
})

test_that("testDigestionEnzymeProtein",{
  f = EmpiricalFormula$new()

  regex_description = ""
  psi_id = ""
  xtandem_id = ""
  comet_id = 0
  omssa_id = 0
  e = DigestionEnzymeProtein$new("testEnzyme", "K", list(), regex_description,
                                      f, f, psi_id, xtandem_id, comet_id, omssa_id)
  expect_true(T)
})

test_that("testMRMAssay",{
  e = MRMAssay$new()
  expect_true(!is.null(e))
})

test_that("testMRMIonSeries",{
  e = MRMIonSeries$new()
  expect_true(!is.null(e))
})

test_that("testPeptideIndexing",{
  e = PeptideIndexing$new()
  expect_true(!is.null(e))
})

test_that("testPeptideProteinResolution",{
  e = PeptideProteinResolution$new(F)
  expect_true(!is.null(e))
})

test_that("testPercolatorOutfile",{
  e = PercolatorOutfile$new()
  expect_true(!is.null(e))
})

test_that("testHiddenMarkovModel",{

  hmm = HiddenMarkovModel$new()
  expect_true(!is.null(hmm))

  expect_equal(hmm$getNumberOfStates(), 0)

  ss = String$new("testState")
  hmm$addNewState(ss)

  expect_true(hmm$getNumberOfStates() == 1)

  e = HMMState$new()
  # hmm$addNewState(e) # Segfault !

  r = hmm$getState(String$new("testState"))
  expect_true(!is.null(r))

})


# test HMMState -----------------------------------------------------------

test_that("HMMState",{
  e = HMMState$new()
  expect_true(!is.null(e))
  e$setName(String$new("somename"))
  expect_true(e$getName() == "somename")
  e$setHidden(T)
  expect_true(e$isHidden())

  pre = HMMState$new()
  pre$setName(String$new("pre"))
  suc = HMMState$new()
  suc$setName(String$new("suc"))

  e$addPredecessorState(pre)
  e$addSuccessorState(suc)

  expect_true(!is.null(e$getPredecessorStates()))
  expect_true(!is.null(e$getSuccessorStates()))
})


# test ProteaseDB ---------------------------------------------------------

test_that("test ProteaseDB",{

  edb = ProteaseDB$new()

  f = EmpiricalFormula$new()
  synonyms = c("dummy", "other")

  expect_true(edb$hasEnzyme(String$new("Trypsin")))

  trypsin = edb$getEnzyme(String$new("Trypsin"))

  names = vector()
  edb$getAllNames(names)
  expect_true("Trypsin" %in% names)

})


# test ElementDB ----------------------------------------------------------

test_that("test ElementDB",{

  edb = ElementDB$new()
  rm(edb)

  # create a second instance of ElementDB without anything bad happening
  edb = ElementDB$new()

  expect_true(edb$hasElement(16))
  edb$hasElement(String$new("O"))

  e = edb$getElement(16)

  expect_true(e$getName() == "Sulfur")
  expect_true( e$getSymbol() == "S")
  expect_false(is.null(e$getIsotopeDistribution()))

  e2 = edb$getElement(String$new("O"))

  expect_true(e2$getName() == "Oxygen")
  expect_true(e2$getSymbol() == "O")
  expect_false(is.null(e$getIsotopeDistribution()))

  # assert e == e2

  #  not yet implemented
  #
  # const Map[ String, Element * ]  getNames() nogil except +
  # const Map[ String, Element * ] getSymbols() nogil except +
  # const Map[unsigned int, Element * ] getAtomicNumbers() nogil except +
})


# test DPosition ----------------------------------------------------------
test_that("Dposition",{

  dp = DPosition1$new()
  dp = DPosition1$new(1.0)
  expect_true( dp[1] == 1.0)

  dp = DPosition2$new()
  dp = DPosition2$new(1.0, 2.0)

  expect_true(dp[1] == 1.0)
  expect_true(dp[2] == 2.0)

})


# test ResidueDB ----------------------------------------------------------

test_that("ResidueDB",{

  rdb = ResidueDB$new()
  rm(rdb)

  # create a second instance of ResidueDB without anything bad happening
  rdb = ResidueDB$new()

  expect_true(rdb$getNumberOfResidues() >= 20)
  expect_true(length(rdb$getResidueSets() ) >= 1)
  el = rdb$getResidues(String$new(rdb$getResidueSets()[[1]]))

  expect_true(length(el) >= 1)

  expect_true(rdb$hasResidue(String$new("Glycine")))
  expect_true(!is.null(rdb$getResidue(String$new("Glycine"))))

  nrr = rdb$getNumberOfResidues()

})


# test ModificationsDB ----------------------------------------------------

test_that("test ModificationsDB",{

  mdb = ModificationsDB$new()
  rm(mdb)

  # create a second instance of ModificationsDB without anything bad happening
  mdb = ModificationsDB$new()

  expect_true(mdb$getNumberOfModifications() > 1)
  m = mdb$getModification(1)

  expect_true(mdb$getNumberOfModifications() > 1)
  m = mdb$getModification(1)
  expect_true(!is.null(m))

  mods = list()
  mdb$searchModifications(mods, String$new("Phosphorylation"), String$new("T"), ResidueModification$TermSpecificity()$ANYWHERE)
  expect_true(length(mods) == 1)

  mods = list()
  mdb$searchModifications(mods, String$new("NIC"), String$new("T"), ResidueModification$TermSpecificity()$N_TERM)
  expect_true(length(mods) == 1)

  mods = list()
  mdb$searchModifications(mods, String$new("NIC"), String$new("T"), ResidueModification$TermSpecificity()$N_TERM)
  expect_true(length(mods) == 1)

  mods = list()
  mdb$searchModifications(mods, String$new("Acetyl"), String$new("T"), ResidueModification$TermSpecificity()$N_TERM)
  expect_true(length(mods) == 1)
  expect_true(mods[[1]]$getFullId() == "Acetyl (N-term)")

  m = mdb$getModification(String$new("Carboxymethyl (C)"), "", ResidueModification$TermSpecificity()$NUMBER_OF_TERM_SPECIFICITY)
  expect_true(m$getFullId() == "Carboxymethyl (C)")

  m = mdb$getModification( String$new("Phosphorylation"), String$new("S"), ResidueModification$TermSpecificity()$ANYWHERE)
  expect_true(m$getId() == "Phospho")

  # get out all mods (there should be many, some known ones as well!)
  mods = list()
  m = mdb$getAllSearchModifications(mods)
  expect_true(length(mods) > 100)

  expect_true("Phospho (S)" %in% mods)
  expect_true("Sulfo (S)" %in% mods)
  expect_false("Phospho" %in% mods)

  # search for specific modifications by mass
  m = mdb$getBestModificationByDiffMonoMass( 80.0, 1.0, "T",  ResidueModification$TermSpecificity()$ANYWHERE)
  expect_true(!is.null(m))
  expect_true(m$getId() == "Phospho")
  expect_true(m$getFullName() == "Phosphorylation")
  expect_true(m$getUniModAccession() == "UniMod:21")

  m = mdb$getBestModificationByDiffMonoMass(80, 100, "T",  ResidueModification$TermSpecificity()$ANYWHERE)
  expect_true(!is.null(m))
  expect_true(m$getId() == "Phospho")
  expect_true(m$getFullName() == "Phosphorylation")
  expect_true(m$getUniModAccession() == "UniMod:21")

  m = mdb$getBestModificationByDiffMonoMass(16, 1.0, "M", ResidueModification$TermSpecificity()$ANYWHERE)
  expect_true(!is.null(m))
  expect_true(m$getId() == "Oxidation")
  expect_true(m$getFullName() == "Oxidation or Hydroxylation")
  expect_true(m$getUniModAccession() == "UniMod:35")

})



# test RNaseDb ------------------------------------------------------------

test_that("test RNaseDB",{

  db = RNaseDB$new()
  names = list()
  db$getAllNames(names)

  e = db$getEnzyme("RNase_T1")
  expect_true(e$getRegEx() == '(?<=G)')
  expect_true( e$getThreePrimeGain() == 'p')

  expect_true( db$hasRegEx('(?<=G)'))
  expect_true(db$hasEnzyme("RNase_T1"))

})


# test RibonucleotideDB ---------------------------------------------------

test_that("test RibonucleotideDB",{

  r = RibonucleotideDB$new()

  uridine = r$getRibonucleotide("U")

  expect_true( uridine$getName() == 'uridine')
  expect_true( uridine$getCode() == 'U')
  expect_true(uridine$getFormula()$toString() == 'C9H12N2O6')
  expect_true( uridine$isModified() == F)
})


# test Ribonucleotide -----------------------------------------------------
test_that("test Ribonucleotide",{
  r = Ribonucleotide$new()

  expect_false(r$isModified())

  r$setHTMLCode("test")
  expect_true(r$getHTMLCode() == "test")

  r$setOrigin("A")
  expect_true(r$getOrigin() == "A")

  r$setNewCode("A")
  expect_true( r$getNewCode() == "A")

})


# test RNaseDigestion -----------------------------------------------------
test_that("test RNaseDigestion",{

  dig = RNaseDigestion$new()
  dig$setEnzyme("RNase_T1")
  expect_true(dig$getEnzymeName() == "RNase_T1")

  oligo = NASequence$fromString("pAUGUCGCAG")

  result = list()
  dig$digest(oligo, result)
  # Pass by reference like effect currently only works for corresponding non-overloaded
  # methods
  expect_true(length(result) == 0)
  result = list()
  dig$digest_0(oligo, result)
  expect_true(length(result) == 3)
})


# test NASequence ---------------------------------------------------------
test_that("test NASequence",{

  oligo = NASequence$fromString("pAUGUCGCAG");

  expect_true(oligo$size() == 9)
  seq_formula = oligo$getFormula()
  expect_equal(seq_formula$toString(),'C86H108N35O64P9')

  oligo_mod = NASequence$fromString("A[m1A][Gm]A")
  seq_formula = oligo_mod$getFormula()
  seq_formula$toString() == 'C42H53N20O23P3'


  expect_true(oligo_mod[2]$isModified())

  charge = 2
  oligo_mod$getMonoWeight(NASequence$NASFragmentType()$WIon, charge)
  oligo_mod$getFormula(NASequence$NASFragmentType()$WIon, charge)
})


# test ExperimentalDesign -------------------------------------------------

test_that("test ExperimentalDesign",{

  f = ExperimentalDesignFile$new()
  fourplex_fractionated_design = ExperimentalDesign$new()
  ed_filename = paste(normalizePath("..", winslash = "/",mustWork = T),fsep = .Platform$file.sep,"ExperimentalDesign_input_2.tsv", sep = "")
  fourplex_fractionated_design = ExperimentalDesignFile$load(ed_filename, F)
  expect_true(fourplex_fractionated_design$getNumberOfSamples() == 8)
  expect_true(fourplex_fractionated_design$getNumberOfFractions() == 3)
  expect_true(fourplex_fractionated_design$getNumberOfLabels() == 4)
  expect_true(fourplex_fractionated_design$getNumberOfMSFiles() == 6)
  expect_true(fourplex_fractionated_design$getNumberOfFractionGroups() == 2)
  expect_true(fourplex_fractionated_design$getSample(1, 1) == 1)
  expect_true(fourplex_fractionated_design$getSample(2, 4) == 8)
  expect_true(fourplex_fractionated_design$isFractionated())
  expect_true(fourplex_fractionated_design$sameNrOfMSFilesPerFraction())

})


# test String -------------------------------------------------------------

test_that("test String",{

  rstr = String$new()
  rstr = String$new("blah")
  expect_true(rstr$toString() == "blah")
  rstr = String$new("blah")
  expect_true(rstr$toString() == "blah")
  rstr = String$new("blah")
  expect_true(rstr$toString() == "blah")
  rstr = String$new(rstr)
  expect_true(rstr$toString() == "blah")
  expect_true(nchar(rstr$toString()) == 4)

  print(rstr)
  print(rstr$toString())
  expect_true(rstr$toString() == "blah")

  rstr = String$new("bläh")
  expect_true(rstr$toString() == "bläh")
  rstr = String$new("bläh")
  expect_true(rstr$toString() == "bläh")
  rstr = String$new(rstr)
  expect_true(rstr$toString() == "bläh")

  print(rstr)
  print(enc2utf8(rstr$toString()))

  expect_true(nchar(rstr$toString()) == 4)
  print(rstr)
  print(enc2utf8(rstr$toString())) # this prints the correct String

  pystr1 = String$new("bläh")
  pystr2 = String$new("bläh")
  expect_true(pystr1 == pystr2)

  pystr1 = String$new("bläh")
  pystr2 = String$new("bläh")
  expect_true(pystr1 == pystr2)

  # Handling of different Unicode Strings:
  # - unicode is natively stored in OpenMS::String
  # - encoded bytesequences for utf8, utf16 and iso8859 can be stored as
  #   char arrays in OpenMS::String (and be accessed using c_str())
  # - encoded bytesequences for anything other than utf8 cannot use
  #   "toString()" as this function expects utf8
  ustr = "bläh"
  pystr = String$new(ustr)
  expect_true(pystr$toString() == "bläh")
  pystr = String$new(enc2utf8(ustr))
  expect_true(pystr$toString() == c("bläh"))


  Encoding(ustr) <- "iso8859-15"

  pystr = String$new(ustr)
  pystr$toString()
  Encoding(ustr) <- "UTF-16"

  pystr = String$new(ustr)
  pystr$toString()

  # should not throw error
  ustr <- "\u{68}\u{65}\u{6c}\u{6c}"
  Encoding(ustr) <- "UTF-16"
  pystr = String$new(ustr)
  didThrow = F

  tryCatch({
    expect_equal(pystr$toString(),"hell")
  }, error = function(e) {
    didThrow <<- T
  })

  expect_false(didThrow)

  # Handling of automatic conversions of String return values
  ustr = "bläh"
  s = MSSpectrum$new()
  s$setNativeID(ustr)
  r = s$getNativeID()
  expect_equal(r, "bläh")

  # IMP: Reticulate handles conversion of unicode strings
  #      from r to python and vice versa.

  x <- "Maur\xC3\xADcio"
  Encoding(x) <- "UTF-8"
  print(x)
  s$setNativeID(x)
  r = s$getNativeID()
  expect_true(r == "Maurício")

  x <- "Maur\xEDcio"
  Encoding(x) <- "latin1"
  print(x)
  s$setNativeID(x)
  r = s$getNativeID()
  expect_true(r == "Maurício")


  s$setNativeID(ustr)
  r = s$getNativeID()
  expect_true(r == "bläh")

  Encoding(ustr) <- "iso8859_15"
  s$setNativeID(ustr)
  r = s$getNativeID()
  print(r)

  t <- "bl\xfch"
  s$setNativeID(t)
  r = s$getNativeID()
  expect_true(r == "blüh")

  s$setNativeID(enc2utf8("abc"))
  r = s$getNativeID()
  expect_true(r == "abc")

  # IMP: something like this will crash the R session (due to segfault)
  # t <- "blüh"
  # Encoding(t) <- "UTF-8"
  #
})

