class SimpleOpenMSSpectraFactory:

    @staticmethod
    def getSpectrumAccessOpenMSPtr(exp):
      is_cached = False

      for i in range(exp.size()):
        for dp in exp[i].getDataProcessing():
          if dp.metaValueExists("cached_data"):
            is_cached = True

      for chrom in exp.getChromatograms():
        for dp in chrom.getDataProcessing():
          if dp.metaValueExists("cached_data"):
            is_cached = True

      if is_cached:
        return SpectrumAccessOpenMSCached( exp.getLoadedFilePath() )
      else:
        return SpectrumAccessOpenMS( exp )



def read_mzML(filepath):
    """
    Convenience function to read mzML files and return an MSExperiment.
    
    This is a simple wrapper around MzMLFile().load() that creates and returns
    an MSExperiment object in a single call. For more advanced options like
    filtering by MS level, RT range, or m/z range, use MzMLFile directly
    with PeakFileOptions.
    
    Parameters
    ----------
    filepath : str
        Path to the mzML file to load
    
    Returns
    -------
    MSExperiment
        The loaded experiment containing spectra and chromatograms
        
    Examples
    --------
    >>> import pyopenms as oms
    >>> exp = oms.read_mzML("data/example.mzML")
    >>> print(f"Loaded {len(exp.getSpectra())} spectra")
    
    See Also
    --------
    MzMLFile : For advanced loading options with filtering
    MSExperiment : The returned data structure
    """
    exp = MSExperiment()
    MzMLFile().load(filepath, exp)
    return exp
