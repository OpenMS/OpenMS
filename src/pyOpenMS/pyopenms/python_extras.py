from . import Constants
from . import Plotting as Plot

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


