

    def estimateType(self, MSSpectrum spec):
        return self.inst.get().estimateType(spec.inst.get().begin(), spec.inst.get().end())
