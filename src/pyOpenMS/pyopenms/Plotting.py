
def plotChromatogram(c):
    import pylab
    x, y = c.get_peaks()
    pylab.plot( x, y )
    pylab.xlabel("Retention time")
    pylab.ylabel("Intensity")
    pylab.show()

def plotSpectrum(s):
    import pylab
    x, y = s.get_peaks()
    pylab.bar( x, y )
    pylab.xlabel("m/z")
    pylab.ylabel("Intensity")
    pylab.show()
