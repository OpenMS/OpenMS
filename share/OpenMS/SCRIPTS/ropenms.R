################### OpenMS in R ###################

#### Some simple scripts how to use OpenMS in R ####

### if not installed:
### - install pyopenms (https://pyopenms.readthedocs.io/en/latest/installation.html)
###   make sure R is using the same python environment as your pyopenms installation
###   eg. reticulate::use_python("/usr/local/miniconda3/envs/py37/bin/python")
###   or  before loading the reticulate library
###       Sys.setenv(RETICULATE_PYTHON = "/usr/local/miniconda3/envs/py37/bin/python")
### - install reticulate (https://rstudio.github.io/reticulate/)
###   eg. install.packages('reticulate')


library("reticulate")
ropenms=import("pyopenms", convert = FALSE)


### load and parse idXML

f="/OpenMS/OpenMS/share/OpenMS/examples/BSA/BSA1_OMSSA.idXML"

idXML=ropenms$IdXMLFile()

pepids=r_to_py(list())
protids=r_to_py(list())

idXML$load(f,protids,pepids)

pepids=py_to_r(pepids)
protids=py_to_r(protids)

pephits=pepids[[1]]$getHits()
pepseq=pephits[[1]]$getSequence()


### load and parse featureXML

f="/OpenMS/OpenMS/share/OpenMS/examples/FRACTIONS/BSA1_F1.featureXML"

featXML=ropenms$FeatureXMLFile()
fmap = ropenms$FeatureMap()

featXML$load(f, fmap)

print(paste0("FeatureID: ", fmap[1]$getUniqueId()))
print(paste0("Charge: ", fmap[1]$getCharge()))
print(paste0("M/z: ", fmap[1]$getMZ()))
print(paste0("RT: ", fmap[1]$getRT()))

### load and parse mzML

f="/OpenMS/OpenMS/share/OpenMS/examples/BSA/BSA1.mzML"

mzML= ropenms$MzMLFile()
msexp = ropenms$MSExperiment()

mzML$load(f,msexp)
spectra = py_to_r(msexp$getSpectra())

#ms1
ms1=sapply(spectra, function(x) x$getMSLevel()==1)
peaks=sapply(spectra[ms1], function(x) cbind(do.call("cbind", x$get_peaks()),x$getRT()))
peaks=do.call("rbind", peaks)
peaks_df=data.frame(peaks)
colnames(peaks_df)=c('MZ','Intensity','RT')
peaks_df$Intensity=log10(peaks_df$Intensity)
ggplot(peaks_df, aes(x=RT, y=MZ) )+geom_point(size=1, aes(colour = Intensity), alpha=0.25) + theme_minimal() + scale_colour_gradient(low = "blue", high = "yellow")

#ms2
ms2=spectra[!ms1][[1]]$get_peaks()
peaks_ms2=do.call("cbind", ms2)
peaks_ms2=data.frame(peaks_ms2)
colnames(peaks_ms2)=c("MZ","Intensity")
ggplot(peaks_ms2, aes(x=MZ, y=Intensity)) +
geom_segment( aes(x=MZ, xend=MZ, y=0, yend=Intensity)) +
geom_segment( aes(x=MZ, xend=MZ, y=0, yend=-Intensity)) + # mirror spectrum possibly useful for synthetic peptide spectra comparison 
theme_minimal()


### Spectrum

spectrum = ropenms$MSSpectrum()
mz = seq(1500, 500, -100)
i = seq(10, 2000, length.out = length(mz))
spectrum$set_peaks(list(mz, i))

# Sort the peaks according to ascending mass-to-charge ratio
spectrum$sortByPosition()

# Iterate using the reticulate::iterate() function
iterate(spectrum, function(x) {print(paste0("M/z :" , x$getMZ(), " Intensity: ", x$getIntensity()))})


# Iterate over spectrum of those peaks
for (i in seq(0,py_to_r(spectrum$size())-1)) {
  print(spectrum[i]$getMZ())
  print(spectrum[i]$getIntensity())
}

# More efficient peak access with get_peaks()
peak_df=do.call("cbind", py_to_r(spectrum$get_peaks()))
apply(peak_df,1,c)

# Access a peak by index
print(c(spectrum[1]$getMZ(), spectrum[1]$getIntensity()))
