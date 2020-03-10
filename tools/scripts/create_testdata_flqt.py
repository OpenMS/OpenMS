
import pyopenms
"""
Producing the test data for TOPP_FeatureLinkerUnlabeledQT_5 and TOPP_FeatureLinkerUnlabeledQT_6
"""

fmaps = [ pyopenms.FeatureMap() for i in range(3)]
pepids = []
pepseq = ["PEPTIDEA", "PEPTIDEK", "PEPTIDER"]
for s in pepseq:
    pepid = pyopenms.PeptideIdentification()
    hit = pyopenms.PeptideHit()
    hit.setSequence(pyopenms.AASequence.fromString(s, True))
    pepid.insertHit(hit)
    pepid.setIdentifier("Protein0")
    pepids.append(pepid)

protid = pyopenms.ProteinIdentification()
protid.setIdentifier("Protein0")
for i,fmap in enumerate(fmaps):
    fmap.setProteinIdentifications( [protid])
    # add 3 features to each map, but with a twist (adding different peptide ids to different maps)
    for k in range(3):
        f = pyopenms.Feature()
        f.setRT(300 + k*100 + i*10)
        f.setMZ(500 + k*0.001 + i*0.01)
        f.setIntensity(500 + i*100)
        f.setMetaValue("sequence", pepseq[ (i+k) % 3]) # easier viewing in TOPPView
        f.setPeptideIdentifications( [pepids[(i+k) % 3]] )
        fmap.push_back(f)
    pyopenms.FeatureXMLFile().store("output_%s.featureXML" % i, fmap)


