import pyopenms

def main():
    import optparse
    # in: file
    in_formats = "mzData,mzXML,mzML,DTA,DTA2D,mgf,featureXML,consensusXML,"\
                 "ms2,fid,tsv,peplist,kroenik,edta".split()
    # out: file
    out_formats = "mzData,mzXML,mzML,DTA2D,mgf,featureXML,"
                  "consensusXML".split(",")

if __name__ == "__main__":
    main()
