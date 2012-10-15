# --------------------------------------------------------------------------
#                   OpenMS -- Open-Source Mass Spectrometry
# --------------------------------------------------------------------------
# Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
# ETH Zurich, and Freie Universitaet Berlin 2002-2012.
#
# This software is released under a three-clause BSD license:
#  * Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#  * Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
#  * Neither the name of any author or any participating institution
#    may be used to endorse or promote products derived from this software
#    without specific prior written permission.
# For a full list of authors, refer to the file AUTHORS.
# --------------------------------------------------------------------------
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
# INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
# OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
# OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# --------------------------------------------------------------------------
# $Maintainer: Mathias Walzer $
# $Authors: Mathias Walzer $
# --------------------------------------------------------------------------

#~ from lxml import etree
import mzQCML
import os
import base64
import sys
from subprocess import Popen, PIPE, STDOUT
import codecs

currentpath = os.getcwd()

print sys.argv

#~ script needs to be called by #python script
#~ datapath = sys.argv[1]
mzmlfile = sys.argv[2]
idxmlfile = sys.argv[3]
precursortsv = sys.argv[4]
tictsv = sys.argv[5]
idtsv = sys.argv[6]
precursorpng = sys.argv[7]
ticpng = sys.argv[8]
idpng = sys.argv[9]
#~ print mzmlfile,idxmlfile,precursortsv,tictsv,idtsv
outfilepath = sys.argv[1] + '/'
#~ mzmlfile = '20100219_SvNa_SA_Ecoli_PP.mzML'
#~ idxmlfile = '20100219_SvNa_SA_Ecoli_OMSSA_pep_ind_fdr_001.idXML'
qcfilename = outfilepath + 'genericwrapper' + '.qcML'

root = mzQCML.MzQualityMLType()
#root.SetQuality = list()
#SetElement = mzQCML.SetQualityAssessmentType()
root.RunQuality = list()
RunElement = mzQCML.RunQualityAssessmentType()

			#~ <cloptions> "%1"/QCExporter.py "%2" "%3" "%4" "%5" "%6" "%7" "%8" "%9" "%10"</cloptions>
			#~ <path>python</path>
			#~ <mappings>
				#~ <mapping id="1" cl="%%scriptpath" />
				#~ <mapping id="2" cl="%TMP" />
				#~ <mapping id="3" cl="%%inML" />
				#~ <mapping id="4" cl="%%inID" />
				#~ <mapping id="5" cl="%%inTSVspec" />
				#~ <mapping id="6" cl="%%inTSVtic" />
				#~ <mapping id="7" cl="%%inTSVacc" />
				#~ <mapping id="8" cl="%%inFIGspec" />
				#~ <mapping id="9" cl="%%inFIGtic" />
				#~ <mapping id="10" cl="%%inFIGacc" />
				#~ <file_post location="%TMP/genericwrapper.qcML" target="out" />
			#~ </mappings>

if tictsv != "":
	#~ rscript = 'ProduceQCFigures_tic.R ' + tictsv + ' ' + outfilepath
	#~ Rcmd = 'Rscript --vanilla %s'
	#~ print Rcmd % rscript
	#~ p = Popen( (Rcmd % rscript), shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).wait()
	#~ file = open(outfilepath+"tic.png", "rb")

	file = open(ticpng , "rb")
	data = file.read()
	file.close()
	byte_arr = base64.b64encode(data)
	atta = mzQCML.AttachmentType()
	atta.set_binary(byte_arr)
	atta.set_name("tic")
	atta.set_cvRef("QCML:1000007") #tic
	#~ atta.set_UnitCvRef("QCML:1000006") #png
	atta.set_accession(os.urandom(9).encode('hex'))
	RunElement.add_Attachment(atta)

if precursortsv != "" and idtsv != "":
	#~ rscript = 'ProduceQCFigures_spec.R ' + precursortsv + ' '+ idtsv + ' ' + outfilepath
	#~ Rcmd = 'Rscript --vanilla %s'
	#~ p = Popen( (Rcmd % rscript), shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).wait()
	#~ file = open(outfilepath+"spec.png", "rb")

	file = open(precursorpng, "rb")
	data = file.read()
	file.close()
	byte_arr = base64.b64encode(data)
	atta = mzQCML.AttachmentType()
	atta.set_binary(byte_arr)
	atta.set_name("spec")
	atta.set_cvRef("QCML:1000008") #precursors
	#~ atta.set_UnitCvRef("QCML:1000006") #png
	atta.set_accession(os.urandom(9).encode('hex'))
	RunElement.add_Attachment(atta)

if precursortsv != "" and idtsv != "":
	#~ rscript = 'ProduceQCFigures_acc.R ' + idtsv + ' ' + outfilepath
	#~ Rcmd = 'Rscript --vanilla %s'
	#~ p = Popen( (Rcmd % rscript), shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True).wait()
	#~ file = open(outfilepath+"accuracy.png", "rb")

	file = open(idpng, "rb")
	data = file.read()
	file.close()
	byte_arr = base64.b64encode(data)
	atta = mzQCML.AttachmentType()
	atta.set_binary(byte_arr)
	atta.set_name("accuracy")
	atta.set_cvRef("QCML:1000009") # accuracy
	#~ atta.set_UnitCvRef("QCML:1000006") #png
	atta.set_accession(os.urandom(9).encode('hex'))
	RunElement.add_Attachment(atta)

kbsize = os.path.getsize(mzmlfile)/1024
qp = mzQCML.CVParamType()
qp.set_value(kbsize)
qp.set_name("size")
qp.set_accession(os.urandom(9).encode('hex'))
qp.set_cvRef('QCML:1000010') # rawfilesize
qp.set_unitCvRef('UO:0000234') #kb
RunElement.add_QualityParameter(qp)


FIcmd = 'FileInfo -in %s'
p = Popen( (FIcmd % mzmlfile), shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
output = p.stdout.read()
output = output.split('\n')

for p in output:
	p = p.strip()
	if p:
		if p.startswith('retention time: '):
			qp = mzQCML.CVParamType()
			qp.set_value(p.split(' ')[-3])
			qp.set_name("rtstart")
			qp.set_accession(os.urandom(9).encode('hex'))
			qp.set_cvRef('QCML:1000001') #interval start
			qp.set_unitCvRef('MS:1000894')
			RunElement.add_QualityParameter(qp)
			qp = mzQCML.CVParamType()
			qp.set_value(p.split(' ')[-1])
			qp.set_name("rtend")
			qp.set_accession(os.urandom(9).encode('hex'))
			qp.set_cvRef('QCML:1000002') #interval end
			qp.set_unitCvRef('MS:1000894')
			RunElement.add_QualityParameter(qp)
			#<QualityParameter accession="String" name="rtstart" value="0.18" cvRef="rtstart"/>
			#<QualityParameter accession="String" name="rtend" value="14399.93" cvRef="rtend"/>

		elif p.startswith('retention time: '):
			qp = mzQCML.CVParamType()
			qp.set_value(p.split(' ')[-3])
			qp.set_name("mzstart")
			qp.set_accession(os.urandom(9).encode('hex'))
			qp.set_unitCvRef('MS:1000040')
			qp.set_cvRef('QCML:1000001') #interval start
			RunElement.add_QualityParameter(qp)
			qp = mzQCML.CVParamType()
			qp.set_value(p.split(' ')[-1])
			qp.set_name("mzend")
			qp.set_accession(os.urandom(9).encode('hex'))
			qp.set_untiCvRef('MS:1000040')
			qp.set_cvRef('QCML:1000002') #interval end
			RunElement.add_QualityParameter(qp)
			#~ <QualityParameter accession="String" name="mzstart" value="87.30" cvRef="mzstart"/>
			#~ <QualityParameter accession="String" name="mzend" value="1999.98" cvRef="mzend"/>

		elif p.startswith('level 1: '):
			qp = mzQCML.CVParamType()
			qp.set_value(p.split(' ')[-3])
			qp.set_name("ms1 peaknumber")
			qp.set_accession(os.urandom(9).encode('hex'))
			qp.set_cvRef('QCML:1000004') # ms1 peaks
			qp.set_unitCvRef('QCML:1000003') #count
			RunElement.add_QualityParameter(qp)
		elif p.startswith('level 2: '):
			qp = mzQCML.CVParamType()
			qp.set_value(p.split(' ')[-1])
			qp.set_name("ms2 peaknumber")
			qp.set_accession(os.urandom(9).encode('hex'))
			qp.set_cvRef('QCML:1000005') # ms2 peaks
			qp.set_unitCvRef('QCML:1000003') #count
			RunElement.add_QualityParameter(qp)
			#~ <QualityParameter accession="String" name="#ms1spectra" value="4536" cvRef="mzstart"/>
			#~ <QualityParameter accession="String" name="#ms2spectra" value="30358" cvRef="mzstart"/>

#run r script get pictures


root.RunQuality.append(RunElement)
#TODO Set qualities

cvl = mzQCML.CVListType()

#~ that might be neccessary due to unicode
f = open(qcfilename, 'w')

xmlHeader = """<?xml version="1.0" encoding="ISO-8859-1"?>
<?xml-stylesheet type="text/xml" href="#stylesheet"?>
<!DOCTYPE catelog [
<!ATTLIST xsl:stylesheet
  id    ID  #REQUIRED>
]>"""

xmlStyle = """
<xsl:stylesheet id="stylesheet" version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
<xsl:template match="/">
	<html>
	<body>
		<h2>The Quality Parameters</h2>
			<table border="1">
				<tr bgcolor="#9acd32">
					<th>Parameter</th>
					<th>Value</th>
				</tr>
				<xsl:for-each select="MzQualityMLType/RunQuality/QualityParameter">
				<tr>
					<td><xsl:value-of select="@name" /></td>
					<td><xsl:value-of select="@value" /></td>
				</tr>
				</xsl:for-each>
			</table><br/>
		<h2>The Quality Plots</h2>
	<xsl:for-each select="MzQualityMLType/RunQuality/Attachment">
		<img>
			<xsl:attribute name="src">
				data:image/png;base64,<xsl:value-of select="binary" />
			</xsl:attribute>
		</img> <br/>
	</xsl:for-each>
	</body>
	</html>
</xsl:template>
</xsl:stylesheet>
"""

#~ h = etree.fromstring(xmlHeader)
#~ root.append(h)
#~ s = etree.fromstring(xmlStyle)
#~ root.append(s)

# c = codecs.lookup ("UTF-8")
# eout = c.streamwriter(f)
# root.export(eout, 0, namespace_='', name_='MzQualityMLType', namespacedef_='')
root.export(f, 0, namespace_='', name_='MzQualityMLType', namespacedef_='')
# eout.write ("\n")

f.close()
xmlstring = xmlHeader
f = open(qcfilename,'r')
L = f.readlines()
L.insert(0,xmlHeader)
L.insert(2,xmlStyle)
f.close()
f = open(qcfilename,'w')
f.write('\n'.join(L))
f.close()

#~ os.chdir(currentpath)

#~ for i in obo:
	#~ cv = mzQCML.CVParamType()
	#~ cv.set_uri("localhost")
	#~ cv.set_fullname("localhost")
	#~ cv.set_id("localhost")
	#~ cvl.add_cv(cv)


#~ def parseOBO(self, oboFile):
        #~ if os.path.exists(oboFile):
            #~ with open(oboFile) as obo:
                #~ collections = {}
                #~ collect = False
                #~ for line in obo:
                    #~ if line.strip() == '[Term]':
                        #~ self.add(collections)
                        #~ collect = True
                        #~ collections = {}
                    #~ else:
                        #~ if line.strip() != '' and collect == True:
                            #~ k = line.find(":")
                            #~ collections[line[:k]] = line[k+1:].strip()
        #~ else:
            #~ print("No obo file version {0} (psi-ms-{0}.obo) found.".format(self.version), file=sys.stderr)
            #~ exit(1)
        #~ return
