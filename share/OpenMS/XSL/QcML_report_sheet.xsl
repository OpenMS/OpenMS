<xsl:stylesheet id="stylesheet" version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:ns="http://www.prime-xs.eu/ms/qcml" xmlns="">

	<xsl:template match="/">
		<html>
      <style><!-- css settings for div elements -->
        H1 {
          font-family:Arial; font-size:20pt
        }
        div.rahmen {  
          background: #FFFAF0;  
          border: 1px solid;  
          padding: 5px;  
          margin: 40px;  
        }
        div.thickline {  
          border-width: 2px;  
          border-style: solid;  
        }  
      </style>
		
			<body>
<!-- Document caption-->
      <center><H1>Run Quality Report</H1></center>
<!-- MetaData Box-->
      <p style="border-color:#000000; border-width:2px; border-style:solid; padding:4px">
        Metadata
        <center><table border="1" rules="none" color="black">
          <tr>       
            <td>Filename</td>
            <td>
              <xsl:for-each select="ns:qcML/ns:runQuality">
                <xsl:apply-templates select="ns:qualityParameter[@accession = 'MS:1000577']"/>
              </xsl:for-each>
            </td>
          </tr>
          <tr>
            <td>Instrument</td>
            <td>
              <xsl:for-each select="ns:qcML/ns:runQuality">
                <xsl:apply-templates select="ns:qualityParameter[@accession = 'MS:1000031']"/>
              </xsl:for-each>
            </td>
          </tr>
          <tr>
            <td>Date</td>
            <td>
              <xsl:for-each select="ns:qcML/ns:runQuality">
                <xsl:apply-templates select="ns:qualityParameter[@accession = 'MS:1000747']"/>
              </xsl:for-each>
            </td>
          </tr>								
        </table></center>			
        <br/>
      </p>
        
<!-- Overview Box-->	
      <p style="border-color:#000000; border-width:2px; border-style:solid; padding:4px">
				Overview Plots
				<table border="1">
          <tr>
            <td>
            TIC
              <xsl:for-each select="ns:qcML/ns:runQuality">
                <xsl:apply-templates select="ns:attachment[@accession = 'MS:1000235']"/>
              </xsl:for-each>
            </td>
            <td>
            Map overview
              <xsl:for-each select="ns:qcML/ns:runQuality">
                <xsl:apply-templates select="ns:attachment[@accession = 'QC:0000055']"/>
              </xsl:for-each>           
            </td>
            <td>
            Mass accuracy
              <xsl:for-each select="ns:qcML/ns:runQuality">
                <xsl:apply-templates select="ns:attachment[@accession = 'QC:0000053']"/>
              </xsl:for-each>
            </td>
          </tr>
        </table>
      </p>
				
<!-- Details Box-->
			<p style="border-color:#000000; border-width:2px; border-style:solid; padding:4px">
				Details
				<table>
          <tr>
            <td rowspan="3" width="200" height="20">Aquisition details</td> <td width="400" height="20">Name</td>	<td width="400" height="20">Value</td>
          </tr>				
          <tr>			
              <xsl:for-each select="ns:qcML/ns:runQuality">
                <xsl:apply-templates select="ns:qualityParameter[@accession = 'QC:0000006']"/>
              </xsl:for-each>
          </tr>				
          <tr>			
              <xsl:for-each select="ns:qcML/ns:runQuality">
                <xsl:apply-templates select="ns:qualityParameter[@accession = 'QC:0000007']"/>
              </xsl:for-each>
          </tr>			
				</table>
				<br/>
				
				<table>
          <tr>
            <td rowspan="3" width="200" height="20">Aquisition ranges</td> <td width="400" height="20"></td> <td width="400" height="20">Minimum</td>	<td width="400" height="20">Maximum</td>
          </tr>				
          <tr>
            <td>RT [s] </td>
              <xsl:for-each select="ns:qcML/ns:runQuality">
                <xsl:apply-templates select="ns:attachment[@accession = 'QC:0000012']"/>
              </xsl:for-each>
          </tr>				
          <tr>	
            <td>MZ [amu] </td>          
              <xsl:for-each select="ns:qcML/ns:runQuality">
                <xsl:apply-templates select="ns:attachment[@accession = 'QC:0000009']"/>
              </xsl:for-each>
          </tr>			
				</table>
				<br/>

				<table>
          <tr>
            <td rowspan="7" width="200" height="20">Identification details</td> <td width="400" height="20">Name</td>	<td width="400" height="20">Value</td>
          </tr>				
          <tr>			
              <xsl:for-each select="ns:qcML/ns:runQuality">
                <xsl:apply-templates select="ns:qualityParameter[@accession = 'QC:0000029']"/>
              </xsl:for-each>
          </tr>				
          <tr>			
              <xsl:for-each select="ns:qcML/ns:runQuality">
                <xsl:apply-templates select="ns:qualityParameter[@accession = 'QC:0000030']"/>
              </xsl:for-each>
          </tr>	          
          <tr>			
              <xsl:for-each select="ns:qcML/ns:runQuality">
                <xsl:apply-templates select="ns:qualityParameter[@accession = 'QC:0000031']"/>
              </xsl:for-each>
          </tr>			
          <tr>			
              <xsl:for-each select="ns:qcML/ns:runQuality">
                <xsl:apply-templates select="ns:qualityParameter[@accession = 'QC:0000032']"/>
              </xsl:for-each>
          </tr>			
          <tr>			
              <xsl:for-each select="ns:qcML/ns:runQuality">
                <xsl:apply-templates select="ns:qualityParameter[@accession = 'QC:0000033']"/>
              </xsl:for-each>
          </tr>			
          <tr>			
              <xsl:for-each select="ns:qcML/ns:runQuality">
                <xsl:apply-templates select="ns:qualityParameter[@accession = 'QC:0000037']"/>
              </xsl:for-each>
          </tr>			
				</table>
				<br/>

				<table>
				<tr>
          <td rowspan="1" width="200" height="30">Database</td>		
          <td>
              <xsl:for-each select="ns:qcML/ns:runQuality">
                <xsl:apply-templates select="ns:attachment[@accession = 'QC:0000026']"/>
              </xsl:for-each>
          </td>
				</tr>
				</table>
				<br/>
				

				<table>
          <tr> 
            <td rowspan="2" width="200" height="20">Feature Finding</td>				
            <td width="400" height="30">Name</td>				
            <td width="400" height="30">Value</td>
          </tr>				
          <tr>			
              <xsl:for-each select="ns:qcML/ns:runQuality">
                <xsl:apply-templates select="ns:qualityParameter[@accession = 'QC:0000046']"/>
              </xsl:for-each>
          </tr>	          
				</table>

				</p>
				
				<p style="border-color:#000000; border-width:2px; border-style:solid; padding:4px">
				<table>
              <xsl:for-each select="ns:qcML/ns:runQuality">
          
                <!-- 
                <xsl:call-template name="qp-attachments" select="ns:attachment[not(@accession = 'QC:0000026') and not(@accession = 'QC:0000009') and not(@accession = 'QC:0000012') and not(@accession = 'QC:0000053') and not(@accession = 'QC:0000055') and not(@accession = 'MS:1000235')]"/>
                <xsl:apply-templates select="ns:qualityParameter[not(@accession = 'QC:0000046') and not(@accession = 'QC:0000037') and not(@accession = 'QC:0000029') and not(@accession = 'QC:0000030') and not(@accession = 'QC:0000031') and not(@accession = 'QC:0000032') and not(@accession = 'QC:0000033') and not(@accession = 'QC:0000037') and not(@accession = 'QC:0000007') and not(@accession = 'QC:0000006') and not(@accession = 'MS:1000747') and not(@accession = 'MS:1000031') and not(@accession = 'MS:1000577')]"/>
                -->
                <xsl:apply-templates select="ns:attachment[not(@accession = 'QC:0000026') and not(@accession = 'QC:0000009') and not(@accession = 'QC:0000012') and not(@accession = 'QC:0000053') and not(@accession = 'QC:0000055') and not(@accession = 'MS:1000235')]"/>
              </xsl:for-each>
        
        <!-- rest of the qps
          <xsl:for-each select="ns:qcML/ns:runQuality">
          <tr>			
              <xsl:apply-templates select="ns:qualityParameter[not(@accession = 'MS:1000577') and not(@accession = 'MS:1000031') and not(@accession = 'MS:1000747' ) and not(@accession = 'QC:0000046')]"/>
          </tr>	 
          </xsl:for-each>
        
				rest of the qps-->
				</table>
				</p>
        
				<xsl:for-each select="ns:qcML/ns:setQuality">
					<xsl:apply-templates/>
				</xsl:for-each>

      </body>
		</html>
	</xsl:template>
  
<!-- for MetaData Box qp template-->
	<xsl:template match="ns:qualityParameter[(@accession = 'MS:1000577') or (@accession = 'MS:1000031') or (@accession = 'MS:1000747' )]">
		<b><xsl:value-of select="@value"/></b>
	</xsl:template>
<!-- for further detail Box qp template-->	
	<xsl:template match="ns:qualityParameter[not(@accession = 'MS:1000577') and not(@accession = 'MS:1000031') and not(@accession = 'MS:1000747' ) and @value]">
		<td><xsl:value-of select="@name"/></td>
		<td><xsl:value-of select="@value"/></td>
	</xsl:template>
<!-- for generic qp name + attachment output if no value in qp-->	
	<xsl:template match="ns:qualityParameter[not(@value)]">
		<td><xsl:value-of select="@name"/></td>
		<td><xsl:call-template name="qp-attachments">
			<xsl:with-param name="qpref" select="@ID"/>
		</xsl:call-template></td>
	</xsl:template>
<!-- called from generic qp template -->	
	<xsl:template name="qp-attachments">
		<xsl:param name="qpref"/>
		<xsl:for-each select="../ns:attachment[@qualityParameterRef=$qpref]"> <xsl:value-of select="@name"/><br/>
			<xsl:choose>
				<xsl:when test="ns:binary">
					<img>
						<xsl:attribute name="src"> data:image/png;base64,<xsl:value-of select="ns:binary"/>
						</xsl:attribute>
					</img>
					<br/>
				</xsl:when>
				<xsl:otherwise>
					<table border="0">
						<tr bgcolor="#B2CCFF">
							<xsl:call-template name="output-header">
								<xsl:with-param name="list"><xsl:value-of select="ns:table/ns:tableColumnTypes"/></xsl:with-param>
							</xsl:call-template>
						</tr>
						<xsl:for-each select="ns:table/ns:tableRowValues">
							<tr>
								<xsl:call-template name="output-row">
									<xsl:with-param name="list"><xsl:value-of select="." /></xsl:with-param>
								</xsl:call-template>
							</tr>
						</xsl:for-each>
					</table><br/>
				</xsl:otherwise>
			</xsl:choose>
		</xsl:for-each>
	</xsl:template>  
<!-- for Overview plots Box template-->	 
	<xsl:template match="ns:attachment[(@accession = 'MS:1000235') or (@accession = 'QC:0000055') or (@accession = 'QC:0000053' )]">
			<xsl:choose>
				<xsl:when test="ns:binary">
					<img style="width: 100%">
						<xsl:attribute name="src"> data:image/png;base64,<xsl:value-of select="ns:binary"/>
						</xsl:attribute>
					</img>
				</xsl:when>
				<xsl:otherwise>
						<xsl:for-each select="ns:table/ns:tableRowValues">
								<xsl:call-template name="output-row">
									<xsl:with-param name="list"><xsl:value-of select="." /></xsl:with-param>
								</xsl:call-template>
  						</xsl:for-each>
 				</xsl:otherwise>
			</xsl:choose>
  </xsl:template>
<!-- for Details Box aquisition range template-->	 
	<xsl:template match="ns:attachment[(@accession = 'QC:0000009') or (@accession = 'QC:0000012' ) or (@accession = 'QC:0000026' )]">
    <xsl:for-each select="ns:table/ns:tableRowValues">
				<xsl:call-template name="output-row">
					<xsl:with-param name="list"><xsl:value-of select="." /></xsl:with-param>
				</xsl:call-template>
  		</xsl:for-each>
  </xsl:template>
<!-- generic attachment template for att as table row-->
	<xsl:template match="ns:attachment">
    <tr>	
      <td><xsl:value-of select="@name"/></td>
      <td>
        <xsl:choose>
          <xsl:when test="ns:binary">
            <img>
              <xsl:attribute name="src"> data:image/png;base64,<xsl:value-of select="ns:binary"/>
              </xsl:attribute>
            </img>
          </xsl:when>
          <xsl:otherwise> subtable </xsl:otherwise> 
        </xsl:choose>   				
      </td>
    </tr>
	</xsl:template>
 
	<xsl:template name="output-header">
		<xsl:param name="list"/>
		<xsl:variable name="newlist" select="concat(normalize-space($list), ' ')"/>
		<xsl:variable name="first" select="substring-before($newlist, ' ')"/>
		<xsl:variable name="remaining" select="substring-after($newlist, ' ')"/>
		<th>
			<xsl:value-of select="$first"/>
		</th>
		<xsl:if test="$remaining">
			<xsl:call-template name="output-header">
				<xsl:with-param name="list" select="$remaining"/>
			</xsl:call-template>
		</xsl:if>
	</xsl:template>
	
  <xsl:template name="output-row">
		<xsl:param name="list"/>
		<xsl:variable name="newlist" select="concat(normalize-space($list), ' ')"/>
		<xsl:variable name="first" select="substring-before($newlist, ' ')"/>
		<xsl:variable name="remaining" select="substring-after($newlist, ' ')"/>
		<td>
			<xsl:value-of select="$first"/>
		</td>
		<xsl:if test="$remaining">
			<xsl:call-template name="output-row">
				<xsl:with-param name="list" select="$remaining"/>
			</xsl:call-template>
		</xsl:if>
	</xsl:template>
</xsl:stylesheet>