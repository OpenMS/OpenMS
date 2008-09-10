<?xml version="1.0"	encoding="utf-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
<xsl:template match="IdXML">
  <html> 
  <body>  
	  <xsl:apply-templates select="//IdentificationRun"/>
  </body> 
  </html>  
</xsl:template>

<xsl:template match="SearchParameters">
    <tr>
  	<td>Database</td>
  	<td align="center">
    		<xsl:value-of select="@db"/> 
  	</td>
  	</tr>
    <tr>
  	<td>Database version</td>
  	<td align="center">
    		<xsl:value-of select="@db_version"/> 
  	</td>
  	</tr>
    <tr>
  	<td>Taxonomy</td>
  	<td align="center">
    		<xsl:value-of select="@taxonomy"/> 
  	</td>
  	</tr>
    <tr>
  	<td>Mass type</td>
  	<td align="center">
    		<xsl:value-of select="@mass_type"/> 
  	</td>
  	</tr>
    <tr>
  	<td>Charges</td>
  	<td align="center">
    		<xsl:value-of select="@charges"/> 
  	</td>
  	</tr>
    <tr>
  	<td>Enzyme</td>
  	<td align="center">
    		<xsl:value-of select="@enzyme"/> 
  	</td>
  	</tr>
    <tr>
  	<td>Missed cleavages</td>
  	<td align="center">
    		<xsl:value-of select="@missed_cleavages"/> 
  	</td>
  	</tr>
    <tr>
  	<td>Precursor peak tolerance</td>
  	<td align="center">
    		<xsl:value-of select="@precursor_peak_tolerance"/> 
  	</td>
  	</tr>
    <tr>
  	<td>Peak mass tolerance</td>
  	<td align="center">
    		<xsl:value-of select="@peak_mass_tolerance"/> 
  	</td>
  	</tr>
</xsl:template>

<xsl:template match="IdentificationRun">
	<table border="1" cellpadding="2" style="border-style:solid; border-collapse:collapse;" width="100%">
    <tr><th>Search parameter</th>
    <th>Value</th></tr>
    <tr>
  	<td>Date</td>
  	<td align="center">
    		<xsl:value-of select="@date"/> 
  	</td>
  	</tr>
    <tr>
  	<td>Search engine</td>
  	<td align="center">
    		<xsl:value-of select="@search_engine"/> 
  	</td>
		</tr>
    <tr>
  	<td>Search engine version</td>
  	<td align="center">
    		<xsl:value-of select="@search_engine_version"/> 
  	</td>
		</tr>
	 	<xsl:variable name="search_parameters_reference">
	 		<xsl:value-of select="@search_parameters_ref"/>
  	</xsl:variable>
	  <xsl:apply-templates select="parent::IdXML/SearchParameters[@id=$search_parameters_reference]"/>
	</table>
	<br/>
	<table border="1" cellpadding="2" style="border-style:solid; border-collapse:collapse;" width="100%">
    <tr><th>Protein</th>
    	<th align="center">Score</th>
    	<th align="center"> Corresponding identified peptides </th>	
    </tr>
	  <xsl:apply-templates select="ProteinIdentification/ProteinHit">
  			<xsl:sort select="@score" data-type="number" order="descending"/>
	  </xsl:apply-templates>
	</table>
	<br/>
	<table border="1" cellpadding="2" width="100%" style="border-style:solid; border-collapse:collapse;">
  	<tr><th><nobr>Unassigned peptides</nobr></th>
   	<th>Score</th>
   	<th><nobr>Charge</nobr></th>
   	<th><nobr>RT</nobr></th>
   	<th><nobr>MZ</nobr></th>
   	<th><nobr>AA before</nobr></th>
   	<th><nobr>AA after</nobr></th>
	  <xsl:if test="PeptideIdentification/PeptideHit[not(@protein_refs)]/UserParam">
			<th width="20%">Additional information</th>
	  </xsl:if>
 		</tr>
		<xsl:choose>
		  <xsl:when test="PeptideIdentification/PeptideHit[not(@protein_refs)]/ancestor::PeptideIdentification[@higher_score_better='true']">
			  <xsl:apply-templates select="PeptideIdentification/PeptideHit[not(@protein_refs)]">
	  				<xsl:sort select="@score" data-type="number" order="descending"/>	
			  </xsl:apply-templates>
  		</xsl:when>
		  <xsl:otherwise>
			  <xsl:apply-templates select="PeptideIdentification/PeptideHit[not(@protein_refs)]">
	  				<xsl:sort select="@score" data-type="number" order="ascending"/>	
			  </xsl:apply-templates>
			</xsl:otherwise>
		</xsl:choose>
	</table>
	<br/><br/>
</xsl:template>

<xsl:template match="ProteinHit">
  <tr>
    <td>
      <xsl:value-of select="@accession"/> 
    </td>
    <td align="center">
      <xsl:value-of select="@score"/> 
    </td>
    <td>
	  	<xsl:variable name="id_string">
	  		<xsl:value-of select="@id"/><xsl:text> </xsl:text> 
 	  	</xsl:variable>
	  	<xsl:variable name="id_string2">
	  		<xsl:value-of select="@id"/>
 	  	</xsl:variable>
    
			<table border="1" cellpadding="2" width="100%" style="border-style:none; border-collapse:collapse;" >
    	<tr><th width="20%">Sequence</th>
	  	<th>Score</th>
	  	<th><nobr>charge</nobr></th>
	   	<th><nobr>RT</nobr></th>
	   	<th><nobr>MZ</nobr></th>
	   	<th><nobr>AA before</nobr></th>
	   	<th><nobr>AA after</nobr></th>
	 		<th width="20%"><nobr>Additional information</nobr></th>
	  	</tr>
			<xsl:choose>
			  <xsl:when test="ancestor::IdentificationRun/PeptideIdentification/PeptideHit[contains(@protein_refs, $id_string) or (contains(@protein_refs, $id_string2) and substring-after(@protein_refs, $id_string2) = '')]/ancestor::PeptideIdentification[@higher_score_better='true']">
		  		<xsl:apply-templates select="ancestor::IdentificationRun/PeptideIdentification/PeptideHit[contains(@protein_refs, $id_string) or (contains(@protein_refs, $id_string2) and substring-after(@protein_refs, $id_string2) = '')]">
		  				<xsl:sort select="@score" data-type="number" order="descending"/>
					</xsl:apply-templates>
				</xsl:when>
				<xsl:otherwise>
		  		<xsl:apply-templates select="ancestor::IdentificationRun/PeptideIdentification/PeptideHit[contains(@protein_refs, $id_string) or (contains(@protein_refs, $id_string2) and substring-after(@protein_refs, $id_string2) = '')]">
		  				<xsl:sort select="@score" data-type="number" order="ascending"/>
					</xsl:apply-templates>
				</xsl:otherwise>
			</xsl:choose>
			</table>
    </td>
  </tr>
</xsl:template>
<xsl:template match="PeptideHit">
	<tr>
  	<td>
    		<xsl:value-of select="@sequence"/> 
  	</td>
  	<td align="center">
    		<xsl:value-of select="@score"/> 
  	</td>
  	<td align="center">
    		<xsl:value-of select="@charge"/> 
  	</td>
  	<td align="center">
    		<xsl:value-of select="parent::PeptideIdentification/@RT"/> 
  	</td>
  	<td align="center">
    		<xsl:value-of select="parent::PeptideIdentification/@MZ"/> 
  	</td>
  	<td align="center">
    		<xsl:value-of select="@aa_before"/> 
  	</td>
  	<td align="center">
    		<xsl:value-of select="@aa_after"/> 
  	</td>
  	<td><nobr>
  		<xsl:if test="UserParam">
 				<xsl:for-each select="UserParam">
  				<xsl:value-of select="@name"/><xsl:text>: </xsl:text><xsl:value-of select="@value"/><br/>
	  		</xsl:for-each>
  	</xsl:if>
  	</nobr></td>
	</tr>
</xsl:template>
</xsl:stylesheet>


