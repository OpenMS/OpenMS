<?xml version="1.0"	encoding="utf-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:template match="consensusXML">
	  <html> 
	  <body>  
		  <xsl:apply-templates select="mapList"/>
		  <br/>
		  <b>Consensus elements:</b>
		  <table border="1" align="center" style="border-style:solid; border-collapse:collapse;" width="100%" cellpadding="2">
				<tr>
					<th><nobr></nobr></th>
					<th><nobr>RT</nobr></th>
					<th><nobr>MZ</nobr></th>
					<th><nobr>Intensity</nobr></th>
					<th><nobr>Quality</nobr></th>
					<th><nobr>Map Id</nobr></th>
					<th><nobr>Element Id</nobr></th>
				</tr>
		  	<xsl:apply-templates select="consensusElementList/consensusElement"/>
		  </table>
	  </body> 
	  </html>  
	</xsl:template>
	
	<xsl:template match="mapList">
		<b>Source map list:</b>
		<table border="1" style="border-style:solid; border-collapse:collapse;" cellpadding="2">
			<tr>
				<th><nobr>Map Id</nobr></th>
				<th><nobr>File name</nobr></th>
				<th><nobr>Label</nobr></th>
				<th><nobr>Size</nobr></th>
			</tr>
		  <xsl:apply-templates select="map"/>
		</table>
	</xsl:template>
	
	<xsl:template match="map">
		<tr>
			<td align="center"><xsl:value-of select="@id"/></td>
			<td align="center"><xsl:value-of select="@name"/></td>
			<td align="center"><xsl:value-of select="@label"/></td>
			<td align="center"><xsl:value-of select="@size"/></td>
		</tr>
	</xsl:template>
	
	<xsl:template match="consensusElement">
		<tr>
			<td align="center"><b><nobr>Consensus element <xsl:value-of select="position()"/></nobr></b></td>
			<td align="center"><b><xsl:value-of select="centroid/@rt"/></b></td>
			<td align="center"><b><xsl:value-of select="centroid/@mz"/></b></td>
			<td align="center"><b><xsl:value-of select="centroid/@it"/></b></td>
			<td align="center"><b><xsl:value-of select="@quality"/></b></td>
			<td align="center"><b></b></td>
			<td align="center"><b></b></td>
		</tr>
		<xsl:apply-templates select="groupedElementList/element">
			<xsl:sort select="@map" data-type="number" order="ascending"/>
		</xsl:apply-templates>
	</xsl:template>

	<xsl:template match="element">
		<tr>
			<td align="center"></td>
			<td align="center"><xsl:value-of select="@rt"/></td>
			<td align="center"><xsl:value-of select="@mz"/></td>
			<td align="center"><xsl:value-of select="@it"/></td>
			<td align="center"></td>
			<td align="center"><xsl:value-of select="@map"/></td>
			<td align="center"><xsl:value-of select="@id"/></td>
		</tr>
	</xsl:template>

</xsl:stylesheet>


