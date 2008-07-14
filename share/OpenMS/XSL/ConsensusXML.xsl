<?xml version="1.0"	encoding="utf-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:template match="consensusXML">
	  <html> 
	  <body>  
		  <xsl:apply-templates select="mapList"/>
		  <table>
				<tr><th></th><th>RT</th><th>MZ</th><th>Intensity</th><th>Map</th><th>id</th></tr>
		  	<xsl:apply-templates select="consensusElementList/consensusElement"/>
		  </table>
	  </body> 
	  </html>  
	</xsl:template>
	
	<xsl:template match="mapList">
		<table>
			<tr><td>id</td><td>name</td><td>label</td><td>size</td></tr>
		  <xsl:apply-templates select="map"/>
		</table>
	</xsl:template>
	
	<xsl:template match="map">
			<tr><td><xsl:value-of select="@id"/></td><td><xsl:value-of select="@name"/></td><td><xsl:value-of select="@label"/></td><td><xsl:value-of select="@size"/></td></tr>
	</xsl:template>
	
	<xsl:template match="consensusElement">
				<tr><td></td><td><xsl:value-of select="centroid/@rt"/></td><td><xsl:value-of select="centroid/@mz"/></td><td><xsl:value-of select="centroid/@it"/></td><td></td><td></td></tr>
			  <xsl:apply-templates select="groupedElementList/element">
		  			<xsl:sort select="@map" data-type="number" order="ascending"/>
			  </xsl:apply-templates>
	</xsl:template>

	<xsl:template match="element">
				<tr><td></td><td><xsl:value-of select="@rt"/></td><td><xsl:value-of select="@mz"/></td><td><xsl:value-of select="@it"/></td><td><xsl:value-of select="@id"/></td><td><xsl:value-of select="@map"/></td></tr>
	</xsl:template>


</xsl:stylesheet>


