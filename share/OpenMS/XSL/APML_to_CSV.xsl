<?xml version="1.0" ?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
<xsl:output method="html"/>



 <xsl:template match="/">
     <html> 
		 <body>
			<xsl:apply-templates select="apml/dataProcessing/software" />
		 <br/>
		 #Format: RT mz Intensity Charge Mass<br/>
       <xsl:apply-templates select="apml/data/peak_lists/peak_list/features/feature">
 				<xsl:sort select="@id" />
       </xsl:apply-templates>
		 </body>
     </html> 
 </xsl:template>

  <xsl:template match="feature">
		<!--<xsl:value-of select="@id" />-->
		<xsl:apply-templates select="coordinate" />

  </xsl:template>

  <xsl:template match="coordinate">
    <xsl:value-of select="@rt" /><xsl:text> </xsl:text>
    <xsl:value-of select="@mz" /><xsl:text> </xsl:text>
    <xsl:value-of select="@intensity" /><xsl:text> </xsl:text>
    <xsl:value-of select="@charge" /><xsl:text> </xsl:text>
    <xsl:value-of select="@mass" /><br/>
  </xsl:template>

		 <xsl:template match="software">
			<xsl:text>#</xsl:text>
		 	<xsl:value-of select="@name" /><xsl:text> (</xsl:text>
		 	<xsl:value-of select="@version" /><xsl:text>) - </xsl:text>
		 	<xsl:value-of select="@type" />
	</xsl:template>



</xsl:stylesheet>
