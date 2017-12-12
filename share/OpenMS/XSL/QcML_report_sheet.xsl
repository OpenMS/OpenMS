<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet id="openms-qc-stylesheet" version="1.1" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:ns="https://github.com/qcML/qcml" xmlns="">
    <xsl:template match="/">
        <html>
            <head>
                <script type="text/javascript" src="http://d3js.org/d3.v3.min.js"></script>
                <style>
 <!-- css settings for div elements -->
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

                path { 4
                        stroke: steelblue;
                        stroke-width: 2;
                        fill: none;
                }

                .axis path,
                .axis line {
                  fill: none;
                  stroke: #000;
                  shape-rendering: crispEdges;
                }

                .line {
                  fill: none;
                  stroke: steelblue;
                  stroke-width: 1.5px;
                }

                .grid .tick {
                        stroke: lightgrey;
                        opacity: 0.7;
                }

                .grid path {
                          stroke-width: 0;
                }

                .dot {
                        stroke: #000;
                        stroke-width: 2;
                        radius: 5;
                }

                div.tooltip {
                        position: absolute;
                        text-align: center;
                        width: 80px;
                        height: 40px;
                        padding: 2px;
                        font: 12px sans-serif;
                        background: #eee;
                        border: 1px solid #ccc;
                        border-radius: 8px;
                }
                </style>
            </head>
            <body>
<!-- Document caption-->
                <xsl:choose>
                        <xsl:when test="ns:qcML/ns:setQuality">
                                <center><H1>Set Quality Report</H1></center>
                        </xsl:when>
                        <xsl:otherwise>
                                <center><H1>Run Quality Report</H1></center>
                        </xsl:otherwise>
                </xsl:choose>

<!-- MetaData Box-->
                <p style="border-color:#000000; border-width:2px; border-style:solid; padding:4px">
                Metadata
                <xsl:for-each select="ns:qcML/ns:runQuality">
                        <center><table border="1" rules="none" color="black">
                                <tr>
                                        <td>Filename</td>
                                <td>
                                        <xsl:apply-templates select="ns:qualityParameter[@accession = 'MS:1000577']"/>
                                </td>
                </tr>
                <tr>
                                <td>Instrument</td>
                                <td>
                                        <xsl:apply-templates select="ns:qualityParameter[@accession = 'MS:1000031']"/>
                                </td>
                </tr>
                <tr>
                                        <td>Date</td>
                                <td>
                                        <xsl:apply-templates select="ns:qualityParameter[@accession = 'MS:1000747']"/>
                                </td>
                </tr>
                </table></center>
                <br/>
                        </xsl:for-each>
                </p>

<!-- Overview Box-->
                <xsl:choose>
                    <xsl:when test="ns:qcML/ns:setQuality">
                        <xsl:variable name="set_setting" select="true()"/>
                        <p style="border-color:#000000; border-width:2px; border-style:solid; padding:4px">
                            Control charts
                        </p>
                        <div id="control-chart"></div>
    
                        <script>
var data_control_chart = [<xsl:for-each select="ns:qcML/ns:runQuality">
    {
            idx: "<xsl:value-of select="position()" />" ,
            MS1count: "<xsl:value-of select="ns:qualityParameter[@accession = 'QC:0000006']/@value"/>" ,
            MS2count: "<xsl:value-of select="ns:qualityParameter[@accession = 'QC:0000007']/@value"/>" ,
            PSMcount: "<xsl:value-of select="ns:qualityParameter[@accession = 'QC:0000029']/@value"/>" ,
            FEATUREcount: "<xsl:value-of select="ns:qualityParameter[@accession = 'QC:0000046']/@value"/>" ,
            idFEATUREcount: "<xsl:value-of select="ns:qualityParameter[@accession = 'QC:0000058']/@value"/>" ,
    },
    </xsl:for-each>];
                        </script>
                        <script>
var margin = {top: 20, right: 100, bottom: 30, left: 100},
    width = 960 - margin.left - margin.right,
    height = 500 - margin.top - margin.bottom;

var color = d3.scale.category10();

data_control_chart.forEach( function(d) {
        console.log(d);
        color.domain(d3.keys(d).filter(function(key) { return key !== "idx"; }));
});
console.log(data_control_chart);
console.log(color.domain());

var svg = d3.select("#control-chart").append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .append("g")
    .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

var div = d3.select("#control-chart").append("div").attr("class", "tooltip")
            .style("opacity", 0);

var control_values = color.domain().map(function(name) {
                return {
                  name: name,
                  values: data_control_chart.map(function(d) {
                        return {idx: d.idx, vals: +d[name]};
                  })
                };
          });
        console.log(control_values);

var line = d3.svg.line()
                //.interpolate("basis")
                .x(function(d) { return x(d.idx); })
                .y(function(d) { return y(d.vals); });

var x = d3.scale.linear()
                .range([0, width]);

var y = d3.scale.linear()
                .range([height, 0]);

function xAxis()  {
                return d3.svg.axis()
                .scale(x)
                .ticks(data_control_chart.length)
                .tickSubdivide(0)
                .orient("bottom")
                }

function yAxis() {
                return d3.svg.axis()
                .scale(y)
                .orient("left")
                }

x.domain(d3.extent(data_control_chart, function(d) { return d.idx; }));

y.domain([
    d3.min(control_values, function(c) { return d3.min(c.values, function(v) { return v.vals; }); }),
    d3.max(control_values, function(c) { return d3.max(c.values, function(v) { return v.vals; }); })
]);

svg.append("g")
      .attr("class", "grid")
      .attr("transform", "translate(0," + height + ")")
      .call(xAxis()
                        .tickSize(-height, 0, 0)
                        .tickFormat("")
                );

svg.append("g")
      .attr("class", "x axis")
           .attr("transform", "translate(0," + height + ")")
      .call(xAxis()	  )
        .append("text")
      .attr("transform", "translate(" + width +",0)")
      .attr("y", -6)
      .attr("x", -6)
      .attr("dx", ".71em")
      .style("text-anchor", "end")
      .text("runs");

svg.append("g")
      .attr("class", "y axis")
      .call(yAxis())
    .append("text")
      .attr("transform", "rotate(-90)")
      .attr("y", 6)
      .attr("dy", ".71em")
      .style("text-anchor", "end")
      .text("counts");

var control_metric = svg.selectAll(".control_metric")
      .data(control_values)
    .enter().append("g")
      .attr("class", "control_metric");

control_metric.append("path")
      .attr("class", "line")
      .attr("d", function(d) { return line(d.values); })
      .style("stroke", function(d) { return color(d.name); });

control_metric.append("text")
      .datum(function(d) { return {name: d.name, value: d.values[d.values.length - 1]}; })
      .attr("transform", function(d) { return "translate(" + x(d.value.idx) + "," + y(d.value.vals) + ")"; })
      .attr("x", 3)
      .attr("dy", ".35em")
      .text(function(d) { return d.name; });

control_metric.append("g").selectAll(".dot")
                .data(function (d) { return d.values })
                .enter().append("circle")
                .attr("stroke", function (d) { return color(this.parentNode.__data__.name) })
                .attr("cx", function (d) { return x(d.idx); })
                .attr("cy", function (d) { return y(d.vals); })
                .attr("r", 2)
                .attr("fill", "white").attr("fill-opacity", .5)
                .attr("stroke-width", 2)
                        .on("mouseover", function (d) {
                                div.transition()
                                        .duration(200)
                                        .style("opacity", .9);
                                div.html(
                                        this.parentNode.__data__.name + "&lt;br&gt;" + d.vals + "&lt;br&gt;" + "&lt;a href='#" + d.idx + "'&gt; run " +d.idx + "&lt;/a&gt;")
                                        .style("left", (d3.event.pageX + 18) + "px").style("top", (d3.event.pageY - 18) + "px").attr('r', 8);
                                d3.select(this).attr('r', 4)
                        })
                        .on("mouseout", function (d) {
                                div.transition().duration(2200).style("opacity", 0)
                                d3.select(this).attr('r', 2);
                        });
                        </script>


                    </xsl:when>
                    <xsl:otherwise>
                        <xsl:variable name="set_setting" select="false()"/>
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
                    </xsl:otherwise>
                </xsl:choose>



<!-- Details Box-->
                <xsl:for-each select="ns:qcML/ns:runQuality">
                        <p style="border-color:#000000; border-width:2px; border-style:solid; padding:4px">
                                <a name="{position()}"> Details <xsl:value-of select="@ID"/> </a>

                <table>
                  <tr>
                    <td rowspan="3" width="200" height="20">Aquisition details</td> <td width="400" height="20">Name</td>	<td width="400" height="20">Value</td>
                  </tr>
                  <tr>
                                        <xsl:apply-templates select="ns:qualityParameter[@accession = 'QC:0000006']"/>
                  </tr>
                  <tr>
                                        <xsl:apply-templates select="ns:qualityParameter[@accession = 'QC:0000007']"/>
                  </tr>
                </table>
                <br/>

                <table>
                  <tr>
                    <td rowspan="3" width="200" height="20">Aquisition ranges</td> <td width="400" height="20"></td> <td width="400" height="20">Minimum</td>	<td width="400" height="20">Maximum</td>
                  </tr>
                  <tr>
                    <td>RT [s] </td>
                                                <xsl:apply-templates select="ns:attachment[@accession = 'QC:0000012']"/>
                  </tr>
                  <tr>
                    <td>MZ [amu] </td>
                        <xsl:apply-templates select="ns:attachment[@accession = 'QC:0000009']"/>
                  </tr>
                </table>
                <br/>

                <table>
                  <tr>
                    <td rowspan="7" width="200" height="20">Identification details</td> <td width="400" height="20">Name</td>	<td width="400" height="20">Value</td>
                  </tr>
                  <tr>
                        <xsl:apply-templates select="ns:qualityParameter[@accession = 'QC:0000029']"/>
                  </tr>
                  <tr>
                        <xsl:apply-templates select="ns:qualityParameter[@accession = 'QC:0000030']"/>
                  </tr>
                  <tr>
                        <xsl:apply-templates select="ns:qualityParameter[@accession = 'QC:0000031']"/>
                  </tr>
                  <tr>
                        <xsl:apply-templates select="ns:qualityParameter[@accession = 'QC:0000032']"/>
                  </tr>
                  <tr>
                        <xsl:apply-templates select="ns:qualityParameter[@accession = 'QC:0000033']"/>
                  </tr>
                  <tr>
                        <xsl:apply-templates select="ns:qualityParameter[@accession = 'QC:0000037']"/>
                  </tr>
                </table>
                <br/>

                <table>
                  <tr>
                      <td rowspan="1" width="200" height="30">Database</td>
                      <td>
                            <xsl:apply-templates select="ns:attachment[@accession = 'QC:0000026']"/>
                      </td>
                  </tr>
                </table>
                <br/>

                <table>
                  <tr>
                    <td rowspan="3" width="200" height="20">Feature Finding</td>
                    <td width="400" height="30">Name</td>
                    <td width="400" height="30">Value</td>
                  </tr>
                  <tr>
                        <xsl:apply-templates select="ns:qualityParameter[@accession = 'QC:0000046']"/>
                  </tr>
                  <tr>
                        <xsl:apply-templates select="ns:qualityParameter[@accession = 'QC:0000058']"/>
                  </tr>
                </table>
                </p>

                <p style="border-color:#000000; border-width:2px; border-style:solid; padding:4px">
                Metric plots <xsl:value-of select="@ID"/>
                <br/>

                <xsl:if test="ns:qcML/ns:setQuality">
                    <table border="1">
                      <tr>
                        <td>
                        TIC
                            <xsl:apply-templates select="ns:attachment[@accession = 'MS:1000235']"/>
                        </td>
                        <td>
                        Map overview
                            <xsl:apply-templates select="ns:attachment[@accession = 'QC:0000055']"/>
                        </td>
                        <td>
                        Mass accuracy
                            <xsl:apply-templates select="ns:attachment[@accession = 'QC:0000053']"/>
                        </td>
                      </tr>
                    </table>
                </xsl:if>

              <table>
                      <xsl:apply-templates select="ns:attachment[
                      not(@accession = 'QC:0000055') and not(@accession = 'QC:0000053') and not(@accession = 'MS:1000235') and
                      not(@accession = 'QC:0000026') and not(@accession = 'QC:0000009') and not(@accession = 'QC:0000012') and
                      not(@accession = 'QC:0000044') and not(@accession = 'QC:0000022') and not(@accession = 'QC:0000038') and
                      not(@accession = 'QC:0000047') and not(@accession = 'QC:0000018')]"/>
              </table>
           </p>
        </xsl:for-each>
                        <!--
                        <xsl:for-each select="ns:qcML/ns:setQuality">
                        <xsl:apply-templates/>
                        </xsl:for-each>
                        -->

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
                    <img style="width: 100%" alt="?missing?">
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
<!-- generic output table header-->
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
<!-- generic output table row-->
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
