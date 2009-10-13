<?php
# -*- mode: C++; tab-width: 2; -*-
# vi: set ts=2:
#
# --------------------------------------------------------------------------
#                   OpenMS Mass Spectrometry Framework
# --------------------------------------------------------------------------
#  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
#
#  This library is free software; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 2.1 of the License, or (at your option) any later version.
#
#  This library is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public
#  License along with this library; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
# --------------------------------------------------------------------------
# $Maintainer: $
# $Authors: Marcel Grunert$
# --------------------------------------------------------------------------

	error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

	include "common_functions.php";

	// increase the memory limit of the particular PHP script
	ini_set('memory_limit', "128M");

	####################################################################################################################################
	########################################################## Usage ###################################################################
	####################################################################################################################################

	function printUsage()
	{
		print "\n";
		print "Evaluate FeatureFinder algorithm using OpenMS TOPP tools (FeatureFinder,\n";
		print "FeatureLinker, PeakPicker, TextExporter) and Gnuplot for data plotting. \n";
		print "Version: 1.0\n";
		print "\n";
		print "Usage:\n";
		print "  php evaluate_ff.php <options>";
		print "\n\n";
		print "Options (mandatory options marked with '*'):\n";
		print "  -in      <file>*  -- Input file (valid formats: 'mzData')\n";
		print "  -log     <file>*  -- Output file for evaluation\n";
		print "  -out_ff  <file>*  -- Output feature list (valid formats: 'featureXML')\n";
		print "  -type    <name>*  -- FeatureFinder algorithm type (default: 'simple')\n";
		print "  -ini_ff  <file>   -- Use the given FeatureFinder TOPP INI file\n";
		print "\n";
		print "Valid FeatureFinder algorithm types:\n";
		foreach ($GLOBALS["types"] as $name => $desc)
		{
				print "  $name  $desc\n";
		}
		print "\n";
		print "Common options:\n";
		print "  -in_an  <file>  -- Annotated feature list (valid formats: 'featureXML')\n";
		print "  -out_co <file>  -- Output consensus feature list (valid formats: 'consensusXML')\n";
    print "  -use_pp <n>     -- Use PeakPicker to find peaks in raw data before FeatureFinder (default: '0')\n";
		print "  -ini_pp <file>  -- Use the given PeakPicker TOPP INI file\n";
		print "  -ini_fl <file>  -- Use the given FeatureLinker TOPP INI file (mandatory if '-in_an' is given)\n";
		print "  -calibrate <n>  -- Parameter to calibrate the mass decimal values (default: '1', i.e. no calibration)\n";
    print "  -eval_ony <n>   -- Evaluation without feature finding (default: '0', i.e. feature finding before evaluation)\n";
		print "  --help          -- Shows this help\n";
		print "\n";
	}

	####################################################################################################################################
	##################################################### Declarations #################################################################
	####################################################################################################################################

	$GLOBALS["types"] = array(
		"simple" => "               -- Evaluate simple FeatureFinder (default)",
		"simplest" => "             -- Evaluate simplest FeatureFinder",
		"picked_peak"	=> "          -- Evaluate picked FeatureFinder",
		"isotope_wavelet"	=> "      -- Evaluate FeatureFinder with IsotopeWavelet",
		"isotope_wavelet_nofit"	=> "-- Evaluate FeatureFinder with IsotopeWavelet without model fitting"
	);
	
	$options = array("-in","-log","-ini_ff", "-out_ff", "-type", "-in_an", "-out_co", "-use_pp", "-ini_pp", "-ini_fl", "-delete", "-calibrate", "-eval_only", "--help");

	####################################################################################################################################
	################################################## Parameter handling ##############################################################
	####################################################################################################################################

	#wrong parameter count
	for ($i=1; $i<count($argv); ++$i)
	{
		#option
		if (beginsWith($argv[$i],"-"))
		{
			#not registered option or flag
			if (!in_array($argv[$i],$options))
			{
				print "\nError: Unregistered option '".$argv[$i]."'!\n\n";
				printUsage();
				exit;
			}
			
			# no argument to an option
			if ($argv[$i]!="--help" && in_array($argv[$i],$options) && ( !isset($argv[$i+1]) || beginsWith($argv[$i+1],"-")))
			{
				print "\nError: No argument to option '".$argv[$i]."'!\n\n";
				printUsage();
				exit;
			}			
		}
	}

	// no parameters
	if ($argc==1 || in_array("--help",$argv))
	{
		if (!in_array("--help",$argv))
		{
			print "\nNo options given. Aborting!\n\n";
		}
		printUsage();
		exit;
	}
  
  // eval_only - evaluation without feature finding
	$eval_only = 0;
	if (in_array("-eval_only",$argv))
  {
      $eval_only = $argv[array_search("-eval_only",$argv)+1];
      
      if ($eval_only!=1 && $eval_only!=0)
      {
          print "\nWrong input format for '-eval_only' (valid formats: '0' or '1') given. Aborting!\n\n";
          printUsage();
          exit;   
      }
  }

	// input file* (raw data)
	$in = "";
	if (in_array("-in",$argv))
	{
		$in = $argv[array_search("-in",$argv)+1];
		
		$end = explode(".",$in);
		if ($end[count($end)-1] != "mzData")
		{
			print "\nWrong input file format (valid formats: 'mzData') given. Aborting!\n\n";
			printUsage();
			exit;
		}
 	}
	else if ($eval_only == 0)
	{
		print "\nNo input file (-in) given. Aborting!\n\n";
		printUsage();
		exit;
	}

	// log file*
	$log = "";
	if (in_array("-log",$argv))
	{
		$log = $argv[array_search("-log",$argv)+1];
	}
	else
	{
		print "\nNo output file (-log) given. Aborting!\n\n";
		printUsage();
		exit;
	}

	// featureXML file (annotated feature list)
	$gold = "";
	if (in_array("-in_an",$argv))
	{
		$gold = $argv[array_search("-in_an",$argv)+1];
		
		$end = explode(".",$gold);
		if ($end[count($end)-1] != "featureXML")
		{
			print "\nWrong annotated feature list format (valid formats: 'featureXML') given. Aborting!\n\n";
			printUsage();
			exit;
		}
	}

	// algorithm type*
	$test = "simple";
	if (in_array("-type",$argv))
	{
		$test = $argv[array_search("-type",$argv)+1];
		if (!in_array($test, array_keys($GLOBALS["types"])))
		{
			print "\nError: Unknown FeatureFinder algorithm '$test'!\n\n";
			printUsage();
			exit;
		}
	}
	$types = $test;

	// out_ff* - FeatureFinder output feature list
	$out_ff = "";
	if (in_array("-out_ff",$argv))
	{
		$out_ff = $argv[array_search("-out_ff",$argv)+1];
		
		$end = explode(".",$out_ff);
		if ($end[count($end)-1] != "featureXML")
		{
			print "\nWrong feature list format (valid formats: 'featureXML') given. Aborting!\n\n";
			printUsage();
			exit;
		}
	}
	else
	{
		print "\nNo output feature list (-out_ff) given. Aborting!\n\n";
		printUsage();
		exit;
	}

	// ini_ff - FeatureFinder TOPP INI file
	$ini_ff = "";
	if (in_array("-ini_ff",$argv))
	{
		$ini_ff = $argv[array_search("-ini_ff",$argv)+1];
	}
  
  // ini_pp - use PeakPicker before FeatureFinder
	$use_pp = 0;
	if (in_array("-use_pp",$argv))
  {
      $use_pp = $argv[array_search("-use_pp",$argv)+1];
      
      if ($use_pp!=1 && $use_pp!=0)
      {
          print "\nWrong input format for '-use_pp' (valid formats: '0' or '1') given. Aborting!\n\n";
		      printUsage();
		      exit;   
      }
  }

	// ini_pp - PeakPicker TOPP INI file
	$ini_pp = "";
	if (in_array("-ini_pp",$argv))
	{
		$ini_pp = $argv[array_search("-ini_pp",$argv)+1];
	}

	// ini_fl - FeatureLinker TOPP INI file
	$ini_fl = "";
	if (in_array("-ini_fl",$argv))
	{
		$ini_fl = $argv[array_search("-ini_fl",$argv)+1];
	}
	else
	{
		if ($gold!="")
		{
			print "\nNo FeatureLinker TOPP INI file (-ini_fl) given (mandatory if '-in_an' is given). Aborting!\n\n";
			printUsage();
			exit;
		}	
	}

	// out_co - FeatureLindeer output consensus feature list
	$out_co = $out_ff.".consensusXML";
	if (in_array("-out_co",$argv))
	{
		$out_co = $argv[array_search("-out_co",$argv)+1];
		
		$end = explode(".",$out_co);
		if ($end[count($end)-1] != "consensusXML")
		{
			print "\nWrong output consensus feature list format (valid formats: 'consensusXML') given. Aborting!\n\n";
			printUsage();
			exit;
		}
	}

	// calibrate - calibrate the mass decimal values
	// theoretical λ = 1.000495 ... see "Analytical model of peptide mass cluster centres with applications"
	$calibrate = 1;
	if (in_array("-calibrate",$argv))
	{
		$calibrate = $argv[array_search("-calibrate",$argv)+1];
	}

	####################################################################################################################################
	########################### Run OpenMS TOPP Tools (FeatureFinder, FeatureLinker, TextExporter) #####################################
	####################################################################################################################################

	// FeatureFinder running time
	$runningTimeFF = 0;
  
	// run PeakPicker
	if ($types=="picked_peak" && $use_pp==1 && $eval_only==0)
	{
		echo "Running PeakPicker ... \n";
		if ($ini_pp!="") passthru("PeakPicker -ini ".$ini_pp." -in ".$in." -out ".$in."_pp.mzData");
		else passthru("PeakPicker -in ".$in." -out ".$in."_pp.mzData");
		$in = $in."_pp.mzData";
	}
	echo "\n";

  	// run FeatureFinder
	if ($eval_only==0)
  {
    echo "Running FeatureFinder ... \n";
    if ($ini_ff != "")
    {
      $starttime = time();
      passthru("FeatureFinder -ini ".$ini_ff." -in ".$in." -out ".$out_ff);
      $runningTimeFF = gmdate("H:i:s", (time() - $starttime));
    }
    else
    {
      $starttime = time();
      passthru("FeatureFinder -type ".$types." -in ".$in." -out ".$out_ff);
      $runningTimeFF = gmdate("H:i:s", (time() - $starttime));
    }
    echo "\n";
  }
    
	// run TextExporter for featureXML
	$out_ff_txt = $out_ff.".txt";
	echo "Running TextExporter for ".$out_ff." ... \n";
	passthru("TextExporter -in ".$out_ff." -out ".$out_ff_txt);

	// run FeatureLinker
	if ($gold!="")
	{
		echo "Running FeatureLinker ... \n";
		if ($ini_fl!="") passthru("FeatureLinker -ini ".$ini_fl." -type unlabeled -in ".$out_ff.",".$gold." -out ".$out_co);
		else passthru("FeatureLinker -type unlabeled -in ".$out_ff.",".$gold." -out ".$out_co);

		// run TextExporter for annotated feature list
		$gold_txt = $out_co.".txt";
		echo "Running TextExporter for ".$gold." ... \n";
		passthru("TextExporter -in ".$gold." -out ".$gold_txt);
	}

	####################################################################################################################################
	################################################## Helper functions  ###############################################################
	####################################################################################################################################

	// Associates the specified value with the specified key in the hash map using for masses.
	function mass_hashmap_put($key, $value)
	{
		global $masses;
		$previous = 0;
		$newValue = 0;

		if (array_key_exists($key, $masses))
		{
			$previous = $masses[$key];
			$newValue = $previous + $value;
		}
		else
		{
			$previous = $value;
			$newValue = $previous;
		}

		$masses[$key] =& $newValue;

		return $previous;
	}

	// Associates the specified value with the specified key in the hash map using for charges.
	function charge_hashmap_put($key, $value)
	{
		global $charges;
		$previous = 0;
		$newValue = 0;

		if (array_key_exists($key, $charges))
		{
			$previous = $charges[$key];
			$newValue = $previous + $value;
		}
		else
		{
			$previous = $value;
			$newValue = $previous;
		}

		$charges[$key] =& $newValue;

		return $previous;
	}

	// Associates the specified value with the specified key in the hash map using for feature pairs.
	function pair_hashmap_put($key, $value)
	{
		global $pairs;
		$previous = 0;
		$newValue = 0;

		if (array_key_exists($key, $pairs))
		{
			$previous = $pairs[$key];
			$newValue = $previous + $value;
		}
		else
		{
			$previous = $value;
			$newValue = $previous;
		}

		$pairs[$key] =& $newValue;

		return $previous;
	}

	// Function to call when the parser encounters a start element.
	function startElement($parser, $name, $attrs)
	{
		global $consensusID, $pair_elements;

		switch(strtolower($name))
		{
	    	case "consensuselement":
	    	{
	    		$consensusID = $attrs['ID'];
	    		break;
	    	}
	    	case "element":
	    	{
				$elem = $consensusID." ".$attrs['MAP']." ".$attrs['ID']." ".$attrs['RT']." ".$attrs['MZ']." ".$attrs['IT'];
				array_push($pair_elements,$elem);
	    		pair_hashmap_put($consensusID,1);
	    		break;
	    	}
		}
	}

	// Function to call when the parser encounters an end element.
	function endElement($parser, $name)
	{
	}

	// XML file parsing
	function parseXmlFile($filename)
	{
		// Defining the XML Parser
		$xml_parser = xml_parser_create();
		// Defining the Element Handlers
		xml_set_element_handler($xml_parser, "startElement", "endElement");

		// Starting the parser
		$file = file($filename);
		foreach($file as $elem)
		{
			xml_parse($xml_parser, $elem);
		}

		// Cleaning up
		xml_parser_free($xml_parser);
	}

	// Charge searching for RT/MZ value
	function getChargeFromFile($lines, $rt, $mz)
	{
		$ch = -1;
		foreach ($lines as $key => $value)
		{
			if ($key > 0)
			{
				$sp = split(" ",$value);
				if ($sp[0]==$rt AND $sp[1]==$mz) $ch = $sp[3];
			}
		}

		return $ch;
	}

	####################################################################################################################################
	################################################# FF evaluation ####################################################################
	####################################################################################################################################

	########################################################################################
	######################### Evaluation for feature list ##################################
	########################################################################################

	// open log file
	$logFile = fopen($log,"w") or die ("Could not open file ".$log);
	
	// some vaiables
 	$charges = array();
 	$mass_keys = array();
	$mass_values = array();
	$max = 0;
	$start = TRUE;
	$start2 = TRUE;
	$min = 0;
	$count = 0;
	$sum = 0;
	$minIntens = 0;
	$sumIntens = 0;
	$maxIntens = 0;

	// open file for mass decimal
	$mdFile = fopen(($out_ff.".massDecimal.txt"),"w") or die ("Could not open file ".$out_ff.".massDecimal.txt");
	fwrite($mdFile, "# charge charge_mass_decimal\n");

	// evaluate feature list
	echo "Read and evaluate feature list (".$out_ff_txt.") ... \n";
	$lines_ff = file($out_ff_txt);
	foreach ($lines_ff as $key => $value)
	{
		if ($key > 0 && trim($key)!="")
		{
			// split ... rt, mz, intensity, charge, overall_quality, rt_quality, mz_quality, rt_start, rt_end
			$sp = split(" ",$value);

			// mean overall_quality
			$sum += $sp[4];

			// mean intensity
			$sumIntens += $sp[2];

			// minimum and maximum mz_quality
			if ($max <= $sp[4] ) $max = $sp[4];
			if ($maxIntens <= $sp[2] ) $maxIntens = $sp[2];
			if ($count==0)
			{
				$min = $sp[4];
				$minIntens = $sp[2];
			}
			if ($min >= $sp[4] ) $min = $sp[4];
			if ($minIntens >= $sp[2] ) $minIntens = $sp[2];

			// monoisotopic mass = (m/z of monoisotopic peak)*(charge) - charge 
			$mass = ($sp[1]*$sp[3])-$sp[3];
			if ($calibrate != 1) $mass = $mass * $calibrate;
			// mass decimal
			$decimal_mass = $mass - (int)$mass;
			
			array_push($mass_keys,$mass);
			array_push($mass_values,number_format($decimal_mass,3));
			$md = number_format($decimal_mass,2);

			// charge
			charge_hashmap_put($sp[3],1);
			fwrite($mdFile, $sp[3]." ".$md."\n");

			// number of features
			$count++;
		}
	}
	fclose($mdFile); // close file for mass decimal

	// number of feature pairs
	$featurePairs = 0;

	########################################################################################
	######## Evaluation for annotated feature list 	in comparison with FF results ##########
	########################################################################################

	// parse annotated feature list and compare with feature list from FeatureFinder
	if ($gold!="")
	{
		// feature pairs (hashmap)
		$pairs = array();

		// consensus elements
		$pair_elements = array();

		// feature attributes
		$intensity_ff = array();
		$intensity_an = array();
		$rt_ff = array();
		$rt_an = array();
		$mz_ff = array();
		$mz_an = array();
		$charge_ff = array();
		$charge_an = array();
		$id_an = array();
		
		// read file (annotated feature list) into array, i.e. line by line 
		$lines_an = file($gold_txt);
		
		// number of features in annotated list
		$featureInAnnotated = 0;
		foreach ($lines_an as $key => $value)
		{
			if ($key > 0) $featureInAnnotated++;
		}

		// parse FeatureLinker output file
		echo "Parse consensus XML file (".$out_co.") ...\n";
		parseXmlFile($out_co);

		// count feature pairs
		$keys = array_keys($pairs);
		for ($i=0; $i<count($pairs); $i++)
		{
			if (($pairs[$keys[$i]]) >= 2) $featurePairs++;
		}

		if ($featurePairs != 0)
		{
			// for all feature pairs
			for ($i=0; $i<count($pair_elements); $i++)
			{
				// split ... sp[0] -> consensusID, sp[1] -> mapID, sp[2] -> id, sp[3] -> rt, sp[4] -> mz, sp[5] -> intensity
				$sp = split(" ", $pair_elements[$i]);

				// if we have a feature pair
				if ($pairs[$sp[0]]>1)
				{
					if ($sp[1]==0)
					{
						array_push($rt_ff, $sp[3]);
						array_push($mz_ff, $sp[4]);
						array_push($intensity_ff, $sp[5]);
						array_push($charge_ff, getChargeFromFile($lines_ff, $sp[3], $sp[4]) );
					}
					else
					{
						array_push($id_an, $sp[2]);
						array_push($rt_an, $sp[3]);
						array_push($mz_an, $sp[4]);
						array_push($intensity_an, $sp[5]);
						array_push($charge_an, getChargeFromFile($lines_an, $sp[3], $sp[4]) );
					}
				}
			}

			// write some statistics to text files
			$intensityFile = fopen(($out_ff.".intensity_distribution.txt"),"w") or die ("Could not open file ".$out_ff.".intensity_distribution.txt");
			$accuracyFile = fopen(($out_ff.".accuracy.txt"),"w") or die ("Could not open file ".$out_ff.".accuracy.txt");
			$chargePredictionFile = fopen(($out_ff.".chargePrediction.txt"),"w") or die ("Could not open file ".$out_ff.".chargePrediction.txt");
			$featureIdFile = fopen(($out_ff.".annotated_pair_ids.txt"),"w") or die ("Could not open file ".$out_ff.".annotated_pair_ids.txt");
																			
			fwrite($intensityFile, "# log2(Intensity_1/Intensity_2) 1/2*log2(Intensity_1*Intensity_2)\n" );
			fwrite($accuracyFile, "# feature_pair mz_abs_error mz_rel_error(%) rt_abs_error rt_rel_error(%)\n" );
			fwrite($chargePredictionFile, "# feature_pair charge_FF charge_annotated\n" );
			fwrite($featureIdFile, "# feature_id_from_annotated_featureXML\n" );
			
			$mean_mz_deviance = 0;
			$mean_mz_rel_error = 0;
			$mean_rt_deviance = 0;
			$mean_rt_rel_error = 0;
			$trueCharges = 0;
			
			for ($i=0; $i<count($charge_ff); $i++)
			{
				// annotated pair IDs
				fwrite($featureIdFile, ($id_an[$i])."\n");
			
				// intensity
				$x_axis = log( $intensity_ff[$i] * $intensity_an[$i] )/log(2) * 0.5;
				$y_axis = log( $intensity_ff[$i] / $intensity_an[$i] )/log(2);
				fwrite($intensityFile, ($x_axis."	".$y_axis)."\n" );

				// mz/rt deviance and relative error
				$mean_mz_deviance += number_format(abs($mz_ff[$i]-$mz_an[$i]),3);
				$mean_mz_rel_error += number_format((abs($mz_ff[$i]-$mz_an[$i])/$mz_an[$i])*100,3);
				$mean_rt_deviance += number_format(abs($rt_ff[$i]-$rt_an[$i]),3);
				$mean_rt_rel_error += number_format((abs($rt_ff[$i]-$rt_an[$i])/$rt_an[$i])*100,3);
														
				$table = ($i+1);
				$table .= "	".number_format(abs($mz_ff[$i]-$mz_an[$i]),3);
				$table .= "	".number_format((abs($mz_ff[$i]-$mz_an[$i])/$mz_an[$i])*100,3);
				$table .= "	".number_format(abs($rt_ff[$i]-$rt_an[$i]),3);
				$table .= "	".number_format((abs($rt_ff[$i]-$rt_an[$i])/$rt_an[$i])*100,3);
				$table .= "\n";
				fwrite($accuracyFile, ($table));
				
				// charge prediction
				fwrite($chargePredictionFile, ($i+1)."	".($charge_ff[$i]."	".$charge_an[$i]."\n"));	
				if ($charge_ff[$i]==$charge_an[$i]) $trueCharges++;
			}
			fclose($intensityFile);
			fclose($accuracyFile);
			fclose($chargePredictionFile);
			fclose($featureIdFile);
			
			// MA plot of feature intensities
			$gnuplotScript = fopen("gnuplot_script","w") or die ("Could not open file gnuplot_script");
			fwrite($gnuplotScript, "set terminal png\n");
			fwrite($gnuplotScript, "set output \"".$out_ff.".intensity_distribution.png\"\n");
			fwrite($gnuplotScript, "unset key\n");
			fwrite($gnuplotScript, "set title 'Feature intensities (".$types." FeatureFinder)'\n");
			fwrite($gnuplotScript, "set xlabel 'log2(I1/I2)'\n");
			fwrite($gnuplotScript, "set ylabel '1/2*log2(I1*I2)'\n");
			fwrite($gnuplotScript, "plot '".$out_ff.".intensity_distribution.txt' ps .5");
			fclose($gnuplotScript);
			passthru("gnuplot ./gnuplot_script");
			passthru("rm -f gnuplot_script");
		
		} // end ... for all feature pairs

	} // end ... parse annotated feature list

	########################################################################################
	############################## Write and plot results ##################################
	########################################################################################

	// output
	fwrite($logFile, "FeatureFinder type: ".$types."\n\n");

	// some output files
	$outCharges = $out_ff.".charges_states.txt";
	$outChargePNG = $out_ff.".charges_states.png";
	$outMz = $out_ff.".mz.txt";
 	$outMzPNG = $out_ff.".mz.png";
	$outMD = $out_ff.".massDecimal.txt";

	// number of features
	fwrite($logFile, "Number of features: ".$count."\n");

	// feature pairs
 	if ($gold!="") 
	{
		// number of feature pairs
		fwrite($logFile, "Feature pairs: ".$featurePairs."/".$featureInAnnotated."\n");
		
		// absolute and realtive errror of mz and rt
		if ($featurePairs!=0)
		{
			fwrite($logFile, "  mean absolute error of mz: ".number_format($mean_mz_deviance/$featurePairs,3)."\n");
			fwrite($logFile, "  mean relative error of mz: ".number_format($mean_mz_rel_error/$featurePairs,3)."\n");
			fwrite($logFile, "  mean absolute error of rt: ".number_format($mean_rt_deviance/$featurePairs,3)."\n");
			fwrite($logFile, "  mean relative error of rt: ".number_format($mean_rt_rel_error/$featurePairs,3)."\n");
			fwrite($logFile, "  true charge prediction: ".$trueCharges."/".$featurePairs."\n");	
		}
	}

	// running time
  if ($eval_only==0) fwrite($logFile, "Running time (hh:mm:ss): ".$runningTimeFF."\n");

	// correlation
	fwrite($logFile, "Correlation:\n");
	fwrite($logFile, "  minimum: ".number_format($min,5)."\n");
	if ($count != 0) fwrite($logFile, "     mean: ".number_format($sum/$count,5)."\n");
	else fwrite($logFile, "     mean: 0\n");
	fwrite($logFile, "  maximum: ".number_format($max,5)."\n");

	// intensity
	fwrite($logFile, "Intensity:\n");
	fwrite($logFile, "  minimum: ".$minIntens."\n");
	if ($count != 0) fwrite($logFile, "     mean: ".number_format($sumIntens/$count,1)."\n");
	else fwrite($logFile, "     mean: 0\n");
	fwrite($logFile, "  maximum: ".$maxIntens."\n");

	// generate scatterplot for mass decimals (output: file and PNG (generated by Gnuplot))
	{
		// open output file for mass decimal
		$massDecimalFile = fopen($outMz,"w") or die ("Could not open file ".$outMz);
		fwrite($massDecimalFile, "# mass_decimal mass\n");
		$keys = $mass_keys;
		for ($i=0; $i<count($keys); $i++)
		{
			fwrite($massDecimalFile, ($mass_values[$i])." ".$keys[$i]."\n");
		}
		fclose($massDecimalFile);

		// run Gnuplot for scatterplot
		$gnuplotScript = fopen("gnuplot_script","w") or die ("Could not open file gnuplot_script");
		fwrite($gnuplotScript, "set terminal png\n");
		fwrite($gnuplotScript, "set output \"".$outMzPNG."\"\n");
		fwrite($gnuplotScript, "unset key\n");
		fwrite($gnuplotScript, "set title 'Mass distribution (".$types." FeatureFinder)'\n");
		fwrite($gnuplotScript, "set xlabel 'Mass decimal (Da)'\n");
		fwrite($gnuplotScript, "set ylabel 'Mass (Da)'\n");
		fwrite($gnuplotScript, "plot '".$outMz."' ps 0.5");
		fclose($gnuplotScript);
		passthru("gnuplot ./gnuplot_script");
		passthru("rm -f gnuplot_script");
	}

	// charge distribution
	{
		fwrite($logFile, "Charges: \n");

		// sort charges
		ksort($charges);

		// open output file for charges
		$chargeDistributionFile = fopen($outCharges,"w") or die ("Could not open file ".$outCharges);
		fwrite($chargeDistributionFile,"# charge charge_frequency\n");
		
		// ... needed by gnuplot axial ranges
		$maxCharge = 0;
		$maxChargeNumber = 0;

		$keys = array_keys($charges);
		for ($i=0; $i<count($charges); $i++)
		{
			$num = $charges[$keys[$i]];
			if ($count != 0) $num2 = $num*100/$count;
			else $num2 = 0;
			fwrite($chargeDistributionFile, $keys[$i]."	".$num."\n");
			fwrite($logFile, "         +".$keys[$i].": ".number_format($num2,0)."% (".$num.")\n");

			if ($keys[$i] > $maxCharge) $maxCharge = $keys[$i];
			if ($num > $maxChargeNumber) $maxChargeNumber = $num;
		}

		// close output file for charges
		fclose($chargeDistributionFile);

		// increase for gnuplot axial ranges
		$maxCharge = $maxCharge +1;
		$maxChargeNumber = $maxChargeNumber +5;

		// generate histogram for charge distribution using Gnuplot
		$gnuplotScript = fopen("gnuplot_script","w") or die ("Could not open file gnuplot_script");
		fwrite($gnuplotScript, "set terminal png\n");
		fwrite($gnuplotScript, "set output \"".$outChargePNG."\"\n");
		fwrite($gnuplotScript, "unset key\n");
		fwrite($gnuplotScript, "set title 'Charge distribution (".$types." FeatureFinder)'\n");
		fwrite($gnuplotScript, "set xlabel 'Charge states'\n");
		fwrite($gnuplotScript, "set ylabel 'Frequency'\n");
		fwrite($gnuplotScript, "set xrange[0:".$maxCharge."]\n");
		fwrite($gnuplotScript, "set xtics 0,1,".$maxCharge."\n");
		fwrite($gnuplotScript, "set yrange [0:".$maxChargeNumber."]\n");
		fwrite($gnuplotScript, "set style fill solid\n");
		fwrite($gnuplotScript, "set style data boxes\n");
		fwrite($gnuplotScript, "set boxwidth 0.2\n");
		fwrite($gnuplotScript, "plot '".$outCharges."'\n");
		fclose($gnuplotScript);
		passthru("gnuplot ./gnuplot_script");
		passthru("rm -f gnuplot_script");
  	}

	// write mass decimal for specific charge
	for ($key=0; $key<count($keys); $key++)
	{
		$mdFileForSpecficCharge = fopen($outMD.".charge".$keys[$key].".txt","w") or die ("Could not open file ".$outMD.".charge".$keys[$key].".txt");
		fwrite($mdFileForSpecficCharge, "# mass_decimal frequency\n");
		$md = 0;
		$masses = array();

		$file = file($outMD);
		foreach ($file AS $line)
		{
			list($c, $md) = explode(' ', $line);
			if ($keys[$key] == $c) mass_hashmap_put(trim($md),1);
		}

    	$maxFrequency = 0;
		$mass_keys = array_keys($masses);
		for ($i=0; $i<count($masses); $i++)
		{
			if ($masses[$mass_keys[$i]] > $maxFrequency) $maxFrequency = $masses[$mass_keys[$i]];
			fwrite($mdFileForSpecficCharge, $mass_keys[$i]."  ".$masses[$mass_keys[$i]]."\n");
		}
		fclose($mdFileForSpecficCharge);

		// generate histogram for specific charge distribution using Gnuplot
		$gnuplotScript = fopen("gnuplot_script","w") or die ("Could not open file gnuplot_script");
		fwrite($gnuplotScript, "set terminal png\n");
		fwrite($gnuplotScript, "set output \"".$outMD.".charge".$keys[$key].".png"."\"\n");
		fwrite($gnuplotScript, "unset key\n");
		fwrite($gnuplotScript, "set title 'Mass distribution of charge +".$keys[$key]."'\n");
		fwrite($gnuplotScript, "set xlabel 'Mass decimal'\n");
		fwrite($gnuplotScript, "set ylabel 'Frequency'\n");
		fwrite($gnuplotScript, "set yrange[0:".($maxFrequency + $maxFrequency*0.1)."]\n");
		fwrite($gnuplotScript, "set xrange[0:1]\n");
		fwrite($gnuplotScript, "set style fill solid\n");
		fwrite($gnuplotScript, "set style data boxes\n");
		fwrite($gnuplotScript, "set boxwidth 0.01\n");
		fwrite($gnuplotScript, "plot '".$outMD.".charge".$keys[$key].".txt"."'\n");
		fclose($gnuplotScript);
		passthru("gnuplot ./gnuplot_script");
		passthru("rm -f gnuplot_script");
	}

	// debug output for log file
	fwrite($logFile, "\n");
	fwrite($logFile, "Input and output files ...\n");
	fwrite($logFile, "  FeatureFinder input file:       ".$in."\n");
	fwrite($logFile, "  FeatureFinder output file:      ".$out_ff."\n");
	fwrite($logFile, "    TextExporter output file:     ".$out_ff_txt."\n");
	if ($out_co!="") fwrite($logFile, "  Output consensus feature list:  ".$out_co."\n");
	if ($ini_pp!="") fwrite($logFile, "  PeakPicker TOPP INI file:       ".$ini_pp."\n");
	
	if ($gold!="")
	{
		fwrite($logFile, "  Annotated feature list:         ".$gold."\n");
		fwrite($logFile, "    TextExporter output file:     ".$gold_txt."\n");
		fwrite($logFile, "  FeatureLinker TOPP INI file:    ".$ini_fl."\n");
	}

	fwrite($logFile, "\n");
  fwrite($logFile, "  Intensity distribution                         ".$out_ff.".intensity_distribution.txt\n");
  fwrite($logFile, "                                                 ".$out_ff.".intensity_distribution.png\n");
	fwrite($logFile, "  Charge distribution:                           ".$outCharges."\n");
	fwrite($logFile, "                                                 ".$outChargePNG."\n");
	fwrite($logFile, "  Masses and mass decimal values:                ".$outMz."\n");
	fwrite($logFile, "                                                 ".$outMzPNG."\n");
	fwrite($logFile, "  Mass decimal values of all charges:            ".$outMD."\n");
	fwrite($logFile, "  Mass decimal distribution of specific charge:  ".$outMD.".charge<*>.txt\n");
	fwrite($logFile, "                                                 ".$outMD.".charge<*>.png\n");
	
	if ($gold!="")
	{
		fwrite($logFile, "  Absolute and relative error of mz and rt:      ".$out_ff.".accuracy.txt\n");
		fwrite($logFile, "  Charge prediction:                             ".$out_ff.".chargePrediction.txt\n");
		fwrite($logFile, "  Matched feature IDs from annotated featureXML: ".$out_ff.".annotated_pair_ids.txt\n");
	}
	
	// close output file
	fclose($logFile);
	echo "Write results to log file (".$log.").\n\n";

  	// command line output
	$handle = fopen ($log, "r") or die ("Could not open file ".$log);
	while ( $content = fgets ($handle, 4096 ))
	{
		echo "$content";
	}
	fclose($handle);

  	####################################################################################################################################
	######################################################## Cleaning up ###############################################################
	####################################################################################################################################

	// delete log files from file system
	passthru("rm -f TOPP.log");
	passthru("rm -f featurefinder.log");

?>
