<?php
//assign gene symbol
$sOrthoGroup = "/data/projects/rcui//macropodus_compare/UPhO/assigngenesymbol/UPhO_orthogroups.genesymbol.txt";
$sGFF = "mhk.full.gff3";
$sOut = "mhk.full.genesymbol.gff3";
$sGFF = "mhk.longest_isoform.gff3";
$sOut = "mhk.longest_isoform.genesymbol.gff3";
$sThisSp = 'MHK';
$sNum2ID = "/data/projects/rcui/macropodus_compare/UPhO/allmacropodusproteins.fa";
$sPredictor = "maker"; //replace column 2 with this value

$hOut = fopen($sOut , 'w');

echo("$sNum2ID \n");

$arrNum2ID  = fnLoadNum2ID($sNum2ID );
$arrID2Num = array_flip($arrNum2ID);

//echo($arrID2Num['NORv12scf1-frag-merged-region2-clust1']);

//print_r($arrID2Num);

//die();

$arrGeneSymbolCounts = array();

$arrGene2RNAMap = array();
$hGFF = fopen($sGFF, 'r');

$arrAvailableRNAIDs = array();
$arrRNAIDsWithSymbols = array();
$arrRNAIDsWithSymbols2 = array();

	while( false !== ($sLn = fgets($hGFF ) ) ) {
		if ($sLn == '#') continue;
		$sLn = trim($sLn);
		if ($sLn == '') continue;

		$arrF = explode("\t" , $sLn);

		if (count($arrF) != 9) continue;
		
		if ($arrF[2] == 'mRNA') {
			$arrAnnot = fnParseFields($arrF[8]);
			if (!array_key_exists("ID" , $arrAnnot) ) continue;
			if (!array_key_exists("Parent" , $arrAnnot) ) continue;
			$arrGene2RNAMap[$arrAnnot['Parent'] ] = $arrAnnot['ID'];
			$sRNAID =str_replace(".", "_", $arrAnnot['ID']);
			$arrAvailableRNAIDs[$sRNAID] = true;

			continue;
		}

	}
fclose($hGFF);

$arrID2Symbol = fnLoadOrthoGroups($sOrthoGroup);

//print_r($arrID2Symbol);


$hGFF = fopen($sGFF, 'r');

	while( false !== ($sLn = fgets($hGFF ) ) ) {
		if ($sLn == '#') {
			fwrite( $hOut , $sLn);
			continue;
		}

		$sLn = trim($sLn);
		if ($sLn == '') continue;

		$arrF = explode("\t" , $sLn);
		$arrF[1] = $sPredictor; //change predictor

		if (count($arrF) != 9) {
			fwrite( $hOut , $sLn."\n");
			continue;
		}
		
		if ($arrF[2] == 'gene' || $arrF[2] == 'mRNA') {
			$arrAnnot = fnParseFields($arrF[8]);
			if (!array_key_exists("ID" , $arrAnnot) ) continue;
			if ($arrF[2] == 'gene') {
				if (!array_key_exists($arrAnnot['ID'] , $arrGene2RNAMap) ) {
					fwrite( $hOut , implode("\t", $arrF)."\n");
					continue;
				}
			}
			$sRNAID = ($arrF[2] == 'gene')? $arrGene2RNAMap[$arrAnnot['ID']] : $arrAnnot['ID'];
			$sRNAID =str_replace(".", "_", $sRNAID);

			if (array_key_exists($sRNAID, $arrID2Num)) {
				
			
			$arrAnnot['spgeneid'] = $arrID2Num[$sRNAID];

			if (!array_key_exists($sRNAID , $arrID2Symbol) ) {
				$arrAnnot['gene'] = 'unknown';
				$arrAnnot['description'] = 'unknown';
			} else {
				$arrAnnot['gene'] = $arrID2Symbol[$sRNAID]['genesymbol'];
				$arrAnnot['description'] = $arrID2Symbol[$sRNAID]['description'];
				$arrAnnot['ensemblorthologs'] = $arrID2Symbol[$sRNAID]['ensembl_orthologs'];


			}
			}

			$arrF[8] = fnArr2Annot($arrAnnot);
			fwrite( $hOut , implode("\t", $arrF)."\n");
			continue;
		} else {
			fwrite( $hOut , implode("\t", $arrF)."\n");
			continue;
		}



	}
fclose($hGFF);



function fnLoadOrthoGroups($sOrthoGroup) {
	global $sThisSp,  $arrNum2ID , $arrGeneSymbolCounts, $arrAvailableRNAIDs, $arrRNAIDsWithSymbols, $arrRNAIDsWithSymbols2;

	$hF = fopen($sOrthoGroup , 'r');
	$arrRet = array(); //key is gene id
	$arrGeneSymbolCounts = array();
	while( false !== ($sLn = fgets($hF) )) {
		$sLn = trim($sLn);
		if ($sLn == '') continue;
		list($sGroupID, $sGeneSymbol, $sDescription, $sGenes) = explode("\t", $sLn); 
		$arrGenes = explode(',' , $sGenes);

		$arrEnsemblOrthologs = array();
		$arrAddedGeneIDs = array();
		foreach( $arrGenes as $sGene) {
			if (substr($sGene,0,7) == '#trees_') continue;
			list($sSp, $sNum, $sGeneID) = explode('|' , $sGene);
			//echo("geneid:".$sGeneID."\n");
 

			if (substr($sGeneID , 0, 3) == 'ENS' ) {
				$arrEnsemblOrthologs[] = preg_replace("/_\d+/" , '', $sGeneID);
			}



			if ($sSp == $sThisSp) {

				if (array_key_exists($sNum, $arrNum2ID )) {
					$sGeneID = $arrNum2ID[$sNum];
					//echo("$sGeneID\n");
				} 

				$sGeneID =str_replace(".", "_", $sGeneID);

				$arrAddedGeneIDs[] = $sGeneID;
				$arrRet[$sGeneID] = array();
				$arrRet[$sGeneID]['internalspeciesgeneid'] = $sNum;
				$arrRet[$sGeneID]['genesymbol'] = $sGeneSymbol;
				$arrRet[$sGeneID]['description'] = $sDescription;
				if (!array_key_exists($sGeneSymbol , $arrGeneSymbolCounts)) {
					$arrGeneSymbolCounts[$sGeneSymbol] = 0;
				}
				if (array_key_exists($sGeneID , $arrAvailableRNAIDs) && (!array_key_exists($sGeneID, $arrRNAIDsWithSymbols) )) {
					$arrGeneSymbolCounts[$sGeneSymbol] += 1;
					$arrRNAIDsWithSymbols[$sGeneID] = true;
				}
			}
		}


		foreach($arrAddedGeneIDs as $sGeneID) {
			$arrRet[$sGeneID]['ensembl_orthologs'] = implode(',' , $arrEnsemblOrthologs);
		}
	}

	$arrGeneSymbolCounter = array();
	foreach($arrRet as $sGeneID => &$arrInfo) {
		if ($arrGeneSymbolCounts[$arrInfo['genesymbol']] > 1) { //
			if (!array_key_exists($arrInfo['genesymbol'] , $arrGeneSymbolCounter) ) {
				$arrGeneSymbolCounter[$arrInfo['genesymbol']] = 0;
			}
			if (array_key_exists($sGeneID , $arrAvailableRNAIDs) && (!array_key_exists($sGeneID, $arrRNAIDsWithSymbols2))) {
				$arrGeneSymbolCounter[$arrInfo['genesymbol']]  += 1;
				$arrRNAIDsWithSymbols2[$sGeneID] = true;
			}
			$arrInfo['genesymbol'] = $arrInfo['genesymbol']." (".$arrGeneSymbolCounter[$arrInfo['genesymbol']]." of ".$arrGeneSymbolCounts[$arrInfo['genesymbol']].")";

		}
	}

	return $arrRet;
}


function fnParseFields($s) {
	$arrF1 = explode(";" , $s);
	$arrRet = array();
	foreach($arrF1 as $sF) {
		$sF = trim($sF);
		if ($sF == '') continue;
		$arrF2 = explode("=" , $sF);
		if (count($arrF2)!=2) {
			echo("parse error: $sF\n");
		}
		$arrRet[trim($arrF2[0]) ] = trim($arrF2[1]);
	}

	return $arrRet;
}

function fnArr2Annot($arr) {
	$s = "";
	
	foreach($arr as $sKey => $sVal) {
		if ($sVal === false) $sVal=0;
		$s .= "$sKey=$sVal;";
	}

	return substr($s, 0, strlen($s)-1);
}

function fnLoadNum2ID($sNum2ID ) {
	global $sThisSp;
	$hF = popen("grep '>".$sThisSp."|' $sNum2ID" , 'r');
	$arrRet = array();
	while(false !== ($sLn = fgets($hF) )) {
		$arrF = explode('|', trim($sLn) );
		if (count($arrF)!=3) continue;
		$arrRet[$arrF[1]] = $arrF[2];
	}

	return $arrRet;
}


?>
