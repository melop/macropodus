<?php
$sRunDir = "rundir";
$bSortScfAsNum = true; // sort scffold name as number

$sScfPrefix = ""; //to be removed from scaffold names, leaving only numbers
$sSp = "Anabas_testudineus";
$sVersion = "104";
$sLongestIsoformProt = "/data/projects/rcui/macropodus_compare/ensembl/anabas/longest_isoform.prot.fa";

/*
$sScfPrefix = ""; //to be removed from scaffold names, leaving only numbers
$sSp = "Betta_splendens";
$sVersion = "104";
$sLongestIsoformProt = "/data/projects/rcui/macropodus_compare/ensembl/betta/longest_isoform.prot.fa";
*/
/*
$bSortScfAsNum = false;
$sScfPrefix = ""; //to be removed from scaffold names, leaving only numbers
$sSp = "Channa_argus";
$sVersion = "104";
$sLongestIsoformProt = "/data/projects/rcui/macropodus_compare/ensembl/channa_argus/longest_isoform.prot.fa";
*/
/*
$sScfPrefix = "mhkscf_"; //to be removed from scaffold names, leaving only numbers
$sSp = "Macropodus_hongkongensis";
$sVersion = "1.0";
$sLongestIsoformProt = "/data/projects/rcui/macropodus_compare/improve_genemodels_UPhO/mhk/longest.addedUPhO.genesymbol.spgeneid.prot.fa";
*/
/*
$sScfPrefix = "mopscf_"; //to be removed from scaffold names, leaving only numbers
$sSp = "Macropodus_opercularis";
$sVersion = "1.0";
$sLongestIsoformProt = "/data/projects/rcui/macropodus_compare/improve_genemodels_UPhO/mop/longest.addedUPhO.genesymbol.spgeneid.prot.fa";
*/
$bSortScfAsNum = false;
$sScfPrefix = ""; //to be removed from scaffold names, leaving only numbers
$sSp = "Gasterosteus_aculeatus";
$sVersion = "104";
$sLongestIsoformProt = "/data2/projects/zwang/macropodus_compare/ensembl/stickleback/longest_isoform.prot.fa";

$bSortScfAsNum = false;
$sScfPrefix = ""; //to be removed from scaffold names, leaving only numbers
$sSp = "Oreochromis_niloticus";
$sVersion = "104";
$sLongestIsoformProt = "/data2/projects/zwang/macropodus_compare/ensembl/tilapia/longest_isoform.prot.fa";


$sWD = "$sRunDir/rawGenomes/$sSp/$sVersion/annotation";
exec("mkdir -p $sWD");

$hGFF = false;
if ($bSortScfAsNum) {
	$hGFF = popen("sort -k1,1n -k4,4n | gzip -c > $sWD/$sSp.gene.gff.gz", 'w');
} else {
	$hGFF = popen("sort -k1,1 -k4,4n | gzip -c > $sWD/$sSp.gene.gff.gz", 'w');
}

$hFas = popen("gzip -c > $sWD/$sSp.pep.fa.gz", 'w');

$hIn = popen("zcat -f $sLongestIsoformProt", 'r');

while(false!== ($sLn = fgets($hIn) )) {
	$sLn = trim($sLn);
	if ($sLn == '') {
		continue;
	}

	if ($sLn[0] == '>') {
		$arrF = explode(' ', substr($sLn, 1));
		$sSeqName = $arrF[0];
		$sSeqName = preg_replace('/[:|]/', '_', $sSeqName);
		$sPos = $arrF[count($arrF)-1];
		preg_match('/([^:]+):([^-]+)-([^(]+)\\((\\S)\\)/', $sPos, $arrM);
		if (count($arrM) != 5) {
			die("Failed to parse header: $sLn\n");
		}

		list($sChr, $nStart, $nEnd, $sOrient) = array_slice($arrM, 1);
		$sChr = str_replace($sScfPrefix, '', $sChr);
		fwrite($hGFF, "$sChr\tfunannotate\tgene\t$nStart\t$nEnd\t\t$sOrient\t\tlocus=$sSeqName\n");
		fwrite($hFas, ">$sSeqName\n");
	} else {
		$sLn = str_replace('*', '', $sLn);
		fwrite($hFas, $sLn."\n");
	}

	
}

?>
