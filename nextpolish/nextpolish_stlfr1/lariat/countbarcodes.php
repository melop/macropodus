<?php
$sRd1Files = "/data/projects/shareddata/mop.stLFR/stLFR/raw_data_2021-07-16/genome/NGS/muscle/128.3G/D2104505A.clean_1.fq.gz"; //separate each file by a space, multiple files are allowed. for example: "test_1_1.fq.gz test_2_1.fq.gz"


$nLn = -1;            // current line number of the input stream.
$arrBarcodes = array();
$nKeptPairs = 0;      // counter for kept pairs
$nDiscardPairs = 0;   // counter for discarded pairs.

$h1 = popen("zcat -f $sRd1Files", 'r');

while(true) {
	$nLn++;										// increase line number counter, line 1 = 0
	$sLn1 = fgets($h1);							// get the next line from the read 1 stream

	if ($sLn1 === false) break;

	$sLn1 = trim($sLn1); // remove trailing "\n"

	if ($sLn1 == '') continue;

	if ($nLn % 4 == 0) { // title line
		//check title lines:
		if ($sLn1[0] != '@' ) { //check if title lines start with @
			die("Error: unexpected title lines on ln $nLn\n$sLn1\n");			
		}
		list($sNamePart , $sCmt) = preg_split('/\s/', $sLn1, 2);
		$arrR1 = explode('/', $sNamePart); // remove the trailing "/1, /2" from read names
		list($s, $sBarcode) = explode('#', $arrR1[0]);

		if ($sBarcode == '0_0_0') {
			$nDiscardPairs++;
			continue;
		}

		if (!array_key_exists($sBarcode, $arrBarcodes)) {
			$arrBarcodes[$sBarcode] = 0;
		}

		$nKeptPairs++;

		$arrBarcodes[$sBarcode]  += 1;
	}
	
	
}

echo("Discarded $nDiscardPairs pairs, kept $nKeptPairs, ". count($arrBarcodes) ." barcodes\n\n");

?>
