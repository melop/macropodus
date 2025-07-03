<?php
$sRd1Files = "/data/projects/shareddata/mop.stLFR/stLFR/raw_data_2021-07-16/genome/NGS/muscle/128.3G/D2104505A.clean_1.fq.gz"; //separate each file by a space, multiple files are allowed. for example: "test_1_1.fq.gz test_2_1.fq.gz"
$sRd2Files = "/data/projects/shareddata/mop.stLFR/stLFR/raw_data_2021-07-16/genome/NGS/muscle/128.3G/D2104505A.clean_2.fq.gz"; //separate each file by a space, multiple files are allowed. for example: "test_1_1.fq.gz test_2_1.fq.gz"

$sOut = "sorted_for_lariat.fq.gz";
$nSortThreads=32;

$sSampleBarcode = "TAAGGCGC"; //just for record keeping, can set to anything

$sTmpOut = "$sOut.tsv";

$nUniqueBarcodeNum =  47378361 ;//obtained by countbarcodes.php;
$nReportEvery = 1000000;


$sBinMaxNum = decbin($nUniqueBarcodeNum);
$nBinLen = strlen($sBinMaxNum);
$nBinLen += ($nBinLen % 2 == 1)? 1:0;

$arrCode2Base = array('00' => 'A', '01' => 'C', '10' => 'G', '11' => 'T' );
//die();

$h1 = popen("zcat -f $sRd1Files", 'r'); // read from the decompression stream of the zcat command, this command will concatenate multiple gzip files into one stream.
$h2 = popen("zcat -f $sRd2Files", 'r');

$hOTmp = fopen($sTmpOut, 'w');

$nLn = -1;            // current line number of the input stream.
$arrLns1 = array();   // buffer for file 1
$nKeptPairs = 0;      // counter for kept pairs
$nDiscardPairs = 0;   // counter for discarded pairs.

$arrBGIto10XBarcode = array(); //BGI 2 10x barcode map
$bSkipThisPair = false;

while(true) {
	$nLn++;										// increase line number counter, line 1 = 0
	$sLn1 = fgets($h1);							// get the next line from the read 1 stream
	$sLn2 = fgets($h2);							// get the next line from the read 2 stream


	if ($sLn1 === false || $sLn2 === false) { // if end of file, sLn1 and sLn2 become false
		if ($sLn1 !== $sLn2) {
			echo("Warning: two read files do not reach the end at the same point\n"); //warn if the two input streams don't reach the end at the same time.
		}
		break;
	}
	
	$sLn1 = trim($sLn1); // remove trailing "\n"
	$sLn2 = trim($sLn2);

	if ($sLn1 == '' || $sLn2 == '') {
		echo("Empty line encountered on line $nLn, exit\n");
		break;
	}
	
	if ($nLn % 4 == 0) { // title line
		//check title lines:
		if ($sLn1[0] != '@' || $sLn2[0] != '@') { //check if title lines start with @
			die("Error: unexpected title lines\n");			
		}
		$arrR1 = explode('/', $sLn1); // remove the trailing "/1, /2" from read names
		$arrR2 = explode('/', $sLn2);
		
		if ($arrR1[0] != $arrR2[0]) { // check if two read names are identical 
			die("Error: on line $nLn, two read names differ. $sRd1Name, $sRd2Name\n"); // if two read names differ, errors out.
		}

		list($sRdName, $sBGIBarcode) = explode("#", $arrR1[0]);
		
		if ($sBGIBarcode == '0_0_0') {
			$bSkipThisPair = true;
			$nDiscardPairs++;
			continue;
		} 

		$bSkipThisPair = false;

		if (!array_key_exists($sBGIBarcode, $arrBGIto10XBarcode)) {

			$s10XBarcode = fnBin2DNA(count($arrBGIto10XBarcode));//compute new barcode
			$s10XBarcode .= "-1"; //append -1 suffix
			$arrBGIto10XBarcode[$sBGIBarcode] = $s10XBarcode;

		}

		$s10XBarcode = $arrBGIto10XBarcode[$sBGIBarcode] ;
		//clear buffers:
		$arrLns1 = array_fill(0, 9, '');
		$arrLns1[0] = $sRdName;
		$arrLns1[5] = $s10XBarcode;
		$arrLns1[6] = str_repeat('I', strlen($s10XBarcode)-2);
		
		continue; //next line;

	}

	if ($bSkipThisPair) {
		continue; //skip this pair
	}
	
	if ($nLn % 4 == 1) { // dna sequence line

		$arrLns1[1] = $sLn1; // perform the cutting, and pushes the cut line into the buffer.
		$arrLns1[3] = $sLn2; // perform the cutting, and pushes the cut line into the buffer.
		
		continue;
		
	}
	
	if ($nLn % 4 == 2) { // + line
		if ($sLn1[0] != '+' || $sLn2[0] != '+') {
			die("+ is not found on expected line.\n"); //check if this line starts with "+", if not, errors out.
			
		}
		
		
		continue;
	}
	
	if ($nLn % 4 == 3) { // quality line
		
	

		$arrLns1[2] = $sLn1;
		$arrLns1[4] = $sLn2;

		$arrLns1[7] = $sSampleBarcode ;
		$arrLns1[8] = str_repeat('I', strlen($sSampleBarcode ));
		
		//write to output streams
		fwrite($hOTmp, implode("\t", $arrLns1) . "\n" );


		$nKeptPairs++; // counter +1

		if ( $nKeptPairs % $nReportEvery == 0 ) {
			echo($nKeptPairs."\n");
		}

	}
	
}

pclose($h1);
pclose($h2);
pclose($hOTmp);

echo("Discarded $nDiscardPairs pairs, trimmed and kept $nKeptPairs\n\n");
echo("sorting ...\n");

$sCmd = "( sort -k6,6 --parallel=$nSortThreads $sTmpOut | ".' tr "\t" "\n" | '. " gzip -c > $sOut ; ) && rm $sTmpOut ";
echo($sCmd);
exec($sCmd);
echo('done\n');

function fnBin2DNA($nNum) {
	global $nBinLen, $arrCode2Base;
	$sBinNum = decbin($nNum);
	$sBinNumPad = str_pad($sBinNum,$nBinLen, '0', STR_PAD_LEFT );

	//echo("$sBinNumPad\n");
	$nStart = $nBinLen -2;
	$sBase = '';
	for($i=0;$i<=$nStart;$i+=2) {
		$sBase .= $arrCode2Base[substr($sBinNumPad, $i, 2)];
	}

	return $sBase;

	
}
?>
