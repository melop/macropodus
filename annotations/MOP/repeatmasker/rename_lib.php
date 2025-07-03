<?php
$sIn="../repeatmodeler/MOP-families.fa";
$sOut = "mop-te-families.fa";
$sSpecies = "Macropodus_opercularis";

$h = fopen($sIn, 'r');
$hO = fopen($sOut, 'w');

while(false !== ($sLn = fgets($h)) ) {
	if ($sLn[0] != '>') {
		fwrite($hO, $sLn);
		continue;
	}

	list($s1, $s2)  =  explode(' ' , $sLn);

	fwrite($hO, "$s1 @$sSpecies\n");
}
?>
