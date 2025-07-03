<?php

$nMaxNum = 41453226 ;//99999999;
$sBinMaxNum = decbin($nMaxNum);
$nLen = strlen($sBinMaxNum);
$nLen += ($nLen % 2 == 1)? 1:0;
$nNum = 2; //41453226 ; //2323321;
$sBinNum = decbin($nNum);
$sBinNumPad = str_pad($sBinNum,$nLen, '0', STR_PAD_LEFT );

echo("$sBinNumPad\n");
$arrCode2Base = array('00' => 'A', '01' => 'C', '10' => 'G', '11' => 'T' );
$nStart = $nLen -2;
$sBase = '';
for($i=0;$i<=$nStart;$i+=2) {
	$sBase .= $arrCode2Base[substr($sBinNumPad, $i, 2)];
}

echo("$sBase\n");

?>
