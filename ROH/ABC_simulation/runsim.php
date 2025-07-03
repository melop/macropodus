<?php
$sCfgTemplate = "admixsimul.cfg.template";
$sRecTemplate = "rec.template";
$arrDatafiles = array("genes.txt", "markers.txt","mutations.txt","null.txt","phenotype.txt");
$sRunDir = "run";
$sOutDir = "out";
$sWinSize = 50000;
$sStepSize = 10000;
$nMaxIndForStats = 20;
$nTotalJobs = intval($argv[1]);
$nThisJob = intval($argv[2]);
$nReps = intval($argv[3]);
$arrGenPrior = array(1.072, 4.09); //priors were obtained by first performing a preliminary run using linear priors, then plotting them in R to find the best region
$arrPopSizePrior = array(0.63, 3.76);
$arrRecPerChr = array(-0.03400595, 0.568054); //in Betta splendens, the rec rate is 6.6 centimorgan (cM)/Mb https://www.science.org/doi/10.1126/sciadv.abm4950, that is 1.849371 Morgan/chr1 in Macropodus, prior is half of that and double that
$bPriorInLog10Scale = true; //priors are given in log10 scale
$sRscript = "Rscript admixem_roh.R";
$sAdmix = "/public/software/admixem/bin/admixemp";
$sWD = getcwd();

exec("mkdir -p $sOutDir/");

$sCfgTemplate = file_get_contents($sCfgTemplate);
$sRecTemplate = file_get_contents($sRecTemplate);

$sOut = "$sOutDir/out.$nThisJob.$nTotalJobs.txt";
$hO = 0;
if (!file_exists($sOut)) {
	$hO = fopen($sOut, 'w');
	fwrite($hO, "gen\tpopsize\trec_per_chr\twinlen_mean\twinlen_sd\tFroh\n");
} else {
	$hO = fopen($sOut, 'a');
}

for($i=0;$i<=$nReps;$i++) {
	$sRepRunDir = "$sRunDir/run.$nThisJob.$nTotalJobs";
	exec("rm -r $sRepRunDir; mkdir -p $sRepRunDir");
	foreach($arrDatafiles as $sFile) {
		exec("ln -sf $sWD/$sFile $sRepRunDir/");
	}

	$nGen = 0;
	$nPopSize = 0;
	$nRec = 0;

	if (!$bPriorInLog10Scale) {
		$nGen = mt_rand($arrGenPrior[0], $arrGenPrior[1]);
		$nPopSize = mt_rand($arrPopSizePrior[0], $arrPopSizePrior[1]);
		$nRec = ($arrRecPerChr[1]-$arrRecPerChr[0])* (mt_rand()/mt_getrandmax()) + $arrRecPerChr[0];
	} else {
		$nGen = round(pow(10, ($arrGenPrior[1]-$arrGenPrior[0])* (mt_rand()/mt_getrandmax()) + $arrGenPrior[0]));
		$nPopSize = round(pow(10, ($arrPopSizePrior[1]-$arrPopSizePrior[0])* (mt_rand()/mt_getrandmax()) + $arrPopSizePrior[0]));
		$nRec = pow(10, ($arrRecPerChr[1]-$arrRecPerChr[0])* (mt_rand()/mt_getrandmax()) + $arrRecPerChr[0]);
	}

	$sCfg = str_replace('%MAXGEN%', $nGen, $sCfgTemplate);
	$sCfg = str_replace('%SAMPLEFREQ%', $nGen, $sCfg);
	$sCfg = str_replace('%POPSIZE%', $nPopSize, $sCfg);

	$sRec = str_replace('%REC%', $nRec, $sRecTemplate);
	file_put_contents("$sRepRunDir/cfg.txt" , $sCfg);
	file_put_contents("$sRepRunDir/rec.txt" , $sRec);

	exec("cd $sRepRunDir; $sAdmix cfg.txt > /dev/null", $arrOut, $nRetCode);

	$sAdmixOutput = "$sRepRunDir/out/Gen".$nGen."_markers.txt";
	if ($nRetCode !== 0 || (!file_exists($sAdmixOutput) )) {
		echo("Warning: admixem error out gen=$nGen popsize=$nPopSize $sRepRunDir\n");
		continue;
	}

	$hW = popen("$sRscript chr1.pos.txt $sAdmixOutput $sWinSize $sStepSize $nMaxIndForStats ", 'r');
	$sLn = fgets( $hW);
	pclose($hW);

	if ($sLn!== false && trim($sLn)!=="") {
		fwrite($hO, "$nGen\t$nPopSize\t$nRec\t$sLn\n");
	}
}


?>
