
arrOutgroupSpp='Anabas_testudineus,Betta_splendens,Danio_rerio,Gadus_morhua,Monopterus_albus,Oryzias_latipes,Takifugu_rubripes,Astyanax_mexicanus,Gasterosteus_aculeatus,Oreochromis_niloticus,Scleropages_formosus,Xiphophorus_maculatus,Lepisosteus_oculatus';
#arrOutgroupSpp='Cheilinus_undulatus';
arrUnusedClades="NA"; #c('Aplocheilidae', 'root'); #mark as unused clades "only for hyphy relax"
arrMarkUnusedCladeChildren='F'; #c(T , F);

#######################################################################
#arrForegroundClades='BTP'; #mark as foreground
#arrForegroundClades='MHK'; #mark as foreground
#arrForegroundClades='Macropodus'; #mark as foreground
#arrForegroundClades='Clade2'; #mark as foreground
arrForegroundClades='MOP'; #mark as foreground
arrMarkForegroundChildren='T';

sTaxonRequirements="min_taxon_requirements_MHK.txt";

nMarkStyle='relax'; #either "codeml" or "relax"
#sOutDIR="Relax_BTP";
#sOutDIR="Relax_MHK";
#sOutDIR="Relax_Macropodus";
#sOutDIR="Relax_Bifurcation";
sOutDIR="Relax_MOP";

mkdir -p $sOutDIR
Rscript labelbranches.R $sOutDIR  $nMarkStyle  $sTaxonRequirements  $arrOutgroupSpp  $arrForegroundClades  $arrMarkForegroundChildren  $arrUnusedClades  $arrMarkUnusedCladeChildren > $sOutDIR/log.txt &

#sOutDIR="Relax_Callopanchax_mrca";
#arrMarkForegroundChildren='F';
#mkdir -p $sOutDIR
#Rscript labelbranches.R $sOutDIR  $nMarkStyle  $sTaxonRequirements  $arrOutgroupSpp  $arrForegroundClades  $arrMarkForegroundChildren  $arrUnusedClades  $arrMarkUnusedCladeChildren > $sOutDIR/log.txt &

#sOutDIR="Codeml_BTP";
#sOutDIR="Codeml_MHK";
#sOutDIR="Codeml_Bifurcation";
sOutDIR="Codeml_MOP";

arrMarkForegroundChildren='F';
nMarkStyle='codeml';
mkdir -p $sOutDIR
Rscript labelbranches.R $sOutDIR  $nMarkStyle  $sTaxonRequirements  $arrOutgroupSpp  $arrForegroundClades  $arrMarkForegroundChildren  $arrUnusedClades  $arrMarkUnusedCladeChildren > $sOutDIR/log.txt &

#sOutDIR="Codeml_Callopanchax_mrca";
#arrMarkForegroundChildren='F';
#mkdir -p $sOutDIR
#Rscript labelbranches.R $sOutDIR  $nMarkStyle  $sTaxonRequirements  $arrOutgroupSpp  $arrForegroundClades  $arrMarkForegroundChildren  $arrUnusedClades  $arrMarkUnusedCladeChildren > $sOutDIR/log.txt &
########################################################################

wait
