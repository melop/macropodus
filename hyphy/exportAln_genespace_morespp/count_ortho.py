Fasta="/data/projects/zwang/macropodus_compare/hyphy/exportAln_genespace_morespp/ret_improved_cleanDNA.fasta"
taxa=["MOP", "MHK", "MOP_YN", "MHK_HN","MSP", "MOC",  "Anabas_testudineus", "Betta_splendens", "Danio_rerio", "Gadus_morhua", "Monopterus_albus", "Oryzias_latipes", "Takifugu_rubripes","Astyanax_mexicanus" ,"Gasterosteus_aculeatus" ,"Lepisosteus_oculatus" ,"Oreochromis_niloticus" ,"Scleropages_formosus" ,"Xiphophorus_maculatus"]
s1=s2=s3=s4=s5=s6=s7=s8=s9=s10=s11=s12=s13=s14=s15=s16=s17=s18=s19=0
with open(Fasta, "r") as f1:
    for line in f1:
        if line[0] == ">":
            if taxa[0] in line:
               s1+=1
            elif taxa[1] in line:
                s2+=1
            elif taxa[2] in line:
                s3+=1
            elif taxa[3] in line:
                s4+=1
            elif taxa[4] in line:
                s5+=1
            elif taxa[5] in line:
                s6+=1
            elif taxa[6] in line:
                s7+=1
            elif taxa[7] in line:
                s8+=1
            elif taxa[8] in line:
                s9+=1
            elif taxa[9] in line:
                s10+=1
            elif taxa[10] in line:
                s11+=1
            elif taxa[11] in line:
                s12+=1
            elif taxa[12] in line:
                s13+=1
            elif taxa[13] in line:
                s14+=1
            elif taxa[14] in line:
                s15+=1
            elif taxa[15] in line:
                s16+=1
            elif taxa[16] in line:
                s17+=1
            elif taxa[17] in line:
                s18+=1
            else:
                s19+=1
        else:
            pass
print(taxa[0],s1,"\n",
    taxa[1],s2,"\n",
    taxa[2],s3,"\n",
    taxa[3],s4,"\n",
    taxa[4],s5,"\n",
    taxa[5],s6,"\n",
    taxa[6],s7,"\n",
    taxa[7],s8,"\n",
    taxa[8],s9,"\n",
    taxa[9],s10,"\n",
    taxa[10],s11,"\n",
    taxa[11],s12,"\n",
    taxa[12],s13,"\n",
    taxa[13],s14,"\n",
    taxa[14],s15,"\n",
    taxa[15],s16,"\n",
    taxa[16],s17,"\n",
    taxa[17],s18,"\n",
    taxa[18],s19)








