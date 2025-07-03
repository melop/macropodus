sPath="/data/projects/zwang/macropodus_compare/hyphy/exportAln_genespace_morespp/genespace.orthogroups.txt"
sList=[]
sEp=s0=s1=s2=s3=s4=s5=s6=s7=s8=s9=s10=s11=s12=s13=s14=s15=sMp=sMEp=sMS=sP=sMp_8=0
with open(sPath, "r") as f1:
    for line in f1:
        if line.startswith("Group"):
            sList=line.split()
            for i in sList[3:18]:
                if i=="NA":
                    sEp+=1
                else:
                    pass
            if sEp==1:
                s1+=1
            elif sEp==2:
                s2+=1
            elif sEp==3:
                s3+=1
            elif sEp==4:
                s4+=1
            elif sEp==5:
                s5+=1
            elif sEp==6:
                s6+=1
            elif sEp==7:
                s7+=1
            elif sEp==8:
                s8+=1
            elif sEp==9:
                s9+=1
            elif sEp==10:
                s10+=1
            elif sEp==11:
                s11+=1
            elif sEp==12:
                s12+=1
            elif sEp==13:
                s13+=1
            elif sEp==14:
                s14+=1
            elif sEp==15:
                s15+=1
            else:
                s0+=1
            sEp=0
            if sList[10]!="NA" and sList[11]!="NA":
                sMp+=1
                for j in sList[3:18]:
                    if j!="NA":
                        sP+=1
                    else:
                        pass
                if sP>=10:
                    sMp_8+=1
                else:
                    pass
                sP=0
            elif sList[10]=="NA" and sList[11]=="NA":
                sMEp+=1
            else:
                sMS+=1

        else:
            pass
    print(s0,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,"Macropodus all:",sMp,"Macropodus all and other 8 spp.:",sMp_8,"Macropodus single:",sMS,"Macropodus empty:",sMEp)
