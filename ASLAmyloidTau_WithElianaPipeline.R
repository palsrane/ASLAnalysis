library(mgcv)
library(plyr)
library(nlme)
library(lme4)
library(ggplot2)
library(plsdepot)
library(stringr)
library(tidymv)
library(olsrr)
library(MatchIt)
library(tidyr)
library(patchwork)
library(ADNIMERGE)
library(grid)
library(cowplot)

rm(list = ls())
`%notin%` <- Negate(`%in%`)
setwd("C:/Projects/NextProj/ASLAmyloidTau") #Directory path

cols <- c("0"="green","1"="darkgreen","2"="blue","3"="red","9"="black") #Colors for CAA values

getASL <- function(){
  studydate= read.csv("MAYOADIRL_MRI_QUALITY_ADNI3_07Dec2023.csv") #MAYOADIRL_MRI_QUALITY_ADNI3_25May2023.csvSERIES_DATE is the acquision date
  studydate$SERIES_Time=sapply(strsplit(studydate$SERIES_DATE," "),`[`,2) #Because the Date has date and time, and the time portion will be saved here
  studydate$SERIES_DATE=sapply(strsplit(studydate$SERIES_DATE," "),`[`,1) #Because the Date has date and time, and the date portion will be saved here
  
  # ASLtmp=read.csv("pijp_asl_withScanDate.csv")
  ASLPerf2=read.csv("GE_ASL_data_passedQC_120623.csv") #GE_ASL_data_passedQC_070523.csv
  ASLPerf2$RIDYear=sapply(strsplit(ASLPerf2$ImageCode,"_"),`[`,4)
  ASLPerf2$RID=sapply(strsplit(sapply(strsplit(ASLPerf2$ImageCode,"_"),`[`,4),"y"),`[`,1)
  ASLPerf2$StudyYear=sapply(strsplit(sapply(strsplit(ASLPerf2$ImageCode,"_"),`[`,4),"y"),`[`,2)
  ASLPerf2$LONI_IMAGE=sapply(strsplit(sapply(strsplit(ASLPerf2$ImageCode,"_"),`[`,5),"i"),`[`,2)
  ASLPerf2=ASLPerf2[ASLPerf2$LONI_IMAGE %in% studydate$LONI_IMAGE[studydate$SERIES_SELECTED=='TRUE'],]
  
  FSQC= read.csv("FSQCResults_09_18_2023.csv")
  FSQC=FSQC[FSQC$LONI_IMAGE %in% studydate$LONI_IMAGE[studydate$SERIES_SELECTED=='TRUE'],]
  
  ASLtmp = ASLPerf2[,c("RID","ProjectID","RIDYear","ImageCode","StudyYear","LONI_IMAGE","CBFQC","SegmentationQC","RegistrationQC")]
  ASLtmp = ASLtmp[!duplicated(ASLtmp$LONI_IMAGE),]
  for (ROI in unique(ASLPerf2$ROI)) {
    ASLtmp= merge(ASLtmp,ASLPerf2[ASLPerf2$ROI==ROI ,c("LONI_IMAGE","Mean")],by="LONI_IMAGE") #& ASLPerf2$Type=="CBF_pvcorr" was removed because that was the only thing that was downloaded from Azure database
    ROI=str_replace_all(ROI,"[[*]]","")
    ROI=str_replace_all(ROI,"[[-]]","_")
    colnames(ASLtmp)[colnames(ASLtmp)=="Mean"]=str_replace_all(ROI,"[[*]]","") # To remove the * from "Right-Thalamus-Proper*", Not sure why the processing has added that *
  }
  
  ASLtmp=merge(ASLtmp,studydate[studydate$LONI_IMAGE %in% ASLtmp$LONI_IMAGE,c("LONI_IMAGE","SERIES_DATE","SERIES_SELECTED","SERIES_COMMENTS")],by="LONI_IMAGE")
  ASLtmp=ASLtmp[order(ASLtmp$RID,ASLtmp$SERIES_DATE),]
  
  ASLtmp<-ASLtmp[ASLtmp$CBFQC==1 & ASLtmp$SegmentationQC==1 & ASLtmp$RegistrationQC==1,] 
  #None of the following were done in this dataset & (ASLtmp$FreeSurfer_OverallQC=="Pass" | (ASLtmp$FreeSurfer_OverallQC=="Partial" & ASLtmp$FS_PartialPass=="TRUE" & ASLtmp$FreeSurfer_TemporalQC=="Pass" & ASLtmp$FreeSurfer_FrontalQC =="Pass" & ASLtmp$FreeSurfer_ParietalQC =="Pass" & ASLtmp$FreeSurfer_InsulaQC =="Pass" & ASLtmp$FreeSurfer_OccipitalQC =="Pass" & ASLtmp$FreeSurfer_BasalGangliaQC =="Pass" & ASLtmp$FreeSurfer_CerebralWMQC =="Pass")),]
  
  FSQC=FSQC[FSQC$RIDYear%in%ASLtmp$RIDYear,]
  FSQC=merge(FSQC,studydate[studydate$LONI_IMAGE %in% FSQC$LONI_IMAGE,c("LONI_IMAGE","SERIES_DATE","SERIES_SELECTED","SERIES_COMMENTS")],by="LONI_IMAGE")
  FSQC = subset(FSQC, select = -c(LONI_IMAGE) ) #Removing LONI_IMAGE of the T1 scan
  
  ASLtmp=merge(ASLtmp,FSQC,by=c("RIDYear","SERIES_DATE"))
  
  ASLtmp=ASLtmp[!duplicated(ASLtmp$RIDYear),]
  # Following was corrected in pipeline by checking for MayoQC passed and paired T1 scans. Instead added the part to remove any duplicates above.
  # # Correct ASL ImageCode and T1 Scancode Combos
  # # ASL ImageCode                 T1 ScanCode
  # # ADNI3_005_S_0610y00_i906093  ADNI3_005_S_0610y00_i906096 => FALSE for SERIES_SELECTED as per MAYODIRL, TRUE = LONI_IMAGE 906097
  # # ADNI3_127_S_5200y00_i871940  ADNI3_127_S_5200y00_i871945 => FALSE for SERIES_SELECTED as per MAYODIRL, TRUE = LONI_IMAGE 871944 
  # # ADNI3_099_S_6175y01_i1189636 ADNI3_099_S_6175y01_i1189640 => TRUE for SERIES_SELECTED as per MAYODIRL, FALSE = ADNI3_099_S_6175y01_i1189639
  # # ADNI3_129_S_6288y00_i980826  ADNI3_129_S_6288y00_i980830 => TRUE for SERIES_SELECTED as per MAYODIRL, FALSE = ADNI3_129_S_6288y00_i980829
  # # ADNI3_129_S_6288y01_i1156423 ADNI3_129_S_6288y01_i1156427 => FALSE for SERIES_SELECTED as per MAYODIRL, TRUE = LONI_IMAGE 1156426
  # # ADNI3_127_S_6436y00_i1021431 ADNI3_127_S_6436y00_i1021434 => TRUE for SERIES_SELECTED as per MAYODIRL, FALSE = ADNI3_127_S_6436y00_i1021435
  # # ADNI3_135_S_6703y00_i1147708 ADNI3_135_S_6703y00_i1147712 => TRUE for SERIES_SELECTED as per MAYODIRL, FALSE ADNI3_135_S_6703y00_i1147711
  # 
  # # #Following have multiple ASL acquisitions per scan. Need to find out which one is correct.
  # #Going to use the later one for now under the assumption that the reson for the second acq was the failure in the first.
  # # ADNI3_127_S_4197y01_i1081533 =acq10@"2018-12-03 14:35:06.0" ADNI3_127_S_4197y01_i1081537 => As per MAYODIRL either ASLs as OK
  # # ADNI3_127_S_4197y01_i1081534 =acq6@"2018-12-03 14:05:12.0" ADNI3_127_S_4197y01_i1081537
  # # ADNI3_126_S_4896y01_i1214810 =acq10@"2019-08-20 16:09:28.0" ADNI3_126_S_4896y01_i1214816 => Selected
  # # ADNI3_126_S_4896y01_i1214812 =acq6@"2019-08-20 15:40:34.0" ADNI3_126_S_4896y01_i1214816 => ADNI3_126_S_4896y01_i1214812 was not selected as per MAYODIRL
  # # **ADNI3_027_S_5169y00_i881751 = FALSE, ADNI3_027_S_5169y00_i881752=TRUE (but not in results, probably failed Reg, Seg, CBF QCs) ADNI3_027_S_5169y00_i881755 (T1)
  # 
  # #May be this can be achieved by choosing the T1 scan closed to ASL scan in time, by keeping the time stamp from studydate
  # ASLtmp = ASLtmp[!((ASLtmp$ImageCode=='ADNI3_005_S_0610y00_i906093' & ASLtmp$ScanCode!='ADNI3_005_S_0610y00_i906097') | #removed ScanCode ADNI3_005_S_0610y00_i906096
  #                     (ASLtmp$ImageCode=='ADNI3_127_S_5200y00_i871940' & ASLtmp$ScanCode!='ADNI3_127_S_5200y00_i871944')| #removed ADNI3_127_S_5200y00_i871945
  #                     (ASLtmp$ImageCode=='ADNI3_099_S_6175y01_i1189636' & ASLtmp$ScanCode!='ADNI3_099_S_6175y01_i1189640')|  
  #                     (ASLtmp$ImageCode=='ADNI3_129_S_6288y00_i980826' & ASLtmp$ScanCode!='ADNI3_129_S_6288y00_i980830')|
  #                     (ASLtmp$ImageCode=='ADNI3_129_S_6288y01_i1156423' & ASLtmp$ScanCode!='ADNI3_129_S_6288y01_i1156426')| #removed ADNI3_129_S_6288y01_i1156427
  #                     (ASLtmp$ImageCode=='ADNI3_127_S_6436y00_i1021431' & ASLtmp$ScanCode!='ADNI3_127_S_6436y00_i1021434')|
  #                     (ASLtmp$ImageCode=='ADNI3_135_S_6703y00_i1147708' & ASLtmp$ScanCode!='ADNI3_135_S_6703y00_i1147712')|
  #                     (ASLtmp$ImageCode=='ADNI3_127_S_4197y01_i1081534')| #selecting ASL image ADNI3_127_S_4197y01_i1081533 instead
  #                     (ASLtmp$ImageCode=='ADNI3_126_S_4896y01_i1214812')),] #selecting ASL image ADNI3_126_S_4896y01_i1214810 instead
  
  
  
  colnames(ASLtmp)[colnames(ASLtmp)=='SERIES_DATE']='ASL_DATE'
  
  ASLtmp$Entorhinal=as.numeric(((ASLtmp$ctx_lh_entorhinal/ASLtmp$ctx_lh_precentral)+(ASLtmp$ctx_rh_entorhinal/ASLtmp$ctx_rh_precentral))/2)
  # Amyloid PET was associated with lower CBF in temporo-parietal regions.
  ASLtmp$Hippocampus=as.numeric(((ASLtmp$Left_Hippocampus/ASLtmp$ctx_lh_precentral)+(ASLtmp$Right_Hippocampus/ASLtmp$ctx_rh_precentral))/2) #hippocampus
  ASLtmp$InferiorTemporal=as.numeric(((ASLtmp$ctx_lh_inferiortemporal/ASLtmp$ctx_lh_precentral)+(ASLtmp$ctx_rh_inferiortemporal/ASLtmp$ctx_rh_precentral))/2) #inferior temporal
  ASLtmp$InferiorParietal=as.numeric(((ASLtmp$ctx_lh_inferiorparietal/ASLtmp$ctx_lh_precentral)+(ASLtmp$ctx_rh_inferiorparietal/ASLtmp$ctx_rh_precentral))/2) #inferior parietal
  ASLtmp$PosteriorCingulate=as.numeric(((ASLtmp$ctx_lh_posteriorcingulate /ASLtmp$ctx_lh_precentral)+(ASLtmp$ctx_rh_posteriorcingulate /ASLtmp$ctx_rh_precentral))/2) #posterior cingulate
  ASLtmp$Precuneus=as.numeric(((ASLtmp$ctx_lh_precuneus/ASLtmp$ctx_lh_precentral)+(ASLtmp$ctx_rh_precuneus/ASLtmp$ctx_rh_precentral))/2) #precuneous
  ASLtmp$MedialOrbitofronal=as.numeric(((ASLtmp$ctx_lh_medialorbitofrontal/ASLtmp$ctx_lh_precentral)+(ASLtmp$ctx_rh_medialorbitofrontal /ASLtmp$ctx_rh_precentral))/2) #medial orbitofronal
  ASLtmp$Pericalcarine=as.numeric(((ASLtmp$ctx_lh_pericalcarine/ASLtmp$ctx_lh_precentral)+(ASLtmp$ctx_rh_pericalcarine/ASLtmp$ctx_rh_precentral))/2) #pericalcarine
  ASLtmp$InsthumasCingulate=as.numeric(((ASLtmp$ctx_lh_isthmuscingulate/ASLtmp$ctx_lh_precentral)+(ASLtmp$ctx_rh_isthmuscingulate/ASLtmp$ctx_rh_precentral))/2)#insthumas cingulate based on Braak4 region
  return(ASLtmp)
}

getAPOE <-function(){
  apoe=read.csv('APOERES_1.csv')
  apoe=read.csv('APOERES_20Sep2023.csv')
  apoe1=read.csv('APOE_ADNI3.csv') #Send by Loni rep, not yet on ADNI Loni website
  
  apoe=apoe[, c('RID','APGEN1','APGEN2')]
  apoe$APOE4=NA
  apoe$APOE4[apoe$APGEN1==4 | apoe$APGEN2==4]=1
  apoe$APOE4[apoe$APGEN1==4 & apoe$APGEN2==4]=2
  apoe$APOE4[apoe$APGEN1!=4 & apoe$APGEN2!=4]=0
  apoe$APOE2=NA
  apoe$APOE2[apoe$APGEN1==2 | apoe$APGEN2==2]=1
  apoe$APOE2[apoe$APGEN1==2 & apoe$APGEN2==2]=2
  apoe$APOE2[apoe$APGEN1!=2 & apoe$APGEN2!=2]=0
  apoe$APOE = 'NA'
  apoe$APOE[apoe$APOE2==0 & apoe$APOE4==0] = 'E3/E3'
  apoe$APOE[apoe$APOE2==1 & apoe$APOE4==0] = 'E3/E2'
  apoe$APOE[apoe$APOE2==0 & apoe$APOE4==1] = 'E3/E4'
  apoe$APOE[apoe$APOE2==2 & apoe$APOE4==0] = 'E2/E2'
  apoe$APOE[apoe$APOE2==0 & apoe$APOE4==1] = 'E3/E4'
  apoe$APOE[apoe$APOE2==1 & apoe$APOE4==1] = 'E2/E4'
  
  
  apoe=apoe[apoe$APOE!='NA',]
  apoe$APOE = factor(apoe$APOE)
  apoe$APOE <- relevel(apoe$APOE, ref='E3/E3')
  apoe=apoe[!duplicated(apoe[c("RID")]),c('RID','APOE')]
  
  apoe1=apoe1[apoe1$RID %notin% apoe$RID,]
  apoe1$APOE=paste0(apoe1$APOE.Value.1,"/",apoe1$APOE.Value.2)
  apoe1$APOE[apoe1$APOE=='E2/E3']='E3/E2'
  apoe1$APOE[apoe1$APOE=='E4/E3']='E3/E4'
  apoe1=apoe1[!duplicated(apoe1[c("RID")]),c('RID','APOE')]
  apoe=rbind(apoe,apoe1)
  return(apoe)
}

getDemog <- function(){
  dem=read.csv('PTDEMOG_3.csv') #as of 5/30/2023
  dem=dem[!is.na(dem$PTGENDER) & !is.na(dem$PTEDUCAT) & dem$PTGENDER!=-4,]
  colnames(dem)[colnames(dem)=='PTGENDER']='Sex'
  dem$Sex=factor(dem$Sex)
  levels(dem$Sex)=c('M','F') #initially 1 and 2
  dem$refDOB=as.Date(paste0(as.character(dem$PTDOBYY),'-',as.character(dem$PTDOBMM),'-01'),format="%Y-%m-%d") #added 01 of birthday month as reference date to calculate  at different tests
  colnames(dem)[colnames(dem)=='USERDATE']='dem.DATE'
  dem=dem[,c('RID','PTEDUCAT','Sex','dem.DATE','refDOB',"PTRACCAT","PTETHCAT")]
  #PTRACCAT: 1=American Indian or Alaskan Native; 2=Asian; 3=Native Hawaiian or Other Pacific Islander;
  #          4=Black or African American; 5=White; 6=More than one race; 7=Unknown
  #PTETHCAT: 1=Hispanic or Latino; 2=Not Hispanic or Latino; 3=Unknown
  
  return(dem)
  
}

getDXSUM <- function(){
  dxsm=read.csv('DXSUM_PDXCONV_ADNIALL_14Jun2023.csv')
  # dxsm$dt=as.Date(dxsm$EXAMDATE,format="%Y-%m-%d")-as.Date(dxsm$USERDATE,format="%Y-%m-%d")
  # dxsm$EXAMDATE[dxsm$EXAMDATE!='' & dxsm$USERDATE!='' & dxsm$dt>0]=dxsm$USERDATE[dxsm$EXAMDATE!='' & dxsm$USERDATE!='' & dxsm$dt>0]
  # dxsm$EXAMDATE[dxsm$EXAMDATE=='' & dxsm$USERDATE!='']=dxsm$USERDATE[dxsm$EXAMDATE=='' & dxsm$USERDATE!='']
  dxsm=dxsm[(dxsm$EXAMDATE!='' & !is.na(dxsm$EXAMDATE))| (dxsm$USERDATE!=''& !is.na(dxsm$USERDATE)),]
  dxsm$EXAMDATE[(dxsm$EXAMDATE=='' | is.na(dxsm$EXAMDATE)) & (dxsm$USERDATE!='' & !is.na(dxsm$USERDATE))]=dxsm$USERDATE[(dxsm$EXAMDATE=='' | is.na(dxsm$EXAMDATE)) & (dxsm$USERDATE!='' & !is.na(dxsm$USERDATE))]
  dxsm$USERDATE[(dxsm$USERDATE=='' | is.na(dxsm$USERDATE))& (dxsm$EXAMDATE!='' & !is.na(dxsm$EXAMDATE))]=dxsm$EXAMDATE[(dxsm$USERDATE=='' | is.na(dxsm$USERDATE))& (dxsm$EXAMDATE!='' & !is.na(dxsm$EXAMDATE))]
  dxsm$dt=as.Date(dxsm$EXAMDATE,format="%Y-%m-%d")-as.Date(dxsm$USERDATE,format="%Y-%m-%d")
  dxsm$EXAMDATE[dxsm$dt>0]=dxsm$USERDATE[dxsm$dt>0]
  
  colnames(dxsm)[colnames(dxsm)=='EXAMDATE']='DX.DATE'
  dxsm=dxsm[!is.na(dxsm$DXCHANGE) | !is.na(dxsm$DIAGNOSIS),]
  dxsm$DXSUM=NA
  dxsm$DXSUM[dxsm$DIAGNOSIS==1 | dxsm$DXCHANGE==1 | dxsm$DXCHANGE==7 | dxsm$DXCHANGE==9]='CN'
  dxsm$DXSUM[dxsm$DIAGNOSIS==2 | dxsm$DXCHANGE==2 | dxsm$DXCHANGE==4 | dxsm$DXCHANGE==8]='MCI'
  dxsm$DXSUM[dxsm$DIAGNOSIS==3 | dxsm$DXCHANGE==3 | dxsm$DXCHANGE==5 | dxsm$DXCHANGE==6]='AD'
  dxsm$DXSUM<-factor(dxsm$DXSUM)
  dxsm=dxsm[!is.na(dxsm$DX.DATE) & !is.na(dxsm$DXSUM),c('RID','DX.DATE','DXSUM')]
  return(dxsm)
}

getCDR <- function(){
  cdr=read.csv('CDR_11Jun2023.csv') 
  cdr=cdr[(cdr$EXAMDATE!='' & !is.na(cdr$EXAMDATE))| (cdr$USERDATE!=''& !is.na(cdr$USERDATE)),]
  cdr$EXAMDATE[(cdr$EXAMDATE=='' | is.na(cdr$EXAMDATE))& (cdr$USERDATE!='' & !is.na(cdr$USERDATE))]=cdr$USERDATE[(cdr$EXAMDATE=='' | is.na(cdr$EXAMDATE))& (cdr$USERDATE!='' & !is.na(cdr$USERDATE))]
  cdr$USERDATE[(cdr$USERDATE=='' | is.na(cdr$USERDATE))& (cdr$EXAMDATE!='' & !is.na(cdr$EXAMDATE))]=cdr$EXAMDATE[(cdr$USERDATE=='' | is.na(cdr$USERDATE))& (cdr$EXAMDATE!='' & !is.na(cdr$EXAMDATE))]
  
  cdr$dt=as.Date(cdr$EXAMDATE,format="%Y-%m-%d")-as.Date(cdr$USERDATE,format="%Y-%m-%d")
  cdr$EXAMDATE[cdr$dt>0]=cdr$USERDATE[cdr$dt>0]
  
  colnames(cdr)[colnames(cdr)=='USERDATE']='CDR.DATE'
  cdr=cdr[cdr$CDGLOBAL!=-1 & !is.na(cdr$CDGLOBAL),]
  cdr$CDSOB = cdr$CDMEMORY + cdr$CDORIENT + cdr$CDJUDGE + cdr$CDCOMMUN + cdr$CDHOME + cdr$CDCARE
  cdr=cdr[,c('RID','CDR.DATE','CDSOB')]
  cdr <- cdr[order(cdr$RID, cdr$CDR.DATE),]
  return(cdr)
  }

getTau <- function(){
  # TauPET_old=read.csv("UCBERKELEYAV1451_PVC_8mm_02_17_23.csv")
  TauPET=read.csv("UCBERKELEY_TAUPVC_6MM_27Sep2023.csv")
  colnames(TauPET)[colnames(TauPET)=='SCANDATE']='TauPET_DATE' #EXAMDATE
  # TauPET$VISCODE2[TauPET$VISCODE2==""]=TauPET$VISCODE[TauPET$VISCODE2==""]  #If VISCODE2 is empty, copying VISCODE in that location. VISCODE2 seems to have more levels, so maintaining that.
  # TauPET=subset(TauPET, select = -c(VISCODE) )
  
  #Braak Region 1 for Tau = Entorhinal SUVR
  TauPET$Braak1Tau=TauPET$CTX_ENTORHINAL_SUVR#/TauPET$INFERIORCEREBELLUM_SUVR #TauPET$INFERIOR_CEREBGM_SUVR
  #Braak Region 3&4
  # Braak 3 parahippocampal, fusiform, lingual, amygdala
  TauPET$Braak3Tau=(TauPET$CTX_PARAHIPPOCAMPAL_SUVR+TauPET$CTX_FUSIFORM_SUVR+TauPET$CTX_LINGUAL_SUVR+TauPET$AMYGDALA_SUVR)/4#*TauPET$INFERIORCEREBELLUM_SUVR) #TauPET$INFERIOR_CEREBGM_SUVR
  # Braak 4 middletemporal, caudalantcing, rostralantcing, postcing, isthmuscing, insula, inferiortemporal, temporalpole
  TauPET$Braak4Tau=(TauPET$CTX_MIDDLETEMPORAL_SUVR+TauPET$CTX_LH_CAUDALANTERIORCINGULATE_SUVR+TauPET$CTX_ROSTRALANTERIORCINGULATE_SUVR+
                      TauPET$CTX_POSTERIORCINGULATE_SUVR+TauPET$CTX_ISTHMUSCINGULATE_SUVR+TauPET$CTX_INSULA_SUVR+
                      TauPET$CTX_INFERIORTEMPORAL_SUVR+TauPET$CTX_TEMPORALPOLE_SUVR)/8#*TauPET$INFERIORCEREBELLUM_SUVR)#TauPET$INFERIOR_CEREBGM_SUVR the old and the new was just 1 for these ref region values
  return(TauPET)
}

getCentiloid <- function(){
  # The new AV45 and FBB files do not have the SUMMARYSUVR_WHOLECEREBNORM, 
  # and instead have the SUMMARY_SUVR column. Assuming they are the same.
  # ROIs tested in Rubinski paper medial-orbitofrontal, precuneus, posterior cingulate, 
  # inferior parietal, inferior temporal and parahippocampal gyrus
  av45=read.csv("UCBERKELEYAV45_8mm_02_17_23.csv")
  colnames(av45)[colnames(av45)=='EXAMDATE']='AV45.DATE'
  av45=av45[,c('RID','AV45.DATE','VISCODE','VISCODE2',"CENTILOIDS","SUMMARYSUVR_COMPOSITE_REFNORM_0.78CUTOFF",
               "SUMMARYSUVR_WHOLECEREBNORM_1.11CUTOFF",
               "SUMMARY_SUVR", "CTX_MEDIALORBITOFRONTAL_SUVR","CTX_PRECUNEUS_SUVR",
               "CTX_POSTERIORCINGULATE_SUVR", "CTX_INFERIORPARIETAL_SUVR",
               "CTX_INFERIORTEMPORAL_SUVR","CTX_PARAHIPPOCAMPAL_SUVR")]
  
  av45[,c("SUMMARY_SUVR", "CTX_MEDIALORBITOFRONTAL_SUVR","CTX_PRECUNEUS_SUVR",
          "CTX_POSTERIORCINGULATE_SUVR", "CTX_INFERIORPARIETAL_SUVR",
          "CTX_INFERIORTEMPORAL_SUVR","CTX_PARAHIPPOCAMPAL_SUVR")]=196.9* av45[,c("SUMMARY_SUVR", "CTX_MEDIALORBITOFRONTAL_SUVR","CTX_PRECUNEUS_SUVR",
                                                                                  "CTX_POSTERIORCINGULATE_SUVR", "CTX_INFERIORPARIETAL_SUVR",
                                                                                  "CTX_INFERIORTEMPORAL_SUVR","CTX_PARAHIPPOCAMPAL_SUVR")]-196.03
  fbb=read.csv("UCBERKELEYFBB_8mm_02_17_23.csv")
  fbb$FBB.DATE=fbb$EXAMDATE
  fbb=fbb[,c('RID','FBB.DATE','VISCODE','VISCODE2',"CENTILOIDS","SUMMARYSUVR_COMPOSITE_REFNORM_0.74CUTOFF","SUMMARYSUVR_WHOLECEREBNORM_1.08CUTOFF",
             "SUMMARY_SUVR", "CTX_MEDIALORBITOFRONTAL_SUVR","CTX_PRECUNEUS_SUVR",
             "CTX_POSTERIORCINGULATE_SUVR", "CTX_INFERIORPARIETAL_SUVR",
             "CTX_INFERIORTEMPORAL_SUVR","CTX_PARAHIPPOCAMPAL_SUVR")]
  fbb[,c("SUMMARY_SUVR", "CTX_MEDIALORBITOFRONTAL_SUVR","CTX_PRECUNEUS_SUVR",
         "CTX_POSTERIORCINGULATE_SUVR", "CTX_INFERIORPARIETAL_SUVR",
         "CTX_INFERIORTEMPORAL_SUVR","CTX_PARAHIPPOCAMPAL_SUVR")]=159.08*fbb[,c("SUMMARY_SUVR", "CTX_MEDIALORBITOFRONTAL_SUVR","CTX_PRECUNEUS_SUVR",
                                                                                "CTX_POSTERIORCINGULATE_SUVR", "CTX_INFERIORPARIETAL_SUVR",
                                                                                "CTX_INFERIORTEMPORAL_SUVR","CTX_PARAHIPPOCAMPAL_SUVR")]-151.65
  
  Centiloid=merge(av45,fbb,by=c('RID','VISCODE','VISCODE2'),all=T)
  Centiloid$VISCODE2[Centiloid$VISCODE2==""]=Centiloid$VISCODE[Centiloid$VISCODE2==""]  #IF VISCODE2 is empty, copying to VISCODE to it. VISCODE2 seems to have more levels, hence maintaining it.
  Centiloid$CENTILOIDS = Centiloid$CENTILOIDS.x
  Centiloid$SUM_Centiloid=Centiloid$SUMMARY_SUVR.x
  Centiloid$MOF_Centiloid=Centiloid$CTX_MEDIALORBITOFRONTAL_SUVR.x
  Centiloid$MOF_Centiloid[is.na(Centiloid$SUM_Centiloid)]= Centiloid$CTX_MEDIALORBITOFRONTAL_SUVR.y[is.na(Centiloid$SUM_Centiloid)]
  Centiloid$PRECU_Centiloid=Centiloid$CTX_PRECUNEUS_SUVR.x
  Centiloid$PRECU_Centiloid[is.na(Centiloid$SUM_Centiloid)]= Centiloid$CTX_PRECUNEUS_SUVR.y[is.na(Centiloid$SUM_Centiloid)]
  Centiloid$PCC_Centiloid=Centiloid$CTX_POSTERIORCINGULATE_SUVR.x
  Centiloid$PCC_Centiloid[is.na(Centiloid$SUM_Centiloid)]= Centiloid$CTX_POSTERIORCINGULATE_SUVR.y[is.na(Centiloid$SUM_Centiloid)]
  Centiloid$INFP_Centiloid=Centiloid$CTX_INFERIORPARIETAL_SUVR.x
  Centiloid$INFP_Centiloid[is.na(Centiloid$SUM_Centiloid)]= Centiloid$CTX_INFERIORPARIETAL_SUVR.y[is.na(Centiloid$SUM_Centiloid)]
  Centiloid$PHC_Centiloid=Centiloid$CTX_PARAHIPPOCAMPAL_SUVR.x
  Centiloid$PHC_Centiloid[is.na(Centiloid$SUM_Centiloid)]= Centiloid$CTX_PARAHIPPOCAMPAL_SUVR.y[is.na(Centiloid$SUM_Centiloid)]
  Centiloid$Cent.DATE=Centiloid$AV45.DATE
  Centiloid$Cent.DATE[is.na(Centiloid$SUM_Centiloid)]=Centiloid$FBB.DATE[is.na(Centiloid$SUM_Centiloid)]
  Centiloid$AmyloidPos=Centiloid$SUMMARYSUVR_WHOLECEREBNORM_1.11CUTOFF
  Centiloid$AmyloidPos[is.na(Centiloid$CENTILOIDS)]=Centiloid$SUMMARYSUVR_WHOLECEREBNORM_1.08CUTOFF[is.na(Centiloid$CENTILOIDS)]
  Centiloid$AmyloidPos_old=Centiloid$SUMMARYSUVR_COMPOSITE_REFNORM_0.78CUTOFF
  Centiloid$AmyloidPos_old[is.na(Centiloid$SUM_Centiloid)]=Centiloid$SUMMARYSUVR_COMPOSITE_REFNORM_0.74CUTOFF[is.na(Centiloid$SUM_Centiloid)]
  
  Centiloid$SUM_Centiloid[is.na(Centiloid$SUM_Centiloid)]= Centiloid$SUMMARY_SUVR.y[is.na(Centiloid$SUM_Centiloid)]
  Centiloid$CENTILOIDS[is.na(Centiloid$CENTILOIDS)]=Centiloid$CENTILOIDS.y[is.na(Centiloid$CENTILOIDS)]
  
  Centiloid = Centiloid[!is.na(Centiloid$SUM_Centiloid),c('RID','Cent.DATE','AmyloidPos','AmyloidPos_old',"CENTILOIDS","SUM_Centiloid","MOF_Centiloid","PRECU_Centiloid","PCC_Centiloid","INFP_Centiloid","PHC_Centiloid")] #"VISCODE2",
  Centiloid <- Centiloid[order(Centiloid$RID, Centiloid$Cent.DATE),]
  Centiloid <- Centiloid[!duplicated(Centiloid[c("RID","Cent.DATE")]),] #Keep all unique centiloid scan visits
  rm(av45,fbb)
  return(Centiloid)
}

plotdatascatter <-function(tmpdd,plt_labs){
  ggp=ggplot(tmpdd, aes(x = Age, y = reg,
                        group = RID, #Switching to Cross sectional analysis
                        colour = AmyloidPos)) +
    geom_line(data=tmpdd[tmpdd$NumOfScans>1,],aes(group=RID)) + #,size=1.1
    # scale_color_manual(values=cols) + scale_x_continuous(limits = c(x_min, x_max)) +
    geom_point(aes(shape=DXSUM))+ #
    plt_labs+
    theme(panel.background = element_rect(fill = "white", colour = "#6D9EC1",size = 2, linetype = "solid"),
          panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour ="black"),
          strip.background = element_blank(),
          strip.placement = "outside")  
  print(paste("Agerange",min(tmpdd$Age),max(tmpdd$Age)))
  return(ggp)
}


#ASL Data----
# ASL_list=read.csv('PCASLResultsFromAzureDatasetQuery_05-11-2023_Separated.csv')
# MAYOADIRL=read.csv("MAYOADIRL_MRI_QUALITY_ADNI3_25May2023.csv")
# colnames(MAYOADIRL)[colnames(MAYOADIRL)=='SERIES_DATE']='ASL_DATE'
# MAYOADIRL=MAYOADIRL[,c("RID","LONI_IMAGE", "PTID","ASL_DATE", "SERIES_TYPE")] #MAYOADIRL$SERIES_TYPE=="ASL" | MAYOADIRL$SERIES_TYPE=="ASL_CBF"
# MAYOADIRL=MAYOADIRL[MAYOADIRL$RID %in% ASL_list$RID,]
# ASL_list1 = merge(ASL_list,MAYOADIRL, by=c("RID","LONI_IMAGE"),all.x = TRUE) #part of ScanCode after i (ADNI3_128_S_0272y00_i1001500) is LONI_IMAGE 
# 
# # #The following Two were not present in the final database which passed the three Reg, Segmentation, and CBF QCs
# # # Manually adding scan date for  based on B:\RawData\OrganizedExternalData\Dicoms\ADNI3\ADNI3_027_S_5170y02\2019-09-06-13-27-38\ADNI3 Basic Human Protoc\acq6__Axial 3D pCASL (Eyes Open)
# # ASL_list1$ASL_DATE[ASL_list1$RID==5170 & ASL_list1$StudyYear==2 & is.na(ASL_list1$ASL_DATE)]='2019-09-06 13:27:38.0'
# # ASL_list1$PTID[ASL_list1$RID==5170 & ASL_list1$StudyYear==2 & is.na(ASL_list1$PTID)]='027_S_5170'
# # ASL_list1$SERIES_TYPE[ASL_list1$RID==5170 & ASL_list1$StudyYear==2 & is.na(ASL_list1$SERIES_TYPE)]='ASL'
# # #/m/RawData/OrganizedExternalData/Dicoms/ADNI3/ADNI3_005_S_6393y00/2018-06-05-14-28-37/Research ADNI3 Basic Hu/acq6__Axial 3D pCASL (Eyes Open)
# # ASL_list1$ASL_DATE[ASL_list1$RID==6393 & ASL_list1$StudyYear==0 & is.na(ASL_list1$ASL_DATE)]='2018-06-05 14:28:37.0'
# # ASL_list1$PTID[ASL_list1$RID==6393 & ASL_list1$StudyYear==0 & is.na(ASL_list1$PTID)]='005_S_6393'
# # ASL_list1$SERIES_TYPE[ASL_list1$RID==6393 & ASL_list1$StudyYear==0 & is.na(ASL_list1$SERIES_TYPE)]='ASL'

# #Tried this alternate method to find scan dates by getting the scan list from Loni using Advanced Search
# # LONI_ASL=read.csv('ASLResultsFromLoni_5_11_2023.csv')
# # colnames(LONI_ASL)[colnames(LONI_ASL)=='Image.ID']='LONI_IMAGE'
# # LONI_ASL=LONI_ASL[,c("Subject.ID","Study.Date","Sex","Age","Description","LONI_IMAGE","Research.Group","Visit")]
# # ASL_list2 = merge(ASL_list,LONI_ASL, by=c("Subject.ID","LONI_IMAGE"),all.x = TRUE) #part of ScanCode after i (ADNI3_128_S_0272y00_i1001500) is LONI_IMAGE 
# ASLPerf=ASL_list1
# ASLPerf$VISCODE2[ASLPerf$VISCODE2==""]=ASLPerf$VISCODE[ASLPerf$VISCODE2==""]
# 

ASLPerf=getASL()


#Adding demographic and diagnostic data info ----
dem=getDemog()
dem=dem[dem$RID %in% unique(ASLPerf$RID),]
dem <- dem[order(dem$RID,dem$dem.DATE),]
dem <- dem[!duplicated(dem[c("RID")]),]
ASLPerf=merge(ASLPerf,dem,by=c('RID'),all.x=T)
ASLPerf$Age=as.numeric(as.Date(ASLPerf$ASL_DATE,format="%Y-%m-%d")-as.Date(ASLPerf$refDOB,format="%Y-%m-%d"))/365 #Age at the time of ASLPerf scan
ASLPerf=ASLPerf[!is.na(ASLPerf$Age),]
rm(dem)
ASLPerf=subset(ASLPerf, select = -c(refDOB,dem.DATE) ) #dt

dxsm=getDXSUM()
dxsm=dxsm[dxsm$RID %in% unique(ASLPerf$RID),]
dxsm <- dxsm[order(dxsm$RID,dxsm$DX.DATE),]
ASLPerf=merge(ASLPerf,dxsm,by=c('RID'),all.x=T)# keeping all ASL records,all.x=T as only keeping RIDs that have DXSUM on them 
ASLPerf$dt=abs(as.Date(ASLPerf$DX.DATE,format="%Y-%m-%d")-as.Date(ASLPerf$ASL_DATE,format="%Y-%m-%d"))
ASLPerf <- ASLPerf[order(ASLPerf$RID,ASLPerf$ASL_DATE,ASLPerf$dt),]
# ASLPerf<-ASLPerf[ASLPerf$dt<=180 & !is.na(ASLPerf$dt),]#limiting to observations which have the other covariates collected within 180 days
ASLPerf <- ASLPerf[!duplicated(ASLPerf[c("RID","ASL_DATE")]),] #Keep the DXSM exam closest in date to the ASLPerf EXAMDATE
ASLPerf = subset(ASLPerf, select = -c(DX.DATE) )
ASLPerf=ASLPerf[!is.na(ASLPerf$DXSUM),]
rm(dxsm)
#RID 6482 does not have DXSUM available, so removed.

APOE = getAPOE()
APOE = APOE[APOE$RID %in% ASLPerf$RID,]
ASLPerf=merge(ASLPerf,APOE,by=c('RID'),all.x=T)# keeping all ASL records,all.x=T as only keeping RIDs that have DXSUM on them 
rm(APOE)
#No APOE available for RIDs 2245 4891 5083 6077 6082 6093 6155 6286 6289 6370 6389 6433 6482 6549 6567 6579 6687 6793 6840

ASLPerf <- ASLPerf[order(ASLPerf$RID,ASLPerf$ASL_DATE),]
#very first ASL record. This will be used to lock in dates for first Tau and Amyloid PET records.
ASLfirst <- ASLPerf[!duplicated(ASLPerf[c("RID")]),] 

ASLPerf$NumOfScans=c(1)
ASLPerf<-ASLPerf[order(ASLPerf$RID,ASLPerf$ASL_DATE),]
for (id in unique(ASLPerf$RID)) 
{
  ASLPerf$NumOfScans[ASLPerf$RID==id]=length(ASLPerf$RID[ASLPerf$RID==id])
}

#Amyloid PET data----
Centiloid = getCentiloid()
# ASLPerfHold=ASLPerf
ASLPerf= merge(ASLPerf,Centiloid[,c("RID","Cent.DATE","CENTILOIDS","SUM_Centiloid","AmyloidPos","AmyloidPos_old")], by=c('RID'), all.x = TRUE)
# RIDs 6077 6093 6155 6286 6567 6579 6793 have no Amyloid PET available on them
ASLPerf$dt=as.Date(ASLPerf$ASL_DATE,format="%Y-%m-%d")-as.Date(ASLPerf$Cent.DATE,format="%Y-%m-%d")
# ASLPerf = ASLPerf[order(ASLPerf$RID,ASLPerf$ASL_DATE,abs(ASLPerf$dt)),]
# ASLPerf = ASLPerf[!duplicated(ASLPerf[,c("RID","ASL_DATE")]),]
ASLPerf$CentiloidCalculated='No' #The column tells us if the centiloid values are calculated or not.

#Keeping Records that have Amyloid data withing 180 days of ASL, or Amyloid data collected 180 days or more before ASL scna showed amyloid positivity, or Amyloid data collected 180 days or more after ASL showed amyloid negativity
# ASLPerf = ASLPerf[abs(ASLPerf$dt)<=180 | (ASLPerf$dt< -180 & ASLPerf$AmyloidPos ==0) | (ASLPerf$dt > 180 & ASLPerf$AmyloidPos==1) | is.na(ASLPerf$dt),] 

#ASL scans that have Centiloid scan within 180 days
tmpASL=ASLPerf[which(abs(ASLPerf$dt)<=180),]
tmpASL = tmpASL[order(tmpASL$RID,tmpASL$ASL_DATE,abs(tmpASL$dt)),]
tmpASL = tmpASL[!duplicated(tmpASL[,c("RID","ASL_DATE")]),]
tmpASL$CentiloidCalculated='No' 

ASLPerf_notintmpASL=ASLPerf[ASLPerf$RIDYear%notin%tmpASL$RIDYear,]

tmpASL=rbind(tmpASL,ASLPerf_notintmpASL[is.na(ASLPerf_notintmpASL$AmyloidPos),]) #9 records that do not have amyloid info were added back

# ASL scans that had centiloid scans that were collected > 180 days apart
Negdt=ASLPerf_notintmpASL[which(ASLPerf_notintmpASL$dt< -180),]
Negdt=Negdt[order(Negdt$RID,Negdt$ASL_DATE,abs(Negdt$dt)),]
Negdt=Negdt[!duplicated(Negdt[,c("RID","ASL_DATE")]),]
Posdt=ASLPerf_notintmpASL[which(ASLPerf_notintmpASL$dt>180),]
Posdt=Posdt[order(Posdt$RID,Posdt$ASL_DATE,abs(Posdt$dt)),]
Posdt=Posdt[!duplicated(Posdt[,c("RID","ASL_DATE")]),]
rm(ASLPerf_notintmpASL)


tmpNegdt=Negdt[which(Negdt$RIDYear%notin%Posdt$RIDYear & Negdt$AmyloidPos==0),] #ASL records that only have Centiloid scan done prior to ASL scan and they were already amyloid positive at that time.
if(length(Negdt$RIDYear[which(Negdt$RIDYear%notin%Posdt$RIDYear & Negdt$AmyloidPos==0)])>0) tmpNegdt[,c("SUM_Centiloid","CENTILOIDS")] <- NA #"MOF_Centiloid","PRECU_Centiloid","PCC_Centiloid","INFP_Centiloid","PHC_Centiloid"
#Can not use specific centiloid levels for such scans, but can use Amyloid positivity info as the subject was Amyloid Negative at a later date
tmpASL=rbind(tmpASL,tmpNegdt)
Negdt=Negsdt[Negdt$RIDYear%notin%tmpNegdt$RIDYear,]
rm(tmpNegdt)

tmpPosdt=Posdt[which(Posdt$RIDYear%notin%Negdt$RIDYear & Posdt$AmyloidPos==1),] #ASL records that only have Centiloid scan done after ASL scan and they were already amyloid positive at that time.
if(length(Posdt$RIDYear[which(Posdt$RIDYear%notin%Negdt$RIDYear & Posdt$AmyloidPos==1)])>0) tmpPosdt[,c("SUM_Centiloid","CENTILOIDS")] <- NA #"MOF_Centiloid","PRECU_Centiloid","PCC_Centiloid","INFP_Centiloid","PHC_Centiloid"
#Can not use specific centiloid levels for such scans, but can use Amyloid positivity info as the subject was Amyloid Positive at an earlier date
tmpASL=rbind(tmpASL,tmpPosdt)
Posdt=Posdt[Posdt$RIDYear%notin%tmpPosdt$RIDYear,]
rm(tmpPosdt)

#Since the following had either an Amyloid scan collected only >180days prior that shows amyloid negativity, or amyloid scan collected only >180 days post ASL that show amyloid positivity,
#Amyloid info can not be associated with those ASL scans
tmpNoAmyData=rbind(Negdt[Negdt$RIDYear%notin%Posdt$RIDYear & Negdt$AmyloidPos==1,],Posdt[Posdt$RIDYear%notin%Negdt$RIDYear & Posdt$AmyloidPos==0,])
tmpNoAmyData[,c("AmyloidPos","SUM_Centiloid","AmyloidPos_old","CENTILOIDS")] <- NA #"MOF_Centiloid","PRECU_Centiloid","PCC_Centiloid","INFP_Centiloid","PHC_Centiloid"
Negdt=Negdt[Negdt$RIDYear %notin%tmpNoAmyData$RIDYear,]
Posdt=Posdt[Posdt$RIDYear%notin%tmpNoAmyData$RIDYear,]

tmpASL=rbind(tmpASL,tmpNoAmyData)
rm(tmpNoAmyData)

#Computing Amyloid levels for ASL scan date, based on amyloid levels that were seen before and after ASL scan
tmpASL=rbind(tmpASL,Negdt)

for(ridyr in unique(Negdt$RIDYear)){
  # tmpASL$dt[tmpASL$RIDYear==ridyr]= print(paste(Negdt$dt[Negdt$RIDYear==ridyr],"to",Posdt$dt[Posdt$RIDYear==ridyr]))
  tmpASL$SUM_Centiloid[tmpASL$RIDYear==ridyr]=Posdt$SUM_Centiloid[Posdt$RIDYear==ridyr]-(Posdt$SUM_Centiloid[Posdt$RIDYear==ridyr]-Negdt$SUM_Centiloid[Negdt$RIDYear==ridyr])*as.numeric(Posdt$dt[Posdt$RIDYear==ridyr])/as.numeric(Posdt$dt[Posdt$RIDYear==ridyr]-Negdt$dt[Negdt$RIDYear==ridyr])
  tmpASL$CENTILOIDS[tmpASL$RIDYear==ridyr]=Posdt$CENTILOIDS[Posdt$RIDYear==ridyr]-(Posdt$CENTILOIDS[Posdt$RIDYear==ridyr]-Negdt$CENTILOIDS[Negdt$RIDYear==ridyr])*as.numeric(Posdt$dt[Posdt$RIDYear==ridyr])/as.numeric(Posdt$dt[Posdt$RIDYear==ridyr]-Negdt$dt[Negdt$RIDYear==ridyr])
  ifelse(tmpASL$CENTILOIDS[tmpASL$RIDYear==ridyr]>=20,tmpASL$AmyloidPos[tmpASL$RIDYear==ridyr]<-1,tmpASL$AmyloidPos[tmpASL$RIDYear==ridyr]<-0) 
  # ifelse(tmpASL$SUM_Centiloid[tmpASL$RIDYear==ridyr]>=20,tmpASL$AmyloidPos[tmpASL$RIDYear==ridyr]<-1,tmpASL$AmyloidPos[tmpASL$RIDYear==ridyr]<-0) 
  # tmpASL$MOF_Centiloid[tmpASL$RIDYear==ridyr]=Posdt$MOF_Centiloid[Posdt$RIDYear==ridyr]-(Posdt$MOF_Centiloid[Posdt$RIDYear==ridyr]-Negdt$MOF_Centiloid[Negdt$RIDYear==ridyr])*as.numeric(Posdt$dt[Posdt$RIDYear==ridyr])/as.numeric(Posdt$dt[Posdt$RIDYear==ridyr]-Negdt$dt[Negdt$RIDYear==ridyr])
  # tmpASL$PRECU_Centiloid[tmpASL$RIDYear==ridyr]=Posdt$PRECU_Centiloid[Posdt$RIDYear==ridyr]-(Posdt$PRECU_Centiloid[Posdt$RIDYear==ridyr]-Negdt$PRECU_Centiloid[Negdt$RIDYear==ridyr])*as.numeric(Posdt$dt[Posdt$RIDYear==ridyr])/as.numeric(Posdt$dt[Posdt$RIDYear==ridyr]-Negdt$dt[Negdt$RIDYear==ridyr])
  # tmpASL$PCC_Centiloid[tmpASL$RIDYear==ridyr]=Posdt$PCC_Centiloid[Posdt$RIDYear==ridyr]-(Posdt$PCC_Centiloid[Posdt$RIDYear==ridyr]-Negdt$PCC_Centiloid[Negdt$RIDYear==ridyr])*as.numeric(Posdt$dt[Posdt$RIDYear==ridyr])/as.numeric(Posdt$dt[Posdt$RIDYear==ridyr]-Negdt$dt[Negdt$RIDYear==ridyr])
  # tmpASL$INFP_Centiloid[tmpASL$RIDYear==ridyr]=Posdt$INFP_Centiloid[Posdt$RIDYear==ridyr]-(Posdt$INFP_Centiloid[Posdt$RIDYear==ridyr]-Negdt$INFP_Centiloid[Negdt$RIDYear==ridyr])*as.numeric(Posdt$dt[Posdt$RIDYear==ridyr])/as.numeric(Posdt$dt[Posdt$RIDYear==ridyr]-Negdt$dt[Negdt$RIDYear==ridyr])
  # tmpASL$PHC_Centiloid[tmpASL$RIDYear==ridyr]=Posdt$PHC_Centiloid[Posdt$RIDYear==ridyr]-(Posdt$PHC_Centiloid[Posdt$RIDYear==ridyr]-Negdt$PHC_Centiloid[Negdt$RIDYear==ridyr])*as.numeric(Posdt$dt[Posdt$RIDYear==ridyr])/as.numeric(Posdt$dt[Posdt$RIDYear==ridyr]-Negdt$dt[Negdt$RIDYear==ridyr])
  tmpASL$CentiloidCalculated[tmpASL$RIDYear==ridyr]='Yes'
  tmpASL$Cent.DATE[tmpASL$RIDYear==ridyr]<-NA
}
#Something is wrong with the way AmyloidPos has been noted for RID 6842 in the FBB table. The summary centiloid values are 0.985 and 0.998, which after conversion to Centiloid units become 5.0438 and 7.11 resp. However, the AmyloidPos value for the second Cent PET has been recorded as 1. It should have been 0 as it is below 20.

tmpASL$NumOfScans=c(1)
tmpASL<-tmpASL[order(tmpASL$RID,tmpASL$ASL_DATE),]
for (id in unique(tmpASL$RID)) {
  tmpASL$NumOfScans[tmpASL$RID==id]=length(tmpASL$RID[tmpASL$RID==id])
}

tmpASL$BaselineAmyloidPos=tmpASL$AmyloidPos

for (id in unique(tmpASL$RID[tmpASL$NumOfScans==3])) {
  print(paste(id,tmpASL$AmyloidPos[tmpASL$RID==id]))
  print(tmpASL$AmyloidPos[tmpASL$RIDYear==min(tmpASL$RIDYear[tmpASL$RID==id])])
  # tmpASL$BaselineAmyloidPos[tmpASL$RID==id]=tmpASL$AmyloidPos[tmpASL$RIDYear==min(tmpASL$RIDYear[tmpASL$RID==id])]
  
}




tmpASL$BaselineAmyloidPos=factor(tmpASL$BaselineAmyloidPos)

colnames(tmpASL)[colnames(tmpASL)=='dt']='Cent.dt'

tmpASL$reg=tmpASL$Entorhinal
plt_labs <- labs(y = 'Entorhinal perfusion',
                 x = 'Age in Years',
                 colour = 'BaselineAmyloidPos')
ggp=plotdatascatter(tmpASL,plt_labs)

ggp 


#Tau PET data ----
TauPET = getTau()
TauPET = TauPET[TauPET$RID %in% tmpASL$RID,]

#Different Tau ROIs as per Clara's work
# The Mesial Temporal: entorhinal, parahippocampus and amygdala; 
# The Meta Temporal region: entorhinal, parahippocampus, amygdala, fusiform, inferior and middle temporal gyri; 
# The Temporo-Parietal: bankssts, cuneus, inferior-superior parietal, inferior-middle-superior temporal, istmuscingulate, lateral occipital, lingal, posterior cingulate, precuneus and superior marginal; 
# The Frontal: caudate middle frontal, precentral, rostral middle frontal, Superior frontal. 
# Tracer	UniversalMask	MesialTemporal	MetaTemporal	TemporoParietal	Frontal
# 18F-FTP	13.63x–15.85	10.42x–12.11	12.95x–15.37	13.75x–15.92	11.61x–13.01

TauPET$MesialTempCTRz = 10.42*(TauPET$CTX_ENTORHINAL_SUVR + TauPET$CTX_PARAHIPPOCAMPAL_SUVR + TauPET$AMYGDALA_SUVR)/3-12.11
TauPET$MetaTempCTRz = 12.95*(TauPET$CTX_ENTORHINAL_SUVR + TauPET$CTX_PARAHIPPOCAMPAL_SUVR + TauPET$AMYGDALA_SUVR + TauPET$CTX_FUSIFORM_SUVR +
                               TauPET$CTX_INFERIORTEMPORAL_SUVR + TauPET$CTX_MIDDLETEMPORAL_SUVR)/6-15.37
TauPET$TemporoParietalCTRz = 13.75*(TauPET$CTX_BANKSSTS_SUVR + TauPET$CTX_CUNEUS_SUVR + TauPET$CTX_INFERIORPARIETAL_SUVR + TauPET$CTX_SUPERIORPARIETAL_SUVR +
                                      TauPET$CTX_INFERIORTEMPORAL_SUVR + TauPET$CTX_MIDDLETEMPORAL_SUVR + TauPET$CTX_SUPERIORTEMPORAL_SUVR +
                                      TauPET$CTX_ISTHMUSCINGULATE_SUVR + TauPET$CTX_LATERALOCCIPITAL_SUVR + TauPET$CTX_LINGUAL_SUVR +
                                      TauPET$CTX_POSTERIORCINGULATE_SUVR + TauPET$CTX_PRECUNEUS_SUVR + TauPET$CTX_SUPRAMARGINAL_SUVR)/13 -15.92
TauPET$FrontalCTRz = 11.61*(TauPET$CTX_CAUDALMIDDLEFRONTAL_SUVR + TauPET$CTX_ROSTRALMIDDLEFRONTAL_SUVR + TauPET$CTX_PRECENTRAL_SUVR + TauPET$CTX_SUPERIORFRONTAL_SUVR)/4-13.01

TauPET=cbind(TauPET,MesialTempTauPos=NA,MetaTempTauPos=NA,TemporoParietalTauPos=NA,FrontalTauPos=NA)

TauPET$MesialTempTauPos[which(TauPET$MesialTempCTRz<2 & !is.na(TauPET$MesialTempCTRz))] = 0
TauPET$MesialTempTauPos[which(TauPET$MesialTempCTRz>2 & !is.na(TauPET$MesialTempCTRz))] = 1

TauPET$MetaTempTauPos[which(TauPET$MetaTempCTRz<2 & !is.na(TauPET$MetaTempCTRz))] = 0
TauPET$MetaTempTauPos[which(TauPET$MetaTempCTRz>2 & !is.na(TauPET$MetaTempCTRz))] = 1

TauPET$TemporoParietalTauPos[which(TauPET$TemporoParietalCTRz<2 & !is.na(TauPET$TemporoParietalCTRz))] = 0
TauPET$TemporoParietalTauPos[which(TauPET$TemporoParietalCTRz>2 & !is.na(TauPET$TemporoParietalCTRz))] = 1

TauPET$FrontalTauPos[which(TauPET$FrontalCTRz<2  & !is.na(TauPET$FrontalCTRz))] = 0
TauPET$FrontalTauPos[which(TauPET$FrontalCTRz>2  & !is.na(TauPET$FrontalCTRz))] = 1

                                      
TauPET = TauPET[,c("RID","TauPET_DATE","Braak1Tau","Braak3Tau","Braak4Tau","MesialTempCTRz","MetaTempCTRz","TemporoParietalCTRz","FrontalCTRz","MesialTempTauPos","MetaTempTauPos","TemporoParietalTauPos","FrontalTauPos")]                               


#Generating thresholds for TauPositivity ----
#As per Ossenkoppele following were used for Tau stratification:
# MTL (unweighted average of bilateral entorhinal cortex and amygdala) 
# unweighted MTL ROI because we intended the entorhinal cortex to have a 
# relatively higher contribution since this region is involved in the earliest 
# stages of tau accumulation yet it is somewhat smaller compared to the amygdala

tmpASLTau=merge(tmpASL, TauPET,by=c("RID"), all.x = TRUE) #Matching by RID and VISCODE2 is working in this case
tmpASLTau$dt=as.Date(tmpASLTau$ASL_DATE,format="%Y-%m-%d")-as.Date(tmpASLTau$TauPET_DATE,format="%Y-%m-%d")
tmpASLTau = tmpASLTau[order(tmpASLTau$RID,tmpASLTau$ASL_DATE,abs(tmpASLTau$dt)),]

tmpTau = tmpASLTau[which(abs(tmpASLTau$dt)<= 180 | is.na(tmpASLTau$dt)),] # | #Tau was recorded within 180 days OR no Tau data
tmpTau = tmpTau[!duplicated(tmpTau[c("RIDYear")]),]

            
Posdt=tmpASLTau[which(tmpASLTau$RIDYear %notin% tmpTau$RIDYear & tmpASLTau$dt>0),] # dt>180
Posdt=Posdt[!duplicated(Posdt$RIDYear),]
Negdt=tmpASLTau[which(tmpASLTau$RIDYear %notin% tmpTau$RIDYear & tmpASLTau$dt<0),] #dt < -180
Negdt=Negdt[!duplicated(Negdt$RIDYear),]


#join back the records where ASL date was >180 days post Tau, and only a single Tau record exists, but Tau was was already positive
tmp = Posdt[which(Posdt$RIDYear %notin% Negdt$RIDYear & Posdt$MesialTempTauPos == 1 & Posdt$TemporoParietalTauPos==1),]
tmp[,c("Braak1Tau","Braak3Tau","Braak4Tau","MesialTempCTRz","MetaTempCTRz","TemporoParietalCTRz","FrontalCTRz")]<-NA
tmpTau=rbind(tmpTau,tmp)
#join back the records where Tau date was >180 days post ASL, and only a single Tau record exists, but Tau was -ve at the later date
tmp = Negdt[which(Negdt$RIDYear %notin% Posdt$RIDYear & Negdt$MesialTempTauPos == 0 & Negdt$TemporoParietalTauPos==0),]
tmp[,c("Braak1Tau","Braak3Tau","Braak4Tau","MesialTempCTRz","MetaTempCTRz","TemporoParietalCTRz","FrontalCTRz")]<-NA
tmpTau=rbind(tmpTau,tmp)

Posdt = Posdt[Posdt$RIDYear %notin% tmpTau$RIDYear,]
Negdt = Negdt[Negdt$RIDYear %notin% tmpTau$RIDYear,]

#Records that Tau info can't be predicted for ASL date because the only Tau collection they have show Tau negativity >180days before ASL or Tau positivity >180 days after ASL
tmp = Negdt[which(Negdt$RIDYear %notin% Posdt$RIDYear),]
tmp[,c("dt","Braak1Tau","Braak3Tau","Braak4Tau","MesialTempCTRz","MetaTempCTRz","TemporoParietalCTRz","FrontalCTRz","MesialTempTauPos","MetaTempTauPos","TemporoParietalTauPos","FrontalTauPos")]<-NA
tmpTau=rbind(tmpTau,tmp)
tmp = Posdt[which(Posdt$RIDYear %notin% Negdt$RIDYear),]
tmp[,c("dt","Braak1Tau","Braak3Tau","Braak4Tau","MesialTempCTRz","MetaTempCTRz","TemporoParietalCTRz","FrontalCTRz","MesialTempTauPos","MetaTempTauPos","TemporoParietalTauPos","FrontalTauPos")]<-NA
tmpTau=rbind(tmpTau,tmp)
rm(tmp)
tmpTau$TauCalculated="No"


Posdt = Posdt[Posdt$RIDYear %notin% tmpTau$RIDYear,]
Negdt = Negdt[Negdt$RIDYear %notin% tmpTau$RIDYear,]
#Computing Tau vales for records that have Tau records before as well as after ASL date but >180 days before or after ASL
# Tau computation is done assuming linear progression of Tau withing the before and after period


UnmatchedComputed=Posdt
UnmatchedComputed$TauCalculated="Yes"
for (RY in UnmatchedComputed$RIDYear) {
  UnmatchedComputed$dt[UnmatchedComputed$RIDYear==RY]=paste(Negdt$dt[Negdt$RIDYear==RY],"to",Posdt$dt[Posdt$RIDYear==RY])
  UnmatchedComputed$MesialTempCTRz[UnmatchedComputed$RIDYear==RY] = Posdt$MesialTempCTRz[Posdt$RIDYear==RY] + (Posdt$MesialTempCTRz[Posdt$RIDYear==RY] - Negdt$MesialTempCTRz[Negdt$RIDYear==RY]) *
    (as.numeric(as.Date(UnmatchedComputed$ASL_DATE[UnmatchedComputed$RIDYear==RY],format="%Y-%m-%d") - as.Date(Posdt$TauPET_DATE[Posdt$RIDYear==RY],format="%Y-%m-%d")))/
    (as.numeric(as.Date(Posdt$TauPET_DATE[Posdt$RIDYear==RY],format="%Y-%m-%d") - as.Date(Negdt$TauPET_DATE[Negdt$RIDYear==RY],format="%Y-%m-%d")))
  UnmatchedComputed$MetaTempCTRz[UnmatchedComputed$RIDYear==RY] = Posdt$MetaTempCTRz[Posdt$RIDYear==RY] + (Posdt$MetaTempCTRz[Posdt$RIDYear==RY] - Negdt$MetaTempCTRz[Negdt$RIDYear==RY]) *
    (as.numeric(as.Date(UnmatchedComputed$ASL_DATE[UnmatchedComputed$RIDYear==RY],format="%Y-%m-%d") - as.Date(Posdt$TauPET_DATE[Posdt$RIDYear==RY],format="%Y-%m-%d")))/
    (as.numeric(as.Date(Posdt$TauPET_DATE[Posdt$RIDYear==RY],format="%Y-%m-%d") - as.Date(Negdt$TauPET_DATE[Negdt$RIDYear==RY],format="%Y-%m-%d")))
  UnmatchedComputed$TemporoParietalCTRz[UnmatchedComputed$RIDYear==RY] = Posdt$TemporoParietalCTRz[Posdt$RIDYear==RY] + (Posdt$TemporoParietalCTRz[Posdt$RIDYear==RY] - Negdt$TemporoParietalCTRz[Negdt$RIDYear==RY]) *
    (as.numeric(as.Date(UnmatchedComputed$ASL_DATE[UnmatchedComputed$RIDYear==RY],format="%Y-%m-%d") - as.Date(Posdt$TauPET_DATE[Posdt$RIDYear==RY],format="%Y-%m-%d")))/
    (as.numeric(as.Date(Posdt$TauPET_DATE[Posdt$RIDYear==RY],format="%Y-%m-%d") - as.Date(Negdt$TauPET_DATE[Negdt$RIDYear==RY],format="%Y-%m-%d")))
  UnmatchedComputed$FrontalCTRz[UnmatchedComputed$RIDYear==RY] = Posdt$FrontalCTRz[Posdt$RIDYear==RY] + (Posdt$FrontalCTRz[Posdt$RIDYear==RY] - Negdt$FrontalCTRz[Negdt$RIDYear==RY]) *
    (as.numeric(as.Date(UnmatchedComputed$ASL_DATE[UnmatchedComputed$RIDYear==RY],format="%Y-%m-%d") - as.Date(Posdt$TauPET_DATE[Posdt$RIDYear==RY],format="%Y-%m-%d")))/
    (as.numeric(as.Date(Posdt$TauPET_DATE[Posdt$RIDYear==RY],format="%Y-%m-%d") - as.Date(Negdt$TauPET_DATE[Negdt$RIDYear==RY],format="%Y-%m-%d")))
  
}

UnmatchedComputed$MesialTempTauPos[which(UnmatchedComputed$MesialTempCTRz<2 & !is.na(UnmatchedComputed$MesialTempCTRz))] = 0
UnmatchedComputed$MesialTempTauPos[which(UnmatchedComputed$MesialTempCTRz>2 & !is.na(UnmatchedComputed$MesialTempCTRz))] = 1

UnmatchedComputed$MetaTempTauPos[which(UnmatchedComputed$MetaTempCTRz<2 & !is.na(UnmatchedComputed$MetaTempCTRz))] = 0
UnmatchedComputed$MetaTempTauPos[which(UnmatchedComputed$MetaTempCTRz>2 & !is.na(UnmatchedComputed$MetaTempCTRz))] = 1

UnmatchedComputed$TemporoParietalTauPos[which(UnmatchedComputed$TemporoParietalCTRz<2 & !is.na(UnmatchedComputed$TemporoParietalCTRz))] = 0
UnmatchedComputed$TemporoParietalTauPos[which(UnmatchedComputed$TemporoParietalCTRz>2 & !is.na(UnmatchedComputed$TemporoParietalCTRz))] = 1

UnmatchedComputed$FrontalTauPos[which(UnmatchedComputed$FrontalCTRz<2  & !is.na(UnmatchedComputed$FrontalCTRz))] = 0
UnmatchedComputed$FrontalTauPos[which(UnmatchedComputed$FrontalCTRz>2  & !is.na(UnmatchedComputed$FrontalCTRz))] = 1

tmpASLTau=rbind(tmpTau,UnmatchedComputed)

rm(UnmatchedComputed,tmpTau,Negdt,Posdt,ASLfirst,ASLPerf)

tmpASLTau$NumOfScans=c(1)
tmpASLTau<-tmpASLTau[order(tmpASLTau$RID,tmpASLTau$ASL_DATE),]
for (id in unique(tmpASLTau$RID)) {
  tmpASLTau$NumOfScans[tmpASLTau$RID==id]=length(tmpASLTau$RID[tmpASLTau$RID==id])
}

tmpASLTaufirst=tmpASLTau[!duplicated(tmpASLTau$RID),]



#Code after this needs to be worked on----
AmyTauPositivity_table=as.data.frame(matrix(nrow=1,ncol = 12))
AmyTauPositivity_table[1,]=c(paste0("Total N=",length(ASLfirst$RID),", N records=",length(tmpASL$RID)),"A-TMes-","A","","","","","","","","")









#Following code needs to be worked on


#i.e. keep centiloid scan within 6 months of Tau date, or centiloid scan that showed subject as amyloid positive 6 to 12 months prior to tau scan, or centiloid scan that showed subject as amyloid negative 6-12 months post Tau scan
#Determining Tau Positivity thresholds as per Ossenkoppele paper
#Tau threshold for MTL (Amygdata + entorhinal)

#Just to test if matching by RID and VISCODE2 gave same results as matching by RID and eleminating by time between scan dates
# temp <-na.omit(refdata[abs(refdata$dt)>180,c("RID","TauPET_DATE","Cent.DATE","AmyloidPos","dt")])
# rm(temp)
dxsm = getDXSUM()
dxsm=dxsm[dxsm$RID %in% unique(refdata$RID),]
refdata=merge(refdata,dxsm[,c('RID',"DX.DATE","DXSUM")],by=c('RID'))# all.x=T as only keeping RIDs that have DXSUM on them 
# ,'VISCODE2', Merging by RID and VISCODE2 was leading to loosing 9 records. So matching DXSUM by EXAMDATEs
refdata$dt=abs(as.Date(refdata$DX.DATE,format="%Y-%m-%d")-as.Date(refdata$TauPET_DATE,format="%Y-%m-%d"))
refdata <- refdata[order(refdata$RID,refdata$TauPET_DATE,refdata$dt),]
refdata<-refdata[refdata$dt<=180 & !is.na(refdata$dt),]#limiting to observations which have the other covariates collected within 180 days
refdata <- refdata[!duplicated(refdata[c("RID","TauPET_DATE")]),] #Keep the DXSM exam closest in date to the refdata EXAMDATE
# refdata=subset(refdata, select = -c(DX.DATE) )
refdata=refdata[!is.na(refdata$DXSUM),]
rm(dxsm)

CDR = getCDR()
CDR = CDR[CDR$RID%in%refdata$RID,]

refdata=merge(refdata,CDR,by=c('RID'))# removed ,all.x=T as only keeping RIDs that have DXSUM on them 
refdata$dt=abs(as.Date(refdata$CDR.DATE,format="%Y-%m-%d")-as.Date(refdata$TauPET_DATE,format="%Y-%m-%d"))
refdata <- refdata[order(refdata$RID,refdata$TauPET_DATE,refdata$dt),]
refdata<-refdata[refdata$dt<=180 & !is.na(refdata$dt),]#limiting to observations which have the other covariates collected within 180 days
refdata <- refdata[!duplicated(refdata[c("RID","TauPET_DATE")]),] #Keep the CDR exam closest in date to the refdata EXAMDATE
refdata=subset(refdata, select = -c(DX.DATE,dt) )
refdata=refdata[!is.na(refdata$DXSUM),]
rm(CDR)

refdataFirst=refdata[!duplicated(refdata[c("RID")]) & refdata$AmyloidPos==0 & refdata$DXSUM=="CN"  & refdata$CDGLOBAL==0,c("RID","AmyloidPos","DXSUM","CDGLOBAL","UWMTL","WNEOT")]#keeping only the first record for each participant for reference value calculation
TauUWMTLThresh=as.numeric(mean(refdataFirst$UWMTL[refdataFirst$AmyloidPos==0 & refdataFirst$DXSUM=="CN"  & refdataFirst$CDGLOBAL==0])+2*sd(refdataFirst$UWMTL[refdataFirst$AmyloidPos==0 & refdataFirst$DXSUM=="CN" & refdataFirst$CDGLOBAL==0]))
TauWNEOTThresh=as.numeric(mean(refdataFirst$WNEOT[refdataFirst$AmyloidPos==0 & refdataFirst$DXSUM=="CN"  & refdataFirst$CDGLOBAL==0])+2*sd(refdataFirst$WNEOT[refdataFirst$AmyloidPos==0 & refdataFirst$DXSUM=="CN" & refdataFirst$CDGLOBAL==0]))




# Merging Tau PET with ASLdata ----
d= TauPET[TauPET$RID %in% ASLfirst$RID,]
d= merge(d[,c("RID","TauPET_DATE")],ASLfirst[,c("RID","ASL_DATE")],by=c("RID"),all.x = T)
d$dt=abs(as.Date(d$ASL_DATE,format="%Y-%m-%d")-as.Date(d$TauPET_DATE,format="%Y-%m-%d"))
d=d[!is.na(d$dt),]
#Tau data collected closest to the first ASL scans
d <- d[order(d$RID,d$ASL_DATE,d$dt),]

d <- d[!duplicated(d[c("RID","ASL_DATE")]),] #THINK 

d <- d[d$dt<=366,]



# # Merging just ASL and Centiloid Data----
# Centiloid <- Centiloid[Centiloid$RID %in% unique(ASLPerf$RID),]
# d=merge(ASLPerf,Centiloid,by=c("RID"), all.x = TRUE)
# d$dt=abs(as.Date(d$ASL_DATE,format="%Y-%m-%d")-as.Date(d$Cent.DATE,format="%Y-%m-%d"))
# #Let's keep only the ADNIMERGE data collected closest to the ASL scans
# d <- d[order(d$RID,d$ASL_DATE,d$dt),]
# d <- d[!duplicated(d[c("RID","ASL_DATE")]),] 
# d <- d[d$dt<180,]



d = merge(d, Centiloid,by=c("RID"), all.x = TRUE)
# #Let's keep only the ADNIMERGE data collected closest to the ASL scans
# d$dt=abs(as.Date(d$ASL_DATE,format="%Y-%m-%d")-as.Date(d$Cent.DATE,format="%Y-%m-%d"))
# d=d[!is.na(d$dt),]
# d <- d[d$dt<=180,]
#Keeping Centiloid scans within 180 days, or -366 to -180 days if the centiloid scan already has Amyloid positivity as 1, 
# or within 180 to 360 days if this later centiloid scand still shows amyloid positivity as 0
d$dt=as.Date(d$ASL_DATE,format="%Y-%m-%d")-as.Date(d$Cent.DATE,format="%Y-%m-%d")
d = d[abs(d$dt)<=180| ((d$dt>=-366 & -180>d$dt & d$AmyloidPos==1) | (d$dt<=366 & d$dt>180 & d$AmyloidPos==0)),] 
d <- d[order(d$RID,d$ASL_DATE,d$dt),]
d <- d[!duplicated(d[c("RID","ASL_DATE")]) &!is.na(d$RID),] 

d$MTLTPos = 0
d$MTLTPos[d$UWMTL>TauUWMTLThresh]=1
d$NEOTTPos=0
d$NEOTTPos[d$WNEOT>TauWNEOTThresh]=1

d$AT="NA"
d$AT[d$AmyloidPos==0 & d$MTLTPos==0 & d$NEOTTPos==0]="A-T-"
d$AT[d$AmyloidPos==1 & d$MTLTPos==0 & d$NEOTTPos==0]="A+T-"
d$AT[d$AmyloidPos==1 & d$MTLTPos==1 & d$NEOTTPos==0]="A+TMTL+TNEOT-"
d$AT[d$AmyloidPos==1 & d$NEOTTPos==1]="A+TNEOT+"                     #changed this from d$AT[d$AmyloidPos==1 & (d$MTLTPos==1 | d$NEOTTPos==1)]="A+TMTLand/orTNEOT+", because there was overlap between A+Tmtl+Tneo- and A+TMTLand/orTNEOT+
d$AT[d$AmyloidPos==0 & (d$MTLTPos==1 | d$NEOTTPos==1)]="A-T+"

length(d$RID[d$AT=="NA"])
length(d$RID[d$AT=="A-T-"])
length(d$RID[d$AT=="A+T-"])
length(d$RID[d$AT=="A+TMTL+TNEOT-"])
length(d$RID[d$AT=="A+TNEOT+"])
length(d$RID[d$AT=="A-T+"])

#If AT grouping is done at the very first record for a participant
d <- d[order(d$RID,d$ASL_DATE),]
dfirst = d[!duplicated(d[c("RID")]),] 

dfirst$AT="NA"
dfirst$AT[dfirst$AmyloidPos==0 & dfirst$MTLTPos==0 & dfirst$NEOTTPos==0]="A-T-"
dfirst$AT[dfirst$AmyloidPos==1 & dfirst$MTLTPos==0 & dfirst$NEOTTPos==0]="A+T-"
dfirst$AT[dfirst$AmyloidPos==1 & dfirst$MTLTPos==1 & dfirst$NEOTTPos==0]="A+TMTL+TNEOT-"
dfirst$AT[dfirst$AmyloidPos==1 & dfirst$NEOTTPos==1]="A+TNEOT+"                     
dfirst$AT[dfirst$AmyloidPos==0 & (dfirst$MTLTPos==1 | dfirst$NEOTTPos==1)]="A-T+"
dfirst=dfirst[,c("RID","AT")]

d=merge(d,dfirst,by=c("RID"))

length(d$RID[d$AT=="NA"])
length(d$RID[d$AT=="A-T-"])
length(d$RID[d$AT=="A+T-"])
length(d$RID[d$AT=="A+TMTL+TNEOT-"])
length(d$RID[d$AT=="A+TNEOT+"])
length(d$RID[d$AT=="A-T+"])

#Entorhinal Perfusion vs Braak1Tau with Amyloidpositivity
d$NumOfScans=c(1)
d<-d[order(d$RID,d$ASL_DATE),]
for (id in unique(d$RID)) 
{
  d$NumOfScans[d$RID==id]=length(d$RID[d$RID==id])
}
plt_labs <- labs(y = 'Mean Entorhinal Perfusion',
                 x = 'Braak I Tau',
                 colour = 'AT')
ggp=ggplot(d, aes(x = Braak1Tau, y = reg1,
                  group = RID, colour = AT)) + #
  geom_line(data=d[d$NumOfScans>1,],aes(group=RID)) +
  geom_point()+ 
  plt_labs

ggp

plt_labs <- labs(y = 'Mean Entorhinal Perfusion',
                 x = 'Age in Years',
                 colour = 'APOE')
ggp=ggplot(d, aes(x = Age, y = reg1,
                  group = RID, colour = APOE)) + #
  geom_line(data=d[d$NumOfScans>1,],aes(group=RID)) +
  geom_point()+ 
  plt_labs

ggp

#Hippocampal Perfusion vs Braak1Tau with Amyloidpositivity
plt_labs <- labs(y = 'Mean Hippocampal Perfusion',
                 x = 'Braak I Tau',
                 colour = 'AmyloidPos')
ggp=ggplot(d, aes(x = Braak1Tau, y = reg2,
                  group = RID, colour = AmyloidPos)) + #
  geom_line(data=d[d$NumOfScans>1,],aes(group=RID)) +
  geom_point()+ 
  plt_labs

ggp


# Combining datasets without drops
ASLDates=ASLPerf[,c("RID","ASL_DATE","VISCODE2","Age")]
TauDates=TauPET[TauPET$RID %in% ASLPerf$RID,c("RID","TauPET_DATE","VISCODE2")]
CentDates=Centiloid[Centiloid$RID %in% ASLPerf$RID,c("RID","Cent.DATE","VISCODE2")]
colnames(ASLDates)[colnames(ASLDates)=='VISCODE2']='ASLVISCODE2'
colnames(ASLDates)[colnames(ASLDates)=='Age']='ASLAge'
colnames(TauDates)[colnames(TauDates)=='VISCODE2']='TauVISCODE2'
colnames(CentDates)[colnames(CentDates)=='VISCODE2']='CentVISCODE2'

dem=getDemog()
dem = dem[,c("RID","refDOB")]
TauDates = merge(TauDates,dem,by=c("RID"), all.x = TRUE)
CentDates = merge(CentDates,dem,by=c("RID"), all.x = TRUE)

TauDates$TauAge=as.numeric(as.Date(TauDates$TauPET_DATE,format="%Y-%m-%d")-as.Date(TauDates$refDOB,format="%Y-%m-%d"))/365
CentDates$CentAge=as.numeric(as.Date(CentDates$Cent.DATE,format="%Y-%m-%d")-as.Date(CentDates$refDOB,format="%Y-%m-%d"))/365

TauDates=subset(TauDates, select = -c(refDOB) ) #dt
CentDates=subset(CentDates, select = -c(refDOB) ) #dt

ASLTauDates= merge(ASLDates,TauDates,by=c("RID")) #Only keeping RIDs that are present in both 
ASLTauDates = ASLTauDates[order(ASLTauDates$RID, ASLTauDates$ASL_DATE),]
ASLTauDates = ASLTauDates[!duplicated(ASLTauDates[,c("RID","ASL_DATE")]),]

ASLTauMatchedDates=ASLTauDates
ASLTauMatchedDates$dt = abs(as.Date(ASLTauMatchedDates$ASL_DATE,format="%Y-%m-%d")-as.Date(ASLTauMatchedDates$TauPET_DATE,format="%Y-%m-%d"))
ASLTauMatchedDates = ASLTauMatchedDates[order(ASLTauMatchedDates$RID,ASLTauMatchedDates$ASL_DATE,ASLTauMatchedDates$dt),]
ASLTauMatchedDates = ASLTauMatchedDates[!duplicated(ASLTauMatchedDates[,c("RID","ASL_DATE")]),]
ASLTauMatchedDates = ASLTauMatchedDates[ASLTauMatchedDates$dt<=180,]
ASLTauMatchedDates$ATmatched="A+T"

ASLPerf=merge(ASLPerf,ASLTauMatchedDates[,c("RID","ASL_DATE","ATmatched")],by=c("RID","ASL_DATE"),all.x = TRUE)
TauPET=merge(TauPET,ASLTauMatchedDates[,c("RID","TauPET_DATE","ATmatched")],by=c("RID","TauPET_DATE"),all.x = TRUE)



ASLCentDates = merge(ASLDates,CentDates,by=c("RID")) #Only keeping RIDs that are present in both 
ASLCentDates = ASLCentDates[order(ASLCentDates$RID, ASLCentDates$ASL_DATE),]
ASLCentDates = ASLCentDates[!duplicated(ASLCentDates[,c("RID","ASL_DATE")]),]

ASLCentMatchedDates=ASLCentDates
ASLCentMatchedDates$dt = abs(as.Date(ASLCentMatchedDates$ASL_DATE,format="%Y-%m-%d")-as.Date(ASLCentMatchedDates$Cent.DATE,format="%Y-%m-%d"))
ASLCentMatchedDates = ASLCentMatchedDates[order(ASLCentMatchedDates$RID,ASLCentMatchedDates$ASL_DATE,ASLCentMatchedDates$dt),]
ASLCentMatchedDates = ASLCentMatchedDates[!duplicated(ASLCentMatchedDates[,c("RID","ASL_DATE")]),]
ASLCentMatchedDates = ASLCentMatchedDates[ASLCentMatchedDates$dt<=180,]
ASLCentMatchedDates$ACmatched="A+C"

ASLPerf=merge(ASLPerf,ASLCentMatchedDates[,c("RID","ASL_DATE","ACmatched")],by=c("RID","ASL_DATE"),all.x = TRUE)
Centiloid=merge(Centiloid,ASLCentMatchedDates[,c("RID","Cent.DATE","ACmatched")],by=c("RID","Cent.DATE"),all.x = TRUE)

#----

#Biomarker data ----
BIOMARKERS=read.csv("UPENNBIOMK_MASTER_FINAL.csv")
colnames(BIOMARKERS)[colnames(BIOMARKERS)=='EXAMDATE']='BIOM_DATE'
BIOMARKERS=BIOMARKERS[BIOMARKERS$BATCH=="UPENNBIOMK9" | BIOMARKERS$BATCH=="UPENNBIOMK10" | BIOMARKERS$BATCH=="UPENNBIOMK11" | BIOMARKERS$BATCH=="UPENNBIOMK12",]
#batch 13 is being reviewed currently and bathes prior to 9 were analyzed differetnly.
#most data available for Abeta42, pTau, and Tau
BIOMARKERS=BIOMARKERS[,c("RID","BIOM_DATE","BATCH","ABETA42","PTAU","TAU")]







#Merging Biomarker data----
d=merge(d,BIOMARKERS,by=c("RID"),all=T)
d$dt=abs(as.Date(d$ASL_DATE,format="%Y-%m-%d")-as.Date(d$BIOM_DATE,format="%Y-%m-%d"))
#Let's keep only the ADNIMERGE data collected closest to the ASL scans
d <- d[order(d$RID,d$ASL_DATE,d$dt),]
d <- d[!duplicated(d[c("RID","ASL_DATE")]),] 
d <- d[d$dt<180,]










#Adding Neuropathology data to the ones available on NOT ENOUGH RECORDS, IGNORE----
d= merge(d,NeuroPathFile[,c("RID","NPAMY")],by=c("RID"),all.x = T)
d$NPAMY[is.na(d$NPAMY)]=9 #Missing NPAMY(CAA rating)=9
d$Age=as.numeric(d$Age)
d$NPAMY<-factor(d$NPAMY) 
#----