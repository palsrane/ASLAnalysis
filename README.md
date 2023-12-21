# ASLAnalysis
This pipeline analyzes the ASL results from ADNI3.
This one works on data that was available on Azure database as of 12/01/2023. 
The Amyloid and Tau info was taken from Loni that was last updated sometime in May 2023. 
So, far the code assocites Amyloid and Tau PET info with the ASL data based on time between the PET scand and ASL scans.
If the PET scans are collected within 180 days of ASL scan, the PET info is directly linked to that ASL scan. 
If there is only one PET scan correlating with an ASL scan, and the amyloid (or Tau) PET was collected >180 days prior and shows amyloid (/Tau) positivity, OR it was collected >180 days later and shows amyloid (/Tau) negativity, the Positivity.Negativity status was retained, but the actual (Centilois/CTRz) levels were changed to NA.
For other single PET scans (>180 later + positivity OR >180 days earlier + negativity) => no use of that PET scan, and all the PET scan parameters were changed to NA
If there are two PET scans available, one before and one after the ASL scan, then those were used to compute the intermediate Amyloid or Tau values at the ASL date, assuming linear progression.

There are 318 ASL records. 
235 records have both amyloid and Tau positivity info (285 have atleast one of those) 
176 baseline ASL records 
140 baseline records have both Tau and Amyloid positivity info (not necessarily amyloid and Tau levels) 

A-T- 
> length(tmpASLTaufirst$RID[which(tmpASLTaufirst$AmyloidPos==0 & tmpASLTaufirst$MesialTempTauPos==0 & tmpASLTaufirst$TemporoParietalTauPos==0)]) 
75 

A+T- 
> length(tmpASLTaufirst$RID[which(tmpASLTaufirst$AmyloidPos==1 & tmpASLTaufirst$MesialTempTauPos==0 & tmpASLTaufirst$TemporoParietalTauPos==0)]) 
21 

A+MTT+TPT- 
> length(tmpASLTaufirst$RID[which(tmpASLTaufirst$AmyloidPos==1 & tmpASLTaufirst$MesialTempTauPos==1 & tmpASLTaufirst$TemporoParietalTauPos==0)]) 
15 

A+T+ 
> length(tmpASLTaufirst$RID[which(tmpASLTaufirst$AmyloidPos==1 & tmpASLTaufirst$MesialTempTauPos==1 & tmpASLTaufirst$TemporoParietalTauPos==1)]) 
24 
