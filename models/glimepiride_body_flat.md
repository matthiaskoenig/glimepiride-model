# model: glimepiride_body
Autogenerated ODE System from SBML with [sbmlutils](https://github.com/matthiaskoenig/sbmlutils).
```
time: [min]
substance: [mmol]
extent: [mmol]
volume: [l]
area: [m^2]
length: [m]
```

## Parameters `p`
```
BW = 75.0  # [kg] body weight [kg]  
COBW = 1.548  # [ml/s/kg] cardiac output per bodyweight [ml/s/kg]  
COHRI = 150.0  # [ml] increase of cardiac output per heartbeat [ml/min*min]  
FQgu = 0.18  # [-] gut fractional tissue blood flow  
FQh = 0.215  # [-] hepatic (venous side) fractional tissue blood flow  
FQki = 0.19  # [-] kidney fractional tissue blood flow  
FQlu = 1.0  # [-] lung fractional tissue blood flow  
FVar = 0.0257  # [l/kg] arterial fractional tissue volume  
FVgu = 0.0171  # [l/kg] gut fractional tissue volume  
FVhv = 0.001  # [l/kg] hepatic venous fractional tissue volume  
FVki = 0.0044  # [l/kg] kidney fractional tissue volume  
FVli = 0.021  # [l/kg] liver fractional tissue volume  
FVlu = 0.0076  # [l/kg] lung fractional tissue volume  
FVpo = 0.001  # [l/kg] portal fractional tissue volume  
FVve = 0.0514  # [l/kg] venous fractional tissue volume  
Fblood = 0.02  # [-] blood fraction of organ volume  
GU__F_gli_abs = 1.0  # [-] fraction absorbed glimepiride  
GU__GLIABS_k = 0.01  # [1/min] rate of glimepiride absorption  
GU__Ka_dis_gli = 2.0  # [1/hr] Ka_dis [1/hr] dissolution glimepiride  
GU__MEXC_k = 0.1  # [1/min] rate feces excretion of M1 and M2  
GU__MREABS_k = 0.05  # [1/min] rate constant for M1 and M2 reabsorption from plasma to intestine  
GU__Mr_gli = 490.616  # [g/mol] Molecular weight glimepiride [g/mole]  
GU__Vstomach = 1.0  # [l] stomach  
GU__f_absorption = 1.0  # [-] scaling factor for absorption rate  
HCT = 0.51  # [-] hematocrit  
HEIGHT = 170.0  # [cm] height [cm]  
HR = 70.0  # [1/min] heart rate [1/min]  
HRrest = 70.0  # [1/min] heart rate [1/min]  
KI__M1EX_k = 0.1  # [1/min] rate of M1 urinary excretion  
KI__M2EX_k = 0.2  # [1/min] rate of M2 urinary excretion  
KI__f_renal_function = 1.0  # [-] parameter for renal function  
LI__GLI2M1_Km_gli = 0.00225  # [mmol/l] Km gli  
LI__GLI2M1_Vmax = 0.01  # [mmol/min/l] rate glimepiride conversion  
LI__GLIIM_k = 0.01  # [1/min] rate glimepiride import  
LI__M12M2_k = 0.001  # [1/min] rate m1 -> m2 conversion  
LI__M1EX_k = 1.0  # [1/min] rate M1 export  
LI__M2EX_k = 1.0  # [1/min] rate M2 export  
LI__f_cyp2c9 = 1.0  # [-] scaling factor CYP2C9 activity  
Mr_gli = 490.616  # [g/mol] Molecular weight gli [g/mole]  
Mr_m1 = 506.62  # [g/mol] Molecular weight m1 [g/mole]  
Mr_m2 = 520.6  # [g/mol] Molecular weight m2 [g/mole]  
Ri_gli = 0.0  # [mg/min] Ri [mg/min] rate of infusion gli  
Ri_m1 = 0.0  # [mg/min] Ri [mg/min] rate of infusion m1  
Ri_m2 = 0.0  # [mg/min] Ri [mg/min] rate of infusion m2  
Vfeces = 1.0  # [l] feces  
Vstomach = 1.0  # [l] stomach  
Vurine = 1.0  # [l] urine  
f_cardiac_function = 1.0  # [-] heart function  
f_cirrhosis = 0.0  # [-] severity of cirrhosis [0, 0.95]  
ti_gli = 10.0  # [s] injection time gli [s]  
ti_m1 = 10.0  # [s] injection time m1 [s]  
ti_m2 = 10.0  # [s] injection time m2 [s]  
```

## Initial conditions `x0`
```
Afeces_gli = 0.0  # [mmol] Glimepiride (feces) in Vfeces  
Afeces_m1 = 0.0  # [mmol] M1 (feces) in Vfeces  
Afeces_m2 = 0.0  # [mmol] M2 (feces) in Vfeces  
Aurine_m1 = 0.0  # [mmol] M1 (urine) in Vurine  
Aurine_m2 = 0.0  # [mmol] M2 (urine) in Vurine  
Car_gli = 0.0  # [mmol/l] glimepiride (arterial blood plasma) in Var  
Car_m1 = 0.0  # [mmol/l] M1 (arterial blood plasma) in Var  
Car_m2 = 0.0  # [mmol/l] M2 (arterial blood plasma) in Var  
Cgu_plasma_gli = 0.0  # [mmol/l] glimepiride (gut plasma) in Vgu_plasma  
Cgu_plasma_m1 = 0.0  # [mmol/l] M1 (gut plasma) in Vgu_plasma  
Cgu_plasma_m2 = 0.0  # [mmol/l] M2 (gut plasma) in Vgu_plasma  
Chv_gli = 0.0  # [mmol/l] glimepiride (hepatic vein plasma) in Vhv  
Chv_m1 = 0.0  # [mmol/l] M1 (hepatic vein plasma) in Vhv  
Chv_m2 = 0.0  # [mmol/l] M2 (hepatic vein plasma) in Vhv  
Cki_plasma_gli = 0.0  # [mmol/l] glimepiride (kidney plasma) in Vki_plasma  
Cki_plasma_m1 = 0.0  # [mmol/l] M1 (kidney plasma) in Vki_plasma  
Cki_plasma_m2 = 0.0  # [mmol/l] M2 (kidney plasma) in Vki_plasma  
Cli_plasma_gli = 0.0  # [mmol/l] glimepiride (liver plasma) in Vli_plasma  
Cli_plasma_m1 = 0.0  # [mmol/l] M1 (liver plasma) in Vli_plasma  
Cli_plasma_m2 = 0.0  # [mmol/l] M2 (liver plasma) in Vli_plasma  
Clu_plasma_gli = 0.0  # [mmol/l] glimepiride (lung plasma) in Vlu_plasma  
Clu_plasma_m1 = 0.0  # [mmol/l] M1 (lung plasma) in Vlu_plasma  
Clu_plasma_m2 = 0.0  # [mmol/l] M2 (lung plasma) in Vlu_plasma  
Cpo_gli = 0.0  # [mmol/l] glimepiride (portal vein plasma) in Vpo  
Cpo_m1 = 0.0  # [mmol/l] M1 (portal vein plasma) in Vpo  
Cpo_m2 = 0.0  # [mmol/l] M2 (portal vein plasma) in Vpo  
Cre_plasma_gli = 0.0  # [mmol/l] glimepiride (rest plasma) in Vre_plasma  
Cre_plasma_m1 = 0.0  # [mmol/l] M1 (rest plasma) in Vre_plasma  
Cre_plasma_m2 = 0.0  # [mmol/l] M2 (rest plasma) in Vre_plasma  
Cve_gli = 0.0  # [mmol/l] glimepiride (venous blood plasma) in Vve  
Cve_m1 = 0.0  # [mmol/l] M1 (venous blood plasma) in Vve  
Cve_m2 = 0.0  # [mmol/l] M2 (venous blood plasma) in Vve  
GU__gli_lumen = 0.0  # [mmol/l] glimepiride (intestinal lumen) in Vgu  
GU__gli_stomach = 0.0  # [mmol] glimepiride (stomach) in GU__Vstomach  
GU__m1_lumen = 0.0  # [mmol/l] M1 (intestinal lumen) in Vgu  
GU__m2_lumen = 0.0  # [mmol/l] M2 (intestinal lumen) in Vgu  
IVDOSE_gli = 0.0  # [mg] IV bolus dose gli [mg]  
IVDOSE_m1 = 0.0  # [mg] IV bolus dose m1 [mg]  
IVDOSE_m2 = 0.0  # [mg] IV bolus dose m2 [mg]  
LI__gli = 0.0  # [mmol/l] glimepiride (liver) in Vli_tissue  
LI__m1 = 0.0  # [mmol/l] M1 (liver) in Vli_tissue  
LI__m2 = 0.0  # [mmol/l] M2 (liver) in Vli_tissue  
PODOSE_gli = 0.0  # [mg] rate of dose dissolution  
cum_dose_gli = 0.0  # [mg] Cumulative dose due to infusion gli  
cum_dose_m1 = 0.0  # [mg] Cumulative dose due to infusion m1  
cum_dose_m2 = 0.0  # [mg] Cumulative dose due to infusion m2  
```

## ODE system
```
# y
Afeces_m1_m2 = Afeces_m1 + Afeces_m2  # [mmol] Sum of M1 and M2 feces  
Aurine_m1_m2 = Aurine_m1 + Aurine_m2  # [mmol] Sum of M1 and M2 urine  
BSA = 0.024265 * (BW / 1)**0.5378 * (HEIGHT / 1)**0.3964  # [m^2] body surface area [m^2]  
CO = f_cardiac_function * BW * COBW + (HR - HRrest) * COHRI / 60  # [ml/s] cardiac output [ml/s]  
Cve_m1_m2 = Cve_m1 + Cve_m2  # [mmol/l] Sum of M1 and M2 plasma  
FQre = 1 - (FQki + FQh)  # [-] rest of body fractional tissue blood flow  
FVre = 1 - (FVgu + FVki + FVli + FVlu + FVve + FVar)  # [l/kg] rest of body fractional tissue volume  
GU__dissolution_gli = (GU__Ka_dis_gli / 60) * PODOSE_gli / GU__Mr_gli  # [mmol/min] dissolution glimepiride  
GU__mtot_feces = Afeces_m1 + Afeces_m2  # [mmol] total fecal metabolites (M1 + M2)  
Ki_gli = (0.693 / ti_gli) * 60  # [1/min] injection rate IV  
Ki_m1 = (0.693 / ti_m1) * 60  # [1/min] injection rate IV  
Ki_m2 = (0.693 / ti_m2) * 60  # [1/min] injection rate IV  
Var = BW * FVar - (FVar / (FVar + FVve)) * BW * Fblood * (1 - FVve - FVar)  # [l] arterial blood  
Vgu = BW * FVgu  # [l] gut  
Vhv = (1 - HCT) * (BW * FVhv - (FVhv / (FVar + FVve + FVpo + FVhv)) * BW * Fblood * (1 - (FVar + FVve + FVpo + FVhv)))  # [l] hepatic venous plasma  
Vki = BW * FVki  # [l] kidney  
Vli = BW * FVli  # [l] liver  
Vlu = BW * FVlu  # [l] lung  
Vpo = (1 - HCT) * (BW * FVpo - (FVpo / (FVar + FVve + FVpo + FVhv)) * BW * Fblood * (1 - (FVar + FVve + FVpo + FVhv)))  # [l] portal plasma  
Vve = BW * FVve - (FVve / (FVar + FVve)) * BW * Fblood * (1 - FVve - FVar)  # [l] venous blood  
Xfeces_m1 = Afeces_m1 * Mr_m1  # [mg] M1 amount (feces)  
Xfeces_m2 = Afeces_m2 * Mr_m2  # [mg] M2 amount (feces)  
Xurine_m1 = Aurine_m1 * Mr_m1  # [mg] M1 amount (urine) [mg]  
Xurine_m2 = Aurine_m2 * Mr_m2  # [mg] M2 amount (urine) [mg]  
f_shunts = f_cirrhosis  # [-] fraction of portal venous blood shunted by the liver  
f_tissue_loss = f_cirrhosis  # [-] fraction of lost parenchymal liver volume  
Aar_gli = Car_gli * Var  # [mmol] glimepiride amount (arterial blood) [mmole]  
Aar_m1 = Car_m1 * Var  # [mmol] M1 amount (arterial blood) [mmole]  
Aar_m2 = Car_m2 * Var  # [mmol] M2 amount (arterial blood) [mmole]  
Ahv_gli = Chv_gli * Vhv  # [mmol] glimepiride amount (hepatic vein) [mmole]  
Ahv_m1 = Chv_m1 * Vhv  # [mmol] M1 amount (hepatic vein) [mmole]  
Ahv_m2 = Chv_m2 * Vhv  # [mmol] M2 amount (hepatic vein) [mmole]  
Apo_gli = Cpo_gli * Vpo  # [mmol] glimepiride amount (portal vein) [mmole]  
Apo_m1 = Cpo_m1 * Vpo  # [mmol] M1 amount (portal vein) [mmole]  
Apo_m2 = Cpo_m2 * Vpo  # [mmol] M2 amount (portal vein) [mmole]  
Ave_gli = Cve_gli * Vve  # [mmol] glimepiride amount (venous blood) [mmole]  
Ave_m1 = Cve_m1 * Vve  # [mmol] M1 amount (venous blood) [mmole]  
Ave_m2 = Cve_m2 * Vve  # [mmol] M2 amount (venous blood) [mmole]  
GU__M1EXC = GU__MEXC_k * GU__m1_lumen * Vgu  # [mmol/min] excretion M1 (feces)  
GU__M1REABS = GU__MREABS_k * Cgu_plasma_m1 * Vgu  # [mmol/min] reabsorption M1 (lumen)  
GU__M2EXC = GU__MEXC_k * GU__m2_lumen * Vgu  # [mmol/min] excretion M2 (feces)  
GU__M2REABS = GU__MREABS_k * Cgu_plasma_m2 * Vgu  # [mmol/min] reabsorption M2 (lumen)  
GU__absorption = GU__f_absorption * GU__GLIABS_k * Vgu * GU__gli_lumen  # [mmol/min] absorption glimepiride  
QC = (CO / 1000) * 60  # [l/min] cardiac output [L/hr]  
Vgu_plasma = Vgu * Fblood * (1 - HCT)  # [l] plasma volume of gut  
Vgu_tissue = Vgu * (1 - Fblood)  # [l] tissue volume of gut  
Vki_plasma = Vki * Fblood * (1 - HCT)  # [l] plasma volume of kidney  
Vki_tissue = Vki * (1 - Fblood)  # [l] tissue volume of kidney  
Vli_plasma = Vli * Fblood * (1 - HCT)  # [l] plasma volume of liver  
Vli_tissue = Vli * (1 - f_tissue_loss) * (1 - Fblood)  # [l] tissue volume of liver  
Vlu_plasma = Vlu * Fblood * (1 - HCT)  # [l] plasma volume of lung  
Vlu_tissue = Vlu * (1 - Fblood)  # [l] tissue volume of lung  
Vre = BW * FVre  # [l] rest of body  
iv_gli = Ki_gli * IVDOSE_gli / Mr_gli  # [mmol/min] iv glimepiride  
iv_m1 = Ki_m1 * IVDOSE_m1 / Mr_m1  # [mmol/min] iv M1  
iv_m2 = Ki_m2 * IVDOSE_m2 / Mr_m2  # [mmol/min] iv M2  
Agu_plasma_gli = Cgu_plasma_gli * Vgu_plasma  # [mmol] glimepiride amount (gut) [mmole]  
Agu_plasma_m1 = Cgu_plasma_m1 * Vgu_plasma  # [mmol] M1 amount (gut) [mmole]  
Agu_plasma_m2 = Cgu_plasma_m2 * Vgu_plasma  # [mmol] M2 amount (gut) [mmole]  
Aki_plasma_gli = Cki_plasma_gli * Vki_plasma  # [mmol] glimepiride amount (kidney) [mmole]  
Aki_plasma_m1 = Cki_plasma_m1 * Vki_plasma  # [mmol] M1 amount (kidney) [mmole]  
Aki_plasma_m2 = Cki_plasma_m2 * Vki_plasma  # [mmol] M2 amount (kidney) [mmole]  
Ali_plasma_gli = Cli_plasma_gli * Vli_plasma  # [mmol] glimepiride amount (liver) [mmole]  
Ali_plasma_m1 = Cli_plasma_m1 * Vli_plasma  # [mmol] M1 amount (liver) [mmole]  
Ali_plasma_m2 = Cli_plasma_m2 * Vli_plasma  # [mmol] M2 amount (liver) [mmole]  
Alu_plasma_gli = Clu_plasma_gli * Vlu_plasma  # [mmol] glimepiride amount (lung) [mmole]  
Alu_plasma_m1 = Clu_plasma_m1 * Vlu_plasma  # [mmol] M1 amount (lung) [mmole]  
Alu_plasma_m2 = Clu_plasma_m2 * Vlu_plasma  # [mmol] M2 amount (lung) [mmole]  
GU__GLIABS = GU__F_gli_abs * GU__absorption  # [mmol/min] absorption glimepiride  
KI__M1EX = KI__f_renal_function * Vki_tissue * KI__M1EX_k * Cki_plasma_m1  # [mmol/min] M1 excretion (M1EX)  
KI__M2EX = KI__f_renal_function * Vki_tissue * KI__M2EX_k * Cki_plasma_m2  # [mmol/min] M2 excretion (M2EX)  
LI__GLI2M1 = LI__f_cyp2c9 * LI__GLI2M1_Vmax * Vli_tissue * LI__gli / (LI__gli + LI__GLI2M1_Km_gli)  # [mmol/min] glimepiride conversion to M1 (GLI2M1) CYP2C9  
LI__GLIIM = LI__GLIIM_k * Vli_tissue * (Cli_plasma_gli - LI__gli)  # [mmol/min] glimepiride import (GLIIM)  
LI__M12M2 = LI__M12M2_k * Vli_tissue * LI__m1  # [mmol/min] M1 conversion to M2 (M12M2)  
LI__M1EX = LI__M1EX_k * Vli_tissue * (LI__m1 - Cli_plasma_m1)  # [mmol/min] M1 export (M1EX)  
LI__M2EX = LI__M2EX_k * Vli_tissue * (LI__m2 - Cli_plasma_m2)  # [mmol/min] M2 export (M2EX)  
Mar_gli = (Aar_gli / Var) * Mr_gli  # [mg/l] glimepiride concentration (arterial blood) [mg/l]  
Mar_m1 = (Aar_m1 / Var) * Mr_m1  # [mg/l] M1 concentration (arterial blood) [mg/l]  
Mar_m2 = (Aar_m2 / Var) * Mr_m2  # [mg/l] M2 concentration (arterial blood) [mg/l]  
Mhv_gli = (Ahv_gli / Vhv) * Mr_gli  # [mg/l] glimepiride concentration (hepatic vein) [mg/l]  
Mhv_m1 = (Ahv_m1 / Vhv) * Mr_m1  # [mg/l] M1 concentration (hepatic vein) [mg/l]  
Mhv_m2 = (Ahv_m2 / Vhv) * Mr_m2  # [mg/l] M2 concentration (hepatic vein) [mg/l]  
Mpo_gli = (Apo_gli / Vpo) * Mr_gli  # [mg/l] glimepiride concentration (portal vein) [mg/l]  
Mpo_m1 = (Apo_m1 / Vpo) * Mr_m1  # [mg/l] M1 concentration (portal vein) [mg/l]  
Mpo_m2 = (Apo_m2 / Vpo) * Mr_m2  # [mg/l] M2 concentration (portal vein) [mg/l]  
Mve_gli = (Ave_gli / Vve) * Mr_gli  # [mg/l] glimepiride concentration (venous blood) [mg/l]  
Mve_m1 = (Ave_m1 / Vve) * Mr_m1  # [mg/l] M1 concentration (venous blood) [mg/l]  
Mve_m2 = (Ave_m2 / Vve) * Mr_m2  # [mg/l] M2 concentration (venous blood) [mg/l]  
Qgu = QC * FQgu  # [l/min] gut blood flow  
Qh = QC * FQh  # [l/min] hepatic (venous side) blood flow  
Qki = QC * FQki  # [l/min] kidney blood flow  
Qlu = QC * FQlu  # [l/min] lung blood flow  
Qre = QC * FQre  # [l/min] rest of body blood flow  
Vre_plasma = Vre * Fblood * (1 - HCT)  # [l] plasma volume of rest  
Vre_tissue = Vre * (1 - Fblood)  # [l] tissue volume of rest  
Xar_gli = Aar_gli * Mr_gli  # [mg] glimepiride amount (arterial blood) [mg]  
Xar_m1 = Aar_m1 * Mr_m1  # [mg] M1 amount (arterial blood) [mg]  
Xar_m2 = Aar_m2 * Mr_m2  # [mg] M2 amount (arterial blood) [mg]  
Xhv_gli = Ahv_gli * Mr_gli  # [mg] glimepiride amount (hepatic vein) [mg]  
Xhv_m1 = Ahv_m1 * Mr_m1  # [mg] M1 amount (hepatic vein) [mg]  
Xhv_m2 = Ahv_m2 * Mr_m2  # [mg] M2 amount (hepatic vein) [mg]  
Xpo_gli = Apo_gli * Mr_gli  # [mg] glimepiride amount (portal vein) [mg]  
Xpo_m1 = Apo_m1 * Mr_m1  # [mg] M1 amount (portal vein) [mg]  
Xpo_m2 = Apo_m2 * Mr_m2  # [mg] M2 amount (portal vein) [mg]  
Xve_gli = Ave_gli * Mr_gli  # [mg] glimepiride amount (venous blood) [mg]  
Xve_m1 = Ave_m1 * Mr_m1  # [mg] M1 amount (venous blood) [mg]  
Xve_m2 = Ave_m2 * Mr_m2  # [mg] M2 amount (venous blood) [mg]  
Are_plasma_gli = Cre_plasma_gli * Vre_plasma  # [mmol] glimepiride amount (rest) [mmole]  
Are_plasma_m1 = Cre_plasma_m1 * Vre_plasma  # [mmol] M1 amount (rest) [mmole]  
Are_plasma_m2 = Cre_plasma_m2 * Vre_plasma  # [mmol] M2 amount (rest) [mmole]  
Flow_ar_gu_gli = Qgu * Car_gli  # [mmol/min] inflow gut glimepiride  
Flow_ar_gu_m1 = Qgu * Car_m1  # [mmol/min] inflow gut M1  
Flow_ar_gu_m2 = Qgu * Car_m2  # [mmol/min] inflow gut M2  
Flow_ar_ki_gli = Qki * Car_gli  # [mmol/min] inflow kidney glimepiride  
Flow_ar_ki_m1 = Qki * Car_m1  # [mmol/min] inflow kidney M1  
Flow_ar_ki_m2 = Qki * Car_m2  # [mmol/min] inflow kidney M2  
Flow_ar_re_gli = Qre * Car_gli  # [mmol/min] inflow rest glimepiride  
Flow_ar_re_m1 = Qre * Car_m1  # [mmol/min] inflow rest M1  
Flow_ar_re_m2 = Qre * Car_m2  # [mmol/min] inflow rest M2  
Flow_gu_po_gli = Qgu * Cgu_plasma_gli  # [mmol/min] outflow gut glimepiride  
Flow_gu_po_m1 = Qgu * Cgu_plasma_m1  # [mmol/min] outflow gut M1  
Flow_gu_po_m2 = Qgu * Cgu_plasma_m2  # [mmol/min] outflow gut M2  
Flow_hv_ve_gli = Qh * Chv_gli  # [mmol/min] outflow hepatic vein glimepiride  
Flow_hv_ve_m1 = Qh * Chv_m1  # [mmol/min] outflow hepatic vein M1  
Flow_hv_ve_m2 = Qh * Chv_m2  # [mmol/min] outflow hepatic vein M2  
Flow_ki_ve_gli = Qki * Cki_plasma_gli  # [mmol/min] outflow kidney glimepiride  
Flow_ki_ve_m1 = Qki * Cki_plasma_m1  # [mmol/min] outflow kidney M1  
Flow_ki_ve_m2 = Qki * Cki_plasma_m2  # [mmol/min] outflow kidney M2  
Flow_lu_ar_gli = Qlu * Clu_plasma_gli  # [mmol/min] outflow lung glimepiride  
Flow_lu_ar_m1 = Qlu * Clu_plasma_m1  # [mmol/min] outflow lung M1  
Flow_lu_ar_m2 = Qlu * Clu_plasma_m2  # [mmol/min] outflow lung M2  
Flow_re_ve_gli = Qre * Cre_plasma_gli  # [mmol/min] outflow rest glimepiride  
Flow_re_ve_m1 = Qre * Cre_plasma_m1  # [mmol/min] outflow rest M1  
Flow_re_ve_m2 = Qre * Cre_plasma_m2  # [mmol/min] outflow rest M2  
Flow_ve_lu_gli = Qlu * Cve_gli  # [mmol/min] inflow lung glimepiride  
Flow_ve_lu_m1 = Qlu * Cve_m1  # [mmol/min] inflow lung M1  
Flow_ve_lu_m2 = Qlu * Cve_m2  # [mmol/min] inflow lung M2  
Mgu_plasma_gli = (Agu_plasma_gli / Vgu_plasma) * Mr_gli  # [mg/l] glimepiride concentration (gut) [mg/l]  
Mgu_plasma_m1 = (Agu_plasma_m1 / Vgu_plasma) * Mr_m1  # [mg/l] M1 concentration (gut) [mg/l]  
Mgu_plasma_m2 = (Agu_plasma_m2 / Vgu_plasma) * Mr_m2  # [mg/l] M2 concentration (gut) [mg/l]  
Mki_plasma_gli = (Aki_plasma_gli / Vki_plasma) * Mr_gli  # [mg/l] glimepiride concentration (kidney) [mg/l]  
Mki_plasma_m1 = (Aki_plasma_m1 / Vki_plasma) * Mr_m1  # [mg/l] M1 concentration (kidney) [mg/l]  
Mki_plasma_m2 = (Aki_plasma_m2 / Vki_plasma) * Mr_m2  # [mg/l] M2 concentration (kidney) [mg/l]  
Mli_plasma_gli = (Ali_plasma_gli / Vli_plasma) * Mr_gli  # [mg/l] glimepiride concentration (liver) [mg/l]  
Mli_plasma_m1 = (Ali_plasma_m1 / Vli_plasma) * Mr_m1  # [mg/l] M1 concentration (liver) [mg/l]  
Mli_plasma_m2 = (Ali_plasma_m2 / Vli_plasma) * Mr_m2  # [mg/l] M2 concentration (liver) [mg/l]  
Mlu_plasma_gli = (Alu_plasma_gli / Vlu_plasma) * Mr_gli  # [mg/l] glimepiride concentration (lung) [mg/l]  
Mlu_plasma_m1 = (Alu_plasma_m1 / Vlu_plasma) * Mr_m1  # [mg/l] M1 concentration (lung) [mg/l]  
Mlu_plasma_m2 = (Alu_plasma_m2 / Vlu_plasma) * Mr_m2  # [mg/l] M2 concentration (lung) [mg/l]  
Qha = Qh - Qgu  # [l/min] hepatic artery blood flow  
Qpo = Qgu  # [l/min] portal blood flow  
Xgu_plasma_gli = Agu_plasma_gli * Mr_gli  # [mg] glimepiride amount (gut) [mg]  
Xgu_plasma_m1 = Agu_plasma_m1 * Mr_m1  # [mg] M1 amount (gut) [mg]  
Xgu_plasma_m2 = Agu_plasma_m2 * Mr_m2  # [mg] M2 amount (gut) [mg]  
Xki_plasma_gli = Aki_plasma_gli * Mr_gli  # [mg] glimepiride amount (kidney) [mg]  
Xki_plasma_m1 = Aki_plasma_m1 * Mr_m1  # [mg] M1 amount (kidney) [mg]  
Xki_plasma_m2 = Aki_plasma_m2 * Mr_m2  # [mg] M2 amount (kidney) [mg]  
Xli_plasma_gli = Ali_plasma_gli * Mr_gli  # [mg] glimepiride amount (liver) [mg]  
Xli_plasma_m1 = Ali_plasma_m1 * Mr_m1  # [mg] M1 amount (liver) [mg]  
Xli_plasma_m2 = Ali_plasma_m2 * Mr_m2  # [mg] M2 amount (liver) [mg]  
Xlu_plasma_gli = Alu_plasma_gli * Mr_gli  # [mg] glimepiride amount (lung) [mg]  
Xlu_plasma_m1 = Alu_plasma_m1 * Mr_m1  # [mg] M1 amount (lung) [mg]  
Xlu_plasma_m2 = Alu_plasma_m2 * Mr_m2  # [mg] M2 amount (lung) [mg]  
Flow_arli_hv_gli = f_shunts * Qha * Car_gli  # [mmol/min] flow arterial shunts  
Flow_arli_hv_m1 = f_shunts * Qha * Car_m1  # [mmol/min] flow arterial shunts  
Flow_arli_hv_m2 = f_shunts * Qha * Car_m2  # [mmol/min] flow arterial shunts  
Flow_arli_li_gli = (1 - f_shunts) * Qha * Car_gli  # [mmol/min] arterial inflow liver glimepiride  
Flow_arli_li_m1 = (1 - f_shunts) * Qha * Car_m1  # [mmol/min] arterial inflow liver M1  
Flow_arli_li_m2 = (1 - f_shunts) * Qha * Car_m2  # [mmol/min] arterial inflow liver M2  
Flow_li_hv_gli = (1 - f_shunts) * (Qpo + Qha) * Cli_plasma_gli  # [mmol/min] outflow liver glimepiride  
Flow_li_hv_m1 = (1 - f_shunts) * (Qpo + Qha) * Cli_plasma_m1  # [mmol/min] outflow liver M1  
Flow_li_hv_m2 = (1 - f_shunts) * (Qpo + Qha) * Cli_plasma_m2  # [mmol/min] outflow liver M2  
Flow_po_hv_gli = f_shunts * Qpo * Cpo_gli  # [mmol/min] flow portal shunts  
Flow_po_hv_m1 = f_shunts * Qpo * Cpo_m1  # [mmol/min] flow portal shunts  
Flow_po_hv_m2 = f_shunts * Qpo * Cpo_m2  # [mmol/min] flow portal shunts  
Flow_po_li_gli = (1 - f_shunts) * Qpo * Cpo_gli  # [mmol/min] outflow po glimepiride  
Flow_po_li_m1 = (1 - f_shunts) * Qpo * Cpo_m1  # [mmol/min] outflow po M1  
Flow_po_li_m2 = (1 - f_shunts) * Qpo * Cpo_m2  # [mmol/min] outflow po M2  
Mre_plasma_gli = (Are_plasma_gli / Vre_plasma) * Mr_gli  # [mg/l] glimepiride concentration (rest) [mg/l]  
Mre_plasma_m1 = (Are_plasma_m1 / Vre_plasma) * Mr_m1  # [mg/l] M1 concentration (rest) [mg/l]  
Mre_plasma_m2 = (Are_plasma_m2 / Vre_plasma) * Mr_m2  # [mg/l] M2 concentration (rest) [mg/l]  
Xre_plasma_gli = Are_plasma_gli * Mr_gli  # [mg] glimepiride amount (rest) [mg]  
Xre_plasma_m1 = Are_plasma_m1 * Mr_m1  # [mg] M1 amount (rest) [mg]  
Xre_plasma_m2 = Are_plasma_m2 * Mr_m2  # [mg] M2 amount (rest) [mg]  

# odes
d Afeces_gli/dt = 0  # [mmol/min] Glimepiride (feces)  
d Afeces_m1/dt = GU__M1EXC  # [mmol/min] M1 (feces)  
d Afeces_m2/dt = GU__M2EXC  # [mmol/min] M2 (feces)  
d Aurine_m1/dt = KI__M1EX  # [mmol/min] M1 (urine)  
d Aurine_m2/dt = KI__M2EX  # [mmol/min] M2 (urine)  
d Car_gli/dt = (-Flow_ar_ki_gli / Var - Flow_arli_li_gli / Var - Flow_arli_hv_gli / Var) + Flow_lu_ar_gli / Var - Flow_ar_gu_gli / Var - Flow_ar_re_gli / Var  # [mmol/l/min] glimepiride (arterial blood plasma)  
d Car_m1/dt = (-Flow_ar_ki_m1 / Var - Flow_arli_li_m1 / Var - Flow_arli_hv_m1 / Var) + Flow_lu_ar_m1 / Var - Flow_ar_gu_m1 / Var - Flow_ar_re_m1 / Var  # [mmol/l/min] M1 (arterial blood plasma)  
d Car_m2/dt = (-Flow_ar_ki_m2 / Var - Flow_arli_li_m2 / Var - Flow_arli_hv_m2 / Var) + Flow_lu_ar_m2 / Var - Flow_ar_gu_m2 / Var - Flow_ar_re_m2 / Var  # [mmol/l/min] M2 (arterial blood plasma)  
d Cgu_plasma_gli/dt = (Flow_ar_gu_gli / Vgu_plasma - Flow_gu_po_gli / Vgu_plasma) + GU__GLIABS / Vgu_plasma  # [mmol/l/min] glimepiride (gut plasma)  
d Cgu_plasma_m1/dt = Flow_ar_gu_m1 / Vgu_plasma - Flow_gu_po_m1 / Vgu_plasma - GU__M1REABS / Vgu_plasma  # [mmol/l/min] M1 (gut plasma)  
d Cgu_plasma_m2/dt = Flow_ar_gu_m2 / Vgu_plasma - Flow_gu_po_m2 / Vgu_plasma - GU__M2REABS / Vgu_plasma  # [mmol/l/min] M2 (gut plasma)  
d Chv_gli/dt = Flow_arli_hv_gli / Vhv + Flow_po_hv_gli / Vhv + Flow_li_hv_gli / Vhv - Flow_hv_ve_gli / Vhv  # [mmol/l/min] glimepiride (hepatic vein plasma)  
d Chv_m1/dt = Flow_arli_hv_m1 / Vhv + Flow_po_hv_m1 / Vhv + Flow_li_hv_m1 / Vhv - Flow_hv_ve_m1 / Vhv  # [mmol/l/min] M1 (hepatic vein plasma)  
d Chv_m2/dt = Flow_arli_hv_m2 / Vhv + Flow_po_hv_m2 / Vhv + Flow_li_hv_m2 / Vhv - Flow_hv_ve_m2 / Vhv  # [mmol/l/min] M2 (hepatic vein plasma)  
d Cki_plasma_gli/dt = Flow_ar_ki_gli / Vki_plasma - Flow_ki_ve_gli / Vki_plasma  # [mmol/l/min] glimepiride (kidney plasma)  
d Cki_plasma_m1/dt = Flow_ar_ki_m1 / Vki_plasma - Flow_ki_ve_m1 / Vki_plasma - KI__M1EX / Vki_plasma  # [mmol/l/min] M1 (kidney plasma)  
d Cki_plasma_m2/dt = Flow_ar_ki_m2 / Vki_plasma - Flow_ki_ve_m2 / Vki_plasma - KI__M2EX / Vki_plasma  # [mmol/l/min] M2 (kidney plasma)  
d Cli_plasma_gli/dt = Flow_arli_li_gli / Vli_plasma + Flow_po_li_gli / Vli_plasma - Flow_li_hv_gli / Vli_plasma - LI__GLIIM / Vli_plasma  # [mmol/l/min] glimepiride (liver plasma)  
d Cli_plasma_m1/dt = (Flow_arli_li_m1 / Vli_plasma + Flow_po_li_m1 / Vli_plasma - Flow_li_hv_m1 / Vli_plasma) + LI__M1EX / Vli_plasma  # [mmol/l/min] M1 (liver plasma)  
d Cli_plasma_m2/dt = (Flow_arli_li_m2 / Vli_plasma + Flow_po_li_m2 / Vli_plasma - Flow_li_hv_m2 / Vli_plasma) + LI__M2EX / Vli_plasma  # [mmol/l/min] M2 (liver plasma)  
d Clu_plasma_gli/dt = Flow_ve_lu_gli / Vlu_plasma - Flow_lu_ar_gli / Vlu_plasma  # [mmol/l/min] glimepiride (lung plasma)  
d Clu_plasma_m1/dt = Flow_ve_lu_m1 / Vlu_plasma - Flow_lu_ar_m1 / Vlu_plasma  # [mmol/l/min] M1 (lung plasma)  
d Clu_plasma_m2/dt = Flow_ve_lu_m2 / Vlu_plasma - Flow_lu_ar_m2 / Vlu_plasma  # [mmol/l/min] M2 (lung plasma)  
d Cpo_gli/dt = (-Flow_po_li_gli / Vpo - Flow_po_hv_gli / Vpo) + Flow_gu_po_gli / Vpo  # [mmol/l/min] glimepiride (portal vein plasma)  
d Cpo_m1/dt = (-Flow_po_li_m1 / Vpo - Flow_po_hv_m1 / Vpo) + Flow_gu_po_m1 / Vpo  # [mmol/l/min] M1 (portal vein plasma)  
d Cpo_m2/dt = (-Flow_po_li_m2 / Vpo - Flow_po_hv_m2 / Vpo) + Flow_gu_po_m2 / Vpo  # [mmol/l/min] M2 (portal vein plasma)  
d Cre_plasma_gli/dt = Flow_ar_re_gli / Vre_plasma - Flow_re_ve_gli / Vre_plasma  # [mmol/l/min] glimepiride (rest plasma)  
d Cre_plasma_m1/dt = Flow_ar_re_m1 / Vre_plasma - Flow_re_ve_m1 / Vre_plasma  # [mmol/l/min] M1 (rest plasma)  
d Cre_plasma_m2/dt = Flow_ar_re_m2 / Vre_plasma - Flow_re_ve_m2 / Vre_plasma  # [mmol/l/min] M2 (rest plasma)  
d Cve_gli/dt = (iv_gli / Vve + Flow_ki_ve_gli / Vve + Flow_hv_ve_gli / Vve - Flow_ve_lu_gli / Vve) + Flow_re_ve_gli / Vve  # [mmol/l/min] glimepiride (venous blood plasma)  
d Cve_m1/dt = (iv_m1 / Vve + Flow_ki_ve_m1 / Vve + Flow_hv_ve_m1 / Vve - Flow_ve_lu_m1 / Vve) + Flow_re_ve_m1 / Vve  # [mmol/l/min] M1 (venous blood plasma)  
d Cve_m2/dt = (iv_m2 / Vve + Flow_ki_ve_m2 / Vve + Flow_hv_ve_m2 / Vve - Flow_ve_lu_m2 / Vve) + Flow_re_ve_m2 / Vve  # [mmol/l/min] M2 (venous blood plasma)  
d GU__gli_lumen/dt = -GU__GLIABS / Vgu + GU__dissolution_gli / Vgu  # [mmol/l/min] glimepiride (intestinal lumen)  
d GU__gli_stomach/dt = 0  # [mmol/min] glimepiride (stomach)  
d GU__m1_lumen/dt = GU__M1REABS / Vgu - GU__M1EXC / Vgu  # [mmol/l/min] M1 (intestinal lumen)  
d GU__m2_lumen/dt = GU__M2REABS / Vgu - GU__M2EXC / Vgu  # [mmol/l/min] M2 (intestinal lumen)  
d IVDOSE_gli/dt = -iv_gli * Mr_gli + Ri_gli  # [mg/min] IV bolus dose gli [mg]  
d IVDOSE_m1/dt = -iv_m1 * Mr_m1 + Ri_m1  # [mg/min] IV bolus dose m1 [mg]  
d IVDOSE_m2/dt = -iv_m2 * Mr_m2 + Ri_m2  # [mg/min] IV bolus dose m2 [mg]  
d LI__gli/dt = LI__GLIIM / Vli_tissue - LI__GLI2M1 / Vli_tissue  # [mmol/l/min] glimepiride (liver)  
d LI__m1/dt = LI__GLI2M1 / Vli_tissue - LI__M1EX / Vli_tissue - LI__M12M2 / Vli_tissue  # [mmol/l/min] M1 (liver)  
d LI__m2/dt = LI__M12M2 / Vli_tissue - LI__M2EX / Vli_tissue  # [mmol/l/min] M2 (liver)  
d PODOSE_gli/dt = -GU__dissolution_gli * GU__Mr_gli  # [mg/min] rate of dose dissolution  
d cum_dose_gli/dt = Ri_gli  # [mg/min] Cumulative dose due to infusion gli  
d cum_dose_m1/dt = Ri_m1  # [mg/min] Cumulative dose due to infusion m1  
d cum_dose_m2/dt = Ri_m2  # [mg/min] Cumulative dose due to infusion m2  
```