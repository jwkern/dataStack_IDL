;Originally the scratch_stack.pro. Modified by Joshua W Kern on 11/27/20 to
;work with COmultiline_fit_JWK.pro

;RESTORE,'/home/jwkern/Research/Exopl/Processed_Data/IRTF/ishell/2021/ABAur_spec_avg_JWK.dat'

;read in CO data
RESTORE, filename='/home/jwkern/Research/Exopl/Software/idl/utils/CO_syn_spec_em/CO_syn_spec_final/CO_molecdat.dat'
RESTORE, filename='/home/jwkern/Research/Exopl/Software/idl/utils/CO_syn_spec_em/CO_syn_spec_final/X12COv10_id.dat'
RESTORE, filename='/home/jwkern/Research/Exopl/Software/idl/utils/CO_syn_spec_em/CO_syn_spec_final/X12COv21_id.dat'
RESTORE, filename='/home/jwkern/Research/Exopl/Software/idl/utils/CO_syn_spec_em/CO_syn_spec_final/X12COv32_id.dat'
RESTORE, filename='/home/jwkern/Research/Exopl/Software/idl/utils/CO_syn_spec_em/CO_syn_spec_final/X12COv43_id.dat'
RESTORE, filename='/home/jwkern/Research/Exopl/Software/idl/utils/CO_syn_spec_em/CO_syn_spec_final/X12COv54_id.dat'
RESTORE, filename='/home/jwkern/Research/Exopl/Software/idl/utils/CO_syn_spec_em/CO_syn_spec_final/X12COv65_id.dat'
RESTORE, filename='/home/jwkern/Research/Exopl/Software/idl/utils/CO_syn_spec_em/CO_syn_spec_final/X12COv76_id.dat'
RESTORE, filename='/home/jwkern/Research/Exopl/Software/idl/utils/CO_syn_spec_em/CO_syn_spec_final/X12COv87_id.dat'
RESTORE, filename='/home/jwkern/Research/Exopl/Software/idl/utils/CO_syn_spec_em/CO_syn_spec_final/X12COv98_id.dat'
RESTORE, filename='/home/jwkern/Research/Exopl/Software/idl/utils/CO_syn_spec_em/CO_syn_spec_final/X13COv10_id.dat'
RESTORE, filename='/home/jwkern/Research/Exopl/Software/idl/utils/CO_syn_spec_em/CO_syn_spec_final/X13COv21_id.dat'
RESTORE, filename='/home/jwkern/Research/Exopl/Software/idl/utils/CO_syn_spec_em/CO_syn_spec_final/X13COv32_id.dat'
RESTORE, filename='/home/jwkern/Research/Exopl/Software/idl/utils/CO_syn_spec_em/CO_syn_spec_final/XC18Ov10_id.dat'



;12CO v1-0
f12_10 = X12CO10(6,*)
Eup12_10 = X12CO10(5,*)
id12_10 = X12CO10_ID(*)

;12CO v2-1
f12_21 = X12CO21(6,*)
Eup12_21 = X12CO21(5,*)
id12_21 = X12CO21_ID(*)

;12CO v3-2
f12_32 = X12CO32(6,*)
Eup12_32 = X12CO32(5,*)
id12_32 = X12CO32_ID(*)

;12CO v4-3
f12_43 = X12CO43(6,*)
Eup12_43 = X12CO43(5,*)
id12_43 = X12CO43_ID(*)

;12CO v5-4
f12_54 = X12CO54(6,*)
Eup12_54 = X12CO54(5,*)
id12_54 = X12CO54_ID(*)

;12CO v6-5
f12_65 = X12CO65(6,*)
Eup12_65 = X12CO65(5,*)
id12_65 = X12CO65_ID(*)

;12CO v7-6
f12_76 = X12CO76(6,*)
Eup12_76 = X12CO76(5,*)
id12_76 = X12CO76_ID(*)

;12CO v8-7
f12_87 = X12CO87(6,*)
Eup12_87 = X12CO87(5,*)
id12_87 = X12CO87_ID(*)

;12CO v9-8 
f12_98 = X12CO98(6,*)
Eup12_98 = X12CO98(5,*)
id12_98 = X12CO98_ID(*)

;13CO v1-0
f13_10 = X13CO10(6,*)
Eup13_10 = X13CO10(5,*)
id13_10 = X13CO10_ID(*)

;13CO v2-1
f13_21 = X13CO21(6,*)
Eup13_21 = X13CO21(5,*)
id13_21 = X13CO21_ID(*)

;13CO v3-2
f13_32 = X13CO32(6,*)
Eup13_32 = X13CO32(5,*)
id13_32 = X13CO32_ID(*)

;C18O v1-0
f18_10 = XC18O10(6,*)
Eup18_10 = XC18O10(5,*)
id18_10 = XC18O10_ID(*)


;Set wvn and amp arrays 
fbig=fbig
rbig=rbig

;Set wavenumber range for spectral lines
lwvn = 1926
hwvn = 2215

;Index a certain range of spectral lines
index=WHERE(f12_10 GT lwvn AND f12_10 LT hwvn,count)
index2=WHERE(fbig GT  lwvn AND fbig LT hwvn,count2)


;Initialize arrays for individual and average lines
r10=FLTARR(601,count)
v10=FLTARR(N_ELEMENTS(fbig(index2)),count)
ravg=FLTARR(N_ELEMENTS(rbig(index2)))
vavg=v10(*,count/2)
vfine=FINDGEN(601) - 300.


;Convert the wvn information into velocity (km/s)
FOR i=0,count-1 DO v10(*,i)=(f12_10(index(i)) - fbig(index2) )*2.9979e5/f12_10(index(i))


;Interpolate the velocity arrays with rbig to create individual profiles of the same size
FOR i=0,count-1 DO r10(*,i)=INTERPOL(rbig(index2),v10(*,i),vfine)


;Why do we redefine count2 and ravg?
count2=601
ravg=FLTARR(601)


;Create the normalized individual line profiles
FOR i=0,count-1 DO BEGIN
   ind_cl=WHERE((vfine GT 75 OR vfine LT 10) AND (r10(*,i) GT 3.3 OR r10(*,i) LT .8),ct)
   if ct NE 0 THEN r10(ind_cl,i)=1./0
   r10(*,i)=r10(*,i)/MEAN(r10(*,i),/NAN)
ENDFOR


;Create, save, and plot average line profile
FOR i=0,count2-1 DO BEGIN
 ind=WHERE(FINITE(r10(i,*)) EQ 1,count3)
 IF count3 NE 0 THEN ravg(i)=TOTAL(r10(i,ind))/count3
ENDFOR
ravg_orig=ravg
CGPLOT,vfine,ravg,PSYM=10,xra=[-100,100]
;CGOPLOT,v,ravg_12,PSYM=10,xra=[-100,100],color='red'

;Make an exit point for the script in case the average profile is bad
Pause


;Initialize include line array
str_ans=''
incl_line=[]
incl_line_rindex=[]
vvar=[]

;Plot the average and individual profiles, check to include line, check to shift the
;line,collect indicies of included lines 
FOR i=0,count-1 DO BEGIN
    str=['CO',id12_10(index(i)),'wvn =',STRSPLIT(STRING(f12_10(index(i))),/EXTRACT),'cm$\up-1$','Eup =',STRSPLIT(STRING(Eup12_10(index(i))),/EXTRACT),'cm$\up-1$']
    CGPLOT,vfine,ravg,PSYM=10,xra=[-100,100],tit=STRJOIN(str,' '),yra=[0,2.5],xtit='Velocity km$\up-1$'
    CGOPLOT,vfine,r10(*,i),PSYM=10,COLOR='red'
   READ,str_ans,PROMPT='Include Line (y/n)?'
   IF str_ans EQ 'n' THEN r10(*,i)=1./0
   IF str_ans NE 'n' THEN BEGIN
      IF N_ELEMENTS(incl_line) EQ 0 THEN incl_line=[index(i)] ELSE incl_line=[incl_line,index(i)]
      IF N_ELEMENTS(incl_line) EQ 0 THEN incl_line_rindex=[i] ELSE incl_line_rindex=[incl_line_rindex,i]
      READ,str_ans,PROMPT='Shift line (y/n)?'
      IF str_ans EQ 'n' THEN BEGIN
              var=0
              vvar=[vvar,var]
      ENDIF
      IF str_ans NE 'n' THEN BEGIN
         redo:
         READ,var,PROMPT='Enter shift in km/s'
         var=FLOAT(var)
            str=['CO',id12_10(index(i)),'wvn =',STRSPLIT(STRING(f12_10(index(i))),/EXTRACT),'cm$\up-1$','Eup =',STRSPLIT(STRING(Eup12_10(index(i))),/EXTRACT),'cm$\up-1$']
            CGPLOT,vfine,ravg,PSYM=10,xra=[-100,100],tit=STRJOIN(str,' '),yra=[0,2.5],xtit='Velocity (km s$\up-1$)'
            CGOPLOT,vfine,INTERPOL(r10(*,i),vfine,vfine-var),color='red',psym=10
            READ,str_ans,PROMPT='Are you happy with the shift (y/n)?'
            IF str_ans EQ 'n' THEN GOTO, redo
            r10(*,i)=INTERPOL(r10(*,i),vfine,vfine-var)
	    vvar=[vvar,var]
     ENDIF
   ENDIF
ENDFOR

count2=601
ravg2=FLTARR(601)

FOR i=0,count2-1 DO BEGIN
	ind=WHERE(FINITE(r10(i,*)) EQ 1,count3)
	IF count3 NE 0 THEN ravg2(i)=TOTAL(r10(i,ind))/count3
	IF ravg2(i) EQ 0 THEN ravg2(i)=!VALUES.F_NAN
ENDFOR

CGPLOT,vfine,ravg2,PSYM=10,xra=[-100,100]

;hiJ: Jlow > 20, lowJ: 4 < Jlow < 20
vvar12_10=vvar
ravg12_10=ravg2         
r12_10=r10
incl_line_12_10=incl_line
incl_line_12_10_rindex=incl_line_rindex

END

