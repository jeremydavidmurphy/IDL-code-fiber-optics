PRO FRDonlyL

; COMPILE: fratioF.pro

; The following lists are needed:
; fiberdir.list
; wave.list
; fiber.list
; fiberD.list

namout = 'CO_bend' ; The output name for the curve-of-growth. This is just a ascii
                       ; plot of ee as a function of radius. These are saved unless
                       ; namout='none'. The fiber name is added, so, if the fiber name was
                       ; S18, the output is
                       ; namout_S18_CP#.txt
                       ; where CP# is the camera position.

plotting = 'on'     ; set to anything other than 'on' to turn off the plotting
frplot = 'on'       ; set to turn on the f-ration plotting

nstep = 7            ; this is the number of camera positions.
cstepsize = 4.0      ; the camera step size, in mm
psize = 0.000024     ; the FLI camera pixel size, in mm. (change to 0.000024
                     ; for binned data)
thresh = 2000.0      ; the threshold above which the code fits an ellipse for 
                     ; determination of the center of the spot. THIS IS NOT 
                     ; THE THRESHHOLD FOR FRATIOF. That level is set
                     ; within fratioF itself and is based on a percentage
                     ; cut of the high values of the frame.
winsize = 3.0        ; set this to 1.0 for apogee data and 6.0 or 7.0 for FLI data.
                     ; (this is used only when plotting is turned on)
slow = 'off'         ; setting this to 'on' will insert pauses at various points.
writebase = 'off'    ; will writeout the median combined base frame that was
                     ; used as the final BL value. This frame has had the dark
                     ; adjusted and subtracted
radstep = 0.5        ; this sets the resolution of the radius step size in 
                     ; fratioF.pro. This step size leads to EE resolution of
                     ; about 0.1% for the 3056x3056 FLI chip size. 0.25 is
                     ; as low as you need. Set to a higher number for faster
                     ; (yet less precise) results. 1.0 gives EE resolution
                     ; of about 0.3%. This leads to doubling up of output
                     ; f-ratio values at certain EEs (which means lower
                     ; doesn't do much for you.
starting_radius = 10 ; This sets the beginning radius for the FRD tests.
                     ; Set this higher if you're just interested in
                     ; the outer radius and lower when wanting to have
                     ; good COGs.

;*************************************************************************
;*************************************************************************
; NOTE: THE SPOT SIZES MUST BE ORDERED FROM LARGEST TO SMALLEST. IF
; NOT, FRATIOF.PRO WILL CRASH.

eenergy = 0.95 ; How far the fratioF routine integrates up to. DON'T GO CHANGING THIS! IT FUCKS OTHER THINGS UP THAT I DON'T WANT TO TRY AND FIX!

cd, current=homedir ;set your current direcory as home.

readcol,'fiberdir.list',f='a',Fdir

wavelist = 'wave.list'
fibL = 'fiber.list'
fibDL = 'fiberD.list'

readcol,silent=1,wavelist,format='i',wavelength
readcol,silent=1,fibL,format='a',frames
readcol,silent=1,fibDL,format='a',Dframes

; A sample frame is read in so that the size of each frame can be
; determined. The data arrays are then determined. THIS ASSUMES ALL
; FRAMES ARE OF THE SAME SIZE (same number of rows and columns)!

cd,Fdir[0]
test = readfits(frames[0],/silent)
cd,homedir

n0 = n_elements(test[*,0]) ;X-size
n1 = n_elements(test[0,*]) ;Y-size
n2 = n_elements(wavelength) ;number of wavelengths (max = 9)
n3 = n0 * n1
nfibers = n_elements(Fdir)

; The arrays for all the data are created
print,''
print,'The data arrays are being created...'
print,''
fiber =  fltarr(n0,n1,nstep) ;only one fiber is held in memory at a time.
darkarr = fltarr(n0,n1,3)
;*********************THE FINAL DATA ARRAY****************************
FINALvalues = fltarr(n2,3,90,nfibers)
;D1: wavelength
;D2: EE, F/#, X2   Note: This is a trimmed down version from trans_FRDL
;D3: 90 steps in EE calculation
;D4: Fibers
;*********************THE FINAL DATA ARRAY****************************

for j=0,nfibers-1 do begin ;a loop through each fiber tested
    cntwave = 0
    cstep = 0
    if (namout ne 'none') then outname = namout+'_'+Fdir[j]+'_W'+strn(wavelength[cntwave]) else outname = namout

    for k=0,n_elements(frames)-1 do begin ;a loop through a single fiber list
        if (frames[k] ne '*') then begin
            print,'Working on frame '+frames[k]
            cd,Fdir[j]
            one = readfits(frames[k],header,/silent)
            for l=0,2 do darkarr[*,*,l] = readfits(Dframes[l],/silent)
            dark = median(darkarr,dim=3)
            cd,homedir

            t1 = median(one[5:105,5:105])
            t2 = median(dark[5:105,5:105])
            adj = t1 / t2
            
            print,'The dark frame is being scaled by x'+strn(adj)
            fiber[*,*,cstep] = one - (dark * adj)

            print,'Wavelength: '+strn(wavelength[cntwave])
            print,'Camera step: '+strn(cstep)
            jump2:
            cstep = cstep + 1
        endif else begin ;once this loops ends, all the frames for a given wavelength and fiber are in the array "fiber"
            
;********************************** THE FRATIO CODE IS STARTED HERE **********************************
            
            cd,Fdir[j] ;move into the calling directory so the COG files get written there.
            FRtemp = fratiof(fiber,cstepsize,psize,winsize,radius=starting_radius,rstep=radstep,EE=eenergy,NAME=outname,plot=frplot)
            cd,homedir

        ;  what's returned here is a 90x5 array.
        ; Row 1: EE
        ; Row 2: output f-ratio
        ; Row 3: chi-square value (the X^2 of the line fit to the radius
        ; estimates)
        ; Row 4: the trimmed chip totals for each camera step(this just goes from 0:cstep)
        ; Row 5: the % of the total light kept (should be the same amount) as
        ; if temp[#,3]/values[?,1,#,?] = temp[#,4]
        ; Row 6: THE TOTAL TRANSISSION

            FINALvalues[cntwave,0,*,j] = FRtemp[*,0] ;EE
            FINALvalues[cntwave,1,*,j] = FRtemp[*,1] ;F/# (at above EE)
            FINALvalues[cntwave,2,*,j] = FRtemp[*,2] ;X^2 fit
            
;********************************** THE FRATIO CODE IS FINISHED HERE **********************************

            cstep = 0
            cntwave = cntwave + 1
        endelse
    endfor
endfor

if (n2 eq 1) then begin
    form1 = '(a12,f11.3)'
    form2 = '(a12,f11.5)'
endif
if (n2 eq 2) then begin
    form1 = '(a12,f11.3,f11.3)'
    form2 = '(a12,f11.5,f11.5)'
endif
if (n2 eq 3) then begin
    form1 = '(a12,f11.3,f11.3,f11.3)'
    form2 = '(a12,f11.5,f11.5,f11.5)'
endif
if (n2 eq 4) then begin
    form1 = '(a12,f11.3,f11.3,f11.3,f11.3)'
    form2 = '(a12,f11.5,f11.5,f11.5,f11.5)'
endif
if (n2 eq 5) then begin
    form1 = '(a12,f11.3,f11.3,f11.3,f11.3,f11.3)'
    form2 = '(a12,f11.5,f11.5,f11.5,f11.5,f11.5)'
endif
if (n2 eq 6) then begin
    form1 = '(a12,f11.3,f11.3,f11.3,f11.3,f11.3,f11.3)'
    form2 = '(a12,f11.5,f11.5,f11.5,f11.5,f11.5,f11.5)'
endif
if (n2 eq 7) then begin
    form1 = '(a12,f11.3,f11.3,f11.3,f11.3,f11.3,f11.3,f11.3)'
    form2 = '(a12,f11.5,f11.5,f11.5,f11.5,f11.5,f11.5,f11.5)'
endif
if (n2 eq 8) then begin
    form1 = '(a12,f11.3,f11.3,f11.3,f11.3,f11.3,f11.3,f11.3,f11.3)'
    form2 = '(a12,f11.5,f11.5,f11.5,f11.5,f11.5,f11.5,f11.5,f11.5)'
endif
if (n2 eq 9) then begin
    form1 = '(a12,f11.3,f11.3,f11.3,f11.3,f11.3,f11.3,f11.3,f11.3,f11.3)'
    form2 = '(a12,f11.5,f11.5,f11.5,f11.5,f11.5,f11.5,f11.5,f11.5,f11.5)'
endif

;ALL the data is written out...
openw,24,'trans_FRD_all.dat'
for j=0,nfibers-1 do begin
    fname = 'Fiber:'+Fdir[j]
    if (j eq 0) then printf,24,'*'
    printf,24,fname 
    printf,24,'wavelength  ',wavelength,format=form1
    printf,24,'*'
    for k=0,89 do begin
        printf,24,'EE          ',FINALvalues[*,0,k,j],format=form2
        printf,24,'Output_F/#  ',FINALvalues[*,1,k,j],format=form2
        printf,24,'Chi_2       ',FINALvalues[*,2,k,j],format=form2
        if (k ne 89) then printf,24,'*'
    endfor
    printf,24,'*************************END_OF_FIBER*************************'
    printf,24,'*'
endfor
free_lun,24

print,''
print,'trans_FRDL has finished...'
print,''

stop 
END

