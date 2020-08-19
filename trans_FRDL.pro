PRO trans_FRDL

; The final pipeline to reduce VIRUS fiber transmission and FRD
; data. The only change from trans_FRD to this version is that the
; lists have changed and go into directories to find files.

; NOTE: This is a reworking of trans_FRD with the idea that the data
; is kept in separate directories and the code pulls the files from
; these directories, as needed. This came about because of the LabVIEW
; limitation of renaming files. Now all the fiber data is named the
; same thing, so the only thing that keeps things apart is the
; directory name. trans_FRDL is now run from a directory where all the
; relevant baselines and fiber data are in subdirectories. 

; The following lists are needed:
; fiberdir.list
; baselinedir.list
; wave.list
; bl.list 
; blD.list
; fiber.list
; fiberD.list

; COMPILE: fratiof.pro

;OUTPUT: The code outputs two things, with the possibility of a 3rd output.

; trans_FRD_all.dat: the primary output is an ascii file that contains
; all the f-ratio values at 90 different EE values, and the
; transmission. USE Ptrans_FRD to plot (and make sense) of this.
; trans_FRD_trim.dat: A trimmed value showing just the final 7 f-ratio
; out values for the last 7 EE values, and the transmission. Same #'s
; as trans_FRD_all.dat, only in a different configuration.
; The other piece saved (if 'namout' below is NOT set to 'none') is
; the COG ascii files. These are just radius (in pixels) vs. EE. It
; saves every step, so they are just good for plotting purposes with
; high spatial resolution.

;*************************************************************************
;*************************************************************************

namout = 'vpaip2' ; The output name for the curve-of-growth. This is just a ascii
                       ; plot of ee as a function of radius. These are saved unless
                       ; namout='none'. The fiber name is added, so, if the fiber name was
                       ; S18, the output is
                       ; namout_S18_CP#.txt
                       ; where CP# is the camera position.

plotting = 'off'     ; set to anything other than 'on' to turn off the plotting
frplot = 'off'       ; set to turn on the f-ration plotting

nstep = 7            ; this is the number of camera positions.
cstepsize = 4.0      ; the camera step size, in mm
psize = 0.000012     ; the FLI camera pixel size, in mm. (change to 0.000024
                     ; for binned data)
thresh = 2000.0      ; the threshold above which the code fits an ellipse for 
                     ; determination of the center of the spot. THIS IS NOT 
                     ; THE THRESHHOLD FOR FRATIOF. That level is set
                     ; within fratioF itself and is based on a percentage
                     ; cut of the high values of the frame.
winsize = 7.0        ; set this to 1.0 for apogee data and 6.0 or 7.0 for FLI data.
                     ; (this is used only when plotting is turned on)
slow = 'off'         ; setting this to 'on' will insert pauses at various points.
writebase = 'off'    ; will writeout the median combined base frame that was
                     ; used as the final BL value. This frame has had the dark
                     ; adjusted and subtracted
radstep = 1.0        ; this sets the resolution of the radius step size in 
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
satCHECK = 'off'     ; Turn this 'on' if you want to monitor (and step over)
                     ; saturated frames. It's time-expensive, so it's suggested
                     ; that you visually inspect your frames.
;*************************************************************************
;*************************************************************************
; NOTE: THE SPOT SIZES MUST BE ORDERED FROM LARGEST TO SMALLEST. IF
; NOT, FRATIOF.PRO WILL CRASH.

; SEE trans_FRD FOR MORE DETAILED NOTES.

;*************************************************************************

eenergy = 0.95 ; How far the fratioF routine integrates up to. DON'T GO CHANGING THIS! IT FUCKS OTHER THINGS UP THAT I DON'T WANT TO TRY AND FIX!

cd, current=homedir ;set your current direcory as home.

free_lun,15
openw,15,'trans_FRD.log'
printf,15,'start time:  ',systime()

readcol,'fiberdir.list',f='a',Fdir
readcol,'baselinedir.list',f='a',BLdir

;these all live in the main directory, and are the actual lists. the
;only 'list-of-lists' are the Fdir and BLdir from above
wavelist = 'wave.list'
blL = 'bl.list' 
blDL = 'blD.list'
fibL = 'fiber.list'
fibDL = 'fiberD.list'

; The lists are read in; NFIBERS: The number of fibers being tested.

readcol,silent=1,wavelist,format='i',wavelength
readcol,silent=1,blL,format='a',BLlist
readcol,silent=1,blDL,format='a',BLDlist
readcol,silent=1,fibL,format='a',frames
readcol,silent=1,fibDL,format='a',Dframes

if (n_elements(frames) ne n_elements(Dframes)) then begin
    print,'Your data list and dark list for '+fiblist[j]+' is not of the same length!'
    print,'Fix this and rerun...'
    stop
endif

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
BL1 = fltarr(n0,n1,n2)
BL2 = fltarr(n0,n1,n2)
BL1D = fltarr(n0,n1,n2)
BL2D = fltarr(n0,n1,n2)
fiber =  fltarr(n0,n1,nstep) ;only one fiber is held in memory at a time.

BL1T = fltarr(n2)
BL2T = fltarr(n2)
fiberT = fltarr(n2,nstep)
baseline = fltarr(n2,5) ;final BL values (BL1_time, BL1_ttl, BL2_time, BL2_ttl, slope)
values = fltarr(n2,5,nstep,nfibers) ;D1:wavelength, D2:various, D3:steps, D4: different fibers
ellipse = fltarr(2,120)
blcent = fltarr(2,n2)
bldia = fltarr(2,n2)

;*********************THE FINAL DATA ARRAY****************************
FINALvalues = fltarr(n2,6,90,nfibers)
;D1: wavelength
;D2: EE, F/#, X2, total trimmed counts (TTC), % of total light kept
; (TTC is only nsteps deep in the D3 direction, with each
; number being an estimate from a different fiber position. the same
; stands for % kept)
;D3: 90 steps in EE calculation
;D4: Fibers
;*********************THE FINAL DATA ARRAY****************************

;colors = intarr(nfibers)
;for j=0,nfibers-1 do colors[j] = 60 + (j * floor(195.0/nfibers))

colors = (findgen(nfibers) * (195.0/nfibers)) + 60.0
colors = floor(colors)

symbols = [6,7,8,9,10,16,17,11,12]
symbols = symbols[0:nstep-1]

if (plotting eq 'on') then begin
    set_plot,'x'
    window,0,retain=2,xsize=round(n0/winsize),ysize=round(n1/winsize)
    device,decomposed=0
    loadct,37
endif

; BASELINE 1 IS READ IN
cntr1 = 0
cntr2 = 0 ;counts the number of fits files (steps over '*'s)
cntwave = 0
for j=0,n_elements(BLlist)-1 do begin ; loop through each wavelength in BL1
    if (BLlist[j] ne '*') then begin
        cntr2 = cntr2 + 1
    endif else begin
        cd,BLdir[0]
        trim = BLlist[cntr1:cntr2-1]
        cntr2 = cntr2 + 1
        cntr1 = cntr2
        ntemp = n_elements(trim)
        temp = fltarr(n0,n1,ntemp)
        time = fltarr(ntemp)
        goodi = 0
        for k=0,ntemp-1 do begin
            one = readfits(trim[k],header,/silent)
            if (satCHECK eq 'on') then begin
                sone = one[bsort(one)]
                if (median(sone[n3-101:*]) gt 65000.0) then begin ;a check for saturated frames
                    print,'SATURATION IN '+trim[k]
                    cd,homedir
                    printf,15,'SATURATION IN '+trim[k]
                    cd,BLdir[0]
                    goto, jump0
                endif
            endif
            goodi = [goodi,k]
            temp[*,*,k] = one
            date = sxpar(header,'DATE-OBS',count=cnt)
            hhh = strsplit(date,'T',/extract)
            hhh = strsplit(hhh[1],':',/extract)
            hour = float(hhh[0])
            min = float(hhh[1])
            sec = float(hhh[2])
            time[k] = hour + min/60. + sec/3600.
            if (k gt 0) then begin ;a correction is made for 24 hr rollever
                t1 = time[0]
                t2 = time[k]
                if (t1 gt t2) then time[k] = time[k] + 24.00
            endif
            jump0:
        endfor
        cd,homedir ;back to calling directory in order to write to the log

        goodi = goodi[1:*]
        temp = temp[*,*,goodi] ;the saturated frames are dumped
        time = time[goodi]
        trim = trim[goodi]
        ntemp = n_elements(trim)

        if (ntemp gt 1) then begin
            BL1[*,*,cntwave] = median(temp,dimension=3,/even)
            baseline[cntwave,0] = median(time,/even)
            print,'BL1 frames to be combined...'
            print,trim
            printf,15,'BL1 frames to be combined...'
            printf,15,trim
        endif else begin
            BL1[*,*,cntwave] = temp
            baseline[cntwave,0] = time
            print,'Single BL1 frame being used...'
            print,trim
            printf,15,'Single BL1 frame being used...'
            printf,15,trim
        endelse

        print,'The time of exposure for BL1 '+strn(wavelength[cntwave])+$
          ' is: '+strn(baseline[cntwave,0])
        printf,15,'The time of exposure for BL1 '+strn(wavelength[cntwave])+$
          ' is: '+strn(baseline[cntwave,0])
        index = where(BL1[*,*,cntwave] ge thresh)
        ellipse[*,*] = fit_ellipse(index,xsize=n0,ysize=n1,$
                                   center=cen,axes=ax)
        if (plotting eq 'on') then begin
            tvimage,bytscl(BL1[*,*,cntwave],Top=!D.Table_Size-3)
            plots,ellipse/winsize,/device,color=fsc_color('yellow'),thick=2
            if (slow eq 'on') then pause
        endif
        print,'The center of the spot is at: '+strn(cen)
        print,'The spot major and minor axes are: '+strn(ax)
        printf,15,'The center of the spot is at: '+strn(cen)
        printf,15,'The spot major and minor axes are: '+strn(ax)
        blcent[*,cntwave] = cen
        bldia[*,cntwave] = ax
        cntwave = cntwave + 1
    endelse
endfor

if (cntwave gt 1) then BL1center = median(blcent,dimension=2,/even) else BL1center = blcent
BL1rad = median(bldia)/2.0
print,'The final baseline 1 center is: '+strn(BL1center)
print,'The final baseline 1 diameter is :'+strn(BL1rad * 2.0)
printf,15,'The final baseline 1 center is: '+strn(BL1center)
printf,15,'The final baseline 1 diameter is :'+strn(BL1rad * 2.0)
printf,15,'********************** END OF BL1 **********************'

; BASELINE 2 IS READ IN
cntr1 = 0
cntr2 = 0
cntwave = 0
for j=0,n_elements(BLlist)-1 do begin
    if (BLlist[j] ne '*') then begin
        cntr2 = cntr2 + 1
    endif else begin
        cd,BLdir[1]
        trim = BLlist[cntr1:cntr2-1]
        cntr2 = cntr2 + 1
        cntr1 = cntr2
        ntemp = n_elements(trim)
        temp = fltarr(n0,n1,ntemp)
        time = fltarr(ntemp)
        goodi = 0
        for k=0,ntemp-1 do begin
            one = readfits(trim[k],header,/silent)
            if (satCHECK eq 'on') then begin
                sone = one[bsort(one)]
                if (median(sone[n3-101:*]) gt 65000.0) then begin ;a check for saturated frames
                    print,'SATURATION IN '+trim[k]
                    cd,homedir
                    printf,15,'SATURATION IN '+trim[k]
                    cd,BLdir[1]
                    goto, jump1
                endif
            endif
            goodi = [goodi,k]
            temp[*,*,k] = one
            date = sxpar(header,'DATE-OBS',count=cnt)
            hhh = strsplit(date,'T',/extract)
            hhh = strsplit(hhh[1],':',/extract)
            hour = float(hhh[0])
            min = float(hhh[1])
            sec = float(hhh[2])
            time[k] = hour + min/60. + sec/3600.
            if (k gt 0) then begin ;a correction is made for 24 hr rollever
                t1 = time[0]
                t2 = time[k]
                if (t1 gt t2) then time[k] = time[k] + 24.00
            endif
            jump1:
        endfor
        cd,homedir

        goodi = goodi[1:*]
        temp = temp[*,*,goodi] ;the saturated frames are dumped
        time = time[goodi]
        trim = trim[goodi]
        ntemp = n_elements(trim)

        if (ntemp gt 1) then begin
            BL2[*,*,cntwave] = median(temp,dimension=3,/even)
            baseline[cntwave,2] = median(time,/even)
            print,'BL2 frames to be combined...'
            print,trim
            printf,15,'BL2 frames to be combined...'
            printf,15,trim
            if (baseline[cntwave,0] gt baseline[cntwave,2]) then baseline[cntwave,2] = baseline[cntwave,2] + 24.00
        endif else begin
            BL2[*,*,cntwave] = temp
            baseline[cntwave,2] = time
            print,'Single BL2 frame being used...'
            print,trim
            printf,15,'Single BL2 frame being used...'
            printf,15,trim
            if (baseline[cntwave,0] gt baseline[cntwave,2]) then baseline[cntwave,2] = baseline[cntwave,2] + 24.00
        endelse
        print,'The time of exposure for BL2 '+strn(wavelength[cntwave])+$
          ' is: '+strn(baseline[cntwave,2])
        printf,15,'The time of exposure for BL2 '+strn(wavelength[cntwave])+$
          ' is: '+strn(baseline[cntwave,2])
        index = where(BL2[*,*,cntwave] ge thresh)
        ellipse[*,*] = fit_ellipse(index,xsize=n0,ysize=n1,$
                                   center=cen,axes=ax)
        if (plotting eq 'on') then begin
            tvimage,bytscl(BL2[*,*,cntwave],Top=!D.Table_Size-3)
            plots,ellipse/winsize,/device,color=fsc_color('yellow'),thick=2
            if (slow eq 'on') then pause
        endif
        print,'The center of the spot is at: '+strn(cen)
        print,'The spot major and minor axes are: '+strn(ax)
        printf,15,'The center of the spot is at: '+strn(cen)
        printf,15,'The spot major and minor axes are: '+strn(ax)
        blcent[*,cntwave] = cen
        bldia[*,cntwave] = ax
        cntwave = cntwave + 1
    endelse
endfor

if (cntwave gt 1) then BL2center = median(blcent,dimension=2,/even) else BL2center = blcent
BL2rad = median(bldia)/2.0
print,'The final baseline 2 center is: '+strn(BL2center)
print,'The final baseline 2 diameter is :'+strn(BL2rad * 2.0)
printf,15,'The final baseline 2 center is: '+strn(BL2center)
printf,15,'The final baseline 2 diameter is :'+strn(BL2rad * 2.0)
printf,15,'********************** END OF BL2 **********************'

; BASELINE 1 DARK IS READ IN
cntr1 = 0
cntr2 = 0 ;counts the number of fits files (steps over '*'s)
cntwave = 0
for j=0,n_elements(BLDlist)-1 do begin
    if (BLDlist[j] ne '*') then begin
        cntr2 = cntr2 + 1
    endif else begin
        trim = BLDlist[cntr1:cntr2-1]
        cntr2 = cntr2 + 1
        cntr1 = cntr2
        ntemp = n_elements(trim)
        temp = fltarr(n0,n1,ntemp)
        time = fltarr(ntemp); To establish the lists, a couple points to work from. Each list must

        cd,BLdir[0]
        for k=0,ntemp-1 do temp[*,*,k] = readfits(trim[k],/silent)
        cd,homedir

        if (ntemp gt 1) then begin
            BL1D[*,*,cntwave] = median(temp,dimension=3,/even)
            print,'BL1D frames to be combined...'
            print,trim
            printf,15,'BL1D frames to be combined...'
            printf,15,trim
        endif else begin
            BL1D[*,*,cntwave] = temp
            print,'Single BL1D frame being used...'
            print,trim
            printf,15,'Single BL1D frame being used...'
            printf,15,trim
        endelse
        cntwave = cntwave + 1
    endelse
endfor

; BASELINE 2 DARK IS READ IN
cntr1 = 0
cntr2 = 0 ;counts the number of fits files (steps over '*'s)
cntwave = 0
for j=0,n_elements(BLDlist)-1 do begin
    if (BLDlist[j] ne '*') then begin
        cntr2 = cntr2 + 1
    endif else begin
        trim = BLDlist[cntr1:cntr2-1]
        cntr2 = cntr2 + 1
        cntr1 = cntr2
        ntemp = n_elements(trim)
        temp = fltarr(n0,n1,ntemp)
        time = fltarr(ntemp)

        cd,BLdir[1]
        for k=0,ntemp-1 do temp[*,*,k] = readfits(trim[k],/silent)
        cd,homedir

        if (ntemp gt 1) then begin
            BL2D[*,*,cntwave] = median(temp,dimension=3,/even)
            print,'BL2D frames to be combined...'
            print,trim
            printf,15,'BL2D frames to be combined...'
            printf,15,trim
        endif else begin
            BL2D[*,*,cntwave] = temp
            print,'Single BL2D frame being used...'
            print,trim
            printf,15,'Single BL2D frame being used...'
            printf,15,trim
        endelse
        cntwave = cntwave + 1
    endelse
endfor

; The dark frames are subtracted from the baseline frames. An
; adjustment to the background is made, if necessary.
for j=0,n2-1 do begin ;a loop through wavelength
    t1 = median(BL1[5:105,5:105,j])
    t2 = median(BL1D[5:105,5:105,j])
    adj = t1 / t2
    print,'Baseline dark #1 '+strn(wavelength[j])+' is being scaled by x'+strn(adj)
    printf,15,'Baseline dark #1 '+strn(wavelength[j])+' is being scaled by x'+strn(adj)
    BL1[*,*,j] = BL1[*,*,j] - (BL1D[*,*,j] * adj)
    r = BL1rad + 50
    c = BL1center
    baseline[j,1] = total(BL1[c[0]-r:c[0]+r,c[1]-r:c[1]+r,j])
    print,strn(baseline[j,1]/total(BL1[*,*,j])*100)+'% of the baseline frame is being kept.'
    printf,15,strn(baseline[j,1]/total(BL1[*,*,j])*100)+'% of the baseline frame is being kept.'
    if (writebase eq 'on') then begin
        fname = 'BL1_'+strn(wavelength[j])+'.fits'
        writefits,fname,BL1[*,*,j]
    endif
endfor

for j=0,n2-1 do begin ;same as above for BL2
    t1 = median(BL2[5:105,5:105,j])
    t2 = median(BL2D[5:105,5:105,j])
    adj = t1 / t2
    print,'Baseline dark #2 '+strn(wavelength[j])+' is being scaled by x'+strn(adj)
    printf,15,'Baseline dark #2 '+strn(wavelength[j])+' is being scaled by x'+strn(adj)
    BL2[*,*,j] = BL2[*,*,j] - (BL2D[*,*,j] * adj)
    r = BL2rad + 50
    c = BL2center
    baseline[j,3] = total(BL2[c[0]-r:c[0]+r,c[1]-r:c[1]+r,j])
    print,strn(baseline[j,3]/total(BL2[*,*,j])*100)+'% of the baseline frame is being kept.'
    printf,15,strn(baseline[j,3]/total(BL2[*,*,j])*100)+'% of the baseline frame is being kept.'
    if (writebase eq 'on') then begin
        fname = 'BL2_'+strn(wavelength[j])+'.fits'
        writefits,fname,BL2[*,*,j]
    endif
endfor

;the slope of the baseline array is calculated...
for j=0,n2-1 do begin
    rise = baseline[j,3] - baseline[j,1]
    run = baseline[j,2] - baseline[j,0]
    if (rise eq 0.0 or run eq 0.0) then begin
        print,'A ZERO FOR THE SLOPE!'
        printf,15,'A ZERO FOR THE SLOPE!'
        baseline[j,4] = 0.0
        endif else baseline[j,4] = rise / run
    print,'The ratio of BL1 to BL2 at '+strn(wavelength[j])+' is '+strn(baseline[j,1]/baseline[j,3])
    printf,15,'The ratio of BL1 to BL2 at '+strn(wavelength[j])+' is '+strn(baseline[j,1]/baseline[j,3])
endfor

;*******************************************************************************
; At this point both the first and last baseline frames have been
; combined (at each wavelength) and had their darks subtracted. They
; are in arrays named BL1 and BL2 with the totals and times in "baseline"
;*******************************************************************************

;space is made by setting the large arrays to single values.
BL1 = 0.0
BL1D = 0.0
BL2 = 0.0
BL2D = 0.0

time = fltarr(nstep)
tempBL = fltarr(nstep)

for j=0,nfibers-1 do begin ;a loop through each fiber tested
    cntwave = 0
    cstep = 0
    printf,15,'*************** NEW FIBER ***************'
    printf,15,'                   '+Fdir[j]
    printf,15,'*************** NEW FIBER ***************'
    if (namout ne 'none') then outname = namout+'_'+Fdir[j]+'_W'+strn(wavelength[cntwave]) else outname = namout

    for k=0,n_elements(frames)-1 do begin ;a loop through a single fiber list
        if (frames[k] ne '*') then begin
            print,'Working on frame '+frames[k]
            printf,15,'Working on frame '+frames[k]
            cd,Fdir[j]
            one = readfits(frames[k],header,/silent)
            dark = readfits(Dframes[k],/silent)
            cd,homedir

            if (satCHECK eq 'on') then begin
                sone = one[bsort(one)]
                if (median(sone[n3-101:*]) gt 65000.0) then begin ;a check for saturated frames
                    print,'SATURATION IN '+frames[k]
                    printf,15,'SATURATION IN '+frames[k]
                    goto, jump2 ;saturated frames are sent into fratioF as zeros and are handled there.
                endif
            endif

            date = sxpar(header,'DATE-OBS',count=cnt)
            if (cnt eq 1) then begin
                hhh = strsplit(date,'T',/extract)
                hhh = strsplit(hhh[1],':',/extract)
                hour = float(hhh[0])
                min = float(hhh[1])
                sec = float(hhh[2])
                values[cntwave,0,cstep,j] = hour + min/60. + sec/3600.
                time[cstep] = hour + min/60. + sec/3600.
                if (cstep gt 0) then begin
                    t1 = time[0]
                    t2 = time[cstep]
                    if (t1 gt t2) then time[cstep] = time[cstep] + 24.00
                endif
            endif else begin
                print,'You have no date written into your headers!'
                stop
            endelse
            
            t1 = median(one[5:105,5:105])
            t2 = median(dark[5:105,5:105])
            adj = t1 / t2
            
            print,'Dark frame '+Dframes[k]+' is being scaled by x'+strn(adj)
            printf,15,'Dark frame '+Dframes[k]+' is being scaled by x'+strn(adj)
            fiber[*,*,cstep] = one - (dark * adj)
            values[cntwave,1,cstep,j] = total(fiber[*,*,cstep])
                                ; calculating proper BL to divide by
            if (values[cntwave,0,cstep,j] lt baseline[cntwave,0]) then $
              alues[cntwave,0,cstep,j] = values[cntwave,0,cstep,j] + 24.00
            tempBL[cstep] = baseline[cntwave,1] + (baseline[cntwave,4] * values[cntwave,0,cstep,j]-baseline[cntwave,0])
            values[cntwave,2,cstep,j] = tempBL[cstep]
            values[cntwave,3,cstep,j] = values[cntwave,1,cstep,j] / values[cntwave,2,cstep,j] ;non-trimmed TRANSMISSION
            print,'Wavelength step: '+strn(cntwave)
            print,'Camera step: '+strn(cstep)
            print,'The non-trimmed transmission for frame '+frames[k]+' is '+strn(values[cntwave,3,cstep,j])
            printf,15,'Wavelength step: '+strn(cntwave)
            printf,15,'Camera step: '+strn(cstep)
            printf,15,'The non-trimmed transmission for frame '+frames[k]+' is '+strn(values[cntwave,3,cstep,j])
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
            FINALvalues[cntwave,3,*,j] = FRtemp[*,3] ;total trimmed counts (one for each camera position)
            FINALvalues[cntwave,4,*,j] = FRtemp[*,4] ;% of chip light kept (trimmed/total)
            temp = [tempBL,fltarr(n_elements(FRtemp[*,3])-n_elements(tempBL))]
            temp[where (temp eq 0.0)] = 1.0
            FINALvalues[cntwave,5,*,j] = FRtemp[*,3] / temp ;FINAL TRANSMISSION VALUE (FROM ADJUSTED BASELINE AND TRIMMED CHIP)      
            
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
    for k=0,nstep-1 do begin
        printf,24,'EE          ',FINALvalues[*,0,k,j],format=form2
        printf,24,'Output_F/#  ',FINALvalues[*,1,k,j],format=form2
        printf,24,'Chi_2       ',FINALvalues[*,2,k,j],format=form2
        printf,24,'%_flux_kept ',FINALvalues[*,4,k,j],format=form2
        printf,24,'transmission',FINALvalues[*,5,k,j],format=form2
        if (k ne 89) then printf,24,'*'
    endfor
    for k=nstep,89 do begin
        printf,24,'EE          ',FINALvalues[*,0,k,j],format=form2
        printf,24,'Output_F/#  ',FINALvalues[*,1,k,j],format=form2
        printf,24,'Chi_2       ',FINALvalues[*,2,k,j],format=form2
;        printf,24,'Zeros       ',FINALvalues[*,4,k,j],format=form2
;        printf,24,'Zeros       ',FINALvalues[*,5,k,j],format=form2
        if (k ne 89) then printf,24,''
    endfor
    printf,24,'*************************END_OF_FIBER*************************'
    printf,24,'*'
endfor
free_lun,24

;The important data is written out...
openw,24,'trans_FRD_trim.dat'
for j=0,nfibers-1 do begin
    fname = 'Fiber: '+Fdir[j]
    printf,24,fname
    printf,24,'wavelength  ',wavelength,format=form1
    printf,24,'*'
    for k=0,nstep-1 do begin
        printf,24,'EE          ',FINALvalues[*,0,90-nstep+k,j],format=form2
        printf,24,'Output_F/#  ',FINALvalues[*,1,90-nstep+k,j],format=form2
        printf,24,'transmission',FINALvalues[*,5,k,j],format=form2
        if(k ne nstep-1) then printf,24,'*'
    endfor
    printf,24,'*************************END OF FIBER*************************'
endfor
free_lun,24

printf,15,'end time: ',systime()
free_lun,15

if (plotting eq 'on') then begin
    print,'The next ENTER deletes all open windows...'
    pause
    while !d.window ne -1 do wdelete, !d.window
endif

print,''
print,'trans_FRDL has finished...'
print,''

stop 
END
