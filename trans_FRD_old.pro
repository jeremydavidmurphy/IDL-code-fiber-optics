PRO trans_FRD

; The final pipeline to reduce VIRUS fiber transmission and FRD data.

; COMPILE: fratioF.pro
;*************************************************************************
;*************************************************************************

outname = 'VP5T-3-18-10' ;the output name that get's subscripted. The final
                     ;name would be "trans_frd_'outname'.dat'

plotting = 'on' ;set to anything else to turn OFF the on-screen plotting
;plotting = 'off'
nstep = 7 ;this is the number of camera positions.
cstepsize = 4.0 ;the camera step size, in mm
psize = 0.000012 ;the FLI camera pixel size, in mm.
thresh = 5000.0 ;the threshold above which the code fits an ellipse for 
                ;determination of the center of the spot. THIS IS NOT 
                ;THE THRESHHOLD FOR FRATIOF. That level is set
                ;within fratioF itself and is based on a percentage
                ;cut of the high values of the frame.
winsize = 7.0 ;set this to 1.0 for apogee data and 6.0 or 7.0 for FLI data
slow = 'off' ;setting this to 'on' will insert pauses at various points.
writebase = 'off' ;will writeout the median combined base frame that was
                 ;used as the final BL value. This frame has had the dark
                 ;adjusted and subtracted
radstep = 0.5 ; this sets the resolution of the radius step size in 
               ; fratioF.pro. This step size leads to EE resolution of
               ; about 0.1% for the 3056x3056 FLI chip size. 0.25 is
               ; as low as you need. Set to a higher number for faster
               ; (yet less precise) results. 1.0 gives EE resolution
               ; of about 0.3%. This leads to doubling up of output
               ; f-ratio values at certain EEs (which means lower
               ; doesn't do much for you.

;*************************************************************************
;*************************************************************************

; This code works from a set of lists, with the fiber lists being a
; list of lists. The format and logic is as follows...

; Description: The code will use the frames in the BL1 and BL2 lists
; (baseline #1 and baseline #2) as the values used in the calculation
; of the transmission. It is assumed the data was taken as
; follows. Baseline #1, several fibers, Baseline #2. The code reads
; the time the frame was taken out of the header and uses this to
; conduct a linear interpolation between the two baseline frames in
; order to determine the best estimate of the baseline to use for a
; give fiber. As the camera/light source transmission values are shown
; to drift <1% over several hours, the use of a linear interpolation
; is a safe one. 

; Some definitions:
; BASELINE: A data frame of the light from the bench BEFORE passing
; through a fiber.
; DARK: The dark frame taken with a given set-up, only with the
; shutter on the bench (NOT CAMERA) closed. This gives an estimate of
; the amount of background light in a frame and gets subtracted from a
; given frame. As there are times when the dark is taken at a
; different time than the frame it's being subtracted from, and some
; evolution of the level of the dark is seen over this time frame,
; there is a scaling that occurs.
; FIBER: This is the light as it comes out of a fiber.

;*************************************************************************
; Files looked for:
; BL1.list
; BL1D.list
; BL2.list
; BL2D.list
; fiber.list
; fiberD.list
; wave.list

; NFIBERS: The number of fibers being tested.

; The Lists:
; BL1L: A list of the baseline frames, in order from blue to red
; (specifically, in the same order as the wavelengths in wave.list),
; and separated by a '*'.
; 
; Each list is then just a listing of the files going into the data
; analysis, and have the same format. So, BL1_350.list could look
; like:
; BL1-001-350.fit
; BL1-002-350.fit
; BL1-003-350.fit
; *
; BL1-001-365.fit
; BL1-002-365.fit
; BL1-003-365.fit
; *
; etc.,
; with the first number corresponding to the baseline taken (these
; could perhaps never be repeated), the second to the number of the
; exposure (generally 5) and the last the wavelength. These are then
; median combined and used as the value for that baseline at that
; wavelength.

; The other lists (BL2 and both BL1D and BL2D) are of the same format,
; with the capital 'D' indicating a dark frame. EX: BL1-001-380D.fit

; THE FIBER LIST IS DIFFERENT. IT IS A LIST OF LISTS.
; They also all need to be of the same length.
; If you're testing 3 fibers, then the fibL is of the form:
; F01.list
; F02.list
; F03.list
; etc...

; The fiber darks are the same, with fibDL of the form:
; F01D.list
; F02D.list
; F03D.list
; etc...

; Each of these lists then has the following format
; F01_005_365.fit
; F01_004_365.fit
; F01_003_365.fit
; F01_002_365.fit
; F01_001_365.fit
; *
; F01_005_380.fit
; F01_004_380.fit
; F01_003_380.fit
; F01_002_380.fit
; F01_001_380.fit
; *
; F01_005_400.fit
; etc...

; NOTE: THE SPOT SIZES MUST BE ORDERED FROM LARGEST TO SMALLEST. IF
; NOT, FRATIOF.PRO WILL CRASH.

; The '*' separate the wavelengths.
; The fiber lists are therefore the lists of all the fits frames, with
; the first block being all the same wavelength, yet at different
; camera positions.

; Another major difference between the baseline frames and the fiber
; frames is that the fiber frames are taken at different camera
; positions and therefore don't get combined as the baseline frames
; do. Rather, each fiber position is used to make an estimate of the
; transmission at that wavelength. Then the transmission values
; (normalized to 1) are median combined to get the final transmission
; curve for the fiber.

; JAN 22, 2010: The code currently takes one hour per fiber with
; plotting off. I don't imagine turning the plotting on would increase
; this time much.

; SETTING UP THE LISTS...
; To establish the lists, a couple points to work from. Each list must
; have wavelengths listed only once. So, if you've done repeat
; measurements at the same wavelength, these must be sent in as
; separate fibers. The wave.list file is what dictates how you run
; things, and it doesn't like (i.e. will crash) repeat wavelengths.

; MODIFIED 03-19-2010: The code now locates and rejects saturated
; frames properly. For the baselines, it just rejects them
; outright. One key bit is that it DOES NOT DO THIS FOR THE BACKGROUND
; FRAMES. As the saturation determination is time-intensive (it sorts
; the entire frame) it is better to do a ds9 check of all the darks
; before you feed them into the code. 
;*************************************************************************

free_lun,15
openw,15,'trans_FRD.log'
printf,15,'start time:  ',systime()

wavelist = 'wave.list'
BL1L = 'BL1.list'
BL2L = 'BL2.list'
BL1DL = 'BL1D.list'
BL2DL = 'BL2D.list'
fibL = 'fiber.list'
fibDL = 'fiberD.list'

satur = 'not'

; The lists are read in
readcol,silent=1,wavelist,format='i',wavelength
readcol,silent=1,BL1L,format='a',BL1list
readcol,silent=1,BL2L,format='a',BL2list
readcol,silent=1,fibL,format='a',fiblist
readcol,silent=1,BL1Dl,format='a',BL1Dlist
readcol,silent=1,BL2DL,format='a',BL2Dlist
readcol,silent=1,fibDL,format='a',fibDlist

; A sample frame is read in so that the size of each frame can be
; determined. The data arrays are then determined. THIS ASSUMES ALL
; FRAMES ARE OF THE SAME SIZE (same number of rows and columns)!

test = readfits(BL1list[0],/silent);***
n0 = n_elements(test[*,0]) ;X-size
n1 = n_elements(test[0,*]) ;Y-size
n2 = n_elements(wavelength) ;number of wavelengths
n3 = n0 * n1
nfibers = n_elements(fiblist)

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
FINALvalues = fltarr(n2,6,30,nfibers)
;D1: wavelength
;D2: EE, F/#, X2, total trimmed counts (TTC), % of total light kept
; (TTC is only nsteps deep in the D3 direction, with each
; number being an estimate from a different fiber position. the same
; stands for % kept)
;D3: 30 steps in EE calculation
;D4: Fibers
;*********************THE FINAL DATA ARRAY****************************

colors = intarr(nfibers)
for j=0,nfibers-1 do colors[j] = 60 + (j * floor(195.0/nfibers))

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
for j=0,n_elements(BL1list)-1 do begin ; loop through each wavelength in BL1
    if (BL1list[j] ne '*') then begin
        cntr2 = cntr2 + 1
    endif else begin
        trim = BL1list[cntr1:cntr2-1]
        cntr2 = cntr2 + 1
        cntr1 = cntr2
        ntemp = n_elements(trim)
        temp = fltarr(n0,n1,ntemp)
        time = fltarr(ntemp)
        goodi = 0
        for k=0,ntemp-1 do begin
            one = readfits(trim[k],header,/silent)
            sone = one[bsort(one)]
            if (median(sone[n3-101:*]) gt 65000.0) then begin ;a check for saturated frames
                print,'SATURATION IN '+trim[k]
                printf,15,'SATURATION IN '+trim[k]
                goto, jump0
            endif
            goodi = [goodi,k]
            temp[*,*,k] = one
            cntr = 0
            out = 'lost'
            repeat begin ;the exp. time is picked out of the header
                hh = strsplit(header[cntr],' ',/extract)
                if (hh[0] eq 'DATE-OBS=') then out = 'found' else $
                  cntr = cntr + 1
            endrep until (out eq 'found')
            hh = hh[1]
            hhh = strsplit(hh,'T',/extract)
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

        print,'The exposure time for BL1 '+strn(wavelength[cntwave])+$
          ' is: '+strn(baseline[cntwave,0])
        printf,15,'The exposure time for BL1 '+strn(wavelength[cntwave])+$
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
    printf,15,''
endfor

BL1center = median(blcent,dimension=2,/even)
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
for j=0,n_elements(BL2list)-1 do begin
    if (BL2list[j] ne '*') then begin
        cntr2 = cntr2 + 1
    endif else begin
        trim = BL2list[cntr1:cntr2-1]
        cntr2 = cntr2 + 1
        cntr1 = cntr2
        ntemp = n_elements(trim)
        temp = fltarr(n0,n1,ntemp)
        time = fltarr(ntemp)
        goodi = 0
        for k=0,ntemp-1 do begin
            one = readfits(trim[k],header,/silent)
            sone = one[bsort(one)]
            if (median(sone[n3-101:*]) gt 65000.0) then begin ;a check for saturated frames
                print,'SATURATION IN '+trim[k]
                printf,15,'SATURATION IN '+trim[k]
                goto, jump1
            endif
            goodi = [goodi,k]
            temp[*,*,k] = one
            cntr = 0
            out = 'lost'
            repeat begin
                hh = strsplit(header[cntr],' ',/extract)
                if (hh[0] eq 'DATE-OBS=') then out = 'found' else $
                  cntr = cntr + 1
            endrep until (out eq 'found')
            hh = hh[1]
            hhh = strsplit(hh,'T',/extract)
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
        print,'The exposure time for BL2 '+strn(wavelength[cntwave])+$
          ' is: '+strn(baseline[cntwave,2])
        printf,15,'The exposure time for BL2 '+strn(wavelength[cntwave])+$
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
    printf,15,''
endfor

BL2center = median(blcent,dimension=2,/even)
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
for j=0,n_elements(BL1Dlist)-1 do begin
    if (BL1Dlist[j] ne '*') then begin
        cntr2 = cntr2 + 1
    endif else begin
        trim = BL1Dlist[cntr1:cntr2-1]
        cntr2 = cntr2 + 1
        cntr1 = cntr2
        ntemp = n_elements(trim)
        temp = fltarr(n0,n1,ntemp)
        time = fltarr(ntemp)

        for k=0,ntemp-1 do temp[*,*,k] = readfits(trim[k],/silent)
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
for j=0,n_elements(BL2Dlist)-1 do begin
    if (BL2Dlist[j] ne '*') then begin
        cntr2 = cntr2 + 1
    endif else begin
        trim = BL2Dlist[cntr1:cntr2-1]
        cntr2 = cntr2 + 1
        cntr1 = cntr2
        ntemp = n_elements(trim)
        temp = fltarr(n0,n1,ntemp)
        time = fltarr(ntemp)

        for k=0,ntemp-1 do temp[*,*,k] = readfits(trim[k],/silent)
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
    t1 = BL1[*,*,j]
    t2 = BL1D[*,*,j]
    t1 = t1[bsort(t1)]
    t2 = t2[bsort(t2)]
    t1m = median(t1[0:1000])
    t2m = median(t2[0:1000])
    adj = t1m / t2m
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
    t1 = BL2[*,*,j]
    t2 = BL2D[*,*,j]
    t1 = t1[bsort(t1)]
    t2 = t2[bsort(t2)]
    t1m = median(t1[0:1000])
    t2m = median(t2[0:1000])
    adj = t1m / t2m
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
    readcol,silent=1,fiblist[j],format='a',frames ;the list of ALL the fits files for a given fiber.
    readcol,silent=1,fibDlist[j],format='a',Dframes
    if (n_elements(frames) ne n_elements(Dframes)) then begin
        print,'Your data list and dark list for '+fiblist[j]+' is not of the same length!'
        print,'Fix this and rerun...'
        stop
    endif
    cntwave = 0
    cstep = 0
    for k=0,n_elements(frames)-1 do begin ;a loop through the a single fiber list
        if (frames[k] ne '*') then begin
            print,'Working on frame '+frames[k]
            printf,15,'Working on frame '+frames[k]
            one = readfits(frames[k],header,/silent)
            sone = one[bsort(one)]

            if (median(sone[n3-101:*]) gt 65000.0) then begin ;a check for saturated frames
                print,'SATURATION IN '+frames[k]
                printf,15,'SATURATION IN '+frames[k]
                goto, jump2 ;saturated frames are sent into fratioF as zeros and are handled there.
            endif

            dark = readfits(Dframes[k],/silent)
            fiber[*,*,cstep] = one
            cntr = 0
            out = 'lost'
            repeat begin
                hh = strsplit(header[cntr],' ',/extract)
                if (hh[0] eq 'DATE-OBS=') then out = 'found' else $
                  cntr = cntr + 1
            endrep until (out eq 'found')
            hh = hh[1]
            hhh = strsplit(hh,'T',/extract)
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

            adj = t1m / t2m
            print,'Fiber frame '+Dframes[k]+' is being scaled by x'+strn(adj)
            printf,15,'Fiber frame '+Dframes[k]+' is being scaled by x'+strn(adj)
            fiber[*,*,cstep] = fiber[*,*,cstep] - (dark * adj)
            values[cntwave,1,cstep,j] = total(fiber[*,*,cstep])
            ; calculating proper BL to divide by
            if (values[cntwave,0,cstep,j] lt baseline[cntwave,0]) then values[cntwave,0,cstep,j] = values[cntwave,0,cstep,j] + 24.00
            tempBL[cstep] = baseline[cntwave,1] + (baseline[cntwave,4] * (values[cntwave,0,cstep,j]-baseline[cntwave,0]))
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

;********************************** THE FRATIO CODE IS RUN HERE **********************************
            
            FRtemp = fratioF(fiber,cstepsize,psize,radius=150.0,rstep=radstep)
; what's returned here is a 30x5 array.
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

;********************************** THE FRATIO CODE IS RUN HERE **********************************

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
    fname = strsplit(fiblist[j],'.',/extract)
    fname = 'Fiber: '+fname[0]
    printf,24,fname
    printf,24,'wavelength  ',wavelength,format=form1
    printf,24,''
    for k=0,cstep-1 do begin
        printf,24,'EE          ',FINALvalues[*,0,k,j],format=form2
        printf,24,'Output F/#  ',FINALvalues[*,1,k,j],format=form2
        printf,24,'Chi 2       ',FINALvalues[*,2,k,j],format=form2
        printf,24,'% flux kept ',FINALvalues[*,4,k,j],format=form2
        printf,24,'transmission',FINALvalues[*,5,k,j],format=form2
        if (k ne 29) then printf,24,''
    endfor
    for k=cstep,29 do begin
        printf,24,'EE          ',FINALvalues[*,0,k,j],format=form2
        printf,24,'Output F/#  ',FINALvalues[*,1,k,j],format=form2
        printf,24,'Chi 2       ',FINALvalues[*,2,k,j],format=form2
        if (k ne 29) then printf,24,''
    endfor
    printf,24,'*************************END OF FIBER*************************'
    printf,24,''
endfor
free_lun,24

;The important data is written out...
openw,24,'trans_FRD_trim.dat'
for j=0,nfibers-1 do begin
    fname = strsplit(fiblist[j],'.',/extract)
    fname = 'Fiber: '+fname[0]
    printf,24,fname
    printf,24,'wavelength  ',wavelength,format=form1
    printf,24,''
    for k=0,nstep-1 do begin
        printf,24,'EE          ',FINALvalues[*,0,30-nstep+k,j],format=form2
        printf,24,'Output F/#  ',FINALvalues[*,1,30-nstep+k,j],format=form2
        printf,24,'transmission',FINALvalues[*,5,k,j],format=form2
        if(k ne nstep-1) then printf,24,''
    endfor
    printf,24,'*************************END OF FIBER*************************'
    printf,24,''
endfor
free_lun,24

printf,15,'end time: '+systime()
free_lun,15

print,'trans_FRD ended successfully'
stop 
END
