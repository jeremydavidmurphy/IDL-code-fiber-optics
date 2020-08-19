PRO lifetime2, list

; This routine is used to do the primary reduction of the lifetime
; test data. Another routine is run after this is complete that plots
; the header parameters against the quality of the frame (determined
; by this routine) to look for correlations.

; It's primary output is the HEADER_1x1_listname.txt file. This
; contains the critical information (the header values) as well as the
; severity of the FRD event, which is defined in this code.

; UPDATED: Running the reductions without an initial baseline seems to
; give the best results, so I am no longer keeping that loop
; updated. 05/06/2011

; Modified on May 13th to include the two humidity readings.

; MODIFIED ON MAY 18, 2011: There was a bug here in that the
; denominator wasn't being used as the subsequent numerator. There's
; also no reason to re-read in the file and header, so now the first
; loop to find the first numerator is run through only once, then
; falls into the denominator loop.

; MODIFIED ON JUNE 1, 2011: A loop was added to include a rolling
; median division. 

;************************************************************************
p2s = 'ps' ;set to 'screen' if you want to see the plots on screen. when set to 'ps', a ps file is generated
plotting1 = 'on' ;this one plots the strips
;plotting1 = 'off'
;fitdiv = 'on' ;turn out to write as fits files the divided frames that are flagged
fitdiv = 'off'
;p2s = 'ps'
c1x1 = 28.0   ;the IDL coords for camera #1 (PL9000)
c1x2 = 1489.0
c1y1 = 11.0
c1y2 = 1488.0
c2x1 = 19.0   ;the IDL coords for camera #2 (Lara's)
c2x2 = 1020.0
c2y1 = 11.0
c2y2 = 1016.0

upperT = 1.003
lowerT = 0.997

;ccdthresh = -32.5 ;if the temp goes above this, you get a flag that pauses the reductions
ccdthresh = -24.5 ;if the temp goes above this, you get a flag that pauses the reductions

cutoff = 80.0 ;This value sets the threshold that a set of divided frames will get flagged.
;If the number of pixel values in the smoothed, divided frame that
;fall above/below the upperT and lowerT thresholds, then these frames
;are flagged and a divided fits file is written out.
cutoff2 = 350.0
;a fudge to get around the partially saturated frames. These numbers
;have been tuned to pick out the bad frames coming from partially a
;opened shutter.
badx1a = 210 & bady1a = 474 & badAtot = 3800.0 
badx2a = 235 & bady2a = 629
badx1b = 210 & bady1b = 1   & badBtot = 2100.0
badx2b = 235 & bady2b = 49
;************************************************************************

readcol,list,f='a',files
n0 = n_elements(files)

free_lun,7
openw,7,'lifetime2_log_'+list+'.txt'

test = readfits(files[0],/silent)
nt1 = n_elements(test[*,0])

if (nt1 gt 1100) then fsize = 'C1' else fsize = 'C2' ;the CCD size is determined.
if (fsize eq 'C1') then begin
    printf,7,'************************************'
    printf,7,'The PL9000 was used for this data...'
    printf,7,'************************************'
    sizex = (c1x2 - c1x1) + 1.0
    sizey = (c1y2 - c1y1) + 1.0
    x1 = c1x1
    x2 = c1x2
    y1 = c1y1
    y2 = c1y2
endif
if (fsize eq 'C2') then begin
    printf,7,'************************************'
    printf,7,'The PL1001E was used for this data...'
    printf,7,'************************************'
    sizex = (c2x2 - c2x1) + 1.0
    sizey = (c2y2 - c2y1) + 1.0
    x1 = c2x1
    x2 = c2x2
    y1 = c2y1
    y2 = c2y2
endif
sl1A = (sizey/2.0)+70.0
sl2A = (sizey/2.0)+100.0
sl1B = (sizey/2.0)-100.0
sl2B = (sizey/2.0)-70.0
n2 = sizex * sizey

cstart = 'first'
cntr = 0
cntrw = 1
cntrw2 = 1
cntrBL = 0
ccdthreshcount = 0

;headARR = strarr(9,n0)
headARR = strarr(11,n0) ;added the two humidity readings on 05/13/2011

free_lun,6
free_lun,16
free_lun,26
openw,6,'flagged_1x1_'+list+'.txt'
openw,16,'BAD_1x1_'+list+'.txt'
openw,26,'HEADER_1x1_'+list+'.txt'
;    form1='(a35,1x,f8.1,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f6.1,1x,f6.1)'
form1='(a35,1x,f8.1,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f6.1,1x,f6.1)'

repeat begin

;*********************************************************************
;The numerator is dutifully sorted out of the mess of data...
;*********************************************************************
    tryagain1:
    file1 = files[cntr]
    top = readfits(file1,header,/silent)
    
    if (cntr eq 0) then begin ; a check to make sure you're on the same track.
        track = sxpar(header,'TRAJDEC',count=cnt)
        printf,7,'****************************************'
        printf,7,'Trajectory '+strn(track)+' is being used!'
        printf,7,'****************************************'
    endif
    
    top = float(top[x1:x2,y1:y2])
    t = sxpar(header,'IMAGETYP')
    t = strsplit(t,' ',/extract)
    t = t[0]
    if (t eq 'DARK') then begin
        print,'Frame ',file1,' is a DARK!'
        printf,7,'Frame ',file1,' is a DARK!'
        cntr = cntr + 1
        goto, tryagain1
    endif
    
    if (fsize eq 'C2') then begin ;now the corrupted files are located and skipped. THERE IS CURRENTLY NO CHECK ON THE PL9000 CAMERA. THIS IS DONE ONLY FOR LARA'S CAMERA.
        if (median(top[badx1a:badx2a,bady1a:bady2a]) gt badAtot) then begin
            print,'Frame ',file1,' is corrupted in the middle!'
            printf,7,'Frame ',file1,' is corrupted in the middle!'
            cntr = cntr + 1
            goto,tryagain1
        endif
        if (median(top[badx1b:badx2b,bady1b:bady2b]) gt badBtot) then begin
            print,'Frame ',file1,' is corrupted along the bottom!'
            printf,7,'Frame ',file1,' is corrupted along the bottom!'
            cntr = cntr + 1
            goto,tryagain1
        endif
    endif
    
    trackcheck = sxpar(header,'TRAJDEC',count=cnt)
    if (trackcheck ne track) then begin
        print,'Whoa! You are mixing trajectories!!!'
        print,'This happens at file'
        print,file1
;        goto,jumpend1
    endif
    
    headARR[0,cntr] = file1
    h0 = sxpar(header,'ITER',count=cnt)     & if (cnt eq 1)then headARR[1,cntr] =  h0
    h0 = sxpar(header,'CCDTEMP',count=cnt)  & if (cnt eq 1)then headARR[2,cntr] =  h0
    h0 = sxpar(header,'X_STRT',count=cnt)   & if (cnt eq 1)then headARR[3,cntr] =  h0
    h0 = sxpar(header,'Y_STRT',count=cnt)   & if (cnt eq 1)then headARR[4,cntr] =  h0
    h0 = sxpar(header,'RHO_STRT',count=cnt) & if (cnt eq 1)then headARR[5,cntr] =  h0
    h0 = sxpar(header,'OB_STRT',count=cnt)  & if (cnt eq 1)then headARR[6,cntr] =  h0
    h0 = sxpar(header,'AMBTEMP1',count=cnt) & if (cnt eq 1)then headARR[7,cntr] =  h0
    h0 = sxpar(header,'AMBTEMP2',count=cnt) & if (cnt eq 1)then headARR[8,cntr] =  h0
    h0 = sxpar(header,'HUMIDTY1',count=cnt) & if (cnt eq 1)then headARR[9,cntr] =  h0
    h0 = sxpar(header,'HUMIDTY2',count=cnt) & if (cnt eq 1)then headARR[10,cntr] = h0
    
    ccdtemp = float(headARR[2,cntr])
    if (ccdtemp gt ccdthresh) then begin
        print,'CCD TEMPERATURE PROBLEMS!'
        print,'CCD temp at ',headARR[2,cntr]
        print,file
        printf,7,'WARNING!!! CCD TEMP ABOVE NOMINAL VALUE!!! ',headARR[2,cntr],' C'
        printf,7,'High CCD temp in frame ',file1
        ccdthreshcount = ccdthreshcount + 1
    endif
    
    printf,26,headARR[*,cntr],-1.0,-1.0,form=form1
    
    cntr = cntr + 1
    
;*********************************************************************
; And now the same goes for the denominator...which is now the primary
; loop run...
;*********************************************************************
    tryagain2:
    if (cntr eq n0-1) then goto, jumpend1
    file2 = files[cntr]
    
    bot = readfits(file2,header,/silent)
    bot = float(bot[x1:x2,y1:y2])
    t = sxpar(header,'IMAGETYP')
    t = strsplit(t,' ',/extract)
    t = t[0]

    trackcheck = sxpar(header,'TRAJDEC',count=cnt)
    if (trackcheck ne track) then begin
        print,'Whoa! You are mixed trajectories!!!'
        print,'This happens at file'
        print,file2
;        goto,jumpend1
    endif
    
    if (t eq 'DARK') then begin
        print,'Frame ',file2,' is a DARK!'
        printf,7,'Frame ',file2,' is a DARK!'
        cntr = cntr + 1
        goto, tryagain2
    endif
    
    if (fsize eq 'C2') then begin
        if (median(bot[badx1a:badx2a,bady1a:bady2a]) gt badAtot) then begin
            print,'Frame ',file2,' is corrupted in the middle!'
            printf,7,'Frame ',file2,' is corrupted in the middle!'
            cntr = cntr + 1
            goto,tryagain2
        endif
        if (median(bot[badx1b:badx2b,bady1b:bady2b]) gt badBtot) then begin
            print,'Frame ',file2,' is corrupted along the bottom!'
            printf,7,'Frame ',file2,' is corrupted along the bottom!'
            cntr = cntr + 1
            goto,tryagain2
        endif
    endif
    
    headARR[0,cntr] = file2
    h0 = sxpar(header,'ITER',count=cnt)     & if (cnt eq 1)then headARR[1,cntr] =  h0
    h0 = sxpar(header,'CCDTEMP',count=cnt)  & if (cnt eq 1)then headARR[2,cntr] =  h0
    h0 = sxpar(header,'X_STRT',count=cnt)   & if (cnt eq 1)then headARR[3,cntr] =  h0
    h0 = sxpar(header,'Y_STRT',count=cnt)   & if (cnt eq 1)then headARR[4,cntr] =  h0
    h0 = sxpar(header,'RHO_STRT',count=cnt) & if (cnt eq 1)then headARR[5,cntr] =  h0
    h0 = sxpar(header,'OB_STRT',count=cnt)  & if (cnt eq 1)then headARR[6,cntr] =  h0
    h0 = sxpar(header,'AMBTEMP1',count=cnt) & if (cnt eq 1)then headARR[7,cntr] =  h0
    h0 = sxpar(header,'AMBTEMP2',count=cnt) & if (cnt eq 1)then headARR[8,cntr] =  h0
    h0 = sxpar(header,'HUMIDTY1',count=cnt) & if (cnt eq 1)then headARR[9,cntr] =  h0
    h0 = sxpar(header,'HUMIDTY2',count=cnt) & if (cnt eq 1)then headARR[10,cntr] = h0
    
    ccdtemp = float(headARR[2,cntr])
    if (ccdtemp gt ccdthresh) then begin
        print,'CCD TEMPERATURE PROBLEMS!'
        print,'CCD temp at ',headARR[2,cntr]
        print,file2
        printf,7,'WARNING!!! CCD TEMP ABOVE NOMINAL VALUE!!! ',headARR[2,cntr],' C'
        printf,7,'High CCD temp in frame ',file2
        ccdthreshcount = ccdthreshcount + 1
    endif
    
;*********************************************************************
; The division and quantification of FRD events occurs here.
;*********************************************************************
    divide = float(top / bot)
    
    stripA = median(divide[*,sl1A:sl2A],dim=2) ;the top strip is selected
    sstripA = smooth(stripA,51,/edge) ;it is smoothed
    stripB = median(divide[*,sl1B:sl2B],dim=2) ;the bottom strip is selected
    sstripB = smooth(stripB,51,/edge) ;it is smoothed
    nstd1 = n_elements(sstripA) / 2.0 ;the index of the middle is selected
    nstd2 = n_elements(sstripA)
    
    shift = 1.0 - median(sstripA,/even) ;it's roughly normalized (to account for changes in the flux from the light source)
    sstripA = sstripA + shift
    shift = 1.0 - median(sstripB,/even) ;it's roughly normalized (to account for changes in the flux from the light source)
    sstripB = sstripB + shift
    
    junk = where(sstripA gt upperT,above)
    junk = where(sstripA lt lowerT,below)
    aPbA = above + below
    junk = where(sstripB gt upperT,above)
    junk = where(sstripB lt lowerT,below)
    aPbB = above + below
    
    printf,26,headARR[*,cntr],aPbA,aPbB,form=form1
    
    print,'The final divided frames are'
    print,file1,' / ',file2
    printf,7,file1,' / ',file2,' ',strn(aPbA),' ',strn(aPbB)
    
    if (aPbA gt cutoff or aPbB gt cutoff) then begin ;if the division falls above the lower cutoff
        print,'File ',file1,' / ',file2,' is flagged!'
        printf,7,'File ',file1,' / ',file2,' is flagged!'
        printf,7,''
        printf,6,file1,' ',file2,' ',strn(aPbA),' ',strn(aPbB)
        if (fitdiv eq 'on') then begin
            com1 = '  '+file1
            com2 = '  '+file2
            sxaddpar,header,'COMMENT',com1
            sxaddpar,header,'COMMENT',com2
            writefits,'divide_'+strn(cntrw)+'.fits',divide,htop
        endif
        cntrw = cntrw + 1
    endif
    
    if (aPbA gt cutoff2 or aPbB gt cutoff2) then begin ;if the division falls above the higher cutoff
        printf,16,file1,' ',file2,' ',strn(aPbA),' ',strn(aPbB)
        cntrw2 = cntrw2 + 1
    endif
    
    if (plotting1 eq 'on') then begin
        if (p2s eq 'screen') then begin
            if (cstart eq 'first') then begin
                set_plot,'x'
                window,2,retain=2,xsize=1100,ysize=300
                device,decomposed=0
                loadct,0
                if (Fsize eq 'C1') then plot,sstripA,yrange=[0.98,1.02],/ystyle,xrange=[0,1500],/xstyle,$
                  title='Strip A: above the center'
                if (Fsize eq 'C2') then plot,sstripA,yrange=[0.98,1.02],/ystyle,xrange=[0,1030],/xstyle,$
                  title='Strip A: above the center'
                loadct,27
                oplot,[0,1500],[upperT,upperT],color=150
                oplot,[0,1500],[lowerT,lowerT],color=150
                
                window,3,retain=2,xsize=1100,ysize=300,ypos=145
                loadct,0
                if (Fsize eq 'C1') then plot,sstripB,yrange=[0.98,1.02],/ystyle,xrange=[0,1500],/xstyle,$
                  title='Strip B: below the center'
                if (Fsize eq 'C2') then plot,sstripB,yrange=[0.98,1.02],/ystyle,xrange=[0,1030],/xstyle,$
                  title='Strip B: below the center'
                loadct,27
                oplot,[0,1500],[upperT,upperT],color=150
                oplot,[0,1500],[lowerT,lowerT],color=150
                
                cstart = 'not first'
            endif else begin
                wset,2
                oplot,sstripA,color=cntr
                wset,3
                oplot,sstripB,color=cntr
            endelse
        endif
        
        if (p2s eq 'ps') then begin ;a postscript file is generated instead
            if (cstart eq 'first') then begin
                SA = sstripA + 0.01
                SB = sstripB - 0.01
                set_plot,'ps'
                device,file='StripDivision1X1_'+list+'.ps',/color,/portrait
                loadct,0
                if (Fsize eq 'C1') then plot,SA,yrange=[0.98,1.02],/ystyle,xrange=[30,1430],/xstyle,$
                  title='Strip A & B: above and below the center',xthick=3,ythick=3,charthick=3,thick=1,$
                  position=[0.12,0.11,0.95,0.92],ytitle='Divided values (offset by +/-0.01)',xtitle='Pixels'
                if (Fsize eq 'C2') then plot,SA,yrange=[0.98,1.02],/ystyle,xrange=[0,1030],/xstyle,$
                  title='Strip A & B: above and below the center',xthick=3,ythick=3,charthick=3,thick=1,$
                  position=[0.12,0.11,0.95,0.92],ytitle='Divided values (offset by +/-0.01)',xtitle='Pixels'
                loadct,27
                cstart = 'not first'
            endif else begin
                SA = sstripA + 0.01
                SB = sstripB - 0.01
                oplot,SA,color=round(float(cntr)/float(n0)*255.0)
                oplot,SB,color=round(float(cntr)/float(n0)*255.0)
            endelse
        endif
    endif
    
    top = bot
    file1 = file2
    cntr = cntr + 1
    
    print,''
    print,'The routine is '+strn((float(cntr)/float(n0))*100.0)+'% complete...'
    print,''

    goto,tryagain2

endrep until (cntr ge n0-1)

jumpend1:

if (plotting1 eq 'on') then begin
    if (p2s eq 'ps') then begin
        loadct,0
        oplot,[0,1500],[upperT+0.01,upperT+0.01],thick=3
        oplot,[0,1500],[lowerT+0.01,lowerT+0.01],thick=3
        oplot,[0,1500],[upperT-0.01,upperT-0.01],thick=3
        oplot,[0,1500],[lowerT-0.01,lowerT-0.01],thick=3
        device,/close_file
        set_plot,'x'
    endif
endif

if (cntrw eq 1) then printf,6,'NO BAD FRAMES FOUND!!!'
if (cntrw2 eq 1) then printf,16,'NO VERY BAD FRAMES FOUND!!!'
free_lun,6
free_lun,16
free_lun,26
free_lun,7

print,'******************************************'
print,'lifetime2 has finished successfully...'
print,'******************************************'

stop
END

;A baseline set of frames, defined by the frames bracketing the two
;"*"s, is set (if there are 2 "*"s)
;istar = where(files eq '*',count)
;if (count gt 0) then begin
;    basefiles = files[istar[0]+1:istar[1]-1]
;    nBL = n_elements(basefiles)
;    baseARR = fltarr(sizex,sizey,nBL)
;    for j=0,nBL-1 do begin
;        print,'Reading in baseline frame ',basefiles[j],'...'
;        frame = readfits(basefiles[j],h,/silent)
;        t = sxpar(h,'IMAGETYP')
;        t = strsplit(t,' ',/extract)
;        t = t[0]
;        if (t eq 'OBJECT') then begin
;            frame = float(frame[x1:x2,y1:y2])
;            baseARR[*,*,j] = frame
;        endif
;    endfor
;    baseline = median(baseARR,dim=3,/even)
;    printf,7,'The files used in the baseline are...'
;    for j=0,nBL-1 do printf,7,basefiles[j]
;    basecheck = 'on'
;endif else basecheck = 'off'

