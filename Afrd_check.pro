PRO Afrd_check, LIMIT=limit, EE=ee, RADIUS=starting_radius, STEP=stepsize

; This code is used to align the test bench and reduce artificial FRD
; due to misalignment of the stage that holds the input head.

; It looks for a list called "frdalign.list" that contains all the
; frames you want to compare.

; LIMIT: The % used to make the pixel cut for the ellipse fitting
; routine. The purpose is to limit the ellipse fitting routine to only
; consider pixel values above a certain value. The default is 10% of
; the max value. (Actually, it gets multiplied by the median of the
; highest 100 pixel values.)
; If the ellipses are not roundish and centered on the spot, try incresing
; this value (~0.30 proves to be good for most situations)

; EE: The encircled energy the routine will continue working out to.

; RADIUS: This is the value of the starting radius. This is set low
; (5) so that you can see the inner region of the fiber profile as
; well as the outer.

; STEP: The radius step size. 0.3 is good for the 512x512 camera and
; 2.0 for the FLI (3056x3056) camera. Reduce this value to increase
; the resolution of the steps.

; There is a switch in the code for "single" or "multiple" images. If
; single then your list (frdalign.list) is A LIST OF FITS FILES. This
; means you've taken ONE IMAGE at each position. If multiple then the
; list (frdalign.list) is A LIST OF LISTS where each list contains
; multiple images at the same alignment. These images will be
; median-combined before the code is run. This can help the routine
; from getting hung up on cosmic rays, although the routine seems to
; be robust enough to avoid this.

;***********************************************************************************************
bgsub = 'on' ;change to 'on' if you want to feed in a single dark to get subtracted.
wscale = 7.0 ;the scale of the chip FLI = 6.0, apogee, try 2.
;***********************************************************************************************

; MODIFIED ON DEC 9, 2009: I have simplified this code by removing
; many of the options. The code now accepts single images (although
; the old code that would median combine several frames is still here,
; just commented out) only. Also, as the FLI files are huge, the data
; isn't held in memory. The colors are also handled internally and
; should accept any number of files. Lastly, an option for a
; postscript file is now in the code.

; MODIFIED ON MAY 20, 2011: The list, "frdalign.list" and, if the
; keyword bgsub is on, the list "bg.list" is read in. The backgrounds
; from bg.list are subtracted from the files in frdalign.list. They
; must therefore be of the same length.

ans=''
titl=''

if (n_elements(limit) eq 0) then limit = 0.30
IF (N_Elements(starting_radius) EQ 0) THEN starting_radius = 5.0
IF (N_Elements(stepsize) EQ 0) THEN stepsize = 2.0 ;This value proves a good one, showing minimal repeat.
IF (N_Elements(ee) EQ 0) THEN ee = 0.985
if (ee ge 1.0) then stop

;print,'Single or multiple images? (s/m)'
;read,ans
;if (ans eq 'm') then writefiles = 'yes' ;change to something else if you don't want to write out the median-combined file.

;if (ans eq 's') then begin
readcol,'frdalign.list',format='a',frames
frames = transpose(frames)
n0 = n_elements(frames)
test = readfits(frames[0])
n1 = n_elements(test[*,0])
n2 = n_elements(test[0,*])
center = fltarr(2,n0)
names = frames

if (bgsub eq 'on') then begin
    readcol,'bg.list',format='a',bglist
    nbg = n_elements(bglist)
    if (nbg ne n0) then begin
        print,'Your lists are not of the same length!'
        print,'Correct this and try again!'
        stop
    endif
endif
;endif

;if (ans eq 'm') then begin
;    readcol,'frdalign.list',format='a',lists
;    n0 = n_elements(lists)
;    names = strarr(n0)
;    for j=0,n0-1 do begin
;        readcol,lists[j],format='a',frames
;        names[j] = lists[j]
;        nf1 = n_elements(frames)
;        test = readfits(frames[0])
;        n1 = n_elements(test[*,0])
;        n2 = n_elements(test[0,*])
;        temp = dblarr(n1,n2,nf1)
;        for k=0,nf1-1 do temp[*,*,k] = readfits(frames[k])
;        if (j eq 0) then data = dblarr(n1,n2,n0)
;        data[*,*,j] = median(temp,dimension=3,/even)
;        if (writefiles eq 'yes') then begin
;            name = strsplit(lists[j],'.',/extract)
;            name = name[0]
;            writefits,name+'F.fits',data[*,*,j]
;        endif
;    endfor
;    if (bgsub eq 'on') then begin
;        readcol,'bg.list',format='a',bgfiles
;        bgarr = dblarr(n1,n2,n_elements(bgfiles))
;        for j=0,n_elements(bgfiles)-1 do bgarr[*,*,j] = readfits(bgfiles[j])
;        bg = median(bgarr,dimension=3,/even)
;        if (writefiles eq 'yes') then writefits,'bg.fits',bg
;    endif
;endif
    
colors = intarr(n0)
cs = floor((255.-60.)/n0)
for j=0,n0-1 do colors[j] = 60 + (j*cs)        

set_plot,'x'
window,1,retain=2,xsize=n1/wscale,ysize=n2/wscale,ypos=50
window,0,retain=2,xsize=n1/wscale,ysize=n2/wscale

diameter = dblarr(n0)
rarr = fltarr(n1*n2,n0)
eearr = rarr
ttl375 = fltarr(n0)

for j=0,n0-1 do begin
    ellipse = fltarr(2,120)
    dd = dblarr(60)
    wset,0
    loadct,37
    oneimage = readfits(frames[j],/silent)
    onebg = readfits(bglist[j],/silent)
    smoothbg = smooth(onebg,7,/edge_truncate)
    oneimage = float(oneimage) - float(smoothbg)
    std = oneimage[bsort(oneimage)]
    cut = limit * median(std[n1*n2-100,*])
    ttl = total(oneimage)
    print,'The total counts on chip '+strn(ttl)
    print,'Considering only pixels greater than '+strn(cut)+' when ellipse fitting.'
    index = where(oneimage gt cut)
    TVImage, BytScl(oneimage, Top=!D.Table_Size-3)
    ellipse[*,*] = Fit_Ellipse(index, XSize=n1, YSize=n2, CENTER=c1, AXES=ax)
    PLOTS, ellipse[*,*]/wscale, Color=FSC_Color('yellow'),thick=2,/device
    center[*,j] = round(c1)

    FOR k=0,59 DO BEGIN
        x1 = ellipse[0,k]
        y1 = ellipse[1,k]
        x2 = ellipse[0,k+59]
        y2 = ellipse[1,k+59]
        dx = double(abs(x1-x2))
        dy = double(abs(y1-y2))
        dd[k] = sqrt(dx^2 + dy^2)
    ENDFOR
    diameter[j]=mean(dd)

    ; THE RADIUS ARRAY IS CREATED
    x2 = indgen(n1-center[0,j])+1
    x1 = reverse(indgen(center[0,j]))
    xline = [x1,x2]
    y2 = indgen(n2-center[1,j])+1
    y1 = reverse(indgen(center[1,j]))
    yline = transpose([y1,y2])
    Xarray = fltarr(n1,n2)
    Yarray = Xarray
    
    FOR i=0,n1-1 DO Xarray(*,i) = xline
    FOR i=0,n2-1 DO Yarray(i,*) = yline
    
    rad = sqrt(Xarray^2 + Yarray^2)
    IF (rad[center[0,j]-1,center[1,j]-1] NE 0) THEN BEGIN
        print,'The radius array is wrong.'
        print,'The central value is not zero.'
        print,'Check your code!'
        stop
    ENDIF
    
    if (j eq 0) then begin
        print,'Generate a postscript plot of the EE curves? (y/n)'
        read,ans
    endif

    r=starting_radius
    radius = r
    enclosed = 0.0
    repeat begin
        radius = [radius,r]
        temp = total(oneimage[where(rad le r)])/ttl
        enclosed = [enclosed,temp]
        print,'Current EE at R='+strn(r)+' pixels is '+strn(temp*100)
        r=r+stepsize
    endrep until temp ge ee

    enclosed = enclosed[1:*]
    radius = radius[1:*]

    i375 = where(radius eq 375.0)
    ttl375[j] = enclosed[i375[0]]

;    rarr[*,j] = radius
;    eearr[*,j] = enclosed

    mask = oneimage
    mask[where(rad gt r)] = 0.0
    wset,1
    loadct,37
    TVImage, BytScl(mask, Top=!D.Table_Size-3)

;The plotting routines start here.
    if (j eq 0) then begin
        window,2,retain=2
        device,decomposed=0
        loadct,0
        plot,radius,enclosed,xtitle='Radius (pixels)',title='Artificial FRD Check',$
          ytitle='Encircled Energy',thick=1,charsize=1.5
        xyouts,0.17,0.9,names[0],charsize=1.5,/normal
    endif
    if (j gt 0) then begin
        wset,2
        loadct,4
        oplot,radius,enclosed,color=colors[j-1],thick=1
        xyouts,0.17,(0.9-(j*0.05)),names[j],color=colors[j-1],charsize=1.5,/normal
    endif
    if (ans eq 'y' and j eq 0) then begin
        set_plot,'ps'
        device,file='Afrd_check.ps',/color
        loadct,0
        plot,radius,enclosed,xtitle='Radius (pixels)',title='Artificial FRD Check',$
          ytitle='Encircled Energy',thick=2,charsize=1.2,xthick=3,ythick=3,charthick=2,$
          /nodata,xrange=[min(radius),max(radius)+25],/xstyle
        loadct,4
        oplot,radius,enclosed,thick=2,color=colors[j]
        xyouts,0.18,0.87,names[j],charsize=1.2,charthick=2,color=colors[j],/normal
        set_plot,'x'
    endif
    if (ans eq 'y' and j gt 0) then begin
        set_plot,'ps'
        loadct,4
        oplot,radius,enclosed,color=colors[j],thick=2
        xyouts,0.18,(0.87-(j*0.05)),names[j],color=colors[j],charsize=1.2,charthick=2,/normal
        set_plot,'x'
    endif
    if (j eq n0-1) then begin
        set_plot,'ps'
        device,/close_file
        set_plot,'x'
    endif
endfor

print,'Next ENTER deletes the frames...'
pause
wdelete & wdelete

loadct,0
plot,ttl375,psym=2,/ynozero

stop

END
