PRO FRDcog, list

; This is, as its name implies, the culmination and simplification of
; the FRD analysis using the curve-of-growth or encircled energy
; analysis. It's a technique to qualify the amount of FRD of a fiber.

; LIMIT: The % used to make the pixel cut for the ellipse fitting
; routine. The purpose is to limit the ellipse fitting routine to only
; consider pixel values above a certain value. The default is 30% of
; the max value. (Actually, it gets multiplied by the median of the
; highest 100 pixel values.)

; The routine will iterate on this step, in steps of 3%, in order to
; get a more robust determination of the center.

; Many of the parameters that were free are now hard-coded below.

; EE: The encircled energy the routine will continue working out to.

; RADIUS: This is the value of the starting radius. This is set low
; (5) so that you can see the inner region of the fiber profile as
; well as the outer.

; STEP: The radius step size. 0.3 is good for the 512x512 camera and
; 2.0 for the FLI (3056x3056) camera. Reduce this value to increase
; the resolution of the steps.

; LIST: This is a list, with two columns, giving the light and dark
; frames, respectively. THE FIRST FILE IN THIS LIST SHOULD BE A
; BASELINE FRAME. ALL COGs WILL BE NORMALIZED TO THIS FRAME.

;**********************************************************************
limit = 0.20 ;percentage of the max to use in the ellipse-fitting
starting_radius = 1.0 ;the radius, in pixels, of where the COG begins
stepsize = 3.0 ;the size of the steps made when counting up EE
ee = 0.98 ;the EE the code will work up to, then stop
plotting = 'on'
wscale = 6.0 ;the window scale. Used for plotting purposes.
slow = 'off'
;**********************************************************************

limitbg = limit * 2.0
readcol,list,format='a,a',frames,darks
n0 = n_elements(frames)
names = strarr(n0)

test = readfits(frames[0],/silent)
n1 = n_elements(test[*,0])
n2 = n_elements(test[0,*])

colors = intarr(n0)
cs = floor(255.0/n0)
for j=0,n0-1 do colors[j] =  (j*cs)        

rarr = fltarr(n1*n2,n0)
eearr = rarr
ttl375 = fltarr(n0)


for j=0,n0-1 do begin           ;a loop through each spot. This is just to locate the center of the spot for the start of the the EE calculation.
    print,'Working on frame '+frames[j]+'...'
    oneimage = readfits(frames[j],/silent)
    s = strsplit(frames[j],'.',/extract)
    names[j] = s[0]
    if (j eq 0) then n1 = n_elements(oneimage[*,0])
    if (j eq 0) then n2 = n_elements(oneimage[0,*])

    if (plotting eq 'on' and j eq 0) then begin
        window,0,retain=2,xsize=n1/wscale,ysize=n2/wscale
        device,decomposed=0
        loadct,27,/silent
        xwin0 = !x.window
        ywin0 = !y.window
    endif else begin
        wset,0
        !x.window = xwin0
        !y.window = ywin0
    endelse
    onebg = readfits(darks[j],/silent)
        print,xwin0,ywin0 & pause
    
                                ; The corners of both the image and
                                ; the background are checked. The
                                ; background is scaled and subtracted.
    scale1 = median([onebg[0:50,0:50],onebg[n1-51:n1-1,0:50],onebg[0:50,n2-51:n2-1],onebg[n1-51:n1-1,n2-51:n2-1]])
    scale2 = median([oneimage[0:50,0:50],oneimage[n1-51:n1-1,0:50],oneimage[0:50,n2-51:n2-1],oneimage[n1-51:n1-1,n2-51:n2-1]])
    scale = float(scale2 / scale1)
    print,scale & pause
    oneimage = float(oneimage) - (float(onebg) * scale)
    frametotal = total(oneimage)
    std = oneimage[bsort(oneimage)]
    peak = median(std[n1*n2-100:*])

    center = fltarr(2,11)
    if (j eq 0) then limitdo = limitbg else limitdo = limit
    for k=0,10 do begin
        cut = (limitdo + (0.03*k)) * peak
        index = where(oneimage gt cut)
        ellipse = Fit_Ellipse(index, XSize=n1, YSize=n2, CENTER=c1)
        center[*,k] = c1
        if (plotting eq 'on') then TVImage, BytScl(oneimage,Top=!D.Table_Size-3)
        if (plotting eq 'on') then plots, ellipse/wscale, Color=FSC_Color('blue'),$
          thick=2,/device
        wait,0.2
    endfor

;    curval,oneimage
    center = median(center,dim=2)
    center = round(center)
    plots,float(center[0]/wscale),float(center[1]/wscale),psym=sym(12),symsize=10
    print,center
    pause

    ; THE RADIUS ARRAY IS CREATED
    x2 = indgen(n1-center[0])+1
    x1 = reverse(indgen(center[0]))
    xline = [x1,x2]
    y2 = indgen(n2-center[1])+1
    y1 = reverse(indgen(center[1]))
    yline = transpose([y1,y2])
    Xarray = fltarr(n1,n2)
    Yarray = Xarray
    
    FOR i=0,n1-1 DO Xarray(*,i) = xline
    FOR i=0,n2-1 DO Yarray(i,*) = yline
    
    rad = sqrt(Xarray^2 + Yarray^2)
    IF (rad[center[0]-1,center[1]-1] NE 0) THEN BEGIN
        print,'The radius array is wrong.'
        print,'The central value is not zero.'
        print,'Check your code!'
        stop
    ENDIF
    
    r = starting_radius
    radius = r
    enclosed = 0.0
    repeat begin
        radius = [radius,r]
        temp = total(oneimage[where(rad le r)])/frametotal
        enclosed = [enclosed,temp]
        print,'Current EE at R='+strn(r)+' pixels is '+strn(temp*100)
        r = r + stepsize
    endrep until temp ge ee

    enclosed = enclosed[1:*]
    radius = radius[1:*]

    r50 = where(enclosed ge 0.50)
    r50 = r50[0]
    
    if (j eq 0) then begin
        r50m = r50
        Sradius = radius
    endif else begin
        Rscale = float(r50m)/float(r50)
        Sradius = radius * Rscale
    endelse

    print,Sradius & pause

    free_lun,5
    openw,5,names[j]+'_EE.txt'
    printf,5,['Rad','SRad','EE'],format='(a3,5x,a4,4x,a2)'
    for l=0,n_elements(radius)-1 do printf,5,strn(radius[l]),' ',strn(Sradius[l]),' ',strn(enclosed[l])
    free_lun,5

    if (plotting eq 'on') then begin

        if (j eq 0) then begin
            window,2,retain=2
            device,decomposed=0
            loadct,0,/silent
            plot,radius,enclosed,xtitle='Radius (pixels)',title='EE v radius',$
              ytitle='Encircled Energy',thick=1,charsize=1.5
            xyouts,0.17,0.9,names[0],charsize=1.0,/normal
            xwin2 = !x.window
            ywin2 = !y.window
        endif else begin
            wset,2
            !x.window = xwin2
            !y.window = ywin2
            loadct,33,/silent
            oplot,Sradius,enclosed,color=colors[j]
            xyouts,0.17,(0.9-(j*0.03)),names[j],color=colors[j-1],charsize=1.0,/normal
        endelse

        if (j eq n0-1) then begin
            pause
            window,2,retain=2,xsize=1000,ysize=650
            device,decomposed=0
            loadct,0,/silent
            plot,radius,enclosed,xtitle='Radius (pixels)',title='EE v radius',$
              ytitle='Encircled Energy',thick=1,charsize=1.5
            loadct,33,/silent
            for k=0,n0-1 do begin
                readcol,names[k]+'_EE.txt',f='f,f,f',rad,srad,ee
                oplot,rad,ee,color=colors[k]
                xyouts,0.17,(0.9-(k*0.03)),names[k],color=colors[k-1],charsize=1.0,/normal
            endfor
            pause
            loadct,0,/silent
            plot,Sradius,enclosed,xtitle='Shifted Radius (pixels)',title='EE v radius',$
              ytitle='Encircled Energy',thick=1,charsize=1.5
            loadct,33,/silent
            for k=0,n0-1 do begin
                readcol,names[k]+'_EE.txt',f='f,f,f',rad,srad,ee
                oplot,srad,ee,color=colors[k]
                xyouts,0.17,(0.9-(k*0.03)),names[k],color=colors[k-1],charsize=1.0,/normal
            endfor
        endif

    endif

    slicex = median(oneimage[*,center[1]-10:center[1]+10],dim=2)
    slicey = median(oneimage[center[0]-10:center[0]+10,*],dim=1)
    
    set_plot,'ps'
    device,file=names[j]+'_profiles.ps',/color
    loadct,0,/silent
    plot,slicex,xthick=2,ythick=2,charthick=2,xtitle='Pixels',ytitle='Cross Sectional Profile',$
      title='Fiber '+names[j],/nodata,yrange=[0,max([slicex,slicey])+200],xrange=[0,n1],$
      /xstyle,/ystyle
    loadct,4,/silent
    oplot,slicex,color=60,thick=2
    oplot,slicey,color=110,thick=2
    device,/close_file
    set_plot,'x'
    
endfor

print,'Next ENTER deletes the frames...'
pause
while !d.window ne -1 do wdelete, !d.window

stop
END
