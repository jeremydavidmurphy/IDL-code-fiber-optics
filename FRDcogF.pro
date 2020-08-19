; This version (sept 17th, 2010) is based on fratioF, and now explores
; the output f-ratio of the center of the fiber spot, as more of a
; curve-of-growth analysis. Certain parameters that were user
; controled (EE, starting radius, etc) in fratioF are now internal.

; This code now samples the EE at every 0.5%, starting at 1% EE and
; terminating at 98% EE.

; This function is a reworking into functional form of the code
; fratio_v8 and is used primarily in trans_FRD.pro.

; The change between fratio_longF and fratioF is that now the ee
; counts are kept at each EE step, not just the last 30. Also, the COG
; can be written out as a text file.

; DATA: A 3-D data array of spots. THEY MUST BE ORDERED LARGEST TO SMALLEST!
; CSTEP: The size of the camera steps, in mm.
; PSIZE: pixel size, in mm 
; EE: Encircled Energy. This the the percentage the code will work up
; to and then stop. The default is 0.98
; RADIUS: THIS IS NO LONGER USED. 
; RSTEPSIZE: The size of the steps, in pixels, taken to reach the EE
; value. This is the pixel step size when growing the EE circle. This
; doesn't need to be an integer as not all the radius values are
; integers. NOTE: THIS IS NOT THE STEPSIZE OF THE CAMERA STEPS! If not
; used, the default is set by the size of the data array fed in.
; COG: If you enter name here the curve-of-growth for the largest fiber
; spot (i.e. the first spot in data) will be written out as a text file
; with that name.

;MODIFIED on Sept 15th, 2010: Includes an option to output the
;curve-of-growth for a given frame.
function FRDcogF, data, cstep, psize, winsize
on_error,2
;******************************************************************
plottingF = 'on'
;plottingF = 'off'
;******************************************************************

starting_radius = 3.0
ee = 0.98
eestart = 0.01
nn4 = 195
Rstepsize = 0.25

nn0 = n_elements(data[*,0,0])
nn1 = n_elements(data[0,*,0])
nn2 = n_elements(data[0,0,*])
nn3 = nn0 * nn1

;*********************************************************************
FRvaluesF = fltarr(nn4,4,nn2)
;discrete EE steps of 0.01
;row 1: radius (pixels)
;row 2: EE
;row3 and 4 are funny in that they are not nn2 values deep
;frvalues[*,2,0] = output f/#
;frvalues[*,2,1] = EE (the mean of each of the step positions
;frvalues[*,2,2] = X^2 goodness of fit (data fit to lines for the
;geometry of the output cone.
;frvalues[0:nn2-1,3,0] = total counts of the trimmed chip
;frvalues[0:nn2-1,3,1] = % of the light KEPT (trimmed/total)
;D3: camera positions
;*********************************************************************

x_axis = indgen(nn0,/long)
tottrim = fltarr(nn2)
ellipse = fltarr(2,120)
center = fltarr(2,nn2) ;The X and Y pixel coords of the center of the spot
diameter = dblarr(nn2)
radius = diameter
dd = dblarr(60)
RA = dblarr(nn0,nn1) ;The radius array
width = dblarr(nn2)

if (plottingF eq 'on') then begin
    set_plot,'x'
endif

igood  = 0
for kk=0,nn2-1 do begin         ;A LOOP THROUGH CAMERA POSITIONS.
    oneimage = data[*,*,kk]
    stdimage = oneimage[bsort(oneimage)]
    if (median(stdimage[nn3-101:*]) gt 65000.0) then begin
        FRvaluesF[*,0,kk] = fltarr(nn4) - 666
        FRvaluesF[*,1,kk] = fltarr(nn4) - 666
        print,'Spot #'+strn(kk+1)+' is SATURATED!!!!'
        goto, jump1F ;steps over saturated frames
    endif
    igood = [igood,kk]
    cut = 0.20 * mean(stdimage[nn3-1000:nn3-10])
    print,'The cutoff radius value is '+strn(cut)
    index = where(oneimage GE cut)
    t2 = total(oneimage(index))
    t3 = total(oneimage)
    percent = t2/t3
    ellipse[*,*] = Fit_Ellipse(index, XSize=nn0, YSize=nn1, CENTER=c1)
    center[*,kk] = round(c1)
    if (plottingF eq 'on') then begin
        loadct,37
        window,0,retain=2,xsize=nn0/winsize,ysize=nn1/winsize
        TVImage, BytScl(oneimage, Top=!D.Table_Size-3)
        plots,ellipse[*,*]/winsize,/Device,Color=FSC_Color('yellow'),thick=2
    endif

; PROFILES OF THE X AND Y SECTIONS ARE GENERATED
    for ll=0,59 do begin
        x1 = ellipse[0,ll]
        y1 = ellipse[1,ll]
        x2 = ellipse[0,ll+59]
        y2 = ellipse[1,ll+59]
        dx = double(abs(x1-x2))
        dy = double(abs(y1-y2))
        dd[ll] = sqrt(dx^2 + dy^2)
    endfor
    
    diameter[kk]=mean(dd)
    radius[kk]=diameter[kk]/2.0

; THE RADIUS ARRAY IS CREATED
    x2 = indgen(nn0-center[0,kk])+1
    x1 = reverse(indgen(center[0,kk]))
    xline = [x1,x2]
    y2 = indgen(nn1-center[1,kk])+1
    y1 = reverse(indgen(center[1,kk]))
    yline = transpose([y1,y2])
    Xarray = fltarr(nn0,nn1)
    Yarray = Xarray
    
    for ii=0,nn0-1 do Xarray(*,ii) = xline
    for ii=0,nn1-1 do Yarray(ii,*) = yline
    
    RA[*,*] = sqrt(Xarray^2 + Yarray^2)
    if (RA[center[0,kk]-1,center[1,kk]-1] ne 0) then begin
        print,'The radius array is wrong.'
        print,'The central value is not zero.'
        print,'Check your code!'
        stop
    endif

; THE FIRST WIDTH IS DETERMINED, SETTING THE RATIO FOR THE REMAINDER
; OF THE ROI'S
    if (kk eq 0) then begin
        w2 = min([center[0,0],nn0-center[0,0],center[1,0],nn1-center[1,0]])
        width[0] = 2*w2
        ix1 = center[0,0]-w2+1.0
        ix2 = center[0,0]+w2-1.0
        iy1 = center[1,0]-w2+1.0
        iy2 = center[1,0]+w2-1.0
        tottrim[0] = total(data[ix1:ix2,iy1:iy2,0])
        tot = total(data[*,*,0])
        print,strn(tottrim[0]/tot)+'% of the light for frame 1 is being kept.'
        FRvaluesF[kk,3,0] = tottrim[0]
        FRvaluesF[kk,3,1] = tottrim[0]/tot
        if (plottingF eq 'on') then begin
            window,2,retain=2,xsize=width[0]/winsize,ysize=width[0]/winsize
            loadct,37
            TVImage, BytScl(data[ix1:ix2,iy1:iy2,0], Top=!D.Table_Size-3)
        endif
    endif
    if (kk gt 0) then begin
        width[kk] = round((diameter[kk]/diameter[0])*width[0])
        w2 = width[kk]/2.0 
        ix1 = center[0,kk]-w2+1
        ix2 = center[0,kk]+w2-1
        iy1 = center[1,kk]-w2+1
        iy2 = center[1,kk]+w2-1
        tottrim[kk] = total(data[ix1:ix2,iy1:iy2,kk])
        tot = total(data[*,*,kk])
        print,strn(tottrim[kk]/tot)+'% of the light for frame '+strn(kk+1)+' is being kept.'
        FRvaluesF[kk,3,0] = tottrim[kk]
        FRvaluesF[kk,3,1] = tottrim[kk]/tot
        if (plottingF eq 'on') then begin
            window,(kk*4)+2,retain=2,xsize=width[kk]/winsize,ysize=width[kk]/winsize
            loadct,37
            TVImage, BytScl(data[ix1:ix2,iy1:iy2,kk], Top=!D.Table_Size-3)
        endif
    endif

    r = starting_radius
    banjo = eestart
    guitar = 0
;    if (n_elements(cogname) ne 0 and kk eq 0) then begin
;        free_lun,5
;        openw,5,cogname
        
;    endif else 
    repeat begin             ;THE ENCIRCLED COUNTS ARE DETERMINED HERE
        enclosed = total(oneimage[where(RA le r)])
        percent = enclosed/tottrim[kk]
;        if (kk eq 0) then printf,5,[r,percent]
        print,'Encircled energy at '+strn(r)+' is '+strn(percent)
        if (percent ge banjo) then begin
            FRvaluesF[guitar,0,kk] = r
            FRvaluesF[guitar,1,kk] = percent
            guitar = guitar + 1
            banjo = banjo + 0.005
        endif
        r = r + Rstepsize
    endrep until percent ge ee
;    if (kk eq 0) then free_lun,5

    if (plottingF eq 'on') then begin
        window,1,retain=2,ypos=50,xsize=400,ysize=400
        plot,FRvaluesF[*,0,kk],FRvaluesF[*,1,kk],psym=1,$
          yrange=[eestart-0.01,percent+0.01],/ystyle,$
          xtitle='Radius (pixels)',ytitle='EE',charsize=1.2
    endif
    print,'Spot #'+strn(kk+1)+' is completed.'

    jump1F:

endfor

; all the EE for a given radius is calculated. the next part uses the
; spacing between camera steps and pixel size to calculate the output f-number.

xx = fltarr(nn2)
for j=0,nn2-1 do xx[j] = j*cstep
r = reverse(xx)
xx = reverse(r[igood])

if (plottingF eq 'on') then begin ;just setting the plotting range...
    rr = FRvaluesF[nn4-1,0,igood] ;radius of each step size, in pixels
    rrr = fltarr(n_elements(rr))
    for k=0,n_elements(rr)-1 do rrr[k] = rr[0,0,k] ;just formatting
    rrr = rrr * psize * 1000 ;radius converted to mm
    rr = reverse(rrr)
    rr = rr - rr[0]
    yymax = max(rr)
    print, yymax
endif


for jj=0,nn4-1 do begin
    rr = FRvaluesF[jj,0,igood] ;radius of each step size, in pixels
    rrr = fltarr(n_elements(rr))
    for k=0,n_elements(rr)-1 do rrr[k] = rr[0,0,k] ;just formatting
    rrr = rrr * psize * 1000.0 ;radius converted to mm
    rr = reverse(rrr)
    rr = rr - rr[0]
    eeout = FRvaluesF[jj,1,igood]
    
    ab = linfit(xx,rr,yfit=yout,chisq=chifit)
    FRvaluesF[jj,2,2] = chifit
    frtop = max(xx)
    frbot = max(yout) * 2.0

    if (plottingF eq 'on') then begin
        xxplot = xx
        yyplot = fltarr(n_elements(xxplot))
        for k=0,n_elements(xxplot)-1 do yyplot[k] = ab[0] + (ab[1] * xxplot[k])
        if (jj eq 0) then begin
            window,3,retain=2,ypos=50
            plot,xx,rr,xrange=[min(xx)-1,max(xx)+1],yrange=[-0.5,yymax],$
              /xstyle,psym=2,xtitle='Camera Position (mm)',ytitle='Spot Radius (mm)',$
              title='F-ratio Calculation'
            oplot,xxplot,yyplot
        endif else begin
            wset,3
            oplot,xx,rr,psym=2
            oplot,xxplot,yyplot
        endelse
    endif

    fratio = frtop / frbot
    print,'The output F-ratio at '+strn(median(eeout,/even))+' is '+strn(fratio)
    FRvaluesF[jj,2,1] = median(eeout,/even)
    FRvaluesF[jj,2,0] = fratio
endfor

if (plottingF eq 'on') then wait,2

if (plottingF eq 'on') then begin
    window,0,retain=2,xpos=200,ypos=200,xsize=800,ysize=700
    device,decomposed=0
    loadct,0
    plot,FRvaluesF[*,2,0],FRvaluesF[*,2,1],title='Output F-ratio Values',$
      xtitle='Output F-ratio',ytitle='Encircled Energy',charsize=1.2,$
      yrange=[FRvaluesF[0,2,1]-0.01,FRvaluesF[nn4-1,2,1]+0.01],/ystyle,$
      xrange=[FRvaluesF[0,2,0]+0.1,FRvaluesF[nn4-1,2,0]-0.1],/xstyle,psym=2
;    print,'The next ENTER deletes all open windows...'
;    pause
    wait,3
    while !d.window ne -1 do wdelete, !d.window
endif

jump2F:
;the output is created (out). row 1 = EE, row 2 = output f-ratio
;row 3 = chisq, row 4 = trimmed chip totals
;row 5 = % of total chip light kept
out = [[FRvaluesF[*,2,1]],[FRvaluesF[*,2,0]],[FRvaluesF[*,2,2]],[FRvaluesF[*,3,0]],[FRvaluesF[*,3,1]]]

return,out
stop
end
