; Version 7 is made for 512x512. This version accepts any size data.

PRO fratio_v8, listolists, EE=ee, RADIUS=starting_radius, STEP=stepsize

; LISTOLISTS: A list of the lists of the files you want to
; reduce. Each set of data should be for an individual fiber and
; wavelength. 

; EE: This is the encircled energy the routine will iterate out
; to. This sets the defination for your final output fratio value.
; RADIUS: The starting radius for your outward iteration. This is
; where the counting begins at the center of your spot.
; STEP: This is step size made when making your radius steps. 

;-------------------------------------------
step = 1.27 ; This is the step size, in centimeters, between your measurements.
difference = 3.0 ;This sets when the sorted stepping used to estia
;-------------------------------------------
IF (N_Elements(starting_radius) EQ 0) THEN starting_radius=10
IF (N_Elements(stepsize) EQ 0) THEN stepsize=0.03 ;This value proves a good one, showing minimal repeat.
IF (N_Elements(ee) EQ 0) THEN ee=0.95 

readcol,listolists,silent=1,format='a',lists
n0 = n_elements(lists)

for j=0,n0-1 do begin ;a loop through the lists
    readcol,lists[j],silent=1,format='a',files
    temp = readfits(files[0],/silent)
    n1 = n_elements(temp[*,0])
    n2 = n_elements(temp[0,*])
    n3 = n_elements(files)
    images = fltarr(n1,n2,n3)
    n4 = n_elements(images[*,*,0])

    std_images = dblarr(n4,n3)
    x_axis = indgen(n1,/long)
    totchip = fltarr(n3)
    tottrim = totchip
    ellipse = fltarr(2,120)
    d = double(step)
    center = fltarr(2,n3) ; the X and Y pixel coords of the center of the spot
    diameter = dblarr(n3)
    radius = diameter
    dd = dblarr(60)
    RA = dblarr(n1,n2,n3) ; The radius array
    ax = fltarr(2) & sax = ax
    width=intarr(n3)
    roitot=dblarr(n3)  ;The totals WITHIN THE ROI (which is the new box)
    data = dblarr(n4,2,n3)

; Two sets of totals are determined in this routine. One is the total
; for the entire chip. The other occurs by sorting, and attempting to
; catch the turnoff (as determined by some set number change in the
; counts...)

; THE DATA IS LOADED AND SORTED FROM LOW TO HIGH.  T1, THE FRAME TOTAL, IS
; ALSO CREATED
    for k=0,n3-1 do begin
        images[*,*,k] = readfits(files[k])
        oneimage = images[*,*,k]
        std_images[*,k] = oneimage[bsort(oneimage)]
        cntr = 0
        repeat begin 
            dif = std_images[cntr+1,k] - std_images[cntr,k]
            cntr = cntr + 1
        endrep until (dif gt difference)
        window,2,retain=2
        loadct,0
        plot,stdimage[*,k]
        loadct,4
        plots,cntr-1,stdimage[cntr-1,k],psym=1
        totchip = total(oneimage)
        tottrim = total(stdimage[cntr-1:*,k])
        print,totchip
        print,tottrim
        print,tottrim/totchip
        print,'The chip total, the trimmed total and their ratio...'
        pause
        loadct,37
        window,0,retain=2,xsize=n1,ysize=n2
        TVImage, BytScl(oneimage, Top=!D.Table_Size-3)
        cut = 0.07 * mean(std_images[n4-10:n4-1,k]) ;an arbitrary cut is made in order that the ellipse fitting routine (used only to locate a center) doesn't lose it's way.
        index = where(oneimage GE cut)
        t2 = total(oneimage(index))
        t3 = total(oneimage)
        percent = t2/t3
        ellipse[*,*] = Fit_Ellipse(index, XSize=n1, YSize=n2, CENTER=c1, AXES=ax)
        plots, ellipse[*,*], /Device, Color=FSC_Color('yellow'),thick=2
        center[*,k] = round(c1)
        
                                ; PROFILES OF THE X AND Y SECTIONS ARE GENERATED        
        FOR l=0,59 DO BEGIN
            x1 = ellipse[0,l]
            y1 = ellipse[1,l]
            x2 = ellipse[0,l+59]
            y2 = ellipse[1,l+59]    highi = where
            
            dx = double(abs(x1-x2))
            dy = double(abs(y1-y2))
            dd[j] = sqrt(dx^2 + dy^2)
        ENDFOR
        
        diameter[k]=mean(dd)
        radius[k]=diameter[k]/2
        
        xgraph = median(oneimage[*,center[1,k]-5:center[1,k]+5],dimension=2)
        ygraph = median(oneimage[center[0,k]-5:center[0,k]+5,*],dimension=2)
        xl = fltarr(500)+(center[0,k]-radius[k])
        xr = fltarr(500)+(center[0,k]+radius[k])
        yl = fltarr(500)+(center[1,k]-radius[k])
        yr = fltarr(500)+(center[1,k]+radius[k])
        yup = indgen(500)*100
        window,2,retain=2,xsize=633,ysize=480
        loadct,0
        plot,xgraph,title='X cross section',ytitle='Counts (ADU)',xtitle='The radius is :'+strn(radius[k]),$
          subtitle='The EE is :'+strn(percent),xrange=[0,n1],xstyle=1,yrange=[0,max(ygraph)],charsize=1.3
        loadct,4
        oplot,xl,yup,color=180,thick=3
        oplot,xr,yup,color=180,thick=3
        loadct,0
        window,3,retain=2,xsize=633,ysize=480
        plot,ygraph,title='Y cross section',ytitle='Counts (ADU)',xtitle='The radius is :'+strn(radius[k]),$
          subtitle='The EE is :'+strn(percent),xrange=[0,n2],xstyle=1,yrange=[0,max(ygraph)],charsize=1.3
        loadct,4
        oplot,yl,yup,color=180,thick=3
        oplot,yr,yup,color=180,thick=3
        print,'Take a gander.'
        pause
        wdelete,2,3
        
                                ; THE RADIUS ARRAY IS CREATED
        x2 = indgen(n1-center[0,k])+1
        x1 = reverse(indgen(center[0,k]))
        xline = [x1,x2]
        y2 = indgen(n2-center[1,k])+1
        y1 = reverse(indgen(center[1,k]))
        yline = transpose([y1,y2])
        Xarray = fltarr(n1,n2) & Yarray = Xarray
        
        FOR i=0,n1-1 DO Xarray(*,i) = xline
        FOR i=0,n2-1 DO Yarray(i,*) = yline
        
        RA[*,*,k] = sqrt(Xarray^2 + Yarray^2)
        IF (RA[center[0,k]-1,center[1,k]-1,k] NE 0) THEN BEGIN
            print,'The radius array is wrong.'
            print,'The central value is not zero.'
            print,'Check your code!'
            stop
        ENDIF

    ENDFOR
    wdelete,0

; THE FIRST WIDTH IS DETERMINED, SETTING THE RATIO FOR THE REMAINDER
; OF THE ROI'S
w2 = min([center[0,0],n1-center[0,0],center[1,0],n2-center[1,0]])
width[0] = 2*w2
ix1 = center[0,0]-w2
ix2 = center[0,0]+w2
iy1 = center[1,0]-w2
iy2 = center[1,0]+w2
ttls[0] = total(images[ix1:ix2,iy1:iy2,0])
window,0,retain=2,xsize=width[0],ysize=width[0]
loadct,37
TVImage, BytScl(images[ix1:ix2,iy1:iy2,0], Top=!D.Table_Size-3)
; THE REMAINDER OF THE WIDTHS AND ENCLOSED TOTALS ARE DETERMINED
FOR k=1,n3-1 DO BEGIN
    w = round((diameter[k]/diameter[0])*width[0])
    width[k] = w
    w2 = w/2 
    ix1 = center[0,k]-w2
    ix2 = center[0,k]+w2
    iy1 = center[1,k]-w2
    iy2 = center[1,k]+w2
    ttls[k] = total(images[ix1:ix2,iy1:iy2,k])
    window,k,retain=2,xsize=width[k],ysize=width[k]
    TVImage, BytScl(images[ix1:ix2,iy1:iy2,k], Top=!D.Table_Size-3)
ENDFOR

pause
wdelete & wdelete & wdelete & wdelete

length = intarr(n3)
n=(ee-0.70)*100
values = dblarr(n+1,3,n3)
fnumber=dblarr(n+1,8)
loadct,0

 FOR k=0,n3-1 DO BEGIN
    oneimage = images[*,*,k]
    rad = RA[*,*,k]
    r=starting_radius
    m=0
    REPEAT BEGIN ;THE ENCIRCLED COUNTS ARE DETERMINED HERE
        r=r+stepsize
        data[m,0,k]=r ;The top row of data is the radius
        enclosed = total(oneimage[Where(rad LE r)])
        percent = enclosed/ttls[k]
        print,'Doing frame:',k+1,percent
        data[m,1,k] = percent ;The bottom row is the encircled energy
        m=m+1
    ENDREP UNTIL percent GE ee
    oneimage[where(rad LE r)] = -1000 ;The encircled spots are set to black.
    loadct,37
    window,2+(k*4),xsize=width[k],ysize=width[k],retain=2
    TVImage, BytScl(oneimage, Top=!D.Table_Size-3)
    
    FOR j=0,n3 DO BEGIN
        v = 0.70 + (j*0.01)
        ink = where(data[*,1,k] GT v)
        ink = ink[0]
        values[j,2,k] = data[ink,1,k]
        values[j,0,k] = data[ink,0,k]
        values[j,1,k] = values[j,0,k]*0.02
        fnumber[j,0] = mean(values[j,2,*])
    ENDFOR

    length[k]=m-1

ENDFOR
print,'Next ENTER deletes the frames...'
pause
wdelete & wdelete & wdelete & wdelete
data = data[0:max(length),*,*]

;THE CORRESPONDING F/# FOR EACH EE VALUE IS DETERMINED BY TAKING AN
;AVERAGE OF ALL 6 POSSIBLE COMPARISONS

fnumber[*,1]=d/(values[*,1,0]*2.0-values[*,1,1]*2.0)
fnumber[*,2]=d/(values[*,1,1]*2.0-values[*,1,2]*2.0)
fnumber[*,3]=d/(values[*,1,2]*2.0-values[*,1,3]*2.0)
fnumber[*,4]=(d+d)/(values[*,1,0]*2.0-values[*,1,2]*2.0)
fnumber[*,5]=(d+d)/(values[*,1,1]*2.0-values[*,1,3]*2.0)
fnumber[*,6]=(d+d+d)/(values[*,1,0]*2.0-values[*,1,3]*2.0)

frange = dblarr(30,4)
FOR j=0,n_elements(fnumber[*,0])-1 DO BEGIN
    temp = fnumber[j,1:6]
    std = sort(temp)
    std = std[1:4] ;Effectively a minmax rejection. This requires 6 values for fnumber.
    fnumber[j,7] = mean(temp[std])
    frange[j,*] = temp[std]
ENDFOR

;THE FOLLOWING IS JUST TO ESTABLISH PLOTTING LIMITS

xulim = max(fnumber[*,7])+0.1
xllim = min(fnumber[*,7])-0.1
lx1=[3.65,3.65]
ly1=[0.6,1.1]
lx2=[fnumber[25,7],fnumber[25,7]]

print,'Enter the title of the plot and output data files:'
read,titl

window,2,retain=2
device,decomposed=0
loadct,0
plot,fnumber[*,7],fnumber[*,0],psym=7,xrange=[xulim,xllim],yrange=[0.6,1.1], $
  title=titl,xstyle=1,$
  xtitle='Output F/#',ytitle='Encircled Energy',charsize=1.3,subtitle='with cover plates'
loadct,4
FOR j=0,n_elements(fnumber[*,0])-1 DO BEGIN
    mnmx = [max(fnumber[j,1:6]),min(fnumber[j,1:6])]
;    mnmx = [max(frange[j,*]),min(frange[j,*])]
    ind = [fnumber[j,0],fnumber[j,0]]
    oplot,mnmx,ind,color=155
ENDFOR
oplot,lx1,ly1,color=60,thick=2
xyouts,3.67,0.63,'INPUT BEAM F/#',orientation=90,color=60,charsize=1.2
oplot,lx2,ly1,color=110,thick=2
xyouts,fnumber[25,7]-0.05,0.63,'95% ENCIRCLED ENERGY',orientation=90,color=110,charsize=2
xyouts,fnumber[25,7]-0.1,0.7,'(f/#: '+strn(fnumber[25,7])+')',orientation=90,color=110,charsize=1.5

print,'Save a hardcopy?'
read,ans
ans='y'

IF (ans EQ 'y') THEN BEGIN
set_plot,'ps'
on3=titl+'.ps'
device,filename=on3,/color
loadct,0
plot,fnumber[*,7],fnumber[*,0],psym=7,xrange=[xulim,xllim],yrange=[0.6,1.1],charthick=2, $
  title=titl,xstyle=1,xtitle='Output F/#',ytitle='Encircled Energy',charsize=1.1,xthick=2,$
  ythick=2,thick=2
loadct,4
FOR j=0,n_elements(fnumber[*,0])-1 DO BEGIN
    mnmx = [max(fnumber[j,1:6]),min(fnumber[j,1:6])]
;    mnmx = [max(frange[j,*]),min(frange[j,*])]
    ind = [fnumber[j,0],fnumber[j,0]]
    oplot,mnmx,ind,color=155,thick=3.0
ENDFOR
oplot,lx1,ly1,color=60,thick=3
xyouts,3.67,0.63,'INPUT BEAM F/3.65',charthick=2.0,orientation=90,color=60,charsize=1.0
oplot,lx2,ly1,color=110,thick=3
xyouts,fnumber[25,7]-0.04,0.63,'95% ENCIRCLED ENERGY',orientation=90,color=110,charsize=1.3,charthick=2.0
xyouts,fnumber[25,7]-0.07,0.7,'(f/#: '+strn(fnumber[25,7])+')',orientation=90,color=110,charthick=2.0

device,/close_file
set_plot,'x'
ENDIF

;on1=titl+'.data'
on2=titl+'.frd'

;openw,10,on1
;printf,10,data
openw,11,on2
printf,11,transpose(fnumber)
free_lun,11
print,'The files are:'
print,transpose(files)
print,'The totals defined as 100% EE:'
print,ttls
print,'The ratio of these values to the total chip values:'
print,ttls/t3
pause
wdelete,2

print,'Word.'

stop

END
