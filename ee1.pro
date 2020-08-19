PRO ee1, frame, Fiber=fiber, EE=ee, RAD=starting_radius, STEP=stepsize

;THIS CODE CONTAINS THE GUTS OF FRATIO.PRO, YET IT'S DESIGNED TO
;CALCULATE AN EE PLOT FOR ONE FIBER, RATHER THAN MANY.  IN THIS WAY
;IT'S IDENTICAL TO FRATIO.PRO, ONLY THE X-AXIS CAN NOT BE CONVERTED TO
;F/# (as there's only 1 image) AND THERE IS NO WAVELENGTH DEPENDENCE.

winsize = 7.0 ;set to 6.0 or 7.0 for FLI camera and 1.0 for 512x512 data
cut = 5000
;*****************************************************************

IF (N_Elements(starting_radius) EQ 0) THEN starting_radius = 150
IF (N_Elements(stepsize) EQ 0) THEN stepsize = 2.0
IF (N_Elements(ee) EQ 0) THEN ee = 0.98 
IF (N_Elements(frame) EQ 0) THEN BEGIN
    frame=''
    print,'Enter a data frame (W/O THE .FIT):'
    read,frame
    frame = frame+'.fit'
ENDIF

image = readfits(frame)
image = float(image)

test = median(image)
if (test gt 300.0) then begin
    dark = ''
    print,'This frame does not appear to be BG subtracted.'
    print,'If you do not want nonsense results,'
    print,'enter the name of a dark frame (W/O .fit)'
    read,dark
    bg = readfits(dark+'.fit')
    bg = float(bg)
endif

n1 = n_elements(image[0,*])
n2 = n_elements(image[*,0])
in = bsort(image)
std = image[in]
ellipse = fltarr(2,120)

ans=''
flname=''
center = fltarr(2)
dd = dblarr(60)
RA = dblarr(n1,n2)
data = fltarr(n1*n2,2)

loadct,37
window,0,retain=2,xsize=n1/winsize,ysize=n2/winsize
TVImage, BytScl(image, Top=!D.Table_Size-3)

index = Where(image GE cut)
    t2 = total(image(index))
    percent = t2/total(image)
    ellipse[*,*] = Fit_Ellipse(index, XSize=n1, YSize=n2, CENTER=c1)
    PLOTS, ellipse/winsize, /Device, Color=FSC_Color('yellow'),thick=2
        center = round(c1)

    ; PROFILES OF THE X AND Y SECTIONS ARE GENERATED        
        FOR j=0,59 DO BEGIN
            x1 = ellipse[0,j]
            y1 = ellipse[1,j]
            x2 = ellipse[0,j+59]
            y2 = ellipse[1,j+59]
            dx = double(abs(x1-x2))
            dy = double(abs(y1-y2))
            dd[j] = sqrt(dx^2 + dy^2)
        ENDFOR

        diameter=mean(dd)
        radius=diameter/2
        diff = max(dd)-min(dd)

;A PLOT OF THE VARIATION IN DD IS GENERATED
loadct,0
window,2,retain=2
plot,dd,psym=2,/ynozero,title='60 different estimates of the diameter',xtitle='Range of Diameter: '+strn(diff),charsize=1.3
print,'The center is: '+strn(center)
print,'Next ENTER deletes these windows...'
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

RA[*,*] = sqrt(Xarray^2 + Yarray^2)
IF (RA[center[0]-1,center[1]-1] NE 0) THEN BEGIN
    print,'The radius array is wrong.'
    print,'The central value is not zero.'
    print,'Check your code!'
    stop
ENDIF


r=starting_radius
m=0
REPEAT BEGIN                 ;THE ENCIRCLED COUNTS ARE DETERMINED HERE
    r=r+stepsize
    data[m,0]=r               ;The top row of data is the radius
    enclosed = total(image[Where(RA LE r)])
    percent = enclosed/total(image)
    print,'Encircled energy at '+strn(r)+' is '+strn(percent)
    data[m,1] = percent        ;The bottom row is the encircled energy
    m=m+1
ENDREP UNTIL percent GE ee

data = data[0:m-1,*]

;THE RESULTS ARE PLOTTED
IF (n_elements(fiber) EQ 0) THEN BEGIN
    t = strsplit(frame,'.',/extract)
    ttl = 'Encircled Energy: '+t[0]
    tout = t[0]
ENDIF
IF (n_elements(fiber) NE 0) THEN BEGIN
    ttl='Encirled Energy: '+fiber
    tout = fiber
ENDIF
window,0,retain=2
plot,data[*,0],data[*,1],title=ttl,xtitle='Radius (pixels)',$
  ytitle='Encircled Energy',charsize=1.3

print,'Save the data and a hardcopy?'
read,ans
IF (ans EQ 'y') THEN BEGIN
    set_plot,'ps'
    device,filename = tout+'.ps'
    plot,data[*,0],data[*,1],title=ttl,xtitle='Radius (pixels)',$
      ytitle='Encircled Energy',charsize=1.3,xthick=3,ythick=3,$
      charthick=3,thick=3
    device,/close_file
    set_plot,'x'
    
    flname=strsplit(flname,'.',/extract)
    on1 = flname[0]+'_rad.txt'
    on2 = flname[0]+'_EEvR.txt'
    openw,5,on1
    printf,5,dd
    openw,6,on2
    printf,6,transpose(data)
    free_lun,5,6
ENDIF

print,'Word.'

END

