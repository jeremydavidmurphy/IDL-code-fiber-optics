FUNCTION eeF, frameL, darkL, srad=SRad, SS=stepsize, EE=eee

; This code is a re-working of eeL, which stands separately from
; trans_FRD, which is what's primarily used to call this function.

; frameL: the list of frames
; darkL: the list of corresponding dark frames
;*****************************************************************
winsize = 8.0 ;set to 6.0 or 7.0 for FLI camera and 1.0 for 512x512 data
cut = 5000
ans = ''
range = 300 ; the value, in pixels, to extend the diameter of the spot
;*****************************************************************

IF (N_Elements(SRad) EQ 0) THEN SRad = 10
IF (N_Elements(stepsize) EQ 0) THEN stepsize = 2.0
IF (N_Elements(eee) EQ 0) THEN eee = 0.98

readcol,frameL,format='a',datalist
readcol,darkL,format='a',darklist
if (n_elements(datalist) ne n_elements(darklist)) then stop

n0 = n_elements(datalist)
test = readfits(datalist[0])
n1 = n_elements(test[*,0])
n2 = n_elements(test[0,*])
names = strarr(n0)
counts = intarr(n0)
colors = intarr(n0)

data = fltarr(floor(n1*n2),2,n0)

for j=0,n0-1 do colors[j] = 60 + (j * floor(195.0/n0))

for j=0,n0-1 do begin ;a loop through each fiber image
    fname = datalist[j]
    image = readfits(fname)
    image = float(image)
    
    dname = darklist[j]
    dark = readfits(dname)
    dark = float(dark)

    t = strsplit(fname,'.',/extract)
   names[j] = t[0]

    t1 = image[bsort(image)]
    t2 = dark[bsort(dark)]
    t1m = median(t1[0:1000])
    t2m = median(t2[0:1000])
    adj = t1m / t2m
    
    print,''
    print,'Dark frame '+dname+' is being scaled by '+strn(adj)
    dark = dark * adj
    
    image = image - dark
    
    ellipse = fltarr(2,120)
    center = fltarr(2)
    dd = dblarr(60)
    RA = dblarr(n1,n2)
    
    loadct,37
    window,3,retain=2,xsize=n1/winsize,ysize=n2/winsize,ypos=50
    TVImage, BytScl(image, Top=!D.Table_Size-3)
    
    index = where(image GE cut)
    ellipse[*,*] = Fit_Ellipse(index, XSize=n1, YSize=n2, CENTER=c1)
    center = round(c1)
    plots,center,psym=5
    plots, ellipse/winsize, /Device, Color=FSC_Color('yellow'),thick=2

    for k=0,59 do begin
        x1 = ellipse[0,k]
        y1 = ellipse[1,k]
        x2 = ellipse[0,k+59]
        y2 = ellipse[1,k+59]
        dx = double(abs(x1-x2))
        dy = double(abs(y1-y2))
        dd[k] = sqrt(dx^2 + dy^2)
    endfor

    diameter = mean(dd)
    radius = diameter/2.0
    print,diameter
    print,radius
    radius = floor(radius)
    print,radius
    diff = max(dd)-min(dd)
    push = radius + range

; THE IMAGE IS TRIMMED SO THE RATIO OF BG TO FIBER IS THE SAME
    trimmed = image[center[0]-push:center[0]+push,center[1]-push:center[1]+push]
    n3 = n_elements(trimmed[*,0])
    if (n3 gt n1) then n3 = n1
    n4 = n_elements(trimmed[0,*])
    if (n4 gt n2) then n4 = n2
    window,7,retain=2,ypos=50,xpos=400,xsize=n3/winsize,ysize=n4/winsize
    TVImage, BytScl(trimmed, Top=!D.Table_Size-3)

    imagettl = total(trimmed)

      ;A PLOT OF THE VARIATION IN DD IS GENERATED
    loadct,0
    window,2,retain=2
    plot,dd,psym=2,/ynozero,title='60 different estimates of the diameter',xtitle='Range of Diameter: '+strn(diff),charsize=1.3
    print,'The center is: '+strn(center)

; THE RADIUS ARRAY IS CREATED
    x2 = indgen(n1-center[0])+1
    x1 = reverse(indgen(center[0]))
    xline = [x1,x2]
    y2 = indgen(n2-center[1])+1
    y1 = reverse(indgen(center[1]))
    yline = transpose([y1,y2])
    Xarray = fltarr(n1,n2)
    Yarray = Xarray

    for i=0,n1-1 do Xarray(*,i) = xline
    for i=0,n2-1 do Yarray(i,*) = yline

    RA[*,*] = sqrt(Xarray^2 + Yarray^2)
    if (RA[center[0]-1,center[1]-1] ne 0) then begin
        print,'The radius array is wrong.'
        print,'The central value is not zero.'
        print,'Check your code!'
        stop
    endif

    r = SRad
    m = 0
    repeat begin             ;THE ENCIRCLED COUNTS ARE DETERMINED HERE
        r = r + stepsize
        data[m,0,j] = r           ;The top row of data is the radius
        enclosed = total(image[where(RA LE r)])
        percent = enclosed/imagettl
        print,'Encircled energy at '+strn(r)+' is '+strn(percent)
        data[m,1,j] = percent    ;The bottom row is the encircled energy
        m = m + 1
    endrep until percent ge eee
    counts[j] = m - 1
;    data = data[0:m-1,*]

    data[*,1,j] = data[*,1,j] * 100.0

endfor

window,0,retain=2
device,decomposed=0
loadct,0
for j=0,n0-1 do begin
    radius = data[0:counts[j],0,j]
    eep = data[0:counts[j],1,j]
    if (j eq 0) then begin
        plot,radius,eep,xtitle='Radius (pixels)',$
          ytitle='Encircled Energy',charsize=1.3,/nodata
    endif
    loadct,4
    oplot,radius,eep,color=colors[j]
    xyouts,0.8,0.5-(j*0.03),names[j],/normal,charsize=1.0,color=colors[j]
endfor

print,'Word.'

stop
END

