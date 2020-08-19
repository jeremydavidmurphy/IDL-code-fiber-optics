PRO fiber_finder, fitsfile

; locates a set number of "peaks" in a fits file. Returns position and
; EE for each.

; This code has the guts of the find.pro routine. I just cut and
; pasted the entire code into the body of this code, then modified it
; to hardcode in some of the free parameters.

;*********************************************************************
plotting1 = 'on'
plotting2 = 'off' ;this one turns on the ellipse-fitting routine, which is VERY slow if you turn on 'slow' below.
slow = 'on' ;turn on to slow the ellipse-fitting routine down for better visual inspection.
fwhm = 22.0 ; the FWHM in pixels...?
print = 'fib_loc.txt'
hmin = 300 ;the detection threshold, in pixel values. This is the counts above the background to consider a peak a peak
sharplim = [0.1,2.0]
roundlim = [-1.0,1.0]
nfibers = 448 ;you'll need to dial this down if you have broken fibers
pitch = 42.0 ;an estimate of the pitch between fibers
hp = floor(pitch/2.0)
fiberR = 16.5 ;the radius to consider the flux for a fiber. 
apR = 20.0 ;the outer radius for aperture photometry
;*********************************************************************

;*********************************************************************
;********* The beginning of the stripped find.pro code ***************
;*********************************************************************

image = readfits(fitsfile,/silent) ;JDM
image = float(image) ;JDM

maxbox = 13                ;Maximum size of convolution box in pixels 

type = size(image)
n_x  = type[1] & n_y = type[2]

;Determine if hardcopy output is desired

radius = 0.637 * FWHM > 2.001     ;Radius is 1.5 sigma
radsq = radius^2
nhalf = fix(radius) < (maxbox-1)/2 ;
nbox = 2*nhalf + 1             ;# of pixels inside of convolution box 
middle = nhalf                  ;Index of central pixel

lastro = n_x - nhalf
lastcl = n_y - nhalf
sigsq = ( fwhm/2.35482 )^2
mask = bytarr( nbox, nbox ) ;Mask identifies valid pixels in convolution box 
c = fltarr( nbox, nbox )   ;c will contain Gaussian convolution kernel

dd = indgen(nbox-1) + 0.5 - middle ;Constants need to compute ROUND
dd2 = dd^2
w = 1. - 0.5*(abs(dd)-0.5) / (middle-.5)   
ir = (nhalf-1) > 1

row2 = (findgen(Nbox)-nhalf)^2

for i = 0, nhalf do begin
    temp = row2 + i^2
    c[0,nhalf-i] = temp         
    c[0,nhalf+i] = temp                           
endfor

mask = fix(c LE radsq) ;MASK is complementary to SKIP in Stetson's Fortran
good = where( mask, pixels) ;Value of c are now equal to distance to center

c = c*mask               
c[good] = exp(-0.5*c[good]/sigsq) ;Make c into a Gaussian kernel
sumc = total(c)
sumcsq = total(c^2) - sumc^2/pixels
sumc = sumc/pixels
c[good] = (c[good] - sumc)/sumcsq
c1 = exp(-.5*row2/sigsq)
sumc1 = total(c1)/nbox
sumc1sq = total(c1^2) - sumc1
c1 = (c1-sumc1)/sumc1sq
sumc = total(w)                 ;Needed for centroid computation

message,'Beginning convolution of image', /INF

h = convol(float(image),c)      ;Convolve image with kernel "c"

h[0:nhalf-1,*] = 0 & h[n_x-nhalf:n_x-1,*] = 0
h[*,0:nhalf-1] = 0 & h[*,n_y-nhalf:n_y-1] = 0

message,'Finished convolution of image', /INF

mask[middle,middle] = 0      ;From now on we exclude the central pixel
pixels = pixels -1      ;so the number of valid pixels is reduced by 1
good = where(mask)         ;"good" identifies position of valid pixels
xx= (good mod nbox) - middle	;x and y coordinate of valid pixels 
yy = fix(good/nbox) - middle    ;relative to the center
offset = yy*n_x + xx

SEARCH:                        ;Threshold dependent search begins here

index = where(h GE hmin, nfound) ;Valid image pixels are greater than hmin
if nfound lt nfibers then begin       ;Any maxima found?
    print,'WARNING! Not enough peaks found!'
    print,'Adjust your threshold value (hmin) and try again!'
    goto,FINISH    
endif

for i= 0L, pixels-1 do begin                             
    stars = where (h[index] GE h[index+offset[i]], nfound)
    if nfound lt nfibers then begin ;Do valid local maxima exist?
        print,'WARNING! Not enough peaks found!'
        print,'Adjust your threshold value (hmin) and try again!'
        goto,FINISH    
    endif
    index = index[stars]
 endfor 
 
 ix = index mod n_x              ;X index of local maxima
 iy = index/n_x                  ;Y index of local maxima
 ngood = N_elements(index)       
 message,strtrim(ngood,2)+' local maxima located above threshold',/INF

 nstar = 0L       	;NSTAR counts all stars meeting selection criteria
 badround = 0L & badsharp=0L  &  badcntrd=0L
 
x = fltarr(ngood) & y = x

;  Loop over star positions; compute statistics

for i = 0L,ngood-1 do begin   
    temp = float(image[ix[i]-nhalf:ix[i]+nhalf,iy[i]-nhalf:iy[i]+nhalf])
    d = h[ix[i],iy[i]]          ;"d" is actual pixel intensity        

;  Compute Sharpness statistic

    sharp1 = (temp[middle,middle] - (total(mask*temp))/pixels)/d
    if ( sharp1 LT sharplim[0] ) or ( sharp1 GT sharplim[1] ) then begin
        badsharp = badsharp + 1
        goto, REJECT            ;Does not meet sharpness criteria
    endif

;   Compute Roundness statistic

    dx = total( total(temp,2)*c1)   
    dy = total( total(temp,1)*c1)
    if (dx LE 0) or (dy LE 0) then begin
        badround = badround + 1
        goto, REJECT            ;Cannot compute roundness
    endif

    around = 2*(dx-dy) / ( dx + dy ) ;Roundness statistic
    if ( around LT roundlim[0] ) or ( around GT roundlim[1] ) then begin
        badround = badround + 1
        goto,REJECT             ;Does not meet roundness criteria
    endif
    
; Find X centroid

    derivat = shift(temp,-1,0) - temp
    derivat = total( derivat[0:nbox-2,middle-ir:middle+ir],2)
    sumd = total(w*derivat)
    sumxd = total(w*dd*derivat)
    sumxsq = total(w*dd2) 

    if ( sumxd GE 0. ) then begin
        badcntrd = badcntrd + 1
        goto,REJECT             ;Cannot compute X centroid
    endif
    
    dx =sumxsq*sumd/(sumc*sumxd)
    if abs(dx) GT nhalf then begin
        badcntrd = badcntrd + 1
        goto,REJECT            ;X centroid too far from local X maxima
    endif

    xcen = ix[i]-dx             ;Convert back to big image coordinates

; Find Y centroid                 

    derivat = shift(temp,0,-1) - temp 
    derivat = total( derivat[middle-ir:middle+ir,0:nbox-2], 1 )
    sumd = total( w*derivat )
    sumxd = total( w*dd*derivat )
    sumxsq = total( w*dd2 )
    if (sumxd GE 0) then begin
        badcntrd = badcntrd + 1
        goto, REJECT  
    endif
    
    dy = sumxsq*sumd/(sumc*sumxd)
    if ( abs(dy) GT nhalf ) then begin
        badcntrd = badcntrd + 1
        goto,REJECT 
    endif
    
    ycen = iy[i] - dy

;  This star has met all selection criteria.  Print out and save results

    x[nstar] = xcen & y[nstar] = ycen
    
   nstar = nstar+1

REJECT: 
  
endfor

nstar = nstar-1             ;NSTAR is now the index of last star found

print,' No. of sources rejected by SHARPNESS criteria',badsharp
print,' No. of sources rejected by ROUNDNESS criteria',badround
print,' No. of sources rejected by CENTROID  criteria',badcntrd

if nstar LT 0 then return       ;Any stars found?

x = x[0:nstar]
y = y[0:nstar] ;these are still in index form, although they are not integers (###.### precision)

nx = n_elements(x)
 
if (plotting1 eq 'on') then begin ;a plot of the convolved field
    window,0,retain=2,xsize=(n_x/2.0),ysize=(n_y/2.0)
    device,decomposed=0
    loadct,0,/silent
    tvscale,c
    pause
    plot,x,y,/nodata,xtitle='Pixels',ytitle='Pixels',title='Centers for '+fitsfile
    for j=0,nx-1 do plots,x[j],y[j],psym=2
    loadct,4,/silent
endif

xout = fltarr(nfibers) & yout = xout
cut = pitch * 0.4
cntr = 0
for j=0,nx-1 do begin ;the fibers are now sorted out to reject duplicate versions of the same fiber found earlier.
    xt = x[j]
    yt = y[j]

    ib = where(xt-cut le x and xt+cut gt x and yt-cut le y and yt+cut gt y,cb)
    if (cb eq 1) then begin
        xt = x[ib]
        yt = y[ib]
    endif else begin
        xt = median(x[ib],/even)
        yt = median(y[ib],/even)
    endelse
    test = where(xout gt xt-cut and xout lt xt+cut and yout gt yt-cut and yout lt yt+cut,count)
    if (count eq 0) then begin
        xout[cntr] = xt
        yout[cntr] = yt
        cntr = cntr + 1
    endif
endfor

if (plotting1 eq 'on') then begin
    for j=0,nfibers-1 do plots,xout[j],yout[j],psym=sym(1),color=150
    print,'CHECK YOUR FIBERS! DID THEY ALL GET FOUND?'
    print,'Next ENTER deletes the plot...'
    pause
    wdelete
endif

fibercenter = fltarr(2,nfibers)

for j=0,nfibers-1 do begin
    chunk = image[xout[j]-hp:xout[j]+hp,yout[j]-hp:yout[j]+hp]
    if (j eq 0) then nn0 = n_elements(chunk[*,0])
    if (j eq 0) then si = nn0-1
    if (j eq 0) then ssi = nn0 * nn0 - 400.0
; The plane-fitting routine is included here. This shouldn't need
; running if the data is well flat-fielded.
    cnr1 = [0,0,median(chunk[0:2,0:2])]
    cnr2 = [si,0,median(chunk[si-2:si,0:2])]
    cnr3 = [0,si,median(chunk[0:2,si-2:si])]
    cnr4 = [si,si,median(chunk[si-2:si,si-2:si])]
    
    cnrval = [cnr1[2],cnr2[2],cnr3[2],cnr4[2]]
    bg = min(cnrval)

;pln = plane(n_x,n_y,z,cnr1[0],cnr1[1],cnr1[2],cnr2[0],cnr2[1],cnr2[2],cnr3[0],cnr3[1],cnr3[2],cnr4[0],cnr4[1],cnr4[2])
;stop
    chunk = chunk - float(bg)

    schunk = chunk[bsort(chunk)]
    cut1 = schunk[ssi]
    ihigh = where(chunk ge cut1)
    chunk[ihigh] = cut1

;    low = median(chunk[0:2,0:2])
;    chunk = chunk - low
;    high = median(chunk[hp-1:hp+1,hp-1:hp+1])
    high = cut1

    if (j eq 0 and plotting2 eq 'on') then begin
        window,0,retain=2,xsize=nn0*10,ysize=nn0*10
        loadct,37,/silent
        window,2,retain=2,xsize=800,ysize=600
        loadct,0,/silent
    endif
    if (plotting2 eq 'on') then begin
        wset,2
        loadct,0,/silent
        surface,chunk,charsize=2.5,title='Fiber # '+strn(j+1)
        wset,0
        loadct,37,/silent
    endif

    center = fltarr(2,11)
    for k=0,10 do begin
        cutoff = (high * (0.8 + (0.01 * float(k))))
        index1 = where(chunk gt cutoff)
        ellipse = Fit_Ellipse(index1, XSize=nn0, YSize=nn0, CENTER=c1)
        center[0,k] = c1[0] & center[1,k] = c1[1]

        if (plotting2 eq 'on') then begin
            loadct,37,/silent
            TVImage, BytScl(chunk, Top=!D.Table_Size-3)
            plots,ellipse*10.0,/Device,Color=FSC_Color('white'),thick=2
            if (slow eq 'on') then wait,0.1
        endif

    endfor
    centerM = median(center,dim=2,/even)
; the [xy]out gives coordinates in the native frame. centerM is in the
; trimmed frame.
    fibercenter[0,j] = xout[j] - hp + centerM[0]
    fibercenter[1,j] = yout[j] - hp + centerM[1]
endfor

if (plotting1 eq 'on') then begin
    window,0,retain=2,xsize=512,ysize=512
    device,decomposed=0
    plot,xout,yout,/nodata
    plots,xout+1.0,yout+1.0,psym=sym(6),symsize=1.5
    loadct,4,/silent
    plots,fibercenter,psym=sym(1),color=150,symsize=1.3
    print,'The next ENTER deletes the plot(s)...'
    pause
    while !d.window ne -1 do wdelete, !d.window
endif

;the region file is created...
free_lun,5
openw,5,'fiber_center.reg'
    printf,5,'# Region file format: DS9 version 4.0'
    printf,5,'# Filename: /home/danu/murphy/fibers/potsdam/position/somename.fits'
    printf,5,'global color=green font="helvetica 10 normal" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source'
    printf,5,'image'
for j=0,nfibers-1 do begin
    line = 'circle('+strn(fibercenter[0,j])+','+strn(fibercenter[1,j])+',11)'
    printf,5,line
endfor
free_lun,5

; Now, the EE calculation is made...
EEarr = fltarr(nfibers)
Xarray = fltarr(pitch,pitch)
Yarray = Xarray

;the radius array for this frame is created
x2 = indgen(hp)+1
x1 = reverse(indgen(hp+1))
xline = [x1,x2]
y2 = indgen(hp)+1
y1 = reverse(indgen(hp+1))
yline = transpose([y1,y2])

FOR i=0,pitch-1 DO Xarray(*,i) = xline
FOR i=0,pitch-1 DO Yarray(i,*) = yline
RADarr = sqrt(Xarray^2 + Yarray^2)

for j=0,nfibers-1 do begin
    x = round(fibercenter[0,j])
    y = round(fibercenter[1,j])
    chunk = image[x-hp:x+hp,y-hp:y+hp]

    ifiber = where(RADarr le fiberR) ;the fiber is identified
    iap = where(RADarr gt fiberR and RADarr le apR) ;the background annulus is identified
    
;First attempt at aperture photometry. The average value of the pixel value
;    n1 = n_elements(ifiber)
;    n2 = n_elements(iap)

    sky = chunk[iap]
    stop

endfor

FINISH:

stop
END
