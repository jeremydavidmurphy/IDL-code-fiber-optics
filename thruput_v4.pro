PRO thruput_v4, fibnum, basenum 

;v3 to v4: the code now allows for any number of wavelengths to to
;read in
wave = [3370,3650,4000,4200,4500,5500,6000]

;2 LISTS ARE REQUIRED IN THE CALLING DIRECTORY (blist and flist)-
;THE CODE EXPECTS FILTER 1 TO BE FIRST, AND THEN REVERSES THE ARRAY SO
;THAT THE TIME OF EXPOSURE IS RUNNING FORWARD (this assumes you took
;the red filters first)

;THE LISTS (DATA AND BASELINE) MUST BE IN THE ORDER THE DATA WAS
;TAKEN. 

nwave = n_elements(wave)

n = fibnum*nwave
flist = strarr(n)
n = basenum*nwave
blist = strarr(n)

base = dblarr(nwave,2,basenum)
basename = strarr(nwave,basenum)
data = dblarr(nwave,3,fibnum)
dataname = strarr(nwave,fibnum)
slope = dblarr(nwave,2,basenum)

openr,5,'blist'
readf,5,blist
openr,6,'flist'
readf,6,flist
free_lun, 5,6


count=-1
; The baseline array is read in
FOR k=0,basenum-1 DO BEGIN
    FOR j=0,nwave-1 DO BEGIN
        count=count+1
        im = readfits(blist[count],h)
        h = h[14];14 for subtracted data, 10 for unsubtracted...
        h = strsplit(h,'T:',/extract)
        hour = fix(h[2])
        min = fix(h[3])/60.
        sec = fix(h[4])/3600.
        base[j,0,k] = hour+min+sec
        base[j,1,k] = total(im)
        basename[j,k] = blist[count]
    ENDFOR
ENDFOR

count=-1

; Now the data array is read in similar fashion
FOR k=0,fibnum-1 DO BEGIN
    FOR j=0,nwave-1 DO BEGIN
        count=count+1
        im = readfits(flist[count],h)
        h = h[14]
        h = strsplit(h,'T:',/extract)
        hour = fix(h[2])
        min = fix(h[3])/60.
        sec = fix(h[4])/3600.
        data[j,0,k] = hour+min+sec
        data[j,1,k] = total(im)
        dataname[j,k] = flist[count]
    ENDFOR
ENDFOR

;print,'The initial baseline times are:'
;print,base[*,0,*] & pause
;THE 24 HOUR CLOCK ROLLOVER IS CORRECTED FOR IN THE BASE AND DATA FRAMES
FOR j=0,nwave-1 DO BEGIN
    FOR k=0,basenum-2 DO BEGIN
        IF (base[j,0,k+1] LT base[j,0,k]) THEN base[j,0,k+1] = base[j,0,k+1]+24.00
    ENDFOR
ENDFOR
print,'The adjusted base times are:'
print,base[*,0,*] & pause

;print,'The initial data times are:'
;print,data[*,0,*] & pause

FOR j=0,nwave-1 DO BEGIN
    FOR k=0,fibnum-2 DO BEGIN
        IF (data[j,0,k+1] LT data[j,0,k]) THEN data[j,0,k+1] = data[j,0,k+1]+24.00
    ENDFOR
ENDFOR
print,'The adjusted data times are:'
print,data[*,0,*] & pause

;A quick plot is generated to check the flatness of the baseline data

window,0,retain=2
device,decomposed=0
plot,wave,base[*,1,0],title='Baseline Stability Check',/ynozero,/nodata
loadct,4
FOR j=0,basenum-1 DO oplot,wave,base[*,1,j],psym=4,color=60+(j*10)
pause

;Now the slope for the linear interpolation is determined
count=0
FOR k=0,fibnum-1 DO BEGIN
    IF (data[0,0,k] GT base[0,0,count]) AND (data[0,0,k] LT base[0,0,count+1]) THEN BEGIN
        FOR j=0,nwave-1 DO BEGIN
            run = base[j,0,count+1] - base[j,0,count]
            rise = base[j,1,count+1] - base[j,1,count]
            slope[j,0,count] = double(rise/run)
            slope[j,1,count] = slope[j,0,count]*(data[j,0,k]-base[j,0,count]) + base[j,1,count]
            data[j,2,k] = data[j,1,k]/slope[j,1,count]
        ENDFOR
    ENDIF ELSE BEGIN
        count = count+1
        FOR j=0,nwave-1 DO BEGIN
            run = base[j,0,count+1] - base[j,0,count]
            rise = base[j,1,count+1] - base[j,1,count]
            slope[j,0,count] = double(rise/run)
            slope[j,1,count] = slope[j,0,count]*(data[j,0,k]-base[j,0,count]) + base[j,1,count]
            data[j,2,k] = data[j,1,k]/slope[j,1,count]
        ENDFOR
    ENDELSE
ENDFOR



window,0,retain=2
device,decomposed=0
loadct,0
plot,wave,reverse(data[*,2,0]),xrange=[3000,6300],xstyle=1,yrange=[0.6,1.1],xtitle='Wavelength',ytitle='Absolute Thruput'
loadct,4
FOR k=0,fibnum-1 DO BEGIN
    oplot,wave,reverse(data[*,2,k]),color=60+(k*5)
ENDFOR

STOP
END
