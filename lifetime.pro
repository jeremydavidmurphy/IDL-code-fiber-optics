; The routine used to reduce the fiber lifetime test data

PRO lifetime,list

; Note: Column IDL index #1525 is bad on the FLI camera

readcol,list,silent=1,f='a',files
istar = where(files eq '*',count)

if (count gt 0) then begin
    baseline = files[0:istar-1]
    files = files[istar+1:*]
endif else begin
    print,'Your list is wrong!'
    print,'Select a set of baseline frames and re-run...'
    stop
endelse

n0 = n_elements(files)
n1 = n_elements(baseline)

OUTarr = strarr(4,n0) ;total1, total2, SD, SDoff
;total1 = total for the frame (pre-division)
;total2 = total for the frame (post-division)
;SD = standard deviation of divided frame
;SDoff = SD[j] / SD[j+1]

basearr = fltarr(1508,1508,n1)

for j=0,n1-1 do begin
    print,'Reading in baseline frame #'+strn(j+1)
    frame = readfits(baseline[j],/silent)
    overscan = frame[0:15,*]
    frame = frame[17:1524,5:1512] - median(overscan)
    basearr[*,*,j] = frame
endfor

baseframe = median(basearr,dim=3)
writefits,'baseline.fits',baseframe

;window,0,retain=2,xsize=812,ysize=512
;device,decomposed=0

free_lun,5
openw,5,'lifetime_'+list+'.txt'

cntr = 1.0
for j = 0, n0-1 do begin
    print,'Reading in data frame #'+strn(j+1)
    frame = readfits(files[j],header,/silent)
    overscan = frame[0:16,*]
    frame = frame[17:1524,5:1512] - median(overscan)
;    print,size(frame)
;    exptime = sxpar(header,'EXPTIME',count=cnt1)
;    temp1   = sxpar(header,'CCDTEMP',count=cnt2)
;    xpos    = sxpar(header,'IFUXPOS',count=cnt3)
;    ypos    = sxpar(header,'IFUYPOS',count=cnt4)
;    zpos    = sxpar(header,'IFUZPOS',count=cnt5)
;    time    = sxpar(header,'IFUTIME',count=cnt6)
;    number  = sxpar(header,'IFUITR',count=cnt7)

    divide = frame / baseframe
    if (j eq 0) then divideold = divide
;    std = divide[bsort(divide)]
;    if (j eq 0) then begin
;        loadct,0
;        plot,std,yrange=[-0.7,1.5],/ystyle
;        loadct,33
;    endif else oplot,std,color=j*2
;    TVImage, BytScl(divide, Top=!D.Table_Size-3)
;    wait,0.3
    t1 = total(frame,/nan)
    OUTarr[0,j] = t1
    t1 = strn(t1)
    t2 = total(divide,/nan)
    OUTarr[1,j] = t2
    t2 = strn(t2)
    std = stddev(divide,/double,/nan)
    OUTarr[2,j] = std
    if (j gt 0) then OUTarr[3,j] = OUTarr[2,j-1] / std
    sstd = strn(std)

    printf,5,[files[j],' ',t1,' ',t2,' ',sstd]
    print,[files[j],' ',t1,' ',t2,' ',sstd]

    if (OUTarr[3,j] gt 1.5) or (OUTarr[3,j] lt 0.6) then begin
        print,'BAD FRAME!!!'
        p1 = strn(j+1)
        outname = 'divide_F'+p1+'.fits'
        writefits,outname,divide
        p1 = strn(j)
        outname = 'divide_F'+p1+'.fits'
        writefits,outname,divideold
        if (cntr eq 1.0) then begin
            badframes = files[j]
        endif else begin
            cntr = cntr + 1.0
            badframes = [badframes,files[j]]
        endelse
    endif
    divideold = divide
endfor

free_lun,5

print,'The frames to check are...'
print,badframes

stop
END
