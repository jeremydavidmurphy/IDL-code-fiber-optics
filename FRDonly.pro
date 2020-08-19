PRO FRDonly, list, DARK=dark

; COMPILE: fratioF.pro

; This code is used to feed fratioF.pro directly. The code accepts a
; list of lists, with each list being a list of frames for a single
; measurement of FRD, with one wavelength and one fiber. 

; This code will read in a dark frame and subtract it on the fly. The
; background will be scaled, based on a comparison of the median of
; the background as compared to the median of a 200 x 200 pixel square
; in the upper right corner.

; If the keyword DARK is used, then the background will be subtracted
; by the code. This is also make a scaling of the BG based on the
; corner of the frame...

cstep = 4.00
psize = 0.000012 * 2
winsize = 3.0

readcol,list,f='a',lists
n0= n_elements(lists)

if (n_elements(dark) ne 0) then begin 
    bg = readfits(dark,/silent)
    s1 = median(bg)
    bgsub = 'yes'
endif else bgsub = 'no'

fratio = fltarr(30,2,n0)
for j=0,n0-1 do begin
    list = lists[j]
    readcol,list,f='a',frames
    n1 = n_elements(frames)
    test = readfits(frames[0],/silent)
    s = size(test)
    fiber = fltarr(s[1],s[2],n1)
    for k=0,n1-1 do begin
        fiber[*,*,k] = readfits(frames[k],/silent)
        if (bgsub eq 'yes') then begin
            s2 = median(fiber[s[1]-200:s[1]-1,s[2]-200:s[2]-1,k]) ;change this if the upper right corner has fiber light in it.
            sbg = bg * (s2/s1)
            fiber[*,*,k] = fiber[*,*,k] - sbg
        endif
        if k eq 0 and j eq 0 then writefits,'test1.fits',fiber[*,*,k]
        if k eq 0 and j eq 1 then writefits,'test2.fits',fiber[*,*,k]
        if k eq 0 and j eq 2 then writefits,'test3.fits',fiber[*,*,k]
    endfor
    out = fratioF(fiber,cstep,psize,winsize,ee=0.98)
    fratio[*,0,j] = out[*,0]
    fratio[*,1,j] = out[*,1]
endfor

form1 = '(f6.4,1x,f6.4,1x,f6.4,1x,f6.4,1x,f6.4,1x,f6.4,1x,f6.4,1x,f6.4,1x,f6.4,1x,f6.4,1x,f6.4,1x,f6.4,1x,f6.4,1x,f6.4,1x,f6.4,1x,f6.4,1x,f6.4,1x,f6.4,1x,f6.4,1x,f6.4,1x,f6.4,1x,f6.4,1x,f6.4,1x,f6.4,1x,f6.4,1x,f6.4,1x,f6.4,1x,f6.4,1x,f6.4,1x,f6.4)'
openw,5,'FRDout.txt'
for j=0,n0-1 do begin
    printf,5,lists[j]
    printf,5,fratio[*,0,j],format=form1
    printf,5,fratio[*,1,j],format=form1
endfor
free_lun,5

stop
END
