FUNCTION circ_appphotF,array,centX,centY,major,R1f,R2f,TYPE=type,WIN=winscale,$
                       PLOT=plot

; This routine is used by fiber_trans.pro to conduct circular aperture
; photometry. When the fibers are very circular in the image, this routine
; will give very similar results to ell_appphotF.pro. (The ellipses become
; circles). A SIGNIFICANT DIFFERENCE BETWEEN THE TWO IS THAT THE TOTAL NUMBER
; OF PIXELS FOR THIS ROUTINE IS CONSTANT. This may matter in the case where the
; background hasn't been well-subtracted.

;ARRAY: The array
;centX: The x-position of the best guess for the center (in pixel
;coords, NOT IDL index)
;centY: The y-position of the best guess for the center
;R1f: The inner radius, in pixels
;R2f: The outer radius, in pixels
;TYPE: enter 'mean' or 'median'- this determines how the background in the
;                                annulus is estimated. The default is mean.
;winscale: the scaling for the window size.

if (n_elements(plot) eq 0) then plottingF = 'off' else plottingF = plot
if (n_elements(winscale) eq 0) then winscale = 4.0
if (n_elements(type) eq 0) then type = 'mean'

array = float(array)
centX = round(centX - 1.0)
centY = round(centY - 1.0)
majorC = round(major)
R1g = R1f + majorC
R2g = R2f + majorC

chunk = array[centX-R2g:centX+R2g,centY-R2g:centY+R2g]
nx = n_elements(chunk[*,0])
ny = n_elements(chunk[0,*])

line2 = indgen(R2g)+1
line1 = reverse(indgen(R2g+1))
line = [line1,line2]
Xarray = fltarr(nx,ny) & Yarray = Xarray

FOR i=0,nx-1 DO Xarray(*,i) = line
FOR i=0,ny-1 DO Yarray(i,*) = line
RADarr = sqrt(Xarray^2 + Yarray^2)

ispot = where(RADarr le R1g)
iap = where(RADarr gt R1g and RADarr le R2g)

if (plottingF eq 'on') then begin
   tchunk = chunk
   tchunk[iap] = 10000.0
   wset,2
   tvscale,tchunk
endif

if (type eq 'median') then M1 = median(chunk[iap],/even)
if (type eq 'mean') then M1 = mean(chunk[iap])
chunkout = chunk[ispot] - M1
out = [total(chunkout),M1]

return,out ;the total, and the background per pixel are subtracted.

END
