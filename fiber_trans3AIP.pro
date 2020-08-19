PRO fiber_trans3, fitslist, APPTYPE=apptype, DARKSWITCH=darkswitch

; This routine is used for the relative transmission and astrometry of VIRUS fibers. It reads
; in a list of files. Each line in the list is the fits image of the fiber
; head, a corresponding dark frame (this is an option) and the list of
; dead fibers (see DEADFIBERS for more information) If there are no
; dead fibers, this list must still exist, but be filled with a -1.

; EX:
; frame1top.fits  dark1a.fits  deadfibers1t.txt
; frame1bot.fits  dark1b.fits  deadfibers1b.txt
; frame2top.fits  dark2a.fits  deadfibers2t.txt
; frame2bot.fits  dark2b.fits  deadfibers2b.txt
  
; This code runs either ell_appphotf or circ_appphotf, depending on
; the switch apptype.

; Note the order of the frames. The first frame is the top-illuminated fibers,
; the second is the bottom, and so on. THIS IS ESSENTIAL TO THE REDUCTIONS TO
; KEEP THE FIBERS IN THEIR PROPER ORDER.

; This code has the guts of the find.pro routine. I just cut and
; pasted the entire code into the body of this code, then modified it
; to hardcode in some of the free parameters. Once it has found all the peaks,
; many multiple times, there's a step to reject duplicates.

; FITSLIST: The name of the list containing all the frames you want to reduce.

; The flat-fielding for this system was built up from 49 frames where a single
; fiber in a VIRUS fiber bundle was iterated over the FOV of the
; reimager. This was done to fully capture the behavior of the F/# that the
; fiber is feeding the system at. These frames, the routines used to reduce
; them, and a description of the process, can be found in
; /work1/jdm/fibers/re-imager/flatfield on Dalby. See the README file.

; APPTYPE: Enter either 'circle' or 'ellipse'. If this keyword isn't used,
; the default is to apply circular aperture photometry. This default
; is probably the best with the new re-imaging bench (the old one
; being the one built of the 2 Nikon lenses and showed significantly
; larger distortion).

; DARKSWITCH: turn this to 'off' and you can have lists without the dark
; frames. The default is to have the dark frames in the fitslist and
; subtracted.

; OVERVIEW:
; This routine is used to determine the relative transmission of VIRUS
; fibers. The re-imaging test bench is used to back illuminate the fibers and
; re-image the head of the fiber bundle.

; DEADFIBERS: This routine used to have the dead fiber numbers as an
; internal parameter. THESE ARE THE ACTUAL DEAD FIBER NUMBERS, NOT THE
; IDL INDEX! But this requires changing the code each time
; for each fiber bundle. It's been changed to read in a text
; file that notes which fibers are dead. Sorry this is rather clunky,
; but necessary. IF NO DEAD FIBERS EXIST, THIS LIST MUST STILL EXIST,
; BUT HAVE A SINGLE LINE OF -1. YOU CAN HAPPILY HAVE ONE LIST, CALLED
; PERHAPS NODEADFIBERS.TXT THAT'S GOT THE -1 AND PLACE THAT IN
; EVERY LIST FOR EACH BUNDLE THAT DOESN'T HAVE A DEAD FIBER.

; OUTPUTS: This code outputs 6 files for each single frame sent
; in. They are described here. The original file name (ex:
; IFU34top.fits) is used to create the names of the output files. The
; names are split on the first '.' so it's advisable to not use a
; period anywhere in the original naming convention. I will refer to
; the "IFU34top" as "name" below
; 1) name_SF.fits: This is just a subtracted and flattened fits
; frame. You can use this to visually inspect the quality of the
; background subtraction (to check for under/over subtraction) and
; whether the flattening was done properly.
; 2) name.reg: This is just the regions file that shows where each of
; the fibers found are located. To see this frame, open name_SF.fits
; in a ds9 window, go into the 'regions' menu and then 'load
; regions'. This is just a visual inspection and will tell you quickly
; if there is any issue with the fiber finding routine. Every fiber
; should have a nice green circle drawn around it.
; 3) name_focus.txt: This file is used only when you are focusing the
; reimaging bench. It's a text file that contains an estimate of the
; radius and standard deviation of the 21 estiamtes of radius made
; during the running of fiber_trans2. (Turn on plotting4 to see a
; visual of the process of mapping out the best focus.) To use these
; frames, create a list of their names and send them into the IDL
; routine bestfocus.pro. See that code for more notes.
; 4) name.txt: This is the primary output of fiber_trans2.pro. It is a
; text file that contains the positions (X and Y) and relative
; transmission for all the fibers in the frame. There are other values
; saved in this text file, but probably won't be used much. They are
; more diagnostics which you can safely ignore.
; 5) name_fluxmap.ps: This is just a postscript file showing a visual
; of the relative flux in each fiber. It's just a hardcopy of the
; values in name.txt
; 6) name_BGsub.ps: This is a hardcopy of the level of background
; subtraction (i.e. the # of counts subtracted from each frame during
; the circular aperture photometry step) subtracted from each
; fiber. This is a file to watch closely as it will tell you about
; whether your background is being subtracted in a reasonable
; fashion. 

; MODIFIED ON JULY 17TH, 2011: There were two major changes made from this
; version from the original (fiber_trans.pro). This version now expects a list
; of 2 to 3 columns. You can skip the background subtraction step (NOT
; ADVISED!) and just have 2 columns. The order needs to be topframe,
; then bottomframe, otherwise the fiber counting and resulting
; astrometry will be off. Both circular and elliptical aperture
; photometry is an option.

; MODIFIED ON AUG 2, 2011: A routine has been added so that this code
; can be used to focus the reimaging bench. A text file, named
; frame_focus.txt, is now created for every fits frame. This file
; quantifies the diameter of the spot at various EE for each fiber,
; then the standard deviation (SD) of that at the various EE
; points. These text files are then reduced by another IDL routine
; called bestfocus.pro. 

; MODIFIED ON OCT 30, 2011: I have cleaned up a few things, and have
; included the flat-fielding. This was reduced with data taken by
; T. Jahn and reduced by Jeremy Murphy in
; /Users/jeremymurphy/fibers/potsdam/reimager/flat_field
; on his computer. A note about the flat-fielding. There's a
; limited range over which we had flat-field data. Outside of this
; range I have set the flat-field values to 10,000. I have done this
; so that the resulting files will have very low values over those
; regions where the flat-fielding is not real. Also note that there
; are 3 different flat fields. They are different only in their range,
; with slight shifts to the left and right. This stems from the fact
; that in the initial data, there were three fibers illuminated in
; each of the frames Thomas took. Therefore, I could generate 3
; different flat field frames. So, if your data is shifted either to
; the left or right of the region of flat field, try another one of
; the flats (flat1, flat2 or flat3) and maybe you'll be in
; luck.

;*********************************************************************
plotting1 = 'on' ;this turns on some general plots showing which fibers have been found and final EE values. This should generally be left on as it will indicate whether all the fibers have been found.
plotting2 = 'off' ;this one turns on the ellipse-fitting routine, which is VERY slow if you turn on 'slow' below. You should run this the first few times you run this code as it will give you a better understanding of the fiber-finding routine.
slow = 'on' ;turn on to slow the ellipse-fitting routine down for better visual inspection.
plotting3 = 'off' ;this turns on or off the aperture photometry plotting routine.
plotting4 = 'off' ;turn this on to see where the diameter estimates for the focus routine are falling. THIS PLOTTING ROUTINE CAN'T BE TURNED ON WITH PLOTTING2, OR PLOTTING3.
plottingscale = 'on' ;turn this on if you want to see a plot showing the fiber centroid spacing (i.e. the radial distance between every fiber).
pauseplot = 'off' ;turn this to 'on' and the plots will be paused rather than just a delay
makeps = 'yes' ;set to 'yes' to generate postscript output of relative fiber flux map
focus = 'yes' ; set to 'yes' to return an estimate of the focus of the frame. You will then use the routine bestfocus.txt to make sense of this output.
focustep = 0.025 ;this value sets the step size for the range of focus. If the horozontal lines exploring the spot diameter in the focus routine are too close, raise this number. If they are too far apart, lower this number.
fwhm = 52.0 ; (22.0 for the photometry bench) the FWHM in pixels...? If all your fibers are being found, this parameter has done it's job. It feeds into the find.pro routine.
hmin = 4000 ;The detection threshold, in pixel values. This is the counts above the background to consider a peak a peak. This is only used for the initial determination of peaks
sharplim = [0.001,5.0]
roundlim = [-5.0,5.0] ;To get the elongated fibers in the corners, the roundness criteria is loosened.
pitch = 109 ;MUST BE AN ODD INTEGER. An estimate of the pitch between fibers. This number is used for the fiber centriod finder and for sorting the fibers into order from 1 to nfibers (35 was the pitch for the Nikon lenses set up)
apertureR = 57.0 ;MUST BE AN ODD INTEGER. This aperture is the outer range when fitting circular apertures.
R1c = 12.0 ;used in the circular aperture photometry routine. The amount added to the major axis to include in the total counts. Also the inner radius for the background estimate
R2c = 18.0 ;used in the circular aperture photometry routine. The outer radius for the estimate of the background.
R1e = 5.0 ;used in the elliptical aperture photometry routine. The amount added to the major axis to include in the total counts. Also the inner radius for the background estimate
R2e = 12.0 ;used in the elliptical aperture photometry routine. The outer radius for the estimate of the background.
;the values above are tuned for the VIRUS fibers
;the values below are tuned for the VP2.7 fiber bundle
;pitch = 79 ;MUST BE AN ODD INTEGER. An estimate of the pitch between fibers. This number is used for the fiber centriod finder and for sorting the fibers into order from 1 to nfibers (35 was the pitch for the Nikon lenses set up)
;apertureR = 41.0 ;MUST BE AN ODD INTEGER. This aperture is the outer range when fitting circular apertures.
;R1c = 12.0 ;used in the circular aperture photometry routine. The amount added to the major axis to include in the total counts. Also the inner radius for the background estimate
;R2c = 18.0 ;used in the circular aperture photometry routine. The outer radius for the estimate of the background.
;R1e = 5.0 ;used in the elliptical aperture photometry routine. The amount added to the major axis to include in the total counts. Also the inner radius for the background estimate
;R2e = 12.0 ;used in the elliptical aperture photometry routine. The outer radius for the estimate of the background.

lop = 0.80 ;The percentage of the peak fiber value to set as the cutoff limit for flat-topping the fibers. As the fibers are not top-hat functions, this sets all the pixel values above this value to the same value, effectively flatten the fibers. This helps in determining the centers of the fibers as the ellipse-fitting routine used (fit_ellipse.pro) is a mass-weighting routine which would be skewed by lop-sided fiber intensities.
cutmagcent = 0.60 ; cutmagcent + (cutstepcent * 10) < 1.00 (This is the lower starting point when running the center-finding routine.)
cutstepcent = 0.02 ;the steps to take in the iteration of the center
cutmagell = 0.25 ;the starting point for the ellipse-fitting routine.
cutstepell = 0.005 ;the step size for the ellipse-fitting routine
wscale = 4.0       ;just scales plotting window sizes. Turn up for larger plots and down for smaller ones.
nf = 224 ;The number of fibers that should exist in the frame. IF THERE ARE DEAD FIBERS, DO NOT REDUCE THIS NUMBER! Just put their # into the 'deadfibers' array above.
flatfield = 'no' ;set to 'yes' to use the flat-fielding, set to 'no' to turn it off (I am dubious of using this...)
;*********************************************************************

if (n_elements(APPTYPE) eq 0) then APPTYPE = 'circle'
if (n_elements(DARKSWITCH) eq 0) then DARKSWITCH = 'on'
hp = floor(pitch/2.0)

if (DARKSWITCH eq 'on') then $
  readcol,fitslist,f='a,a,a',silent=1,fitsfile,darkfile,deadfile
if (DARKSWITCH eq 'off') then $
  readcol,fitslist,f='a,a',silent=1,fitsfile,deadfile

nframes = n_elements(fitsfile)

if (flatfield eq 'yes') then begin
;   flat = readfits('/Users/jeremymurphy/fibers/potsdam/reimager/flat_field/flatMASTER.fits')
   flat = readfits('/work1/jdm/flat_field/flatMASTER.fits')
   flat = float(flat)
endif

counter = [0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1]
;Use this counter when you are running both the top and bottom parts
;of the IFU. (This just handles fiber-counting = the code then knows
;when fiber #225 is really #225 and not #1.)
;counter = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
;Use this counter when you want to run each reduction
;independently. This simply keeps the fiber counting correct. 

for all = 0,nframes-1 do begin    ;a loop through each individual fits frame
   print,''
   print,'Working on frame '+fitsfile[all]+'...'
   print,''
   
   readcol,deadfile[all],silent=1,f='i',deadfibers
   
   if (deadfibers[0] ne -1) then begin
       idead = deadfibers - 1.0
       ndead = n_elements(deadfibers)
       nfibers = nf - ndead
   endif else begin
       nfibers = nf
       ndead = 0
   endelse

   EEarr = fltarr(nf) ;the EE for each fiber. dead fibers are left as 0.0
   BG = fltarr(nf) ;the background subtracted during aperture photometry
   fiberarray = fltarr(6,nf) ;X_center, Y_center, PA, Major Axis, Minor Axis, R1
; fiberarray holds the positions and components of the ellipse found by the
; fit_ellipse routine. These components are used later when conducting the
; ellipse-aperture photometry routine. (R1 is the radial distance from fiber #1)
   
   image = readfits(fitsfile[all],/silent,header)
   image = float(image)
   if (DARKSWITCH eq 'on') then begin
       dark = readfits(darkfile[all],/silent,/noscale)
       dark = float(dark)
       image = image - dark
       if (max(image[100:200,100:200]) gt 30000.0) then image = image - 32768.0 ;a fudge to capture some funny-busniess which depends on the system you're running this on
   endif

   if (flatfield eq 'yes') then image = image / flat
   tempname = strsplit(fitsfile[all],'.',/extract)
   nameout1 = tempname[0]+'_SF.fits'
   nameout2 = tempname[0]+'.txt'
   nameout3 = tempname[0]+'_fluxmap.ps'
   nameout4 = tempname[0]+'_focus.txt'
   nameout5 = tempname[0]+'_BGsub.ps'
   nameout6 = tempname[0]+'IFUscale.txt'
;   sxaddpar,header,'COMMENT','  Flattened via the file built_flatBOX.fits'

   writefits,nameout1,image,header ;the dark subtracted and flattened fits image is written out.

;*********************************************************************
;********* The beginning of the stripped find.pro code ***************
;*********************************************************************

   maxbox = 31                  ;Maximum size of convolution box in pixels 
   type = size(image)
   n_x  = type[1] & n_y = type[2]
   radius = 0.637 * FWHM > 2.001 ;Radius is 1.5 sigma
   radsq = radius^20
   nhalf = fix(radius) < (maxbox-1)/2 ;
   nbox = 2*nhalf + 1                 ;# of pixels inside of convolution box 
   middle = nhalf                     ;Index of central pixel
   lastro = n_x - nhalf
   lastcl = n_y - nhalf
   sigsq = ( fwhm/2.35482 )^2
   mask = bytarr( nbox, nbox )  ;Mask identifies valid pixels in convolution box 
   c = fltarr( nbox, nbox )     ;c will contain Gaussian convolution kernel
   
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
   
   mask = fix(c LE radsq)       ;MASK is complementary to SKIP in Stetson's Fortran
   good = where( mask, pixels)  ;Value of c are now equal to distance to center
   
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
   sumc = total(w)              ;Needed for centroid computation
   
   h = convol(float(image),c)   ;Convolve image with kernel "c"
   
   h[0:nhalf-1,*] = 0 & h[n_x-nhalf:n_x-1,*] = 0
   h[*,0:nhalf-1] = 0 & h[*,n_y-nhalf:n_y-1] = 0
   
   mask[middle,middle] = 0      ;From now on we exclude the central pixel
   pixels = pixels -1           ;so the number of valid pixels is reduced by 1
   good = where(mask)           ;"good" identifies position of valid pixels
   xx= (good mod nbox) - middle	;x and y coordinate of valid pixels 
   yy = fix(good/nbox) - middle ;relative to the center
   offset = yy*n_x + xx
   
   SEARCH:                      ;Threshold dependent search begins here
   
   index = where(h GE hmin, nfound) ;Valid image pixels are greater than hmin
   if nfound lt nfibers then begin  ;Any maxima found?
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
   
   ix = index mod n_x           ;X index of local maxima
   iy = index/n_x               ;Y index of local maxima
   ngood = N_elements(index)       
   message,strtrim(ngood,2)+' local maxima located above threshold',/INF
   
   nstar = 0L                   ;NSTAR counts all stars meeting selection criteria
   badround = 0L & badsharp=0L  &  badcntrd=0L
   
   x = fltarr(ngood) & y = x
   
;  Loop over star positions; compute statistics
   for i = 0L,ngood-1 do begin   
      temp = float(image[ix[i]-nhalf:ix[i]+nhalf,iy[i]-nhalf:iy[i]+nhalf])
      d = h[ix[i],iy[i]]        ;"d" is actual pixel intensity        
      
;  Compute Sharpness statistic
      sharp1 = (temp[middle,middle] - (total(mask*temp))/pixels)/d
      if ( sharp1 LT sharplim[0] ) or ( sharp1 GT sharplim[1] ) then begin
         badsharp = badsharp + 1
         goto, REJECT           ;Does not meet sharpness criteria
      endif
      
;   Compute Roundness statistic
      dx = total( total(temp,2)*c1)   
      dy = total( total(temp,1)*c1)
      if (dx LE 0) or (dy LE 0) then begin
         badround = badround + 1
         goto, REJECT           ;Cannot compute roundness
      endif
      
      around = 2*(dx-dy) / ( dx + dy ) ;Roundness statistic
      if ( around LT roundlim[0] ) or ( around GT roundlim[1] ) then begin
         badround = badround + 1
         goto,REJECT            ;Does not meet roundness criteria
      endif

; Find X centroid
      derivat = shift(temp,-1,0) - temp
      derivat = total( derivat[0:nbox-2,middle-ir:middle+ir],2)
      sumd = total(w*derivat)
      sumxd = total(w*dd*derivat)
      sumxsq = total(w*dd2) 
      
      if ( sumxd GE 0. ) then begin
         badcntrd = badcntrd + 1
         goto,REJECT            ;Cannot compute X centroid
      endif
      
      dx =sumxsq*sumd/(sumc*sumxd)
      if abs(dx) GT nhalf then begin
         badcntrd = badcntrd + 1
         goto,REJECT            ;X centroid too far from local X maxima
      endif
      
      xcen = ix[i]-dx           ;Convert back to big image coordinates
      
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
   
   nstar = nstar-1              ;NSTAR is now the index of last star found
   
   print,' No. of sources rejected by SHARPNESS criteria',badsharp
   print,' No. of sources rejected by ROUNDNESS criteria',badround
   print,' No. of sources rejected by CENTROID  criteria',badcntrd
   
   if nstar LT 0 then return    ;Any stars found?
   
   x = x[0:nstar]
   y = y[0:nstar]               ;these are still in index form, although they are not integers (###.### precision)
   
   nx = n_elements(x)
   
;*********************************************************************
;*************** The end of the stripped find.pro code ***************
;*********************************************************************

   if (plotting1 eq 'on') then begin ;a plot of ALL FIBERS found by the find.pro routine
      window,0,retain=2,xsize=(n_x/wscale),ysize=(n_y/wscale)
      device,decomposed=0
      loadct,0,/silent
      plot,x,y,/nodata,xtitle='Pixels',ytitle='Pixels',title='Centers for '+fitsfile[all]
      for j=0,nx-1 do plots,x[j],y[j],psym=2
      loadct,4,/silent
      wait,1.0
   endif
   
;Now, the fibers found more than once are thrown out.
   xout = fltarr(nfibers) & yout = xout
   cut = pitch * 0.4       ;used to pick off fibers found more than once by find.pro
   
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
      print,'Did all your fibers get found???'
      print,'(if not, change some parameters and try again!)'
      print,'Next ENTER deletes the plot...'
      if (pauseplot eq 'on') then pause else wait,1.0
      wdelete
   endif
  
;the fibers are sorted from 1 to nfibers (the dead fibers are not in this array yet...)
   iy = reverse(bsort(yout))
   ytemp = yout[iy]
   xtemp = xout[iy]
   youtS = 0.0
   xoutS = 0.0
   for i=0,nfibers-2 do begin
      iy1 = i
      if (i eq 0) then iy1l = iy1
      iy2 = i+1
      diff = abs(ytemp[iy1] - ytemp[iy2])
      if (diff gt pitch*0.5) or (i eq nfibers-2) then begin
         iy1 = iy1l
         iy1l = iy2
         if (i lt nfibers-2) then iy2 = iy2-1
         x1 = xtemp[iy1:iy2] ;the x values for a row of fibers is identified
         y1 = ytemp[iy1:iy2]
         ix = bsort(x1) ;they are sorted, from low to high, which places them in their proper fiber order
         youtS = [youtS,y1[ix]]
         xoutS = [xoutS,x1[ix]]
      endif
   endfor
   youtS = youtS[1:*] ;the leading zero's are dropped
   xoutS = xoutS[1:*] ;these ordered arrays now become the primary array so that all that follows is in order.
   
   if (plotting1 eq 'on') then begin
      window,0,retain=2,xsize=512,ysize=512
      plot,xoutS,youtS,/nodata,title='Order of fibers',/isotropic
      for j=0,nfibers-1 do begin
         plots,xoutS[j],youtS[j],psym=sym(1)
         wait,0.01
      endfor
  endif
   
   if (ndead gt 0) then begin ;the dead fibers are re-inserted back into the xoutS and youtS arrays
      for j=0,ndead-1 do begin
         idead1 = idead[j]
         px1 = xoutS[0:idead1-1]
         px2 = xoutS[idead1:*]
         xoutS = [px1,-1.0,px2]
         py1 = youtS[0:idead1-1]
         py2 = youtS[idead1:*]
         youtS = [py1,-1.0,py2]
      endfor
   endif

;the radius array for this frame is created
   x2 = indgen(hp)+1
   x1 = reverse(indgen(hp+1))
   xline = [x1,x2]
   y2 = indgen(hp)+1
   y1 = reverse(indgen(hp+1))
   yline = transpose([y1,y2])
   
   Xarray = fltarr(pitch,pitch)
   Yarray = Xarray
   
   FOR i=0,pitch-1 DO Xarray(*,i) = xline
   FOR i=0,pitch-1 DO Yarray(i,*) = yline
   RADarr = sqrt(Xarray^2 + Yarray^2)
   iouter = where(RADarr gt apertureR-2.0)
   
;the first part of the .reg file is generated
   name1 = strsplit(fitsfile[all],'.',/extract)
   free_lun,55
   openw,55,name1[0]+'.reg'
   printf,55,'# Region file format: DS9 version 4.0'
   printf,55,'# Filename: /home/danu/murphy/fibers/potsdam/position/'+fitsfile[all]
   printf,55,'global color=green font="helvetica 10 normal" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source'
   printf,55,'image'

   deadcntr = 0
   pcntr = 0.0
   for j=0,nf-1 do begin ;a loop through each of the fibers
      l = j + deadcntr
       print,'Working on fiber #'+strn(j+1)
;      print,l
      if ((ndead gt 0) and (ndead ne deadcntr)) then if (j eq idead[deadcntr]) then begin
          print,'Fiber #'+strn(idead[deadcntr]+1.0)+' is totally dead!!!'
         deadcntr = deadcntr + 1
         goto, deadjump
      endif
      x1 = uint(round(xoutS[j])-hp)
      x2 = uint(round(xoutS[j])+hp)
      y1 = uint(round(youtS[j])-hp)
      y2 = uint(round(youtS[j])+hp)
      chunk = image[x1:x2,y1:y2]
      chunkorig = chunk
      chunk[iouter] = 0.0 ;the values outside of the outer aperture are set to zero for the fiber location. This helps the fit_ellipse routine not get fooled by light from a neighboring fiber.
      if (j eq 0) then begin
         nn0 = n_elements(chunk[*,0])
         nn1 = n_elements(chunk[0,*])
         ssi = floor((nn0*nn1*lop)) ;the threshold for the top-hat of the fiber profile is set.
      endif
      schunk = chunk[bsort(chunk)]
      cut1 = schunk[ssi]
      ihigh = where(chunk ge cut1)
      chunk[ihigh] = cut1 ;the top is loped off the fiber to improve the center-finding routines.
      
      if (j eq 0 and plotting2 eq 'on') then begin
         window,0,retain=2,xsize=nn0*wscale,ysize=nn0*wscale
         loadct,37,/silent
         window,1,retain=2,xsize=nn0*wscale,ysize=nn0*wscale
         loadct,37,/silent
         window,6,retain=2,xsize=800,ysize=600,xpos=550
         loadct,0,/silent
      endif
      if (plotting2 eq 'on') then begin
         wset,6
         loadct,0,/silent
         surface,chunk,charsize=2.5,title='Fiber # '+strn(l+1)
         wset,0
         loadct,37,/silent
         wset,1
         loadct,37,/silent
      endif
      
      center = fltarr(2,11)
      pa = fltarr(11)
      axes = fltarr(2,11)
      for k=0,10 do begin       ;an iteration on the fit_ellipse routine.
         cutoff1 = (cut1 * (cutmagcent + (cutstepcent * float(k))))
         index1 = where(chunk gt cutoff1)
         cutoff2 = (cut1 * (cutmagell + (cutstepell * float(k))))
         index2 = where(chunk gt cutoff2)
         circle  = Fit_Ellipse(index1, XSize=nn0, YSize=nn0, CENTER=c1)
         ellipse = Fit_Ellipse(index2, XSize=nn0, YSize=nn0, orient=p1,axes=a1)
         center[0,k] = c1[0] & center[1,k] = c1[1]
         pa[k] = p1
         axes[0,k] = a1[0] & axes[1,k] = a1[1]
         if (plotting2 eq 'on') then begin
            wset,0
            TVImage, BytScl(chunk, Top=!D.Table_Size-3)
            plots,circle*wscale,/Device,Color=FSC_Color('white'),thick=2
            if (k eq 0) then begin 
               wset,1
               TVImage, BytScl(chunk, Top=!D.Table_Size-3)
               plots,ellipse*wscale,/Device,Color=FSC_Color('white'),thick=2
            endif
            if (slow eq 'on') then wait,0.2
         endif
      endfor

      centerM = median(center,dim=2,/even)
      paM = median(pa,/even)
      axesM = median(axes,dim=2,/even)
; the [xy]out gives coordinates in the native frame. centerM is in the
; trimmed frame.
      fiberarray[0,j] = float(x1) + float(centerM[0]) + 1.0
      fiberarray[1,j] = float(y1) + float(centerM[1]) + 1.0
      fiberarray[2,j] = paM
      fiberarray[3,j] = axesM[0]
      fiberarray[4,j] = axesM[1]
      xx = fiberarray[0,j]
      yy = fiberarray[1,j]

      if (counter[all] eq 0) then begin
         Xfirst = fiberarray[0,0]
         Yfirst = fiberarray[1,0]
         fiberarray[5,j] = sqrt((xx-Xfirst)^2 + (yy-Yfirst)^2)
      endif
      if (counter[all] eq 1) then begin
         fiberarray[5,j] = sqrt((xx-Xfirst)^2 + (yy-Yfirst)^2)
      endif

;the ds9 coordinates are written
;      line = 'circle('+strn(fiberarray[0,l])+','+strn(fiberarray[1,l])+',10)' ;use for the nikon lens re-imaging set up
;      line = 'circle('+strn(fiberarray[0,j])+','+strn(fiberarray[1,j])+',32)' ;use for the large lens imaging set up
      line = 'circle('+strn(fiberarray[0,j])+','+strn(fiberarray[1,j])+',32) # text={'+strn(j+1)+'}' ;this line adds the fiber numbers to each circle.
      printf,55,line

;the estimate of focus is made, if called for
      if (focus eq 'yes') then begin
          schunk = chunkorig[bsort(chunkorig)]
          cut2 = median(schunk[n_elements(schunk)-1500,*],/even)
          diameter = fltarr(31)
          diaARR = fltarr(60)
          for k=1,31 do begin
            cutoff3 = cut2 - (cut2 * (focustep * k))
            index3 = where(chunkorig gt cutoff3)
            ellipse = Fit_Ellipse(index3, XSize=nn0, YSize=nn0,CENTER=c1)
            if (plotting4 eq 'on' and j eq pcntr) then begin ;plots the cross-sections and where the diameter line falls
               if (k eq 1) then begin
                  xslice = median(chunkorig[*,c1[1]-3:c1[1]+3],dim=2)
                  yslice = median(chunkorig[c1[0]-3:c1[0]+3,*],dim=1)
                  window,5,retain=2,xsize=500,ysize=500,xpos=430,ypos=50
                  plot,xslice,title='X Profile'
                  oplot,[0,n_elements(xslice)],[cutoff3,cutoff3],thick=1
                  window,15,retain=2,xsize=500,ysize=500,xpos=950,ypos=50
                  plot,yslice,title='Y Profile'
                  oplot,[0,n_elements(yslice)],[cutoff3,cutoff3],thick=1
               endif else begin
                  wset,5
                  oplot,[0,n_elements(xslice)],[cutoff3,cutoff3],thick=1
                  wset,15
                  oplot,[0,n_elements(yslice)],[cutoff3,cutoff3],thick=1
                  wait,0.01
               endelse
            endif

            for ll=0,59 do begin ;makes an estimate of the diameter of the spot, based on the ellipse fitting.
               x1 = ellipse[0,ll]
               y1 = ellipse[1,ll]
               x2 = ellipse[0,ll+59]
               y2 = ellipse[1,ll+59]
               dx = double(abs(x1-x2))
               dy = double(abs(y1-y2))
               diaARR[ll] = sqrt(dx^2 + dy^2)
            endfor
            diameter[k-1] = median(diaARR,/even)
         endfor

        MEDdia = median(diameter)
        SDdia = stddev(diameter)
        if (j eq 0) then begin
            free_lun,56
            openw,56,nameout4
        endif 
        printf,56,MEDdia,' ',SDdia
        if (j eq pcntr) then pcntr = pcntr + 10.0
     endif

;now the EE calculation is made, using either ell_appphotF.pro or circ_appphotF.pro
      xc = fiberarray[0,j]
      yc = fiberarray[1,j]
      pa = fiberarray[2,j]
      ma = fiberarray[3,j]
      mi = fiberarray[4,j]
      ma = ma/2.0
      mi = mi/2.0
      if (APPTYPE eq 'ellipse') then out = ell_appphotf(image,xc,yc,pa+90.0,ma,mi,R1e,R2e,TYPE='median',WIN=wscale,PLOT=plotting3)
      if (APPTYPE eq 'circle') then out = circ_appphotf(image,xc,yc,ma,R1c,R2c,TYPE='mean',WIN=wscale,PLOT=plotting3)
      EEarr[j] = out[0]
      BG[j] = out[1]

      deadjump: ;jumps over dead fibers

   endfor ;the end of the loop over each fiber found
   free_lun,55 ;closing the ds9 region file
   free_lun,56 ;closing the focus estimate file

   if (plotting1 eq 'on') then begin
      window,0,retain=2,xsize=512,ysize=512
      device,decomposed=0
      loadct,0,/silent
      plot,xoutS,youtS,/nodata
      plots,xoutS+1.0,youtS+1.0,psym=sym(6),symsize=1.5
      loadct,4,/silent
      plots,fiberarray[0,*],fiberarray[1,*],psym=sym(1),color=150,symsize=1.3
      print,'The next ENTER deletes the plot(s)...'
      if (pauseplot eq 'on') then pause else wait,1.0
      while !d.window ne -1 do wdelete, !d.window
   endif
   
   EEarrN = EEarr / max(EEarr)
   igood = where(EEarrN ne 0.0)
   ibad = where(EEarrN eq 0.0,countbad)

;added on April 8, 2013: this part of the routine calculates the
;average radial spacing between all the fibers in the image.

   radius = 0.0
   angle = 0.0
   for j=0,nf-1 do begin
      x1 = fiberarray[0,j]
      y1 = fiberarray[1,j]
      dx = fiberarray[0,*] - x1
      dy = fiberarray[1,*] - y1
      for k=0,n_elements(dx)-1 do begin
         if (dx[k] lt 90 and dx[k] gt 0 and abs(dy[k]) lt 30) then begin
            arc = atan(dy[k],dx[k]) ;the angle, in radians, between two fibers
            angle = [angle,arc]
         endif
      endfor
      rad = sqrt(dx^2 + dy^2)
      ikeep = where(rad lt 90.0)
      radius = [radius,rad[ikeep]]
   endfor

   ikeep = where(radius ne 0.0)
   radius = radius[ikeep]
   angle = angle[1:*] * (180.0/!pi)
   rotate = median(angle)
   
   m = moment(radius)
   sd = sqrt(m[1])
   m = [m,sd]
   scale = 340.0 / m[0]
   xreal = scale * fiberarray[0,*] ;the x-position, converted to microns
   yreal = scale * fiberarray[1,*] ;the y-position, converted to microns
   rreal = scale * fiberarray[5,*] ;the radial spacing between fibers, in microns
   centerx = ((xreal[123]-xreal[122])/2.0) + xreal[122]
   centery = ((yreal[123]-yreal[122])/2.0) + yreal[122]
   offsetx = centerx - xreal
   offsety = centery - yreal
   radoff = sqrt(offsetx^2 + offsety^2)

  if (plottingscale eq 'on') then begin
      window,0,retain=2
      loadct,0,/silent
      plot,radius,/ynozero,ytitle='Fiber Spacing (pixels)',psym=1
      if (pauseplot eq 'on') then pause else wait,3.0
      plot,angle,ytitle='Angle of rotation (degrees)',psym=1
      oplot,[0,300],[rotate,rotate]
      xyouts,10,rotate+(0.1*rotate),'The median of the angles...',charsize=2
      if (pauseplot eq 'on') then pause else wait,3.0
   endif
  
; the data is written out to a text file
   free_lun,5
   openw,5,nameout6
   printf,5,'The calculated rotation, in degrees, is '+strn(rotate)+'.'
   printf,5,'The calculated offsets between fibers (via the MOMENT IDL routine):'
   printf,5,'Mean (pixel) Variance   Skewness      Kurtosis      S.D.'
   printf,5,strn(m)
   form0 = '(a5,3x,a10,2x,a10,3x,a8)'
   printf,5,['Fib #','X_off (um)','Y_off (um)','Rad (um)'],format=form0
   form1 = '(i3,5x,f10.4,2x,f10.4,2x,f10.4)'
   for k=0,nf-1 do begin
      if (counter[all] eq 0) then piece = [k+1,offsetx[k],offsety[k],radoff[k]]
      if (counter[all] eq 1) then piece = [k+nf+1,offsetx[k],offsety[k],radoff[k]]
      printf,5,piece,format=form1
   endfor
   free_lun,5

   openw,5,nameout2
   form0 = '(a5,4x,a5,6x,a5,6x,a5,6x,a4,6x,a2,8x,a7,3x,a7)'
   printf,5,['Fib #','X_pos','Y_pos','Rad','Flux','PA','R_major','R_minor'],format=form0
   form1 = '(i3,5x,f9.4,2x,f9.4,2x,f9.4,2x,f6.4,2x,f10.5,2x,f8.4,2x,f8.4)'
   for k=0,nf-1 do begin
      if (counter[all] eq 0) then piece = [k+1,fiberarray[0,k],fiberarray[1,k],fiberarray[5,k],EEarrN[k],fiberarray[2,k],fiberarray[3,k],fiberarray[4,k]]
      if (counter[all] eq 1) then piece = [k+nf+1,fiberarray[0,k],fiberarray[1,k],fiberarray[5,k],EEarrN[k],fiberarray[2,k],fiberarray[3,k],fiberarray[4,k]]
      printf,5,piece,format=form1
   endfor
   free_lun,5
   
   fcolors = EEarrN * 255.0
   fcolors = round(fcolors)
   fcolorsBG = round((BG/max(BG)) * 255.0)
   
   if (plotting1 eq 'on') then begin
      window,0,retain=2,xsize=512,ysize=512
      device,decomposed=0
      loadct,0,/silent
      plot,fiberarray[0,igood],fiberarray[1,igood],/nodata,title='Relative Fiber Flux',$
           xrange=[min(fiberarray[0,igood])-50,max(fiberarray[0,igood])+50],/xstyle,$
           yrange=[min(fiberarray[1,igood])-50,max(fiberarray[1,igood])+50],/ystyle,/isotropic
      loadct,33,/silent
      for k=0,nf-1 do if (fiberarray[0,k] gt 0) then plots,fiberarray[0,k],fiberarray[1,k],psym=2,color=fcolors[k]
      print,'The minimum fiber value is: '+strn(min(EEarrN[igood]))
      print,'The maximum fiber value is: '+strn(max(EEarrN[igood]))
      window,3,retain=2,xsize=1024,ysize=1024
      device,decomposed=0
      loadct,0,/silent
      plot,fiberarray[0,igood],fiberarray[1,igood],psym=sym(1),title='Relative Fiber Flux',$
           xrange=[min(fiberarray[0,igood])-30,max(fiberarray[0,igood])+50],/xstyle,$
           yrange=[min(fiberarray[1,igood])-30,max(fiberarray[1,igood])+50],/ystyle,/nodata,/isotropic
      for k=0,nf-1 do xyouts,fiberarray[0,k],fiberarray[1,k],' '+strn((EEarrN[k]*100.0),length=4),$
                                  charsize=1.5,orientation=30,charthick=1.0
      loadct,33,/silent
      for k=0,nf-1 do plots,fiberarray[0,k],fiberarray[1,k],psym=sym(6),color=fcolors[k]
      
      if (makeps eq 'yes') then begin
         set_plot,'ps'
         device,file=nameout3,/color ;the fiber transmission .ps file is created.
         loadct,0,/silent
         plot,fiberarray[0,igood],fiberarray[1,igood],psym=sym(1),title='Relative Fiber Flux: '+fitsfile[all],$
              xrange=[min(fiberarray[0,igood])-30,max(fiberarray[0,igood])+70],/xstyle,$
              yrange=[min(fiberarray[1,igood])-30,max(fiberarray[1,igood])+70],/ystyle,$
              xthick=2,ythick=2,charthick=2,/isotropic,/nodata
         for k=0,nf-1 do xyouts,fiberarray[0,k],fiberarray[1,k],' '+strn((EEarrN[k]*100.0),length=4),$
                                     charsize=0.6,orientation=30,charthick=1.5
         loadct,33,/silent
         for k=0,nf-1 do plots,fiberarray[0,k],fiberarray[1,k],psym=sym(1),color=fcolors[k],symsize=0.5
         device,/close_file

         device,file=nameout5,/color ;the fiber subtracted background is created.
         loadct,0,/silent
         plot,fiberarray[0,igood],fiberarray[1,igood],psym=sym(1),title='Subtracted Background: '+fitsfile[all],$
              xrange=[min(fiberarray[0,igood])-30,max(fiberarray[0,igood])+70],/xstyle,$
              yrange=[min(fiberarray[1,igood])-30,max(fiberarray[1,igood])+70],/ystyle,$
              xthick=2,ythick=2,charthick=2,/isotropic,/nodata
         for k=0,nf-1 do xyouts,fiberarray[0,k],fiberarray[1,k],' '+strn(round(BG[k])),$
                                     charsize=0.6,orientation=30,charthick=1.5
         loadct,33,/silent
         for k=0,nf-1 do plots,fiberarray[0,k],fiberarray[1,k],psym=sym(1),color=fcolorsBG[k],symsize=0.5
         device,/close_file
         set_plot,'x'
      endif
      
      if (pauseplot eq 'on') then pause else wait,3.0

      loadct,0,/silent
      plot,fiberarray[0,igood],fiberarray[1,igood],psym=sym(1),title='Subtracted Background',$
           xrange=[min(fiberarray[0,igood])-20,max(fiberarray[0,igood])+40],/xstyle,$
           yrange=[min(fiberarray[1,igood])-20,max(fiberarray[1,igood])+40],/ystyle,/nodata,/iso
      for k=0,nf-1 do xyouts,fiberarray[0,k],fiberarray[1,k],' '+strn(round(BG[k])),$
                                  charsize=1.5,orientation=30,charthick=1.0
      loadct,33,/silent
      for k=0,nf-1 do plots,fiberarray[0,k],fiberarray[1,k],psym=sym(6),color=fcolorsBG[k]
      wset,0
      loadct,0,/silent
      plot,fiberarray[0,igood],fiberarray[1,igood],/nodata,title='Subtracted Background',$
           xrange=[min(fiberarray[0,igood])-20,max(fiberarray[0,igood])+20],/xstyle,$
           yrange=[min(fiberarray[1,igood])-20,max(fiberarray[1,igood])+20],/ystyle,/iso
      loadct,33,/silent
      for k=0,nf-1 do plots,fiberarray[0,k],fiberarray[1,k],psym=2,color=fcolorsBG[k]
      
      if (pauseplot eq 'on') then pause else wait,1.0

   endif
FINISH:

   wait,3.0
   while !d.window ne -1 do wdelete, !d.window

endfor ;the end of the main loop that runs over different images   

if (plotting1 eq 'on' or plotting2 eq 'on') then begin
   if (pauseplot eq 'on') then begin
      pause 
      while !d.window ne -1 do wdelete, !d.window
   endif else begin
      wait,1.0
      while !d.window ne -1 do wdelete, !d.window
   endelse
endif

print,''
print,'Fiber Finder has completed successfully!'
print,''

stop
END
