;----------------------------------------------------------------------------------
pro jvr_omni_flux, tele_flux, pitch_angles, j_omni, FILLVAL = FillVal, $
   TRAPEZOID = trapezoid

   ;+
   ; Calculates omnidirectional *average* of differential number flux.
   ; Input pitch angles are 'folded' or 'flipped' (180-0 transformed to 0-90 deg), 
   ; then sorted prior to calculating the integral.  
   ; 
   ; Missing values handled following the Moments ATBD:  if one value (per energy)
   ; is missing, (1) replace with nearest neighbor if at lowest or highest folded
   ; pitch angle, otherwise (2) interpolate in pitch angle.
   ; 
   ; In the future could consider allowing two or more fill values.
   ; 
   ; :inputs:
   ;     tele_flux: differential fluxes by telescope and energy, n x 9 x 5 array
   ;     pitch_angles: pitch angles by telescope, 9 x n array
   ;
   ; :keywords:
   ;     FILLVAL: Fill value.  -99999.0 assumed if not set
   ;
   ; :outputs:
   ;     j_omni: omnidirectionally-averaged differential number flux by energy channel, n x 5 array
   ;
   ; :history:
   ;     jvr_omni_flux created October 5, 2013 from omni_flux
   ;     omni_flux created October 19, 2012 from omni_flux_p
   ;     omni_flux_p originally created October 26, 2010 with an update March 29, 2011
   ;     October 5, 2013: modified to estimate flux if one telescope is fill per energy channel
   ;     July 17, 2015: added TRAPEZOID keyword and an alternative (classic trapezoidal) 
   ;                    integration method for comparison and testing.  Using sin^2(alpha) functions
   ;                    at different pitch angle sampling (10, 1 and 0.1 deg), the two methods are
   ;                    equivalent at 1 and 0.1 deg sampling and both slightly underestimate the omnidirectional
   ;                    average at 10 deg sampling (baseline: 1.8%, trapezoid: 2.2% for n = 6, less for
   ;                    smaller n).
   ;
   ; :examples:
   ;-

   compile_opt idl2, strictarrsubs

   if ~keyword_set(FillVal) then FillVal = -99999. 
   if ~keyword_set(trapezoid) then trapezoid = 0
   
   tele_flux = double(tele_flux)
   pitch_angles = double(pitch_angles)
   
   SizeFluxes = size(tele_flux)
   n_pts = SizeFluxes[1]
   n_pitch = SizeFluxes[2]
   n_chans = SizeFluxes[3]

   SizePAs = size(pitch_angles)
   if SizePAs[2] ne n_pts or SizePAs[1] ne n_pitch then begin
      print, 'Mismatch between size of pitch angle and flux arrays'
      stop
   endif
   n_pts = SizePAs[2]
   j_omni = fltarr(n_pts, n_chans)

; Replace fill values with NaN
   WhFluxFill = where(tele_flux eq FillVal, CtFlux)
   if CtFlux gt 0 then tele_flux[WhFluxFill] = !values.f_nan
   WhPAFill = where(pitch_angles eq FillVal, CtPAs)
   if CtPAs gt 0 then pitch_angles[WhPAFill] = !values.f_nan
         
; Fold pitch angles
   pitch_angles_fold = pitch_angles
   WhFold = where(pitch_angles gt 90.0, Count)
   if Count gt 0 then begin
      pitch_angles_fold[WhFold] = 180.0d0- pitch_angles[WhFold] 
   endif
   
   for i = 0, n_pts-1 do begin
; Sort pitch angles
      pasort = sort(pitch_angles_fold[*, i])
      pa_fold_sort = pitch_angles_fold[pasort, i]
      pitch_angle_low = 0.5d0*(pa_fold_sort[0:n_pitch-2] + pa_fold_sort[1:n_pitch-1])
      pitch_angle_high = [pitch_angle_low, 90.0d0]*!dtor
      pitch_angle_low = [0.0d0, pitch_angle_low]*!dtor
; Sort fluxes and replace fill values
      j_sort = fltarr(n_pitch, n_chans)
      for k = 0, n_chans-1 do begin
         j_sort[*, k] = tele_flux[i, pasort, k]
         WhFil = where(~finite(j_sort[*,k]), CtFil)
         if CtFil eq 1 then begin ; Replace only if one fill value.
            if WhFil eq 0 then j_sort[WhFil,k] = j_sort[WhFil+1,k] $
            else if WhFil eq n_pitch-1 then j_sort[WhFil,k] = j_sort[WhFil-1,k] $ 
            else j_sort[WhFil,k] = j_sort[WhFil-1,k] + (j_sort[WhFil+1,k]-j_sort[WhFil-1,k]) / $
               (pa_fold_sort[WhFil+1]-pa_fold_sort[WhFil-1]) * (pa_fold_sort[WhFil]-pa_fold_sort[WhFil-1])
         endif
      endfor
; for trapezoidal calculation: duplicate 0 and 90 deg points
      pa_fold_complete = [0.0d0, pa_fold_sort, 90.0d0]*!dtor
      j_trapez = [j_sort[0, *], j_sort, j_sort[n_pitch-1, *]]
; Sum flux over pitch angle
      flux = fltarr(n_chans)*0.0
      
      if (trapezoid) then begin
         for pa_ind = 0, n_pitch do begin
            pa_term = (-cos(pa_fold_complete[pa_ind+1]) + cos(pa_fold_complete[pa_ind]))
            flux_term = 0.5d0*(j_trapez[pa_ind, *] + j_trapez[pa_ind+1, *])
            flux = flux + pa_term*flux_term
         endfor
      endif else begin
         for pa_ind = 0, n_pitch-1 do begin
            pa_term = (-cos(pitch_angle_high[pa_ind]) + cos(pitch_angle_low[pa_ind]))
            bin_flux = j_sort[pa_ind, *] * pa_term ; a vector
            flux = flux + bin_flux
         endfor         
      endelse


; omni number flux -- The average over pitch angle results in a coefficient of 1/2.
; Since the integral is from 0 to 90 deg, it is necessary to multiply by 2, 
; which results in a unity coefficient
      j_omni[i,*] = flux
      
   endfor
   
; Replace NaNs with fill value
   WhNan = where(~finite(j_omni), CtNan)
   if CtNan gt 0 then j_omni[WhNan] = FillVal

end