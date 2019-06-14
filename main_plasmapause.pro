pro main_plasmapause
;;; plot plasmapause-location-related panels.
thm_init
computer = 'I:'
@folders

;;;; select how many plasma pause crossings to plot
;n_cross_plot = 2 ;; original
n_cross_plot = 1 ;; now

t_onset = '2013 9 30 1 20 16'
probe = 'a'

if n_cross_plot eq 1 then begin
	trange = ['2013 9 29 18 5', '2013 9 30 2 55']
endif else begin
	trange = ['2013 9 29 17 5', '2013 9 30 13 59']
endelse

color_cross = 40 ;; color for plasmapause crossing

;;; load locations of THA
load_bin_data, trange = trange, probes = probe, datatype = 'pos', datafolder = pos_folder
get_data, 'th'+probe+'_state_pos', t, pos
store_data, 'th'+probe+'_state_pos_RE', data={x:t, y:pos/RE}

split_vec, 'th'+probe+'_state_pos_RE'
options, 'th'+probe+'_state_pos_RE_x', ytitle = 'X [R!dE!n]', ysubtitle = ''
options, 'th'+probe+'_state_pos_RE_y', ytitle = 'Y [R!dE!n]', ysubtitle = ''
options, 'th'+probe+'_state_pos_RE_z', ytitle = 'Z [R!dE!n]', ysubtitle = ''
qtts_labl = ['th'+probe+'_state_pos_RE_x', 'th'+probe+'_state_pos_RE_y', 'th'+probe+'_state_pos_RE_z']


;;;;;;;;;;;;;;;;;;;;;; Make data plot ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;;;; load potential-inferred density of THA
;thm_scpot2dens_opt_n, sc=probe, datatype_esa='peer', trange = trange, level = 1 ;; now: this only works for sc = 'a'
;get_data, 'th'+probe+'_peer_density_npot', t, n_pot
;;i_bad = where(n_pot gt 100000, n_bad)
;i_bad = where(n_pot gt 100000 or (t gt time_double('13 9 29 23 58') and t lt time_double('13 9 30 0 5')), n_bad) ;; add new bad data in new data product.
;if n_bad gt 0 then begin
;	n_pot[i_bad] = !values.f_nan
;	store_data, 'th'+probe+'_peer_density_npot', data = {x:t, y:n_pot}
;endif
;options, 'th'+probe+'_peer_density_npot', ytitle = 'n!s!de!r!uscpot!n', ysubtitle = '!c[cm!u-3!n]', labels = '', title = 'THEMIS-P5'
;
;if n_cross_plot eq 1 then begin
;	xsize =  4.63
;	options, 'th'+probe+'_peer_density_npot', ytickname = ['10!u-2', '10!u-1', '1', '10', '10!u2', '10!u3', '10!u4']
;endif else begin
;	xsize =  8.63 ;; minimum number to get more tickmarks.
;end
;ysize = 3
;
;popen, pic_folder+'/tha_ppause'
;print_options,xsize=xsize, ysize=ysize
;tplot, 'th'+probe+'_peer_density_npot', trange = trange, var_label = reverse(qtts_labl)
;timebar, t_onset, line = 1
;if n_cross_plot eq 1 then begin
;	add_pointer, 0.47, 0.72, 'plasmapause', dir = 'south', length = 0.05, charsize = 0.9, color = color_cross
;	xyouts, 0.7, 0.82, '(a)'
;endif else begin
;	add_pointer, 0.3, 0.64, 'plasmapause', dir = 'south', length = 0.05, charsize = 0.9, color = color_cross
;	add_pointer, 0.785, 0.67, 'plasmapause', dir = 'south', length = 0.05, charsize = 0.9, color = color_cross
;	xyouts, 0.15, 0.82, '(a)'
;endelse
;pclose
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;; Make cartoon plot ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;cross_range1 = ['2013 9 29 21', '2013 9 29 23 10']
cross_range1 = ['2013-09-29/21:00:00', '2013-09-30/01:20:00']
cross1 = ['2013-09-29/21:40:00', '2013-09-29/22:27:00']
cross_range2 = ['2013-09-30/11:00:00', '2013-09-30/12:10:00']
cross2 = ['2013-09-30/11:24:00', '2013-09-30/11:43:00']

xrange = [0., -11.]
yrange = [9., -1.5]
kp = double(1.+2./3.) ;; looked up from kp list
xsize = 5.3
ysize = 4

;;;;;;;;; trace location to equator, out put is 'th?_efoot'
;; T89
ttrace2equator, 'th'+probe+'_state_pos_RE', newname='th'+probe+'_efoot', external_model='t89', par=kp
time_clip, 'th'+probe+'_efoot', cross_range1[0], cross_range1[1], newname = 'th'+probe+'_efoot_range1'
time_clip, 'th'+probe+'_efoot', cross_range2[0], cross_range2[1], newname = 'th'+probe+'_efoot_range2'
time_clip, 'th'+probe+'_efoot', cross1[0], cross1[1], newname = 'th'+probe+'_efoot_cross1'
time_clip, 'th'+probe+'_efoot', cross2[0], cross2[1], newname = 'th'+probe+'_efoot_cross2'
get_data, 'th'+probe+'_efoot_range1', t_range1, pos_range1
get_data, 'th'+probe+'_efoot_range2', t_range2, pos_range2
get_data, 'th'+probe+'_efoot_cross1', t_cross1, pos_cross1
get_data, 'th'+probe+'_efoot_cross2', t_cross2, pos_cross2
if n_cross_plot gt 1 then begin
	locations = transpose([pos_range1[0,*], pos_range1[-1,*], pos_range2[0,*], pos_range2[-1,*]])
endif else begin
	locations = transpose([pos_range1[0,*], pos_range1[-1,*]])
endelse

popen, pic_folder+'/tha_ppause_pos'
print_options,xsize=xsize, ysize=ysize

;; first draw the background
probe_locations, 0.5*total(time_double(trange)), probes = replicate(probe, n_elements(locations[0,*])), locations = locations, planes_plot = 'XY', title = 'Plasmapause Encounters', xrange = xrange, yrange = yrange, /earth, /geo, size_sc = 0.5, reference_axis = 'X'
;; sc trajectory
oplot, pos_range1[*,0], pos_range1[*,1], color = 1
if n_cross_plot gt 1 then begin
	oplot, pos_range2[*,0], pos_range2[*,1], color = 1
endif
;; crossing range
thick_cross = 10
oplot, pos_cross1[*,0], pos_cross1[*,1], color = color_cross, thick = thick_cross
if n_cross_plot gt 1 then begin
	oplot, pos_cross2[*,0], pos_cross2[*,1], color = color_cross, thick = thick_cross
endif

;; write the times
charsize_time = 0.7
xyouts, pos_range1[0,0], pos_range1[0,1]-0.3, strmid(cross_range1[0], 6, 10), charsize = charsize_time, color = 1, align = 0.8
xyouts, pos_range1[-1,0], pos_range1[-1,1]-0.3, strmid(cross_range1[1], 6, 10), charsize = charsize_time, color = 1, align = 0.2
if n_cross_plot gt 1 then begin
	xyouts, pos_range2[0,0], pos_range2[0,1]+0.5, strmid(cross_range2[0], 6, 10), charsize = charsize_time, color = 1, align = 0.3
	xyouts, pos_range2[-1,0], pos_range2[-1,1]-0.3, strmid(cross_range2[1], 6, 10), charsize = charsize_time, color = 1, align = 0.7
endif

;; write abc
xyouts, !x.crange[0]+0.1*(!x.crange[1]-!x.crange[0]), !y.crange[0]+0.1*(!y.crange[1]-!y.crange[0]), '(b)'

pclose

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

stop
end
