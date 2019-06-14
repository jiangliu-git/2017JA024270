pro main_plot_solarwind
@folders
;;; plot solarwind data (WIND or ARTEMIS, ARTEMIS has no velocity for event 1784)
;; if want to plot OMNI or ARTEMIS data, use main_plot_events_find. (OMNI has only 1-min resolution)
trange = ['2013 9 30 1 00', '2013 9 30 1 30']

;;;;;;;;;;; load wind with solarwind_load (automatically moved to bowshock nose point) ;;;;;;;;;;;;
;trange_load = time_double(trange)+[-2*24*3600., 3600.]
;solarwind_load, swdata, dstout, trange_load, /wind
;store_data, 'wind_bz', data = {x:swdata[*,0], y:swdata[*,2]}
;store_data, 'wind_pdyn', data = {x:swdata[*,0], y:swdata[*,1]}
;tplot, ['wind_bz', 'wind_pdyn'], trange = trange, title = 'WIND data from solarwind_load'
;timebar, 0, /databar, varname = 'wind_bz'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;; load ARTEMIS data ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; plot parameters
xsize = 6
ysize = 4.8
;; for event 1784, only thc has data
sc = 'c'

;;; labels (abc)
abc = '('+['a', 'b', 'c']+')'
xabc = replicate(0.2, n_elements(abc))
yabc = [0.9, 0.6, 0.35]

trange_load = time_double(trange)+[-2*3600., 2*3600.]
thm_load_state, trange=trange_load, probe = sc, datatype = 'pos', coord = 'gsm', suffix = '_gsm'
thm_load_fit, trange=trange_load, probe = sc, datatype = 'fgs', coord = 'gsm', suffix = '_gsm'
thm_load_fit, trange=trange_load, probe = sc, datatype = 'fgs', coord = 'dsl', suffix = '_dsl'
thm_load_esansst2, trange=trange_load, probe = sc, /fill ;; this is used for 3-sec resolution density

;;;;; this part is to load ~10 min velocity data
;thm_load_esa_pot,probe=sc, trange=trange_load
;combined_i = thm_part_combine(probe=sc, trange=trange_load, $
;	esa_datatype='peif', sst_datatype='psif', $
;	orig_esa=iesa, orig_sst=isst, sst_sun_bins=isst_mask) 
;thm_part_products, dist_array=combined_i, outputs='moments', sc_pot_name = 'th'+sc+'_esa_pot', mag_name = 'th'+sc+'_fgs_dsl'
;thm_cotrans,'th'+sc+'_ptiff_velocity', 'th'+sc+'_ptiff_velocity_gsm', in_coord='dsl', out_coord='gsm' 
;;;; load onboard moments for velocity
thm_load_mom, trange=trange_load, probe = sc, coord = 'gsm', suffix = '_gsm'

;;; compute dynamic pressure
split_vec, 'th'+sc+'_ptim_velocity_gsm'
tinterpol_mxn,'th'+sc+'_ptix_density','th'+sc+'_ptim_velocity_gsm_x',newname='th'+sc+'_ptix_density_int'
get_data, 'th'+sc+'_ptim_velocity_gsm_x', t, vx
get_data, 'th'+sc+'_ptix_density_int', t, ni
store_data, 'th'+sc+'_Pdyn', data={x:t, y:ni2npa*ni*vx^2}

;; get satellite locations
get_data, 'th'+sc+'_state_pos_gsm', t, pos
store_data, 'th'+sc+'_state_pos_RE', data = {x:t, y:pos/6371.}
split_vec, 'th'+sc+'_state_pos_RE'
options, 'th'+sc+'_state_pos_RE_x', ytitle = 'X [R!dE!n]', ysubtitle = ''
options, 'th'+sc+'_state_pos_RE_y', ytitle = 'Y [R!dE!n]', ysubtitle = ''
options, 'th'+sc+'_state_pos_RE_z', ytitle = 'Z [R!dE!n]', ysubtitle = ''

;;; compute subsolar point to get time shift
time_clip, 'th'+sc+'_state_pos_RE_x', trange[0], trange[1]
time_clip, 'th'+sc+'_fgs_gsm', trange[0], trange[1]
time_clip, 'th'+sc+'_ptim_velocity_gsm', trange[0], trange[1]
time_clip, 'th'+sc+'_Pdyn', trange[0], trange[1]
get_data, 'th'+sc+'_state_pos_RE_x', t, xall
x_sc = mean(xall, /nan) 
get_data, 'th'+sc+'_fgs_gsm_tclip', t, ball
Bz = mean(ball[*,2], /nan)
get_data, 'th'+sc+'_ptim_velocity_gsm_tclip', t, vall
Vx = mean(vall[*,0], /nan)
get_data, 'th'+sc+'_Pdyn_tclip', t, Pdynall
Pdyn = mean(Pdynall, /nan)
r0 = (10.22+1.29*tanh(0.184*(Bz+8.14)))*Pdyn^(-1./6.6) ;; subsolar point
time_shift = (x_sc-r0)*6371./(-Vx)

;;;; shift the quantities to be plotted
;tnames_plot = 'th'+sc+'_'+['fgs_gsm', 'ptim_velocity_gsm', 'ptim_density', 'Pdyn']
tnames_plot = 'th'+sc+'_'+['fgs_gsm', 'ptim_velocity_gsm', 'Pdyn']
for i = 0, n_elements(tnames_plot)-1 do begin
	get_data, tnames_plot[i], t, data
	store_data, tnames_plot[i]+'_shift', data = {x:t-time_shift, y:data}
endfor

;; set labels
options, 'th'+sc+'_ptix_density', ylog=0
options, 'th'+sc+'_Pdyn*', ytitle = 'P!ddyn', thick = l_thick, ysubtitle = '!c[nPa]', ylog = 1
options, 'th'+sc+'_fgs_gsm*', ytitle = 'B', thick = l_thick, labels = 'B!d'+['x','y','z'], labflag = 1, colors = [2,4,6], ysubtitle = '!c[nT]'
options, 'th'+sc+'_ptim_velocity_gsm*', ytitle = 'V!di', thick = l_thick, labels = 'V!d'+['x','y','z'], labflag = 1, colors = [2,4,6], ysubtitle = '!c[km/s]'

popen, pic_folder+'/th'+sc
print_options,xsize=xsize, ysize=ysize
;tplot, tnames_plot, trange = trange, title = 'Solar Wind from P2 (time-shifted)', var_label = ['th'+sc+'_state_pos_RE_z', 'th'+sc+'_state_pos_RE_y', 'th'+sc+'_state_pos_RE_x']
tplot, tnames_plot+'_shift', trange = trange, title = 'Solar Wind from P2 (time-shifted)', var_label = ['th'+sc+'_state_pos_RE_z', 'th'+sc+'_state_pos_RE_y', 'th'+sc+'_state_pos_RE_x']
;;; put onset time
timebar, '2013-09-30/01:20:16', line = 1
timebar_mass, 0, /databar, varname = 'th'+sc+['_fgs_gsm_shift', '_ptim_velocity_gsm_shift'], line = 1
xyouts, xabc, yabc, abc, /normal
pclose

;makepng, pic_folder+'/thc_sw'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

stop
end
