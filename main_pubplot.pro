pro main_pubplot
;;;;; plot all data plots (Themis, RBSP, geotail, GBOs). For solar wind use main_plot_solarwind
computer = 'I:'
@folders
thm_init
del_data, '*'

plot_type = 'publication' ;; for publication
;plot_type = 'hidedate' ;; for publication
;plot_type = 'view'
;plot_type = 'long_injection'
;plot_type = 'long_plasmapause'

;;;;;;;; the kp during this event
kp = double(1.+2./3.) ;; looked up from kp list

;;;;;;;; Define time ranges and and onset time
date = '2013 9 30 '
rtrange = date+['0 55', '1 5']
t_onset = date+'1 20 16' ;; change by clicking
if strmatch(plot_type, '*long*') then begin
	;trange_show = date+['1 0 1', '8 39 30'] ;; long range
	trange_show = date+['0 0 1', '8 39 30'] ;; long range
endif else begin
	trange_show = date+['1 0 1', '1 39 30'] ;; regular range
	;trange_show = date+['0 30 1', '1 59 30'] ;; longer range
	;trange_show = date+['1 0 1', '2 29 30'] ;; longer range
endelse

;;;;;;;; values for band pass filter
;;;;; Pi2: values from Saito, 1969, T = 40 - 150 s
mins_filt_lowf = 150./60. ;; minutes for high-pass_filter
mins_filt_highf = 40./60. ;; minutes for low-pass_filter
;;;; Pc5
;mins_filt_lowf = 600./60. ;; minutes for high-pass_filter
;mins_filt_highf = 150./60. ;; minutes for low-pass_filter

;;; correct the altitude-caused waves in RBSP
mins_atitude = 3.5
;;; trend for the electric field to see whether convection happened.
mins_filt_trend = 3.
trange_load = time_double(trange_show)+[-100, 100]
;;; period range for power spectra
powerspec_range = [600., 4.] ;; in seconds

;;;;;;;;; load kyoto AL
load_bin_data, trange = trange_load, datatype = 'kyoto_al', datafolder = kyoto_al_folder
load_bin_data, trange = trange_load, datatype = 'pseudo_al', datafolder = al_folder
options, 'kyoto_al', ytitle = 'Kyoto AL', ysubtitle = '[nT]', thick = l_thick
options, 'thg_idx_al', ytitle = 'THM AL', ysubtitle = '[nT]', thick = l_thick

;; plot settings
if strcmp(plot_type, 'hidedate') then begin
	tplot_options,'vtitle',''
endif



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; THEMIS D, E, A ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;; plot settings

;;; Whehter to plot spectra of probe a to be lines.
if_spec_lines_a = 1

;;; fill or not fill plasma moments
fill_type = 'fill'
;fill_type = 'nofill'

;;; size in the order of a, d, e
size_x = [6,6,6]
size_y = [10,7,7]
x_abc = 0.2

;trange_show  = date+['1 17', '1 30'] ;; change trange_show to zoom in
trange_load = time_double(trange_show)+[-100, 100]

probes  = ['a', 'e', 'd'] ;; e must be before 

for i = 0, n_elements(probes)-1 do begin
	del_data, '*' ;; for faster performence, delete all saved tvnames in each loop.
	sc = probes[i]

	;if ~strcmp(sc, 'a') then continue ;; plot P5 only (for Figure 4)
	;if strcmp(sc, 'a') then continue ;; skip P5
	if ~strcmp(sc, 'd') then continue

	;;; choose specific i_mask and e_mask for SST loading
	case sc of
	'a': begin
		i_mask = [8, 24, 40, 53, 54, 55, 56, 57, 58, 59]
		e_mask = [8, 24, 40, 55, 56, 57, 58, 59]
		end
	'd': begin
		i_mask = [8, 24, 40, 55, 56]
		e_mask = [8, 40, 56]
		end
	'e': begin
		i_mask = [40, 55, 56]
		e_mask = [40, 56]
		end
	else: begin
		i_mask = 0
		e_mask = 0
		end
	endcase

	;;; load data
	;del_data, '*' ;; for faster load

	thm_load_esansst2, trange=trange_load, probe = sc, i_mask = i_mask, e_mask = e_mask, /fill
	;;;;; to keep no fill moments, use the following
	if strcmp(fill_type, 'nofill') then begin
		thm_load_esansst2, trange=trange_load, probe = sc, i_mask = i_mask, e_mask = e_mask
		copy_data, 'th'+sc+'_ptirf_eflux_energy', 'th'+sc+'_ptix_en_eflux'
		copy_data, 'th'+sc+'_pterf_eflux_energy', 'th'+sc+'_ptex_en_eflux'
	endif

	;;; load pitch angle distribution for tha
	if strcmp(sc, 'a') then begin
		datatypes = ['peir', 'peef', 'psif', 'psef'] ;; peer has no angular resolution
		energy_ranges_esa = [[30, 100.], [100., 1000], [1000, 5000], [5000, 10000], [10000, 17000], [17000, 25000]] ;; in eV
		energy_ranges_sst = [[30., 60.], [60, 100], [100, 150.], [150., 1000.]]*1000. ;; in eV
		for i_type = 0, n_elements(datatypes)-1 do begin
			datatype_this = datatypes[i_type]
			thm_load_state, probe=sc, coord='gei', /get_support, trange=trange_load
			thm_load_fit, probe=sc, coord='dsl', trange=trange_load
			thm_part_load, probe=sc, trange=trange_load, datatype=datatype_this
			if strcmp_or(datatype_this, ['peir', 'peef']) then begin
				energy_ranges = energy_ranges_esa
			endif else begin
				energy_ranges = energy_ranges_sst
			endelse

			for i_energy = 0, n_elements(energy_ranges[0,*])-1 do begin
				erange_this = energy_ranges[*,i_energy]
				thm_part_products, probe=sc, datatype=datatypes[i_type], trange=trange_load, outputs='pa', energy = erange_this, suffix = '_'+strcompress(string(round(0.5*total(erange_this))), /remove)
			endfor
		endfor
	endif

	;;; load this last because thm_load_esansst2 will load efs too
	load_efi_data, trange=trange_load, probe=sc, datatype = 'efs', /tclip, e_folder = efs_folder, b_folder = fgs_folder, rtrange = rtrange, angle_good = 15.

	;;;; get the spin axis direction in gsm for diagnose and cartoon use. Spin direction for THA in GSM: [0.99, 0.09, 0.15]
	;get_data, 'th'+sc+'_fgs_gsm', t, nouse
	;store_data, 'th'+sc+'_spinaxis_dsl', data = {x:t, y:rebin(transpose([0.,0.,1.]), n_elements(t), 3)}
	;thm_cotrans, 'th'+sc+'_spinaxis_dsl', 'th'+sc+'_spinaxis_gsm', in_coord = 'dsl', out_coord = 'gsm'
	;tplot, ['th'+sc+'_spinaxis_dsl', 'th'+sc+'_spinaxis_gsm']

	calc,'"th?_state_pos_gsm_RE" = "th?_state_pos_gsm"/6374.4'
	split_vec, 'th'+sc+'_Pth'
	split_vec, 'th'+sc+'_state_pos_gsm_RE'
	;;; store the thermal + magnetic variable
	;; compute Pbx
	tinterpol_mxn, 'th'+sc+'_fgs_gsm', 'th'+sc+'_Pth', newname='th'+sc+'_fgs_gsm_int'
	get_data, 'th'+sc+'_fgs_gsm_int', t, bdata
	Pbx = nTesla2_to_nPa*bdata[*,0]^2
	Pbxy = Pbx+nTesla2_to_nPa*bdata[*,1]^2
	get_data, 'th'+sc+'_Pall', t, pdata
	store_data, 'th'+sc+'_Pbth', data = {x:t, y:pdata[*,0:1]} ;; Pb, Pth
	store_data, 'th'+sc+'_Pbxyth', data = {x:t, y:[[Pbxy], [pdata[*,0:1]]]} ;; Pbxy, Pb, Pth

	;;; compute field-aligned B
	cotrans2fac, 'th'+sc+'_fgs_gsm', b_variable = 'th'+sc+'_fgs_gsm', pos = 'th'+sc+'_state_pos_gsm', smooth_time = 6.*60., newname = 'th'+sc+'_fgs_fac'

	;;; compute spin-plane field-alighed E and B
	cotrans2fac, 'th'+sc+'_efs_dot0_dsl', b_variable = 'th'+sc+'_fgs_dsl', smooth_time = 6.*60., newname = 'th'+sc+'_efs_dot0_sfac', /sc_spin
	cotrans2fac, 'th'+sc+'_fgs_dsl', b_variable = 'th'+sc+'_fgs_dsl', smooth_time = 6.*60., newname = 'th'+sc+'_fgs_sfac', /sc_spin

	;;; compute B strength
	get_data, 'th'+sc+'_fgs_gsm', t, bdata
	store_data, 'th'+sc+'_fgs_b', data={x:t, y:sqrt(total(bdata^2, 2))}

	;;; compute alven and sound speeds
	;; Alfven
	tinterpol_mxn, 'th'+sc+'_ptix_density', 'th'+sc+'_fgs_b', newname='th'+sc+'_ptix_density_int'
	get_data, 'th'+sc+'_fgs_b', t, data_b
	get_data, 'th'+sc+'_ptix_density_int', t, data_n
	v_alfven = 21.8*data_b/sqrt(data_n)
	store_data, 'th'+sc+'_alfven', data = {x:t, y:v_alfven}
	;; sound
	tinterpol_mxn, 'th'+sc+'_Pall', 'th'+sc+'_fgs_b', newname='th'+sc+'_Pall_int'
	get_data, 'th'+sc+'_Pall_int', t, data_pall
	v_sound = sqrt(data_pall[*,1]/(mass_proton*data_n)*1e-15)*1e-3 ;; in km/s
	;; store
	store_data, 'th'+sc+'_speeds', data = {x:t, y:[[v_alfven], [v_sound]]}

	;;; filter
	tbandpass_filter, 'th'+sc+'_fgs_gsm', mins_filt_lowf*60., mins_filt_highf*60.
	tbandpass_filter, 'th'+sc+'_fgs_fac', mins_filt_lowf*60., mins_filt_highf*60.
	tbandpass_filter, 'th'+sc+'_fgs_sfac', mins_filt_lowf*60., mins_filt_highf*60.
	tbandpass_filter, 'th'+sc+'_fgs_b', mins_filt_lowf*60., mins_filt_highf*60.
	tbandpass_filter, 'th'+sc+'_Pth_z', mins_filt_lowf*60., mins_filt_highf*60.
	tbandpass_filter, 'th'+sc+'_Pbth', mins_filt_lowf*60., mins_filt_highf*60.
	tbandpass_filter, 'th'+sc+'_efs_gsm', mins_filt_lowf*60., mins_filt_highf*60.
	tbandpass_filter, 'th'+sc+'_efs_dot0_dsl', mins_filt_lowf*60., mins_filt_highf*60.
	tbandpass_filter, 'th'+sc+'_efs_dot0_sfac', mins_filt_lowf*60., mins_filt_highf*60.
	;; E trend
	tlowpass_filter, 'th'+sc+'_efs_dot0_sfac', mins_filt_trend*60.
	tlowpass_filter, 'th'+sc+'_efs_gsm', mins_filt_trend*60.

	;;; compute poynting vector
	tinterpol_mxn,'th'+sc+'_efs_dot0_sfac_bpfilt','th'+sc+'_fgs_sfac_bpfilt',newname='th'+sc+'_efs_dot0_sfac_bpfilt_int'
	get_data, 'th'+sc+'_efs_dot0_sfac_bpfilt_int', t, e_wave
	get_data, 'th'+sc+'_fgs_sfac_bpfilt', t, b_wave
	poynting = Re^2*1e-6*transpose(crossp_long(transpose(e_wave), transpose(b_wave), n_dir = poynting_dir)) ;; Unit: W/RE^2
	store_data, 'th'+sc+'_poynting_sfac', data={x:t, y:poynting}, dlimits = {data_att:{coord_sys:'sfac'}}
	store_data, 'th'+sc+'_poynting_sfac_dir', data={x:t, y:transpose(poynting_dir)}, dlimits = {data_att:{coord_sys:'sfac'}}

	;;; original quantities labels
	options, 'th'+sc+'_state_pos_gsm_RE_x', ytitle = 'X [R!dE]', ysubtitle = ''
	options, 'th'+sc+'_state_pos_gsm_RE_y', ytitle = 'Y [R!dE]', ysubtitle = ''
	options, 'th'+sc+'_state_pos_gsm_RE_z', ytitle = 'Z [R!dE]', ysubtitle = ''
	options, 'th'+sc+'_fgs_gsm*',colors=[2,4,6],labels=['B!dx','B!dy','B!dz'], ytitle = 'B', ysubtitle = '!c[nT]', labflag = 1, thick = l_thick
	options, 'th'+sc+'_fgs_fac*',colors=[2,4,6],labels=["B!dr'",'B!d'+phi_letter+"'",'B!db'], ytitle = 'B!dFAC', ysubtitle = '!c[nT]', labflag = 1, thick = l_thick
	options, 'th'+sc+'_fgs_b*', ytitle = '|B|', ysubtitle = '!c[nT]', thick = l_thick
	options, 'th'+sc+'_efs_dot0_dsl*',colors=[2,4,6],labels=["E!dx'","E!dy'",'E!dax'], ytitle = 'E!dDSL', ysubtitle = '!c[mV/m]', labflag = 1, thick = l_thick
	options, 'th'+sc+'_efs_dot0_sfac*',colors=[2,4,6],labels=['E!d-ax','E!dsp'+perp_sign,"E!db'"], ytitle = 'E!dsFAC', ysubtitle = '!c[mV/m]', labflag = 1, thick = l_thick
	options, 'th'+sc+'_efs_gsm*',colors=[2,4,6],labels=['E!dx','E!dy','E!dz'], ytitle = 'E', ysubtitle = '!c[mV/m]', labflag = 1, thick = l_thick
	options, 'th'+sc+'_ptix_vperp_gsm*',colors=[2,4,6],labels=['V!dx','V!dy','V!dz'], ytitle = 'V!di'+perp_sign+'!n', ysubtitle = '!c[km/s]', labflag = 1, thick = l_thick
	options, 'th'+sc+'_ptix_velocity_gsm*',colors=[2,4,6],labels=['V!dx','V!dy','V!dz'], ytitle = ' V!di!n', ysubtitle = '!c[km/s]', labflag = 1, thick = l_thick
	options, 'th'+sc+'_poynting_dir_sfac*',colors=[2,4,6],labels=["S!d-ax",'S!dsp'+perp_sign,"S!db'"], ytitle = 'S!dsFAC!n dir', ysubtitle = '', labflag = 1, thick = l_thick
	options, 'th'+sc+'_poynting_sfac*',colors=[2,4,6],labels=["S!d-ax",'S!dsp'+perp_sign,"S!db'"], ytitle = 'S!dsFAC!n', ysubtitle = '!c[W/R!s!u2!r!dE!n]', labflag = 1, thick = l_thick
	options, 'th'+sc+'_Pth_z*', ytitle = 'P!dth!n', ysubtitle = '!c[nPa]', color = 0, labels = '', thick = l_thick
	options, 'th'+sc+'_Pbth*', ytitle = 'P', ysubtitle = '!c[nPa]', colors = [1,4], labels = ['P!dB', 'P!dth'], thick = l_thick, labflag = 1
	options, 'th'+sc+'_Pbxyth*', ytitle = 'P', ysubtitle = '!c[nPa]', colors = [2,1,4], labels = ['P!dBxy', 'P!dB', 'P!dth'], thick = l_thick, labflag = 1
	options, 'th'+sc+'_speeds*', ytitle = 'V', ysubtitle = '!c[km/s]', colors = [2,4], labels = ['V!dAlfven', 'V!dsound'], thick = l_thick, labflag = 1


	;;;;; adjust spectra
	if strcmp(sc, 'a') and if_spec_lines_a then begin
		;;;; remove some bins (also photo electrons) so the energy curves looks better.
		get_data, 'th'+sc+'_ptix_en_eflux', t_ptix, flux_ptix, en_ptix
		get_data, 'th'+sc+'_ptex_en_eflux', t_ptex, flux_ptex, en_ptex
		n_ch_i = n_elements(flux_ptix[0,*])
		n_ch_e = n_elements(flux_ptex[0,*])

		i_bins_keep_i = [25:32]
		i_bins_photo_e = [0:2]
		i_bins_keep_e = [14:42]
		i_bins_bad_e = [28:34]

		inc_i = 2
		inc_e = 2
		i_bins_remove_i = [1,2,linspace(3,i_bins_keep_i[0]-2,inc=inc_i), linspace(i_bins_keep_i[-1]+1,n_ch_i-1,inc=inc_i)]
		i_bins_remove_e = [i_bins_photo_e, linspace(i_bins_photo_e[-1]+2,i_bins_keep_e[0]-2,inc=inc_e), i_bins_bad_e, linspace(i_bins_keep_e[-1]+1,n_ch_e-1,inc=inc_e)]

		flux_ptix[*,i_bins_remove_i] = !values.f_nan
		flux_ptex[*,i_bins_remove_e] = !values.f_nan
		store_data, 'th'+sc+'_ptix_en_eflux', data = {x:t_ptix, y:flux_ptix, v:en_ptix}
		store_data, 'th'+sc+'_ptex_en_eflux', data = {x:t_ptex, y:flux_ptex, v:en_ptex}


		;;;;;;;; set labels and limits
		;;;;; generate labels
		;;;;;;;;; for diagnose
		;lbl_ptix = strcompress(string(indgen(n_ch_i)))
		;lbl_ptex = strcompress(string(indgen(n_ch_e)))

		;;;;;;;; for real plot
		;;;; re-store the data: remove the NaN rows.
		i_finite_i = where(finite(mean(flux_ptix, dim = 1, /nan)), n_finite_i)
		en_ptix = en_ptix[i_finite_i]
		store_data, 'th'+sc+'_ptix_en_eflux', data = {x:t_ptix, y:flux_ptix[*,i_finite_i], v:en_ptix}
		i_finite_e = where(finite(mean(flux_ptex, dim = 1, /nan)), n_finite_e)
		en_ptex = en_ptex[i_finite_e]
		i_finite_e = where(finite(mean(flux_ptex, dim = 1, /nan)), n_finite_e)
		store_data, 'th'+sc+'_ptex_en_eflux', data = {x:t_ptex, y:flux_ptex[*,i_finite_e], v:en_ptex}

		lbl_ptix = strcompress(string(fix(en_ptix, type = 13)), /remove)
		lbl_ptex = strcompress(string(fix(en_ptex, type = 13)), /remove)

		;;;; reduce the number of labels to half
		lbl_ptix[linspace(n_elements(lbl_ptix)-2, 0, inc = -2)] = ''
		lbl_ptex[linspace(n_elements(lbl_ptex)-2, 0, inc = -2)] = ''

		options, 'th'+sc+'_ptix_en_eflux', ytitle = 'p!u+!n', ysubtitle = '!c!c', ztitle = '', spec = 0, labels = lbl_ptix, labsize = 0.45, labflag = -1
		options, 'th'+sc+'_ptex_en_eflux', ytitle = 'e!u'+minus_sign+'!n', ysubtitle = '!c!c', ztitle = '', spec = 0, labels = lbl_ptex, labsize = 0.45, labflag = -1
		ylim, 'th'+sc+'_ptix_en_eflux', 1.05, 1.8e6
		ylim, 'th'+sc+'_ptex_en_eflux', 1.1, 1.8e8
	endif else begin
		;;;;; d and e
		options, 'th'+sc+'_ptix_en_eflux', ytitle = 'p!u+!n', ysubtitle = '!c[eV]', ztitle = ''
		options, 'th'+sc+'_ptex_en_eflux', ytitle = 'e!u'+minus_sign+'!n', ysubtitle = '!c[eV]', ztitle = ''
	endelse

	;;; filtered data
	options, 'th'+sc+'_fgs_gsm_bpfilt', ytitle = 'B!ufilt!n'
	options, 'th'+sc+'_fgs_fac_bpfilt', ytitle = 'B!s!ufilt!r!dFAC!n'
	options, 'th'+sc+'_fgs_b_bpfilt', ytitle = '|B|!ufilt!n'
	options, 'th'+sc+'_Pth_z_bpfilt', ytitle = 'P!dth!n', ylog=0
	options, 'th'+sc+'_Pbth_bpfilt', ytitle = 'P!ufilt', ylog=0
	options, 'th'+sc+'_efs_gsm_bpfilt', ytitle = 'E!ufilt!n'
	options, 'th'+sc+'_efs_dot0_dsl_bpfilt', ytitle = 'E!s!dDSL!r!ufilt!n'
	options, 'th'+sc+'_efs_dot0_sfac_bpfilt', ytitle = 'E!s!dsFAC!r!ufilt!n'
	options, 'th'+sc+'_efs_dot0_sfac_lpfilt', ytitle = 'E!s!dsFAC!r!utrend!n'

	;;; tplot names
	if strcmp(sc, 'a') then begin
		;; THA showing no difference between beta// and betaperp, do not try again.
		tplotnames_line = 'th'+sc+'_'+['fgs_gsm', 'fgs_fac_bpfilt', 'efs_dot0_sfac_bpfilt', 'poynting_sfac', 'Pbxyth', 'Pbth_bpfilt', 'ptix_vperp_gsm']
	endif else begin
		tplotnames_line = 'th'+sc+'_'+['fgs_gsm', 'efs_dot0_dsl', 'ptix_vperp_gsm', 'Pth_z']
	endelse
	tnames_spec = 'th'+sc+'_'+['ptix_en_eflux', 'ptex_en_eflux']
	tplotnames = [tplotnames_line, tnames_spec]
	;labelnames = 'th'+sc+'_'+['state_pos_gsm_RE_z', 'state_pos_gsm_RE_y', 'state_pos_gsm_RE_x'] ;; no need for this because already have location plot

	;;; for publication only: options for two panels
	if strcmp(plot_type, 'publication') then begin
		;;; options for different probes
		if strcmp(sc, 'a') then begin
			abc = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i']
			y_abc = [0.955, 0.84, 0.74, 0.64, 0.54, 0.44, 0.34, 0.18, 0.09]
			if if_spec_lines_a then y_abc[-2] = 0.24 & y_abc[-1] = 0.14
			ylim, 'th'+sc+'_fgs_gsm', -29.9, 29.9
			ylim, 'th'+sc+'_efs_dot0_sfac_bpfilt', -1.49, 1.49
			ylim, 'th'+sc+'_Pbth_bpfilt', -0.059, 0.0599
			zlim, 'th'+sc+'_ptix_en_eflux', 1.1, 9.9e5
		endif
		if strcmp(sc, 'd') then begin
			options, tplotnames, ytitle = '', ysubtitle = '', ytickname = replicate(' ', 50)
			tplot_options,'vtitle',''
			ylim, 'th'+sc+'_'+'Pth_z', 0.4, 1.
			abc = ['b', 'd', 'f', 'h', 'j', 'l']
			y_abc = [0.93, 0.77, 0.63, 0.495, 0.29, 0.16]
		endif
		if strcmp(sc, 'e') then begin
			options, tplotnames, labels = ''
			options, tnames_spec, no_color_scale = 1
			abc = ['a', 'c', 'e', 'g', 'i', 'k']
			y_abc = [0.93, 0.77, 0.63, 0.495, 0.29, 0.16]
		endif
		;; set y limits to be the same
		if ~strcmp(sc, 'a') then begin
			ylim, 'th'+sc+'_fgs_gsm', -15, 48.
			ylim, 'th'+sc+'_efs_gsm', -14, 38.
			ylim, 'th'+sc+'_efs_dot0_dsl', -14, 38.
			if strcmp(fill_type, 'nofill') then begin
				;;; ranges for no fill
				ylim, 'th'+sc+'_Pth_z', 0.12, 1.0
				ylim, 'th'+sc+'_ptix_vperp_gsm', -340, 645
			endif else begin
				;;; ranges for fill
				;ylim, 'th'+sc+'_Pth_z', 0.22, 1.7 ;; older computation
				;ylim, 'th'+sc+'_ptix_vperp_gsm', -399.9, 900 ;; older computation
				ylim, 'th'+sc+'_Pth_z', 0.14, 1.05 ;; new computation
				ylim, 'th'+sc+'_ptix_vperp_gsm', -350, 700 ;; new computation
			endelse
			zlim, 'th'+sc+'_ptix_en_eflux', 1.01, 4e6
			zlim, 'th'+sc+'_ptex_en_eflux', 1.01, 1e8
		endif
	endif

	;;;;;; test plot: to do timing.
	;;tplot, tplotnames[[0, 3, 4, 5]], trange = ['2013 9 30 1 20', trange_show[1]]
	;tplot, tplotnames[[0, 3]], trange = ['2013 9 30 1 20', trange_show[1]]
	;;tplot, tplotnames[[0,3]], trange = trange_show
	;ctime, startandend
	;pm, time_string(startandend)
	;stop


	;;; plot
	popen, pic_folder+'/th'+sc
	print_options,xsize=size_x[i], ysize=size_y[i]
	tplot, tplotnames, trange = trange_show, title = thm_probe_color(sc, /number, /long);, var_label = labelnames
	timebar, t_onset, line = 1
	timebar_mass, 0, /databar, varname = tplotnames_line, line = 1

	;; separation between ESA and SST
	if ~(strcmp(sc,'a') and if_spec_lines_a) then begin
		timebar_mass, 23000., /databar, varname = tnames_spec, line = 2
	endif

	;;;;;;;;;; write eflux unit
	;; for D and E
	if (i eq n_elements(probes)-1) or ~strcmp(plot_type, 'publication') then begin
		xyouts, 0.93, 0.24, 'Energy Flux [eV/cm!u2!n'+dot_sign+'s'+dot_sign+'sr'+dot_sign+'eV]', orientation = 90, align = 0.5
	endif

	;; for A
	if strcmp(sc,'a') then begin
		if if_spec_lines_a then x_fluxlabel = 0.07 else x_fluxlabel = 0.93
		xyouts, x_fluxlabel, 0.17, 'Energy Flux!c!c[eV/cm!u2!n'+dot_sign+'s'+dot_sign+'sr'+dot_sign+'eV]', orientation = 90, align = 0.5
		if if_spec_lines_a then xyouts, 0.885, 0.17, 'Energy [eV]', orientation = -90, align = 0.5, charsize = 0.6
	endif


	;;;;;;;;; put different markers on the thd and the plots.
	;;;; All locations in normal coordinates
	if strcmp_or(sc, ['d', 'e']) then begin
		case sc of
		'd': begin
			;;; bars indicating the particle flux signature of DFBs
			intervals_flux_all = [['2013-09-30/01:21:18', '2013-09-30/01:22:44'], $
				['2013-09-30/01:23:26', '2013-09-30/01:24:02'], $
				['2013-09-30/01:24:18', '2013-09-30/01:24:54'], $
				['2013-09-30/01:25:08', '2013-09-30/01:25:42'], $
				['2013-09-30/01:26:13', '2013-09-30/01:27:40'], $
				['2013-09-30/01:31:40', '2013-09-30/01:32:46'], $
				['2013-09-30/01:34:36', '2013-09-30/01:35:24'], $
				['2013-09-30/01:36:03', '2013-09-30/01:37:16'], $
				['2013-09-30/01:37:34', '2013-09-30/01:38:20']]
			intervals_flux_i = time_double(intervals_flux_all[*, [0, 4, 5, 7, 8]])
			intervals_flux_i[0,1] = intervals_flux_i[0,1]-3.
			intervals_flux_e = intervals_flux_all
			v_line_i = 4e5
			v_line_e = 3e5

			;;;; DFB ranges
			tranges_dfb = [['2013-09-30/01:21:28', '2013-09-30/01:22:52'], $
				['2013-09-30/01:23:22', '2013-09-30/01:24:10'], $ 
				['2013-09-30/01:24:24', '2013-09-30/01:24:50'], $ 
				['2013-09-30/01:25:18', '2013-09-30/01:25:38'], $ 
				['2013-09-30/01:26:18', '2013-09-30/01:26:48'], $ 
				['2013-09-30/01:27:08', '2013-09-30/01:27:37'], $ 
				['2013-09-30/01:27:58', '2013-09-30/01:28:10'], $ 
				['2013-09-30/01:28:30', '2013-09-30/01:29:12'], $ 
				['2013-09-30/01:31:42', '2013-09-30/01:33:06'], $ 
				['2013-09-30/01:34:44', '2013-09-30/01:35:24'], $ 
				['2013-09-30/01:36:12', '2013-09-30/01:37:22']]

			;;;; bars indicating the ULF wave
			trange_ulf = ['2013-09-30/01:15:22', '2013-09-30/01:19:40']
			vb_ulf = 27

			end

		'e': begin
			;;; bars indicating the particle flux signature of DFBs
			intervals_flux_all = [['2013-09-30/01:20:52', '2013-09-30/01:22:28'], $
				['2013-09-30/01:24:00', '2013-09-30/01:25:08'], $
				['2013-09-30/01:26:06', '2013-09-30/01:27:39'], $
				['2013-09-30/01:28:22', '2013-09-30/01:29:00'], $
				['2013-09-30/01:31:42', '2013-09-30/01:32:48']]
			intervals_flux_i = intervals_flux_all[*, 0:-1]
			intervals_flux_e = intervals_flux_all
			v_line_i = 3e5
			v_line_e = 2e5
			
			;;;; DFB ranges
			tranges_dfb = [['2013-09-30/01:20:16', '2013-09-30/01:21:07'], $
				['2013-09-30/01:21:38', '2013-09-30/01:22:34'], $
				['2013-09-30/01:23:12', '2013-09-30/01:23:40'], $
				['2013-09-30/01:24:06', '2013-09-30/01:25:08'], $
				['2013-09-30/01:26:13', '2013-09-30/01:27:40'], $
				['2013-09-30/01:27:50', '2013-09-30/01:28:02'], $
				['2013-09-30/01:28:21', '2013-09-30/01:28:42'], $
				['2013-09-30/01:29:10', '2013-09-30/01:29:42'], $
				['2013-09-30/01:30:16', '2013-09-30/01:30:30'], $
				['2013-09-30/01:31:44', '2013-09-30/01:32:27'], $
				['2013-09-30/01:32:38', '2013-09-30/01:32:58'], $
				['2013-09-30/01:34:26', '2013-09-30/01:35:10'], $
				['2013-09-30/01:35:50', '2013-09-30/01:36:20']]

			;;;; bars indicating the ULF wave
			trange_ulf = ['2013-09-30/01:15:20', '2013-09-30/01:19:56']
			vb_ulf = 11

			end
		endcase

		;;; draw on i panels
		for i_int = 0, n_elements(intervals_flux_i[0,*])-1 do begin
			draw_in_tpanel, intervals_flux_i[*,i_int], [v_line_i, v_line_i], varname = 'th'+sc+'_ptix_en_eflux', color = 1, thick = 5
		endfor
		
		;;; draw on e panels
		for i_int = 0, n_elements(intervals_flux_e[0,*])-1 do begin
			draw_in_tpanel, intervals_flux_e[*,i_int], [v_line_e, v_line_e], varname = 'th'+sc+'_ptex_en_eflux', color = 1, thick = 5
		endfor

		;;; draw the ULF range on magnetic field panel
		draw_in_tpanel, trange_ulf, [vb_ulf, vb_ulf], varname = tplotnames[0], thick = 0.5
		draw_in_tpanel, [trange_ulf[0], trange_ulf[0]], vb_ulf+[-3,3], varname = tplotnames[0], thick = 0.5
		draw_in_tpanel, [trange_ulf[1], trange_ulf[1]], vb_ulf+[-3,3], varname = tplotnames[0], thick = 0.5
		text_in_tpanel, mean(time_double(trange_ulf)), vb_ulf+7, 'ULF!coscillations', varname = tplotnames[0], charsize = 0.5, align = 0.5

		;;;;; mark the DFBs
		;;; bars for easy looking
		;timebar, tranges_dfb[0,*], varname = tplotnames[0], thick = 0.2
		;timebar, tranges_dfb[1,*], varname = tplotnames[0], line = 2, thick = 0.2
		;;; shadows
		shadow_tranges_tpanel, tranges = tranges_dfb, varname = tplotnames[0], color = 5

	endif ;; if d or e.


	;;;;;;;;; write abc
	if strcmp(plot_type, 'publication') then begin
		xyouts, replicate(x_abc, n_elements(abc)), y_abc, '('+abc+')', /normal
	endif

	pclose

	;;;; plot pitch angle distribution
	;if strcmp(sc, 'a') then begin
	;	tvnames_line = 'th'+sc+'_'+['fgs_gsm', 'fgs_fac_bpfilt', 'efs_dot0_sfac_bpfilt', 'Pbth']
	;	tplot, [tvnames_line, 'th'+sc+'_peir_eflux_pa_*', 'th'+sc+'_psif_eflux_pa_*'], title = 'THA, ion PA'
	;	timebar, t_onset, line = 1
	;	makepng, pic_folder+'/th'+sc+'_pa_i'
	;	tplot, [tvnames_line, 'th'+sc+'_peef_eflux_pa_*', 'th'+sc+'_psef_eflux_pa_*'], title = 'THA, electron PA'
	;	timebar, t_onset, line = 1
	;	makepng, pic_folder+'/th'+sc+'_pa_e'
	;endif

	
	;;;;; compute wave velocity and output on screen.
	trange_v = date+['1 10', '1 17']
	time_clip, 'th'+sc+'_fgs_gsm', trange_v[0], trange_v[1], newname = 'th'+sc+'_fgs_gsm_cv'
	time_clip, 'th'+sc+'_Pth_z', trange_v[0], trange_v[1], newname = 'th'+sc+'_Pth_z_cv'
	time_clip, 'th'+sc+'_ptix_density', trange_v[0], trange_v[1], newname = 'th'+sc+'_ptix_density_cv'
	get_data, 'th'+sc+'_fgs_gsm_cv', t_nouse, b_data
	get_data, 'th'+sc+'_Pth_z_cv', t_nouse, Pth_data
	get_data, 'th'+sc+'_ptix_density_cv', t_nouse, ni_data
	b_ttl = sqrt(total(b_data^2, 2))
	b_mean = mean(b_ttl, /nan) ;; nT
	Pth_mean = mean(Pth_data, /nan) ;; nPa
	ni_mean = mean(ni_data, /nan) ;; /cc
	;;; compute velocities using these values
	gma = 5./3 ;; adiabatic constant gamma.
	va = 21.8*b_mean/sqrt(ni_mean) ;; Alfven speed
	cs = 1e-3*sqrt(gma*Pth_mean*1e-9/(mass_proton*ni_mean*1e6)) ;; sound speed
	v_fast = sqrt(va^2+cs^2)
	print, 'Ahead of the DFB on TH'+strupcase(sc)+':'
	print, 'VA: '+string(va)
	print, 'cs: '+string(cs)
	print, 'V_fast: '+string(v_fast)

endfor
;stop
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



;;;;;;;;;;;;;;;;;;;;;;;;; RBSP, a, b ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;; revise the plot range to hide later oscillations on rbsp-B.
trange_show = date+['1 0 1', '1 39']

;;;;;;;; checked RB-A/B all flux line plots (hope p, e, he, o; magneis p, e; rept p, e); no injection, no need to check again.
;;; plot size settings
size_x_rb = 5
size_y_rb = 7.9
size_x_rb_spec = 5
size_y_rb_spec = 5
size_x_discuss = 5
size_y_discuss = 4
x_abc = 0.247

if_spec_lines = 1

figure_plot = 'data' ;; data figure, original (Figure 7 in paper)
;figure_plot = 'data_normal' ;; data figure, with normalized spectra
;figure_plot = 'allchannels' ;; data figure, with normalized spectra
;figure_plot = 'allchannels_pa' ;; data figure, with normalized spectra
;figure_plot = 'discussion' ;; discussion figure of wave mode. (THIS is no-go)
;figure_plot = 'powerspec' ;; powerspectra of magnetic fields for the supporting material (added figure in paper)
;figure_plot = 'plasmapause' ;; check plasmapause location


;;; load settings
probes = ['r', 's']
;probes = 'r'

;;; time range for MVA to get wave propagation direction
trange_mva = date+['1 20', '1 30']

if strcmp(figure_plot, 'discussion') then probes = 'r'
;energy_min = 200. ;; in eV, the energy below this will be removed, all particles
energy_min = [0, 200., 0, 0] ;; in eV, the energy below this will be removed, remove electron only

;;; change trange_show for the discussion plot
;trange_show  = date+['1 17', '1 30']


for i = 0, n_elements(probes)-1 do begin
	;del_data, '*' ;; for faster load
	sc = probes[i]
	case sc of 
	'r': sc_rb = 'a'
	's': sc_rb = 'b'
	endcase

	;;; choose the variable for electric field
	ename = 'efw_esvy_mgse_vxb_removed_coro_removed_spinfit'
	;case sc of 
	;'r': ename = 'efw_esvy_mgse_vxb_removed_coro_removed_spinfit'
	;'s': ename = 'efw_esvy_gsm'
	;endcase

	;;; load data
	rbsp_load, trange=trange_load, probe=sc, datatype = 'fgs', /tclip, rbsp_folder = rbsp_folder ;; both pos and fgs data are load
	rbsp_load, trange=trange_load, probe=sc, datatype = 'fgl', /tclip, rbsp_folder = rbsp_folder
	copy_data, 'th'+sc+'_state_pos', 'th'+sc+'_state_pos_gsm'
	rbsp_load, trange=trange_load, probe=sc, datatype = 'efw_density', /tclip, rbsp_folder = rbsp_folder, level = 3
	rbsp_load, trange=trange_load, probe=sc, datatype = 'rept', /tclip, rbsp_folder = rbsp_folder, level = 2
	rbsp_load, trange=trange_load, probe=sc, datatype = 'mageis', /tclip, rbsp_folder = rbsp_folder, level = 2
	rbsp_load, trange=trange_load, probe=sc, datatype = 'hope_sa', /tclip, rbsp_folder = rbsp_folder, level = 2, reduce_connect = 'algebra'

	;;;;;; Also load Xiaojia's data for electron
	restore, 'rbsp'+sc_rb+'_fesa_com3.sav'
	fesa = transpose(fesa_out)
	store_data, 'th'+sc+'_espec_combined_xj', data = {x:fesa[*,0], y:fesa[*,1:99], v:fesa[*,100:198]}
	time_clip, 'th'+sc+'_espec_combined_xj', trange_load[0], trange_load[1]

	;;;; PA data
	if strmatch(figure_plot, '*pa') then begin
		rbsp_load, trange=trange_load, probe=sc, datatype = 'mageis', rbsp_folder = rbsp_folder, level = 3
		rbsp_load, trange=trange_load, probe=sc, datatype = 'hope_pa', rbsp_folder = rbsp_folder, level = 3
		rbsp_load, trange=trange_load, probe=sc, datatype = 'rept', rbsp_folder = rbsp_folder, level = 3
	endif

	if strmatch_or(figure_plot, ['data*', 'powerspec']) then begin
		;;;; efield data (slowest)
		rbsp_load_efw, trange=trange_load, probe=sc, /tclip, bad_determine = 'bonnell', rtrange = rtrange

		;;; compute spin-plane field-alighed B and E
		;;; For cartoon reference: Spin axis of rbsp-a in GSM: [0.9, 0.21, 0.27]
		cotrans2fac, 'th'+sc+'_fgs_mgse', b_variable = 'th'+sc+'_fgs_mgse', smooth_time = 6.*60., newname = 'th'+sc+'_fgs_sfac', /sc_spin
		cotrans2fac, 'th'+sc+'_'+ename, b_variable = 'th'+sc+'_fgs_mgse', smooth_time = 6.*60., newname = 'th'+sc+'_'+ename+'_sfac', /sc_spin
	endif

	calc,'"th'+sc+'_state_pos_gsm_RE" = "th'+sc+'_state_pos_gsm"/6374.4'

	;;; compute total field
	get_data, 'th'+sc+'_fgs_gsm', t, bdata, dlimits = dl
	store_data, 'th'+sc+'_fgs_b', data={x:t, y:sqrt(total(bdata^2, 2))}

	;;; transform field to field aligned coordinate system
	cotrans2fac, 'th'+sc+'_fgs_gsm', b_variable = 'th'+sc+'_fgs_gsm', pos = 'th'+sc+'_state_pos_gsm', smooth_time = 6.*60., newname = 'th'+sc+'_fgs_fac'
	cotrans2fac, 'th'+sc+'_fgl_gsm', b_variable = 'th'+sc+'_fgl_gsm', pos = 'th'+sc+'_state_pos_gsm', smooth_time = 6.*60., newname = 'th'+sc+'_fgl_fac'


	;;; combine spectrums
	combine_spec, 'th'+sc+'_mageis_pspec_tclip', 'th'+sc+'_hope_sa_pspec_tclip', newname = 'th'+sc+'_pspec_combinedlow_tclip', /eV2keV_2nd, sep = sep_p_low
	combine_spec, 'th'+sc+'_mageis_espec_tclip', 'th'+sc+'_hope_sa_espec_tclip', newname = 'th'+sc+'_espec_combinedlow_tclip', /eV2keV_2nd, sep = sep_e_low
	combine_spec, 'th'+sc+'_pspec_combinedlow_tclip', 'th'+sc+'_rept_pspec_tclip', newname = 'th'+sc+'_pspec_combined_4moment_tclip', /MeV2keV_2nd, sep = sep_p_high
	combine_spec, 'th'+sc+'_espec_combinedlow_tclip', 'th'+sc+'_rept_espec_tclip', newname = 'th'+sc+'_espec_combined_tclip', /MeV2keV_2nd, sep = sep_e_high
	get_data, 'th'+sc+'_rept_pspec_tclip', t_reptp, data_reptp, v_reptp
	store_data, 'th'+sc+'_rept_pspec_keV_tclip', t_reptp, data_reptp/1000., v_reptp*1000.
	if if_spec_lines then begin
		copy_data, 'th'+sc+'_pspec_combined_4moment_tclip', 'th'+sc+'_pspec_combined_tclip'
	endif else begin
		store_data, 'th'+sc+'_pspec_combined_tclip', data = 'th'+sc+['_pspec_combinedlow_tclip', '_rept_pspec_keV_tclip']
	endelse
	;;;; all population for pressures
	compute_moments, 'th'+sc+['_pspec_combined_4moment_tclip', '_espec_combined_tclip', '_hope_sa_hespec_tclip', '_hope_sa_ospec_tclip'], intypes = 'flux', inunits_energy = ['keV', 'keV', 'eV', 'eV'], inunits_flux = 'keV', energy_min = energy_min, particle_types = ['p', 'e', 'He+', 'O+'], /combine, newname_combine = 'all', rbsp_folder = rbsp_folder

	;;;;;;;;;;;; make flux ratios ;;;;;;;;
	tvnames_ratios = 'th'+sc+['_pspec_combined_4moment_tclip', '_espec_combined_xj_tclip']
	for k = 0, n_elements(tvnames_ratios)-1 do begin
		get_data, tvnames_ratios[k], tp, yp, vp
		yp_start = transpose(mean(yp[0:17,*], dim = 1, /nan)) ;; 3 min average
		yp_divide = rebin(yp_start, n_elements(tp), n_elements(yp_start))
		store_data, tvnames_ratios[k]+'_normalize', data={x:tp, y:yp/yp_divide, v:vp}
		options, tvnames_ratios[k]+'_normalize', ylog = 1, spec = 1, zlog = 1, ysubtitle = '!c[keV]', ztitle = ''
	endfor

	;;; compute alven and sound speeds
	tinterpol_mxn, 'th'+sc+'_efw_density', 'th'+sc+'_fgs_gsm', newname='th'+sc+'_efw_density_int'
	get_data, 'th'+sc+'_fgs_gsm', t, data_b
	get_data, 'th'+sc+'_efw_density_int', t, data_n
	v_alfven = 21.8*sqrt(total(data_b^2, 2))/sqrt(data_n)
	store_data, 'th'+sc+'_alfven', data = {x:t, y:v_alfven}

	;;; locations
	split_vec, 'th'+sc+'_state_pos_gsm_RE'

	;;;;;; Make wave power spectrum
	if strcmp(figure_plot, 'powerspec') then begin
		;;; B field
		split_vec, 'th'+sc+'_fgs_fac'
		thm_powerspec, 'th'+sc+'_fgs_fac_x', period_range = powerspec_range
		thm_powerspec, 'th'+sc+'_fgs_fac_y', period_range = powerspec_range
		thm_powerspec, 'th'+sc+'_fgs_fac_z', period_range = powerspec_range
		;;; E field
		split_vec, 'th'+sc+'_'+ename+'_sfac'
		thm_powerspec, 'th'+sc+'_'+ename+'_sfac_y', period_range = powerspec_range
		thm_powerspec, 'th'+sc+'_'+ename+'_sfac_z', period_range = powerspec_range
	endif

	;;;;;; filter
	tbandpass_filter, 'th'+sc+'_fgs_gsm', mins_filt_lowf*60., mins_filt_highf*60.
	tbandpass_filter, 'th'+sc+'_fgs_fac', mins_filt_lowf*60., mins_filt_highf*60.
	;;;; for matlab filter: out put data
	;dataout4matlab, 'th'+sc+'_fgs_fac'
	;dataout4matlab, 'th'+sc+'_fgs_fac_bpfilt'

	tbandpass_filter, 'th'+sc+'_fgs_sfac', mins_filt_lowf*60., mins_filt_highf*60.
	tbandpass_filter, 'th'+sc+'_fgl_gsm', mins_filt_lowf*60., mins_filt_highf*60.
	tbandpass_filter, 'th'+sc+'_fgl_fac', mins_filt_lowf*60., mins_filt_highf*60.
	tbandpass_filter, 'th'+sc+'_fgs_b', mins_filt_lowf*60., mins_filt_highf*60.
	tbandpass_filter, 'th'+sc+'_efw_density', mins_filt_lowf*60., mins_filt_highf*60.
	if strmatch(figure_plot, 'data*') then begin
		tbandpass_filter, 'th'+sc+'_all_Pth', mins_filt_lowf*60., mins_filt_highf*60.
		tbandpass_filter, 'th'+sc+'_pspec_combined_density', mins_filt_lowf*60., mins_filt_highf*60.
		tbandpass_filter, 'th'+sc+'_espec_combined_density', mins_filt_lowf*60., mins_filt_highf*60.
		tbandpass_filter, 'th'+sc+'_all_density', mins_filt_lowf*60., mins_filt_highf*60.
		tbandpass_filter, 'th'+sc+'_'+ename, mins_filt_lowf*60., mins_filt_highf*60.
		tbandpass_filter, 'th'+sc+'_'+ename+'_sfac', mins_filt_lowf*60., mins_filt_highf*60.
		tbandpass_filter, 'th'+sc+'_efw_esvy_gsm', mins_filt_lowf*60., mins_filt_highf*60.
	endif
	tlowpass_filter, 'th'+sc+'_'+ename+'_sfac', mins_filt_trend*60.

	;; for RB-B filter again because trend is too strong
	if strcmp(sc,'s') then begin
		thigh_pass_filter, 'th'+sc+'_fgs_fac_bpfilt', mins_filt_lowf*60., newname = 'th'+sc+'_fgs_fac_bpfilt' ;; using tbandpass_filter will give a different result because of backward=0
		thigh_pass_filter, 'th'+sc+'_fgs_sfac_bpfilt', mins_filt_lowf*60., newname = 'th'+sc+'_fgs_sfac_bpfilt'
	endif

	;;; remove the atitude trend from data, fgs_fac (xy components only)
	split_vec, 'th'+sc+'_fgs_fac_bpfilt'
	get_data, 'th'+sc+'_fgs_fac_bpfilt_z', t, bfac_z_filt
	tbandpass_filter, 'th'+sc+'_fgs_fac_bpfilt', mins_atitude*60., 0., newname = 'th'+sc+'_fgs_fac_bpfilt_lsplow', /lsp
	if strcmp(sc,'s') then begin
		;; again for ths
		tbandpass_filter, 'th'+sc+'_fgs_fac_bpfilt_lsplow', mins_atitude*60., 0., newname = 'th'+sc+'_fgs_fac_bpfilt_lsplow', /lsp
	endif
	split_vec, 'th'+sc+'_fgs_fac_bpfilt_lsplow'
	get_data, 'th'+sc+'_fgs_fac_bpfilt_lsplow_x', t, bfac_x_filt_lsplow
	get_data, 'th'+sc+'_fgs_fac_bpfilt_lsplow_y', t, bfac_y_filt_lsplow
	;; complete the filtered Bfac three components.
	store_data, 'th'+sc+'_fgs_fac_bpfilt', data = {x:t, y:[[bfac_x_filt_lsplow], [bfac_y_filt_lsplow], [bfac_z_filt]]}
	;;;;; repeat the above for fgs_bsfac? No need! it does not change anything.

	;;;; get wave propagation velocity by applying MVA on the filtered data
	;time_clip, 'th'+sc+'_fgs_fac_bpfilt', trange_mva[0], trange_mva[1], newname = 'th'+sc+'_fgs_fac_bpfilt_4mva'
	;get_data, 'th'+sc+'_fgs_fac_bpfilt_4mva', t_mva, b_4mva
	;minvar, transpose(b_4mva), lmn, lambdas2 = lambdas
	;print, 'Minimum varience direction in the FAC system: '
	;print, lmn[*,2]
	;print, 'lambda_m/lambda_n:'
	;print, lambdas[1]/lambdas[2]
	;stop

	;;; compute wave poynting flux in the spin axis direction
	if strmatch(figure_plot, 'data*') then begin
		tinterpol_mxn, 'th'+sc+'_'+ename+'_sfac_bpfilt', 'th'+sc+'_fgs_sfac_bpfilt', newname='th'+sc+'_'+ename+'_sfac_bpfilt_int'
		get_data, 'th'+sc+'_'+ename+'_sfac_bpfilt_int', t, e_wave
		get_data, 'th'+sc+'_fgs_sfac_bpfilt', t, b_wave
		poynting = Re^2*1e-6*transpose(crossp_long(transpose(e_wave), transpose(b_wave), n_dir = poynting_dir)) ;; unit: W/RE^2
		store_data, 'th'+sc+'_poynting_sfac', data={x:t, y:poynting}, dlimits = {data_att:{coord_sys:'sfac'}}
	endif

	;;; original quantities labels
	options, 'th'+sc+'_state_pos_gsm_RE_x', ytitle = 'X [R!dE!n]', ysubtitle = ''
	options, 'th'+sc+'_state_pos_gsm_RE_y', ytitle = 'Y [R!dE!n]', ysubtitle = ''
	options, 'th'+sc+'_state_pos_gsm_RE_z', ytitle = 'Z [R!dE!n]', ysubtitle = ''
	options, 'th'+sc+'_fg?_gsm*',colors=[2,4,6],labels=['B!dx','B!dy','B!dz'], ytitle = 'B', ysubtitle = '!c[nT]', labflag = 1, thick = l_thick
	options, 'th'+sc+'_fg?_fac*',colors=[2,4,6],labels=["B!dr'",'B!d'+phi_letter+"'",'B!db'], ytitle = 'B!dFAC', ysubtitle = '!c[nT]', labflag = 1, thick = l_thick
	options, 'th'+sc+'_fg?_b*', ytitle = '|B|', ysubtitle = '!c[nT]', thick = l_thick
	options, 'th'+sc+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit*',colors=[2,4,6],labels=['E!dx','E!dy','E!dz'], ytitle = 'E!dmGSE!n', ysubtitle = '!c[mV/m]', labflag = 1, thick = l_thick
	options, 'th'+sc+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_sfac*',colors=[2,4,6],labels=['E!d-ax','E!dsp'+perp_sign,"E!db'"], ytitle = 'E!dsFAC!n', ysubtitle = '!c[mV/m]', labflag = 1, thick = l_thick
	options, 'th'+sc+'_efw_esvy_gsm*',colors=[2,4,6],labels=['E!dx','E!dy','E!dz'], ytitle = 'E', ysubtitle = '!c[mV/m]', labflag = 1, thick = l_thick
	options, 'th'+sc+'_poynting_sfac', ytitle = 'S!dsFAC!n',colors=[2,4,6],labels=['S!d-ax','S!dsp'+perp_sign,"S!db'"], ysubtitle = '!c[W/R!s!u2!r!dE!n]', labflag = 1, thick = l_thick
	options, 'th'+sc+'_efw_density*', ytitle = 'n!defw!n', ysubtitle = '!c[/cc]', thick = l_thick
	options, 'th'+sc+'_alfven*', ytitle = 'V!dAlfven!n', ysubtitle = '!c[km/s]', thick = l_thick
	options, 'th'+sc+'_all_Pth*', ytitle = 'P!dth!n', ysubtitle = '!c[nPa]', thick = l_thick
	options, 'th'+sc+'_pspec_combined_density*', ytitle = 'n!di!n', ysubtitle = '!c[/cc]', thick = l_thick
	options, 'th'+sc+'_espec_combined_density*', ytitle = 'n!de!n ', ysubtitle = '!c[/cc]', thick = l_thick
	options, 'th'+sc+'_all_density*', ytitle = 'n!dall!n', ysubtitle = '!c[/cc]', thick = l_thick
	if ~strcmp(figure_plot, 'allchannels') then begin
		options, 'th'+sc+'_pspec_combined_tclip', ytitle = 'p!u+!n', ysubtitle = '!c[keV]', ztitle = '', spec=1, ylog=1, zlog=1
		options, 'th'+sc+'_espec_combined_tclip', ytitle = 'e!u'+minus_sign+'!n', ysubtitle = '!c[keV]', ztitle = '', spec=1, ylog=1, zlog=1
		options, 'th'+sc+'_espec_combined_xj_tclip', ytitle = 'e!u'+minus_sign+'!n', ysubtitle = '!c[keV]', ztitle = '', spec=1, ylog=1, zlog=1
		options, 'th'+sc+'_pspec_combined_4moment_tclip_normalize', ytitle = 'p!u+!n Ratio', ysubtitle = '!c[keV]'
		options, 'th'+sc+'_espec_combined_xj_tclip_normalize', ytitle = 'e!u'+minus_sign+'!n Ratio', ysubtitle = '!c[keV]'
		options, 'th'+sc+'_rept_?spec*', spec=1, ylog = 1, zlog=1
		options, 'th'+sc+'_*_pa_*spec_*', spec=1, zlog=1
		options, 'th'+sc+'_*_pa_*spec_*', spec=1, zlog=1
		;;;;; limit the range of specs, in keV
		;;; hope+mageis+rept combined
		ylim, 'th'+sc+'_pspec_combined_tclip', 0.24, 2e5
		ylim, 'th'+sc+'_espec_combined_tclip', 0.6, 20000.
		ylim, 'th'+sc+'_pspec_combined_4moment_tclip*', 0.24, 2e5
		ylim, 'th'+sc+'_espec_combined_xj_tclip*', 0.6, 20000.
		zlim,'th'+sc+'_pspec_combined_tclip', 1.05e-6, 0.999e6 
		zlim,'th'+sc+'_espec_combined_tclip', 8e-5, 1.05e6 
		zlim,'th'+sc+'_espec_combined_xj_tclip', 8e-5, 1.05e6 
	endif

	;;; power spectra
	if strcmp(figure_plot, 'powerspec') then begin
		;yrange_spec = [0.001, 0.05]
		yrange_spec = [1/150., 1/40.] ;; show Pi2 band only
		;;; magnetic field
		time_clip, 'th'+sc+'_fg?_fac_x_wv_pow', trange_load[0], trange_load[1]
		time_clip, 'th'+sc+'_fg?_fac_y_wv_pow', trange_load[0], trange_load[1]
		time_clip, 'th'+sc+'_fg?_fac_z_wv_pow', trange_load[0], trange_load[1]
		zrange = [0.0011, 50]
		options, 'th'+sc+'_fg?_fac_x_wv_pow*', ytitle = "!af!nB!dr'", ysubtitle = '!c[Hz]', zrange = zrange, yrange = yrange_spec, ylog = 0, ztitle = ''
		options, 'th'+sc+'_fg?_fac_y_wv_pow*', ytitle = '!af!nB!d'+phi_letter+"'", ysubtitle = '!c[Hz]', zrange = zrange, yrange = yrange_spec, ylog = 0, ztitle = ''
		options, 'th'+sc+'_fg?_fac_z_wv_pow*', ytitle = '!af!nB!db', ysubtitle = '!c[Hz]', zrange = zrange, yrange = yrange_spec, ylog = 0, ztitle = ''
		;;; electric field
		;zrange = [0.001, 2000]
		options, 'th'+sc+'_'+ename+'_sfac_y_wv_pow*', ytitle = '!af!nE!dsp'+perp_sign, ysubtitle = '!c[Hz]', zlog = 1, ylog = 1, yrange = yrange_spec
		options, 'th'+sc+'_'+ename+'_sfac_z_wv_pow*', ytitle = "!af!nE!db'", ysubtitle = '!c[Hz]', zlog = 1, ylog = 1, yrange = yrange_spec
	endif

	;;; Spec data, special treatments if line spec is set
	if if_spec_lines then begin
		get_data, 'th'+sc+'_pspec_combined_tclip', t_ptix, flux_ptix, en_ptix
		get_data, 'th'+sc+'_espec_combined_xj_tclip', t_ptex, flux_ptex, en_ptex
		n_ch_i = n_elements(flux_ptix[0,*])
		n_ch_e = n_elements(flux_ptex[0,*])

		;;;;;;;;; clean offset data points for the electron spec. (caused by mode change)
		for i_ch = 0, n_ch_e-1 do begin
			;;; find the major energy (actual energy) of this channel.
			ens_this = en_ptex[*,i_ch]
			en_vs_this = ens_this[uniq(ens_this, sort(ens_this))]
			if n_elements(en_vs_this) gt 0 then begin
				n_major = 0 ;; to be changed
				i_major = 100 ;; to be changed, the i in en_vs_this which gives the actuall energy of this channel.
				for i_env = 0, n_elements(en_vs_this)-1 do begin
					i_match = where(ens_this eq en_vs_this[i_env], n_match)
					if n_match gt n_major then begin
						n_major = n_match
						i_major = i_env
					endif
				endfor
				en_maj = en_vs_this[i_major]

				;;; find chanels that are not the major energies
				i_correct = where(ens_this eq en_maj, n_correct)
				i_off = where(ens_this ne en_maj, n_off)
				if n_correct gt 0 and n_off gt 0 then begin
					flux_ptex[i_off, i_ch] = interpol(flux_ptex[i_correct, i_ch], t_ptex[i_correct], t_ptex[i_off], /nan)
					;flux_ptex[i_off, i_ch] = !values.f_nan
					en_ptex[i_off, i_ch] = en_maj
				endif
			endif
		endfor


		;;;;;;; multiply values to the REPT bins.
		i_rept = where(en_ptix[0,*] gt 1e4, n_rept)
		if n_rept gt 0 then begin
			if strcmp(sc, 'r') then factor_rept = 100 else factor_rept = 10000
			flux_ptix[*, i_rept] = flux_ptix[*, i_rept]*factor_rept
		endif

		;;;;;;; remove bins to make plot cleanner
		;;;;; generate labels
		;;;; for diagnose
		lbl_ptix = strcompress(string(indgen(n_ch_i)))
		lbl_ptex = strcompress(string(indgen(n_ch_e)))

		case sc of
		'r': begin
			i_bins_keep_i = [linspace(3, 26, inc = 3), linspace(27, 48, inc = 6), linspace(49, 83, inc = 3),linspace(84, 91, inc = 2), linspace(93, 98, inc = 3), 99, 100]
			i_bins_keep_e = [linspace(0, 20, inc = 3), linspace(21, 48, inc = 6), linspace(49, 89, inc = 3), 90, 96]
			end
		's': begin
			i_bins_keep_i = [[16:33], linspace(34, 48, inc = 2), linspace(49, 71, inc = 3), linspace(72, 76, inc = 2), linspace(77, 92, inc = 3), 93, linspace(94, 98, inc = 2), 99, 100]
			i_bins_keep_e = [linspace(0, 11, inc = 4), linspace(12, 34, inc = 5), linspace(35, 72, inc = 6), [73:77], linspace(78, 94, inc = 3), linspace(95, n_ch_e-1, inc = 3)]
			end
		endcase

		flux_ptix = flux_ptix[*, i_bins_keep_i]
		en_ptix = en_ptix[*, i_bins_keep_i]
		lbl_ptix = lbl_ptix[i_bins_keep_i]
		flux_ptex = flux_ptex[*, i_bins_keep_e]
		en_ptex = en_ptex[*, i_bins_keep_e]
		lbl_ptex = lbl_ptex[i_bins_keep_e]


		;;;;;; store the data again.
		store_data, 'th'+sc+'_pspec_combined_tclip', data = {x:t_ptix, y:flux_ptix, v:en_ptix}
		store_data, 'th'+sc+'_espec_combined_xj_tclip', data = {x:t_ptex, y:flux_ptex, v:en_ptex}
		n_ch_i = n_elements(en_ptix[0,*])
		n_ch_e = n_elements(en_ptex[0,*])

		;;;;;;;;;;;;;;;;;; set labels and limits ;;;;;;;;;;;;;;;;

		;;;;; generate labels for real plot
		;;; ion
		lbl_ptix = strcompress(string(fix(en_ptix[0,*], type = 13)), /remove)+' keV'
		i_eV = where(en_ptix[0,*] lt 1., n_eV)
		i_MeV = where(en_ptix[0,*] ge 1000., n_MeV)
		lbl_ptix[i_eV] = strcompress(string(fix(en_ptix[0,i_eV]*1000, type = 13)), /remove)+' eV'
		lbl_ptix[i_MeV] = strcompress(string(fix(en_ptix[0,i_MeV]/1000, type = 13)), /remove)+' MeV'
		i_rept = where(en_ptix[0,*] gt 1e4, n_rept_i)
		if n_rept gt 0 then lbl_ptix[i_rept] = lbl_ptix[i_rept]+'!u*'

		;;; electron
		lbl_ptex = strcompress(string(fix(en_ptex[0,*], type = 13)), /remove)+' keV'
		i_eV = where(en_ptex[0,*] lt 1., n_eV)
		i_MeV = where(en_ptex[0,*] ge 1000., n_MeV)
		lbl_ptex[i_eV] = strcompress(string(fix(en_ptex[0,i_eV]*1000, type = 13)), /remove)+' eV'
		lbl_ptex[i_MeV] = strcompress(string(fix(en_ptex[0,i_MeV]/1000, type = 13)), /remove)+' MeV'
		i_rept = where(en_ptex[0,*] gt 1e4, n_rept_i)
		if n_rept gt 0 then lbl_ptex[i_rept] = lbl_ptex[i_rept]+'!u*'

		;;; reduce the number of labels for better presentation
		if strcmp(sc, 'r') then inc_i = 2 else inc_i = 3
		i_keep_i = linspace(1, n_ch_i-1, inc = inc_i)
		for i_ch = 0, n_elements(lbl_ptix)-1 do begin
			k_match = where(i_keep_i eq i_ch, n_match)
			if n_match lt 1 then lbl_ptix[i_ch] = ''
		endfor
		i_remove_e = linspace(1, n_ch_e-1, inc = 2)
		lbl_ptex[i_remove_e] = ''

		options, 'th'+sc+'_pspec_combined_tclip', ytitle = 'p!u+!n', ysubtitle = '!c', ztitle = '', spec = 0, labels = lbl_ptix, labsize = 0.4, labflag = -1
		options, 'th'+sc+'_espec_combined_xj_tclip', ytitle = 'e!u'+minus_sign+'!n', ysubtitle = '!c', ztitle = '', spec = 0, labels = lbl_ptex, labsize = 0.4, labflag = -1

		case sc of
		'r': begin
			ylim, 'th'+sc+'_pspec_combined_tclip', 1.1e-4, 5e7
			ylim, 'th'+sc+'_espec_combined_xj_tclip', 1.1e-4, 1.8e9
			end
		's': begin
			ylim, 'th'+sc+'_pspec_combined_tclip', 1.1e-2, 1e6
			ylim, 'th'+sc+'_espec_combined_xj_tclip', 1.1e-4, 1.2e9
			end
		endcase

		options, 'th'+sc+'_pspec_combined_tclip', spec = 0
		options, 'th'+sc+'_espec_combined_xj_tclip', spec = 0
	end

	;;; filtered data
	options, 'th'+sc+'_fg?_gsm_bpfilt', ytitle = 'B!ufilt!n'
	options, 'th'+sc+'_fg?_fac_bpfilt', ytitle = 'B!s!ufilt!r!dFAC!n'
	options, 'th'+sc+'_fg?_b_bpfilt', ytitle = '|B|!ufilt!n'
	options, 'th'+sc+'_all_Pth_bpfilt', ytitle = 'P!s!dth!r!ufilt!n'
	options, 'th'+sc+'_all_density_bpfilt', ytitle = 'n!ufilt!n'
	options, 'th'+sc+'_pspec_combined_density_bpfilt', ytitle = 'n!s!di!r!ufilt!n'
	options, 'th'+sc+'_espec_combined_density_bpfilt', ytitle = 'n!s!de!r!ufilt!n'
	options, 'th'+sc+'_efw_esvy_gsm_bpfilt', ytitle = 'E!ufilt!n'
	options, 'th'+sc+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_bpfilt', ytitle = 'E!s!dmGSE!r!ufilt!n'
	options, 'th'+sc+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_bpfilt_sFAC', ytitle = 'E!s!dsFAC!r!ufilt!n'
	options, 'th'+sc+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit_lpfilt_sFAC', ytitle = 'E!s!dsFAC!r!utrend!n'
	options, 'th'+sc+'_efw_density_bpfilt', ytitle = 'n!s!defw!r!ufilt!n'



	;;;;;;;;;;;;;;;; data picture plot ;;;;;;;;;;;;;;;;;;;;;;;
	if strmatch(figure_plot, 'data*') then begin
		;;;; tplot names
		;tnames_line = 'th'+sc+'_'+['fgs_gsm', 'fgs_fac_bpfilt', 'fgs_sfac_bpfilt', ename+'_sfac_bpfilt', 'poynting_sfac', 'efw_density', 'efw_density_bpfilt', 'alfven'] ;; all quantities
		tnames_line = 'th'+sc+'_'+['fgs_gsm', 'fgs_fac_bpfilt', ename+'_sfac_bpfilt', 'poynting_sfac', 'alfven'] ;; publication quantities
		if strcmp(figure_plot, 'data_normal') then begin
			tnames_spec = 'th'+sc+'_'+['pspec_combined_4moment_tclip_normalize', 'espec_combined_xj_tclip_normalize'] ;; increase ratio
		endif else begin
			tnames_spec = 'th'+sc+'_'+['pspec_combined_tclip', 'espec_combined_xj_tclip'] ;; original
		endelse

		;;;; plot only injection data, longer range, mageis only
		if strcmp(plot_type, 'long_injection') then begin
			tnames_line = 'th'+sc+'_fgs_gsm'
			tnames_spec = 'th'+sc+'_'+['mageis_pspec', 'mageis_espec']
			options, tnames_spec, spec = 0, ylog = 1
			size_y_rb = 5.
		endif

		labelnames = 'th'+sc+'_'+['state_pos_gsm_RE_z', 'state_pos_gsm_RE_y', 'state_pos_gsm_RE_x']
		tplotnames = [tnames_line, tnames_spec]


		;;; for publication only: options for two panels and limits
		if strcmp(plot_type, 'publication') then begin
			if strcmp(sc, 'r') then begin
				ylim, 'th'+sc+'_fgs_gsm', -90, 170.
				ylim, 'th'+sc+'_'+ename+'_sfac_bpfilt', -0.14, 0.099
				ylim, 'th'+sc+'_poynting_sfac', -1.9, 0.2
				if ~strcmp(figure_plot, 'data_normal') then options, tnames_spec, no_color_scale = 1
			endif else begin
				ylim, 'th'+sc+'_fgs_gsm', -499., 1999.
				ylim, 'th'+sc+'_fgs_fac_bpfilt', -0.7, 0.599
				ylim, 'th'+sc+'_'+ename+'_sfac_bpfilt', -0.25, 0.31
				ylim, 'th'+sc+'_poynting_sfac', -7.5, 3.5
			endelse

			if i eq 0 then begin
				abc = ['a', 'b', 'c', 'd', 'e', 'f', 'g']
			endif else begin
				abc = ['h', 'i', 'j', 'k', 'l', 'm', 'n']
			endelse

			if i ne 0 then begin
				options, tplotnames, ytitle = '', ysubtitle = '';, ytickname = replicate(' ', 50)
				if ~if_spec_lines then options, tnames_spec, ytickname = replicate(' ', 50)
				tplot_options,'vtitle',''
			endif

			if i ne n_elements(probes)-1 or i eq 0 then begin
				options, tnames_line, labels = ''
				;options, tnames_spec, no_color_scale = 1
			endif

			;;; abc labels
			case n_elements(abc) of
			5: y_abc = [0.885, 0.743, 0.57, 0.397, 0.15]
			7: begin
				if if_spec_lines then begin
					y_abc = [0.915, 0.81, 0.685, 0.49, 0.415, 0.275, 0.16]
				endif else begin
					y_abc = [0.915, 0.81, 0.685, 0.49, 0.415, 0.295, 0.16]
				endelse
				end
			endcase
			abc = '('+abc+')'
		endif

		;;; make plot
		popen, pic_folder+'/th'+sc
		print_options,xsize=size_x_rb,ysize=size_y_rb
		tplot, tplotnames, trange = trange_show, title = thm_probe_color(sc, /number, /long);, var_label = labelnames
		timebar, t_onset, line = 1
		timebar_mass, 0, /databar, varname = tnames_line, line = 1
		if ~if_spec_lines then begin
			timebar_mass, [sep_p_low, sep_p_high], /databar, varname = tnames_spec[0], line = 2
			timebar_mass, [sep_e_low, sep_e_high], /databar, varname = tnames_spec[1], line = 2
		endif
		if strcmp(plot_type, 'publication') then begin
			xyouts, replicate(x_abc, n_elements(abc)), y_abc, abc, /normal

			case n_elements(abc) of
			5: y_flux = 0.28
			7: y_flux = 0.215
			endcase

			if i eq 0 and if_spec_lines then begin
				if ~strcmp(figure_plot, 'data_normal') then xyouts, 0.115, y_flux, 'Flux [/cm!u2!n'+dot_sign+'s'+dot_sign+'sr'+dot_sign+'keV]', orientation = 90, align = 0.5
			endif

			if i eq n_elements(probes)-1 then begin
				if if_spec_lines then begin
					xyouts, 0.87, y_flux, 'Energy', orientation = 90, align = 0.5, charsize = 0.6
				endif else begin
					if ~strcmp(figure_plot, 'data_normal') then xyouts, 0.917, y_flux, 'Flux [/cm!u2!n'+dot_sign+'s'+dot_sign+'sr'+dot_sign+'keV]', orientation = 90, align = 0.5
				endelse
			endif
		endif
		pclose
	endif


	if strcmp(figure_plot, 'allchannels') then begin
		;;;;;;; plot all channels to find enhancement channels
		split_vec, 'th'+sc+'_mageis_pspec_tclip'
		split_vec, 'th'+sc+'_hope_sa_pspec_tclip'
		split_vec, 'th'+sc+'_rept_pspec_tclip'
		split_vec, 'th'+sc+'_espec_combined_xj_tclip'

		;;;; plot the channels
		dtype_plot = ['mageis_pspec_tclip', 'hope_sa_pspec_tclip', 'espec_combined_xj_tclip']
		for i_type = 0, n_elements(dtype_plot)-1 do begin
			dtype = dtype_plot[i_type]
			case dtype of
			'mageis_pspec_tclip': n_max = 2
			'hope_sa_pspec_tclip': n_max = 6
			'espec_combined_xj_tclip': n_max = 9
			endcase
			for k = 0, n_max do begin
				suf_file = strcompress(string(k), /remove)
				suf_next = strcompress(string(k+1), /remove)
				if k eq 0 then suf = '' else suf = suf_file
				tnames_series = 'th'+sc+'_'+dtype+'_'+suf+'?'
				if (k eq n_max) and (~strcmp(dtype, 'espec_combined_xj')) then begin
					tnames_plot = [tnames_series, 'th'+sc+'_'+dtype+'_'+suf_next+'?']
				endif else begin
					tnames_plot = tnames_series
				endelse
				filter_suf = ''

				names = tnames(tnames_plot, n_names)
				for i_var = 0, n_names-1 do begin
					;;; filter the data for RBSP-B. Comment if do not want
					if strcmp(sc, 's') then begin
						thigh_pass_filter, names[i_var], mins_filt_lowf*60., newname = names[i_var]
						thigh_pass_filter, names[i_var], mins_filt_lowf*60., newname = names[i_var] ;; do it twice
						thigh_pass_filter, names[i_var], mins_filt_lowf*60., newname = names[i_var] ;; do it thrice
						thigh_pass_filter, names[i_var], mins_filt_lowf*60., newname = names[i_var] ;; do it four times
						filter_suf = ' high-pass filterred'
					endif
					get_data, names[i_var], t_temp, data_temp, v_this
					options, names[i_var], labels = strcompress(string(fix(v_this[0])), /remove)+'keV'
				endfor

				options, tnames_plot, ylog = 0, spec = 0;, ytitle = dtype, ysubtitle = '[/cm!u2!n'+dot_sign+'s'+dot_sign+'sr'+dot_sign+'keV]'
				tplot, ['th'+sc+'_fgs_fac_bpfilt', tnames_plot], trange = trange_show, title = 'RBSP-'+sc_rb+filter_suf
				timebar, t_onset, line = 1
				makepng, pic_folder+'/th'+sc+'_'+dtype+'_split_'+suf_file
			endfor
		endfor

		filter_suf = ''
		names = tnames('th'+sc+'_rept_pspec_tclip_?', n_names)
		for i_var = 0, n_names-1 do begin
			;;; filter the data for RBSP-B. Comment if do not want
			if strcmp(sc, 's') then begin
				thigh_pass_filter, names[i_var], mins_filt_lowf*60., newname = names[i_var]
				thigh_pass_filter, names[i_var], mins_filt_lowf*60., newname = names[i_var] ;; do it twice
				filter_suf = ' high-pass filterred'
			endif
			get_data, names[i_var], t_temp, data_temp, v_this
			options, names[i_var], labels = strcompress(string(fix(v_this[0])), /remove)+'MeV'
		endfor
		options, 'th'+sc+'_rept_pspec_tclip_?', ylog = 0, spec = 0
		tplot, 'th'+sc+'_'+['fgs_fac_bpfilt', 'rept_pspec_tclip_?'], trange = trange_show, title = 'RBSP-'+sc_rb+filter_suf
		timebar, t_onset, line = 1
		makepng, pic_folder+'/th'+sc+'_rept_pspec_tclip_split'
	endif


	;;;;;;; All picth angles
	if strcmp(figure_plot, 'allchannels_pa') then begin
		tnames_pa = tnames('*_pa_*_*', n_names) ;; all
		;tnames_pa = tnames('th?_mageis_pa_espec_*', n_names) ;; mageis electron only
		for i_name = 0, n_names-1 do begin
			if i_name gt 0 then del_data, tnames_pa[i_name-1]+'*'
			tname_pa = tnames_pa[i_name]
			get_data, tname_pa, t, flux, angles

			;;; correct data for RBSP-B, electron mageis
			if strmatch(tname_pa, 'ths_mageis_pa_espec_*') then begin
				en_str = strmid(tname_pa, 20)
				en_num = float(en_str)
				if (en_num gt 30) and (en_num lt 240) then begin
					i_bad = where((t gt time_double('2013 9 30 1 26')) and (t lt time_double('2013 9 30 1 30')), n_bad)
					if n_bad gt 0 then begin
						flux[i_bad,*] = !values.f_nan
						store_data, tname_pa, data = {x:t, y:flux, v:angles}
					endif
				endif
			endif

			split_vec, tname_pa
			tnames_this = tnames(tname_pa+'_*', n_names_this)
			for j_name = 0, n_names_this-1 do begin
				filter_suf = ''

				;;; filter the data for RBSP-B. Comment if do not want
				if strcmp(sc, 's') then begin
					thigh_pass_filter, tnames_this[j_name], mins_filt_lowf*60., newname = tnames_this[j_name]
					thigh_pass_filter, tnames_this[j_name], mins_filt_lowf*60., newname = tnames_this[j_name] ;; do it twice
					filter_suf = ' high-pass filterred'
				endif


				if size(angles, /n_dim) gt 1 then stop ;; in case angles is a matrix
				angle_str = strcompress(string(fix(angles[j_name])), /remove)
				options, tnames_this[j_name], labels = angle_str, ytitle = 'Flux '+angle_str+'!uo!n', spec = 0, ylog = 0
			endfor

			if n_names_this gt 15 then begin
				n_each = 10
				n_sep = floor(n_names_this/n_each)
				for i_plot = 0, n_sep do begin
					if i_plot eq n_sep then begin
						i_end = -1
					endif else begin
						i_end = (i_plot+1)*n_each
					endelse
					tnames_plot = tnames_this[i_plot*n_each:i_end]
					tplot, ['th'+sc+'_fgs_fac_bpfilt', tnames_plot], trange = trange_show, title = tname_pa+filter_suf+' group '+strcompress(string(i_plot))
					timebar, t_onset, line = 1
					makepng, pic_folder+'/'+tname_pa+'_grp'+strcompress(string(i_plot), /remove)
				endfor
			endif else begin
				tplot, ['th'+sc+'_fgs_fac_bpfilt', tnames_this], trange = trange_show, title = tname_pa+filter_suf
				timebar, t_onset, line = 1
				makepng, pic_folder+'/'+tname_pa
			endelse
		endfor
	endif


	;;;;;;;;;;;;;;; wave power spectra plot ;;;;;;;;;;;;;;;;;;
	if strcmp(figure_plot, 'powerspec') then begin
		;tplotnames = 'th'+sc+['_fgs_fac_bpfilt', '_fgs_fac_'+['x', 'y', 'z']+'_wv_pow_tclip', '_'+ename+'_sfac_'+['y', 'z']+'_wv_pow']
		tplotnames = 'th'+sc+['_fgs_fac_bpfilt', '_fgs_fac_'+['x', 'y', 'z']+'_wv_pow_tclip']

		;;; for publication only: options for two panels and limits
		if strcmp(plot_type, 'publication') then begin
			ylim, 'th'+sc+'_fgs_fac_bpfilt', -0.64, 0.64
			if strcmp(sc, 'r') then begin
				ylim, 'th'+sc+'_poynting_sfac', -1.9, 0.2
				options, tplotnames[1:*], no_color_scale = 1
			endif

			if i eq 0 then begin
				abc = ['a', 'b', 'c', 'd']
			endif else begin
				abc = ['e', 'f', 'g', 'h']
			endelse

			if i ne 0 then begin
				options, tplotnames, ytitle = '', ysubtitle = '', ytickname = replicate(' ', 50)
				options, tnames_spec, ytickname = replicate(' ', 50)
				tplot_options,'vtitle',''
			endif

			if i ne n_elements(probes)-1 then begin
				options, tplotnames, labels = ''
				;options, tnames_spec, no_color_scale = 1
			endif

			;;; abc labels
			case n_elements(abc) of
			4: y_abc = [0.895, 0.68, 0.48, 0.28]
			5: y_abc = [0.885, 0.743, 0.57, 0.397, 0.15]
			7: y_abc = [0.915, 0.81, 0.685, 0.49, 0.415, 0.295, 0.16]
			endcase
			abc = '('+abc+')'
		endif

		ztitle_b = 'Power [nT!u2!n/Hz]'
		ztitle_e = 'Power [mV!u2!n/m!u!n'+dot_sign+'Hz]'

		popen, pic_folder+'/th'+sc+'_powerspec'
		print_options,xsize=size_x_rb_spec, ysize=size_y_rb_spec
		tplot, tplotnames, trange = trange_show, title = thm_probe_color(sc, /number, /long)
		timebar_mass, 0, varname = tplotnames[0], /databar, line = 1
		;timebar_mass, 1/[40., 150.], varname = tplotnames[1:*], color = 1, /databar, thick = 1.5 ;; the Pi2 range
		timebar, t_onset, line = 1, color = 0
		;;; put the z title
		if i eq n_elements(probes)-1 then xyouts, 0.95, 0.44, ztitle_b, align = 0.5, orient = 90
		;;; draw the abc labels
		xyouts, x_abc, y_abc[0], abc[0], /normal
		xyouts, replicate(x_abc, n_elements(abc[1:*])), y_abc[1:*], abc[1:*], /normal, color = 255

		pclose
	endif


	;;;;;;;;;;;;;;; discussion picture plot ;;;;;;;;;;;;;;;;;;
	if strcmp(figure_plot, 'discussion') then begin
		tplotnames = 'th'+sc+'_'+['fgs_b_bpfilt', 'efw_density_bpfilt']
		abc = '('+['a', 'b']+')'
		y_abc = [0.6, 0.3]
		;;; make plot
		popen, pic_folder+'/th'+sc+'_discuss'
		print_options,xsize=size_x_discuss,ysize=size_y_discuss
		tplot, tplotnames, trange = trange_show, title = thm_probe_color(sc, /number, /long)+' Wave Mode'
		timebar, t_onset, line = 1
		timebar_mass, 0, /databar, varname = tplotnames, line = 1
		if strcmp(plot_type, 'publication') then begin
			xyouts, replicate(x_abc, n_elements(abc)), y_abc, abc, /normal
		endif
		pclose
	endif


	;;;;;;;;;;;;;;;; examine plasma pause plot ;;;;;;;;;;;;;;;;;;
	if strcmp(figure_plot, 'plasmapause') then begin
		ttrace2equator, 'th'+sc+'_state_pos_gsm_RE', newname='th'+sc+'_efoot', external_model='t89', par=kp
		tinterpol_mxn, 'th'+sc+'_efoot', 'th'+sc+'_efw_density', newname='th'+sc+'_efoot_int'
		get_data, 'th'+sc+'_efw_density', t, data_density
		get_data, 'th'+sc+'_efoot_int', t, data_efoot
		L_shell = sqrt(total(data_efoot[*,0:1]^2, 2))
		store_data, 'th'+sc+'_L', data = {x:t, y:L_shell}
		;; compute the sheeley and moldwind densities
		n_sheeley = 124*(3./L_shell)^4
		n_moldwin = 10*(6.6/L_shell)^4
		store_data, 'th'+sc+'_efw_density_compare', data = {x:t, y:[[data_density], [n_sheeley], [n_moldwin]]}
		options, 'th'+sc+'_efw_density_compare', colors = [0, 2, 6], labels = ['n!dEFW!n', 'Sheeley', 'Moldwin'], ylog = 1

		tplotnames = 'th'+sc+'_efw_density_compare'
		labelnames = 'th'+sc+'_'+['state_pos_gsm_RE_z', 'state_pos_gsm_RE_y', 'state_pos_gsm_RE_x', 'L']
		;;; make plot
		popen, pic_folder+'/th'+sc+'_plasmapause'
		print_options,xsize=6,ysize=3
		tplot, tplotnames, trange = trange_show, title = thm_probe_color(sc, /number, /long)+' Density', var_label = labelnames
		timebar, t_onset, line = 1
		timebar_mass, 0, /databar, varname = tplotnames, line = 1
		pclose
	endif

	;stop
endfor
;stop
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; GOES 13 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;figure_plot = 'data' ;; data figure, original
;;figure_plot = 'data_wpi' ;; data figure, show images with wave-particle interactions.
;;figure_plot = 'data_normal' ;; data figure, with normalized spectra
;;figure_plot = 'flux_all' ;; print the lines of all flux channels (including epead)
;;figure_plot = 'allchannels_pa' ;; all flux channels with pitch angle
;figure_plot = 'allchannels_gyro' ;; all flux channels with gyro phase
;
;size_x_goes = 5.4
;bandpass_lsp = 0 ;; set 1 to use lsp method, 0 to use smooth method
;
;;avg_1m_flux = 1
;
;abc = ['a', 'b', 'c', 'd', 'e']
;x_abc = replicate(0.21, n_elements(abc))
;
;case n_elements(abc) of
;4: y_abc = [0.9, 0.7, 0.465, 0.25]
;5: y_abc = [0.9, 0.74, 0.58, 0.39, 0.24]
;endcase
;
;if strcmp(figure_plot, 'allchannels_pa') then pa = 1 else pa = 0 ;; select whether to load pitch angel data
;if strcmp(figure_plot, 'allchannels_gyro') then gyro = 1 else gyro = 0 ;; select whether to load pitch angel data
;
;sc = '13'
;goes_load_plasma, trange = trange_load, probes = sc, pa = pa, gyro = gyro, avg_1m = avg_1m_flux
;goes_load_data, trange = trange_load, datatype = 'fgm', probes = sc ;; load fgm after flux because flux may load 1 m data to overwrite
;
;;;;; transform position also to gsm
;cotrans, 'g'+sc+'_pos_gei', 'g'+sc+'_pos_gse', /gei2gse
;cotrans, 'g'+sc+'_pos_gse', 'g'+sc+'_pos_gsm', /gse2gsm
;
;;;;; rotate the FGM data from ENP coordinates to GEI coordinates
;enp_matrix_make, 'g'+sc+'_pos_gei'
;tvector_rotate, 'g'+sc+'_pos_gei_enp_mat', 'g'+sc+'_H_enp_1', /invert
;cotrans, 'g'+sc+'_H_enp_1_rot', 'g'+sc+'_H_gse', /gei2gse
;cotrans, 'g'+sc+'_H_gse', 'g'+sc+'_H_gsm', /gse2gsm
;cotrans2fac, 'g'+sc+'_H_gsm', b_variable = 'g'+sc+'_H_gsm', pos = 'g'+sc+'_pos_gsm', smooth_time = 6.*60., newname = 'g'+sc+'_H_fac'
;
;;;;; create partial pressure
;compute_moments, 'g'+sc+['_protflux', '_elecflux'], intypes = 'flux', inunits_energy = 'keV', inunits_flux = 'keV', particle_types = ['p', 'e'], /combine, newname_combine = 'all' ;; there is no corrected data for this event, so do not use this.
;;; plot density, temperature
;tinterpol_mxn, 'g'+sc+'_all_Pth', 'g'+sc+'_H_gsm', newname='g'+sc+'_all_Pth_int'
;get_data, 'g'+sc+'_all_Pth_int', t, pdata
;get_data, 'g'+sc+'_H_gsm', t, bdata
;Pbx = nTesla2_to_nPa*bdata[*,0]^2
;Pb = nTesla2_to_nPa*total(bdata^2, 2)
;store_data, 'g'+sc+'_Pall', data = {x:t, y:[[Pb], [pdata]]}
;store_data, 'g'+sc+'_Pb', data = {x:t, y:Pb}
;
;;;;;;;;;;;;; make flux ratios ;;;;;;;;
;tvnames_ratios = 'g'+sc+['_protflux', '_elecflux']
;for k = 0, n_elements(tvnames_ratios)-1 do begin
;	time_clip, tvnames_ratios[k], time_double(trange_show[0]), time_double(trange_show[0])+3*60. ;; 3 min average
;	get_data, tvnames_ratios[k]+'_tclip', tst, yst, vst
;	yp_start = transpose(mean(yst, dim = 1, /nan)) 
;
;	get_data, tvnames_ratios[k], tp, yp, vp
;	yp_divide = rebin(yp_start, n_elements(tp), n_elements(yp_start))
;	store_data, tvnames_ratios[k]+'_normalize', data={x:tp, y:yp/yp_divide, v:vp}
;	options, tvnames_ratios[k]+'_normalize', ylog = 1, spec = 1, zlog = 1, ysubtitle = '!c[keV]', ztitle = ''
;endfor
;
;;;;;;;;;;;;; make selected flux channels ;;;;
;;;; proton
;get_data, 'g'+sc+'_protflux', tp, datap, vp
;rate_p = [1,9,70]
;subtract_p = 4e4
;rate_p_str = ' '+strcompress(string(rate_p),/remove)+'x'
;rate_p_str[0] = ''
;multiply_p = rebin(transpose(rate_p), n_elements(tp), n_elements(rate_p))
;store_data, 'g'+sc+'_protflux_less', data = {x:tp, y:datap[*,0:2]*multiply_p-subtract_p, v:vp[*,0:2]}
;options, 'g'+sc+'_protflux_less', ytitle = 'p!u+!n Flux!uUncor', ysubtitle = '!c', ylog = 0, labels = '!d'+transpose(strcompress(string(fix(vp[0,0:2])), /remove))+' keV', labflag = -1, thick = l_thick
;ylim, 'g'+sc+'_protflux_less', 4.01e4-subtract_p, 8.9e4-subtract_p
;;;; electron
;get_data, 'g'+sc+'_elecflux', te, datae, ve
;rate_e = [1,1,2.6,4.9,17,0.54e2]
;subtract_e = 0.8e5
;rate_e_str = ' '+strcompress(string(rate_e),/remove)+'x'
;rate_e_str[0:1] = ['','']
;multiply_e = rebin(transpose(rate_e), n_elements(te), n_elements(rate_e))
;store_data, 'g'+sc+'_elecflux_less', data = {x:te, y:datae[*,0:4]*multiply_e-subtract_e, v:ve[*,0:4]} ;; do not plot 600 keV because it has some confusing peaks.
;options, 'g'+sc+'_elecflux_less', ytitle = 'e!u'+minus_sign+'!n Flux!uUncor', ysubtitle = '!c', ylog = 0, labels = '!d'+strcompress(string(fix(ve[0,0:4])), /remove)+' keV', labflag = -1, thick = l_thick
;ylim, 'g'+sc+'_elecflux_less', 8.5e4-subtract_e, 1.59e5-subtract_e
;
;
;;;;;;;;;;; band pass filter fields ;;;;;;;;;;
;tbandpass_filter, 'g'+sc+'_H_gsm', mins_filt_lowf*60., mins_filt_highf*60., lsp = bandpass_lsp
;tbandpass_filter, 'g'+sc+'_H_fac', mins_filt_lowf*60., mins_filt_highf*60., lsp = bandpass_lsp
;tbandpass_filter, 'g'+sc+'_Pall', mins_filt_lowf*60., mins_filt_highf*60., lsp = bandpass_lsp
;tbandpass_filter, 'g'+sc+'_Pb', mins_filt_lowf*60., mins_filt_highf*60., lsp = bandpass_lsp
;;;; corret quiet time value when lsp is set
;if keyword_set(bandpass_lsp) then begin
;	;;; H
;	get_data, 'g'+sc+'_H_fac_bpfilt', t_fgs, data_fgs, limit = limit, dlimit = dlimit
;	time_clip, 'g'+sc+'_H_fac_bpfilt', date+'1 5', date+'1 15', newname = 'g'+sc+'_H_fac_bpfilt_pre'
;	get_data, 'g'+sc+'_H_fac_bpfilt_pre', t_pre, data_pre
;	bb_pre = mean(data_pre[*,2], /nan)
;	data_fgs[*,2] = data_fgs[*,2]-bb_pre
;	store_data, 'g'+sc+'_H_fac_bpfilt', data = {x:t_fgs, y:data_fgs}, limit = limit, dlimit = dlimit
;	;;; P
;	get_data, 'g'+sc+'_Pall_bpfilt', t_fgs, data_fgs, limit = limit, dlimit = dlimit
;	time_clip, 'g'+sc+'_Pall_bpfilt', date+'1 5', date+'1 15', newname = 'g'+sc+'_Pall_bpfilt_pre'
;	get_data, 'g'+sc+'_Pall_bpfilt_pre', t_pre, data_pre
;	value0 = mean(data_pre[*,0], /nan)
;	value1 = mean(data_pre[*,1], /nan)
;	data_fgs[*,0] = data_fgs[*,0]-value0
;	data_fgs[*,1] = data_fgs[*,1]-value1
;	store_data, 'g'+sc+'_Pall_bpfilt', data = {x:t_fgs, y:data_fgs}, limit = limit, dlimit = dlimit
;endif
;
;
;;;;; write labels
;options, 'g'+sc+'_H_gsm*', ytitle = 'B', ysubtitle = '!c[nT]', colors=[2,4,6], labels=['B!dx','B!dy','B!dz'], labflag = 1, thick = l_thick
;options, 'g'+sc+'_H_fac*', colors=[2,4,6], labels=["B!dr'",'B!d'+phi_letter+"'",'B!db'], ytitle = 'B!dFAC', ysubtitle = '!c[nT]', labflag = 1, thick = l_thick
;options, 'g'+sc+'_protflux', ytitle = 'p!u+!n Energy', ysubtitle = '!c[keV]', spec=1, ztitle = ''
;options, 'g'+sc+'_elecflux', ytitle = 'e!u'+minus_sign+'!n Energy', ysubtitle = '!c[keV]', spec=1, ztitle = ''
;options, 'g'+sc+'_protflux_normalize', ytitle = 'p!u+!n Ratio'
;options, 'g'+sc+'_elecflux_normalize', ytitle = 'e!u'+minus_sign+'!n Ratio'
;options, 'g'+sc+'_Pall*', ytitle = 'P', ysubtitle = '!c[nPa]', labels = ['P!dB', 'P!s!dth!r!upartial'], labflag = 1, colors = [1,4], thick = l_thick
;options, 'g'+sc+'_Pb*', ytitle = 'P!db', ysubtitle = '!c[nPa]', thick = l_thick
;split_vec, 'g'+sc+'_H_gsm'
;split_vec, 'g'+sc+'_protflux'
;split_vec, 'g'+sc+'_elecflux'
;options, 'g'+sc+'_H_gsm_bpfilt*', ytitle = 'B!ufilt'
;options, 'g'+sc+'_H_fac_bpfilt*', ytitle = 'B!s!dFAC!r!ufilt'
;options, 'g'+sc+'_Pall_bpfilt*', ytitle = 'P!ufilt'
;options, 'g'+sc+'_Pb_bpfilt*', ytitle = 'P!s!db!r!ufilt', ysubtitle = '!c[nPa]', thick = l_thick
;options, 'g'+sc+'_????flux_?', spec = 0, ylog = 0
;
;ylim, 'g'+sc+'_H_gsm', -45., 85.
;ylim, 'g'+sc+'_Pall_bpfilt', -0.055, 0.059
;ylim, 'g'+sc+'_Pb_bpfilt', -0.055, 0.059
;
;;tplotnames_line = 'g'+sc+'_'+['H_gsm', 'H_fac_bpfilt', 'H_gsm_x']
;tplotnames_line = 'g'+sc+'_'+['H_gsm', 'H_fac_bpfilt', 'Pb_bpfilt']
;
;case figure_plot of
;'data': begin
;	tplotnames_spec = 'g'+sc+'_'+['protflux', 'elecflux', 'protflux_less', 'elecflux_less'] ;tplotnames_spec = 'g'+sc+'_'+['protflux', 'elecflux'] (original)
;	size_y_goes = 8.2
;	end
;'data_normal': begin
;	tplotnames_spec = 'g'+sc+'_'+['protflux', 'elecflux']+'_normalize'
;	size_y_goes = 8.2
;	end
;'data_wpi': begin
;	tplotnames_spec = 'g'+sc+'_'+['protflux_less', 'elecflux_less']
;	size_y_goes = 6.2
;	end
;else: tplotnames_spec = '' ;; undefined
;end
;;tplotnames_spec = 'g'+sc+'_protflux_'+['0', '1', '2']
;;tplotnames_spec = 'g'+sc+'_elecflux_'+['0', '1', '2', '3', '4', '5']
;
;
;if strcmp(plot_type, 'long_injection') then begin
;	tplotnames_line = 'g'+sc+'_H_gsm'
;	options, 'g'+sc+'_protflux', ylog = 1, spec = 0
;	options, 'g'+sc+'_elecflux', ylog = 1, spec = 0
;endif else begin
;	;;; comment these for line plot
;	ylim, 'g'+sc+'_protflux', 80., 5e5
;	ylim, 'g'+sc+'_elecflux', 38., 4500
;	ylim, 'g'+sc+'_protflux_normalize', 80., 5e5
;	ylim, 'g'+sc+'_elecflux_normalize', 38., 4500
;	;;;;; use line plot
;	;options, 'g'+sc+'_protflux', ylog = 1, spec = 0
;	;options, 'g'+sc+'_elecflux', ylog = 1, spec = 0
;endelse
;
;
;;;;;;;;; make plot ;;;;;;;;;
;
;case 1 of
;
;
;;;;;;;;; data plots ('data' is for publication plot)
;strmatch(figure_plot, 'data*'): begin
;	tplotnames = [tplotnames_line, tplotnames_spec]
;
;	popen, pic_folder+'/g'+sc
;	print_options, xsize=size_x_goes, ysize=size_y_goes
;	tplot, tplotnames, trange = trange_show, title = thm_probe_color(sc, /number, /long); , var_label = labelnames
;	timebar, t_onset, line = 1
;	timebar_mass, 0, /databar, varname = tplotnames_line, line = 1
;
;	case figure_plot of
;	'data_wpi': begin
;		y_abc[3] = 0.33
;		xyouts, x_abc, y_abc, '('+abc+')', /normal
;		xyouts, 0.085, 0.27, '[/cm!u2!n'+dot_sign+'s'+dot_sign+'sr'+dot_sign+'keV]', align = 0.5, orient = 90.
;		end
;	else: begin
;		xyouts, 0.93, 0.3, 'Flux!uUncor!n [/cm!u2!n'+dot_sign+'s'+dot_sign+'sr'+dot_sign+'keV]', align = 0.5, orient = 90.
;		xyouts, x_abc[0:-3], y_abc[0:-3], '('+abc[0:-3]+')', /normal
;		xyouts, x_abc[-2:-1], y_abc[-2:-1], '('+abc[-2:-1]+')', /normal, color = 255
;		end
;	endcase
;
;	pclose
;	end
;
;
;;;;;;;;; All ion channels
;strmatch(figure_plot, 'flux_all'): begin
;	trange_plot = trange_show
;
;	split_vec, 'g'+sc+'_protflux', names_out = names_protfluxes
;	split_vec, 'g'+sc+'_elecflux', names_out = names_elecfluxes
;	names_fluxes = [names_protfluxes, names_elecfluxes]
;	options, names_fluxes, spec = 0
;	for j = 0, n_elements(names_fluxes)-1 do begin
;		;; find the y limit to use
;		time_clip, names_fluxes[j], trange_plot[0], trange_plot[1]
;		get_data, names_fluxes[j]+'_tclip', t_nouse, data_this
;		data_min = min(data_this, /nan)
;		data_max = max(data_this, /nan)
;		data_span = data_max-data_min
;		ylim, names_fluxes[j], data_min-0.07*data_span, data_max+0.07*data_span
;	endfor
;
;	tplot, names_protfluxes, trange = trange_plot
;	timebar, t_onset, line = 1
;	makepng, pic_folder+'/'+'g'+sc+'_allprotflux'
;	tplot, names_elecfluxes, trange = trange_plot
;	timebar, t_onset, line = 1
;	makepng, pic_folder+'/'+'g'+sc+'_allelecflux'
;	end
;
;
;;;;;;;;; Splitted pitch angles
;strmatch(figure_plot, 'allchannels_pa'): begin
;	tnames_pa = tnames('g'+sc+'_mag?d_*_pa*', n_pa_name)
;	trange_plot = time_double(trange_show)+[8*60., 0]
;	pic_folder_pa = pic_folder+'/goes_pa/'
;	for i = 0, n_pa_name-1 do begin
;		get_data, tnames_pa[i], t_nouse, pa_nouse, pa_angles
;		split_vec, tnames_pa[i], names_out = pa_names_this
;		for j = 0, n_elements(pa_names_this)-1 do begin
;			options, pa_names_this[j], labels = strcompress(string(fix(pa_angles[j])), /remove), labflag = 1, spec = 0
;			;; find the y limit to use
;			time_clip, pa_names_this[j], trange_plot[0], trange_plot[1]
;			get_data, pa_names_this[j]+'_tclip', t_nouse, data_this
;			data_min = min(data_this, /nan)
;			data_max = max(data_this, /nan)
;			data_span = data_max-data_min
;			ylim, pa_names_this[j], data_min-0.07*data_span, data_max+0.07*data_span
;		endfor
;		tplot, pa_names_this, trange = trange_plot
;		timebar, t_onset, line = 1
;		if keyword_set(avg_1m_flux) then pic_folder_pa_use = pic_folder_pa+'/1min' else pic_folder_pa_use = pic_folder_pa+'/uncor'
;		makepng, pic_folder_pa_use+'/'+tnames_pa[i]+'_allpa'
;	endfor
;	end
;
;strmatch(figure_plot, 'allchannels_gyro'): begin
;	tnames_gyro = tnames('g'+sc+'_mag?d_*_gyro*', n_gyro_name)
;	trange_plot = time_double(trange_show)+[8*60., 0]
;	pic_folder_gyro = pic_folder+'/goes_gyro/'
;	for i = 0, n_gyro_name-1 do begin
;		get_data, tnames_gyro[i], t_nouse, gyro_nouse, gyro_angles
;		split_vec, tnames_gyro[i], names_out = gyro_names_this
;		for j = 0, n_elements(gyro_names_this)-1 do begin
;			options, gyro_names_this[j], labels = strcompress(string(fix(gyro_angles[j])), /remove), labflag = 1, spec = 0
;			;; find the y limit to use
;			time_clip, gyro_names_this[j], trange_plot[0], trange_plot[1]
;			get_data, gyro_names_this[j]+'_tclip', t_nouse, data_this
;			data_min = min(data_this, /nan)
;			data_max = max(data_this, /nan)
;			data_span = data_max-data_min
;			ylim, gyro_names_this[j], data_min-0.07*data_span, data_max+0.07*data_span
;		endfor
;		tplot, gyro_names_this, trange = trange_plot
;		timebar, t_onset, line = 1
;		if keyword_set(avg_1m_flux) then pic_folder_gyro_use = pic_folder_gyro+'/1min' else pic_folder_gyro_use = pic_folder_gyro+'/uncor'
;		makepng, pic_folder_gyro_use+'/'+tnames_gyro[i]+'_allgyro'
;	endfor
;	end
;
;else: message, 'figure_plot selected wrong!'
;endcase
;
;;stop
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



;;;;;;;;;;;;;;;;;;;;;;;;;;;; GBO ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; select the plot type
;figure_plot = 'gbo_only'
;;figure_plot = 'gbo_sc_sep' ;; if choose this, all the previous must be uncommented.
;;figure_plot = 'gbo_sc_tog' ;; if choose this, all the previous must be uncommented. GBO and SC components together in the same panel.
;
;;;; plot options
;sort_suf = '_slat'
;
;if strcmp(figure_plot, 'gbo_only') then begin
;	xsize = 7
;	ysize = 9
;	x_abc = 0.18
;	y_abc = [0.935, 0.783, 0.633, 0.483, 0.333, 0.183]
;	abc = '('+['b', 'c', 'd', 'e', 'f', 'g']+')'
;endif
;
;if strcmp(figure_plot, 'gbo_sc_sep') then begin
;	xsize = 7
;	ysize = 4.5
;endif
;
;if strcmp(figure_plot, 'gbo_sc_tog') then begin
;	xsize = 7
;	;ysize = 5.5
;	ysize = 7.5
;endif
;
;;;;; footprints of satellites (lat, long, in degrees)
;sc_foots = [[60.643441, -73.026010], $ ;; tha
;      [61.644975, -51.019944], $ ;; thd
;      [61.295331, -62.145895], $ ;; the
;      [55.139415, -99.229428], $ ;; thr
;      [46.006514, -69.522629], $ ;; ths
;      [56.605112, -79.751387]] ;; g13
;
;if strmatch(figure_plot, 'gbo_sc*') then begin
;	sc_mags = ['tha_fgs_fac_bpfilt', '', '', 'thr_fgs_fac_bpfilt', 'ths_fgs_fac_bpfilt', 'g13_H_fac_bpfilt']
;	sc_mags_original = ['tha_fgs_fac', '', '', 'thr_fgs_fac', 'ths_fgs_fac', 'g13_H_fac']
;	sc_elecs = ['tha_efs_dot0_sfac_bpfilt', '', '', 'thr_'+ename+'_sfac_bpfilt', 'ths_'+ename+'_sfac_bpfilt', '']
;endif
;
;;gbos = ['atha', 'ccnv', 'fcc', 'frn', 'fsmi', 'gua', 'hon', 'larg', 'mea', 'pine', 'tpas', 'tuc', 'ukia', 'vic', 'whit', 'ykc', 'bou', 'bfe', 'dob', 'don', 'drby', 'frd', 'gbay', 'glyn', 'jck', 'kapu', 'kar', 'lrv', 'new', 'ott', 'pina', 'roe', 'rvk', 'satx', 'sjg', 'sol', 'stj', 'swno', 'tpas', 'tuc', 'wrth', 'dmh', 'kuv', 'lyr', 'svs', 'thl', 'umq', 'upn', 'abk', 'amk', 'amd', 'and','bjn', 'fcc', 'ghb', 'gill', 'iqa', 'kuuj', 'naq', 'nor', 'nrsq', 'sco', 'skt', 'snkq', 'sor', 'tro'] ;; removed fhb.
;thm_gmag_stations, gbos, locations_ll, magnetic = locations_ll_mag, midnight = midnight_gbos ;; locations_geo: latitude, longitude. This will make gbos into upper case
;;; return to lower case
;gbos = strlowcase(gbos)
;
;;;; remove station "fhb" (low resolution data)
;;stations_1min = ['fhb', 'ghb']
;;i_highres = where(~strcmp(gbos, 'fhb') and ~strcmp(gbos, 'ghb') and ~strcmp(gbos, 'nrsq') and ~strcmp(gbos, 'naq') and ~strcmp(gbos, 'skt'), n_highres)
;;if n_highres gt 0 then begin
;;	gbos = gbos(i_highres)
;;	locations_ll = locations_ll[*, i_highres]
;;	locations_ll_mag = locations_ll_mag[*, i_highres]
;;	midnight_gbos = midnight_gbos(i_highres)
;;endif
;
;mlat_gbos = locations_ll_mag[0,*]
;
;t_onset_str = time_string(t_onset)
;midnight_str_full = '2000/01/01/'+strcompress(midnight_gbos, /remove)
;t_pos_r_str_full = '2000/01/01/'+strmid(t_onset_str, 11)
;mlt_gbos = (time_double(t_pos_r_str_full)-time_double(midnight_str_full))/3600. ;; t_pos_r is the time of first loaded quantity
;i_negative = where(mlt_gbos lt 0, n_negative)
;if n_negative gt 0 then mlt_gbos[i_negative] = mlt_gbos[i_negative]+24.
;
;
;;;; select stations that are closest to the sc footprints
;geo_scfoots = latlong2geo(sc_foots)
;geo_gbos = latlong2geo(locations_ll)
;
;i_gbos_use = intarr(n_elements(sc_foots[0,*]))
;for i = 0, n_elements(sc_foots[0,*])-1 do begin
;	distance = sqrt(total((geo_gbos-rebin(geo_scfoots[*,i], n_elements(geo_scfoots[*,i]), n_elements(gbos)))^2, 1))
;	no_use = min(distance, i_close, /nan)
;	i_gbos_use[i] = i_close
;endfor
;
;gbos_use = gbos[i_gbos_use]
;for i = 0, n_elements(i_gbos_use)-1 do begin
;	i_use = i_gbos_use[i]
;	stat = gbos[i_use]
;	;;; determine the magnetic latitude and MLT
;	;; MLat
;	mlat_stat = mlat_gbos[i_use]
;	;; MLT
;	mlt_stat = mlt_gbos[i_use]
;	;;; load data
;	thm_load_pi2, site = stat, trange = trange_load, secs_highpass = mins_filt_lowf*60., secs_lowpass = mins_filt_highf*60.
;
;	;;;;; special treatment for certain stations
;
;	;;; for KUUJ, remove bad point
;	if strcmp(stat, 'kuuj') then begin
;		get_data, 'thg_mag_'+stat+'_pi2', t, data
;		i_bad = where((t gt time_double('2013 9 30 1 8 42')) and (t lt time_double('2013 9 30 1 10 28')), n_bad)
;		if n_bad gt 0 then begin
;			data[i_bad,*] = !values.f_nan
;			store_data, 'thg_mag_'+stat+'_pi2', data = {x:t, y:data}
;		endif
;	endif
;
;	;;; for DRBY, also plot RBSP-B's compressional magnetic field.
;	if strcmp(stat, 'drby') then begin
;		get_data, 'thg_mag_'+stat+'_pi2', t, data
;		tinterpol_mxn, 'ths_fgs_fac_bpfilt', 'thg_mag_'+stat+'_pi2', newname='ths_fgs_fac_bpfilt_int'
;		get_data, 'ths_fgs_fac_bpfilt_int', t, bfac
;		store_data, 'thg_mag_'+stat+'_pi2', data = {x:t, y:[[data], [bfac[*,2]]]}
;		options, 'thg_mag_'+stat+'_pi2', ytitle = 'DRBY & RB-B!c', ysubtitle = '[nT]', colors = [2,4,6,1], labels = ['H!uDRBY', 'D!uDRBY', 'Z!uDRBY', 'B!s!db!r!uRB-B'], thick = l_thick
;	endif
;
;
;	;options, 'thg_mag_'+stat+'_pi2', ytitle = strupcase(stat)+'!cmlat='+strmid(strcompress(string(mlat_stat), /remove), 0, 4)+'!cMLT='+strmid(strcompress(string(mlt_stat), /remove), 0, 4), ysubtitle = '[nT]', thick = l_thick ;; more information
;
;	if ~strcmp(stat, 'drby') then options, 'thg_mag_'+stat+'_pi2', ytitle = strupcase(stat)+'!c', ysubtitle = '[nT]', thick = l_thick ;; less information, uncomment this if the DRBY part is uncommented.
;	;options, 'thg_mag_'+stat+'_pi2', ytitle = strupcase(stat)+'!c', ysubtitle = '[nT]', thick = l_thick ;; less information
;endfor
;
;tvnames_plot_gbo = 'thg_mag_'+gbos_use+'_pi2'
;tvnames_plot_gbo_original = 'thg_mag_'+gbos_use
;;;; sort the tplot names
;case sort_suf of
;'_slat': i_sort = reverse(sort(mlat_gbos[i_gbos_use]))
;'_smlt': i_sort = sort(mlt_gbos[i_gbos_use])
;endcase
;tvnames_plot_gbo = tvnames_plot_gbo[i_sort]
;tvnames_plot_gbo_original = tvnames_plot_gbo_original[i_sort]
;;
;if strcmp(figure_plot, 'gbo_only') then begin
;	;;;; plot
;	popen, pic_folder+'/gbo'+sort_suf
;	print_options,xsize=xsize, ysize=ysize
;	tplot, tvnames_plot_gbo, trange = trange_show, title = 'GBO Magnetic Field Pi2'
;	timebar, t_onset, line = 1
;	timebar, 0, line = 2, /databar, varname = tvnames_plot_gbo
;	xyouts, replicate(x_abc, n_elements(abc)), y_abc, abc, /normal, charsize = 1.1
;	pclose
;endif
;
;if strmatch(figure_plot, 'gbo_sc*') then begin
;	sc_mags = sc_mags[i_sort]
;	sc_mags_original = sc_mags_original[i_sort]
;	sc_elecs = sc_elecs[i_sort]
;	for i = 0, n_elements(i_sort)-1 do begin
;		if ~strcmp(sc_mags[i], '') then begin
;			thx = strmid(sc_mags[i], 0, 3)
;			if strcmp(thx, 'thr') then thx = 'rbspa'
;			if strcmp(thx, 'ths') then thx = 'rbspb'
;			stat = strmid(tvnames_plot_gbo[i], 8, 4)
;
;			;;;; separate panels
;			if strcmp(figure_plot, 'gbo_sc_sep') then begin
;				if ~strcmp(sc_elecs[i], '') then begin
;					tvnames_plot_this = [tvnames_plot_gbo[i], sc_mags[i], sc_elecs[i]]
;				endif else begin
;					tvnames_plot_this = [tvnames_plot_gbo[i], sc_mags[i]]
;				endelse
;			endif
;
;			;;;; same panels 
;			if strcmp(figure_plot, 'gbo_sc_tog') then begin
;				options, tvnames_plot_gbo_original[i]+'*', colors = [0, 3, 1]
;				split_vec, tvnames_plot_gbo[i]
;				split_vec, tvnames_plot_gbo_original[i]
;				split_vec, sc_mags[i]
;				split_vec, sc_mags_original[i]
;				case 1 of 
;				strcmp(thx, 'tha'): begin
;					store_data, 'combine1', data = [tvnames_plot_gbo[i]+'_x', sc_mags[i]+'_z']
;					store_data, 'combine2', data = [tvnames_plot_gbo[i]+'_z', sc_mags[i]+'_z']
;					store_data, 'combine3', data = [tvnames_plot_gbo[i]+'_y', sc_mags[i]+'_x']
;					store_data, 'combine4', data = [tvnames_plot_gbo[i]+'_x', sc_mags[i]+'_y']
;					tvnames_plot_this = [tvnames_plot_gbo_original[i], sc_mags_original[i], 'combine3', 'combine4', 'combine1', 'combine2']
;					end
;				strcmp(thx, 'g13'): begin
;					store_data, 'combine1', data = [tvnames_plot_gbo[i]+'_x', sc_mags[i]+'_x']
;					store_data, 'combine2', data = [tvnames_plot_gbo[i]+'_z', sc_mags[i]+'_y']
;					store_data, 'combine3', data = [tvnames_plot_gbo[i]+'_y', sc_mags[i]+'_x']
;					store_data, 'combine4', data = [tvnames_plot_gbo[i]+'_x', sc_mags[i]+'_y']
;					store_data, 'g13xy', data = [sc_mags_original[i]+'_x', sc_mags_original[i]+'_y']
;					options, 'g13xy', ytitle = 'B [nT]', ysubtitle = '', labflag = 1
;					options, sc_mags_original[i]+'_z', ytitle = 'B [nT]', ysubtitle = ''
;					tvnames_plot_this = [tvnames_plot_gbo_original[i], 'g13xy', sc_mags_original[i]+'_z', 'combine3', 'combine4', 'combine1', 'combine2']
;					end
;				strcmp(thx, 'rbspa'): begin
;					store_data, 'combine1', data = [tvnames_plot_gbo[i]+'_y', sc_mags[i]+'_z']
;					store_data, 'combine3', data = [tvnames_plot_gbo[i]+'_y', sc_mags[i]+'_x']
;					store_data, 'combine4', data = [tvnames_plot_gbo[i]+'_x', sc_mags[i]+'_y']
;					tvnames_plot_this = [tvnames_plot_gbo_original[i], sc_mags_original[i]+'_z', 'combine3', 'combine4', 'combine1']
;					end
;				strcmp(thx, 'rbspb'): begin
;					store_data, 'combine1', data = [tvnames_plot_gbo[i]+'_x', sc_mags[i]+'_z']
;					store_data, 'combine2', data = [tvnames_plot_gbo[i]+'_y', sc_mags[i]+'_z']
;					store_data, 'combine3', data = [tvnames_plot_gbo[i]+'_y', sc_mags[i]+'_x']
;					store_data, 'combine4', data = [tvnames_plot_gbo[i]+'_x', sc_mags[i]+'_y']
;					tvnames_plot_this = [tvnames_plot_gbo_original[i], sc_mags_original[i]+'_z', 'combine3', 'combine4', 'combine1', 'combine2']
;					end
;				endcase
;
;				options, 'combine*', ytitle = 'B [nT]', ysubtitle = '', labflag = 1
;				options, tvnames_plot_gbo_original[i]+'*_x', labels = stat+'_H'
;				options, tvnames_plot_gbo_original[i]+'*_y', labels = stat+'_D'
;				options, tvnames_plot_gbo_original[i]+'*_z', labels = stat+'_Z'
;				options, sc_mags_original[i]+'*_x', labels = thx+"_r'", colors = 2
;				options, sc_mags_original[i]+'*_y', labels = thx+'_'+phi_letter+"'", colors = 4
;				options, sc_mags_original[i]+'*_z', labels = thx+'_b', colors = 6
;			endif
;			;;; plot and save image
;			popen, pic_folder+'/gbo_'+thx+strmid(figure_plot, 6, 4)
;			print_options,xsize=xsize, ysize=ysize
;			tplot, tvnames_plot_this, trange = trange_show, title = ''
;			timebar, t_onset, line = 1
;			timebar, 0, line = 2, /databar, varname = tvnames_plot_this
;			pclose
;		endif
;	endfor
;endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;




;;;;;;;;;;;;;;;;;;;;;;;;;;; Cartoon for 1-1 corelation ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;trange_show  = date+['1 14', '1 30 59']
;
;rbsp_load, trange=trange_load, probe=['r', 's'], datatype = 'fgs', /tclip, rbsp_folder = rbsp_folder ;; both pos and fgs data are load
;load_bin_data, trange=trange_load, probe='e', datatype = 'fgs', /tclip, datafolder = fgs_folder
;
;tbandpass_filter, 'thr_fgs_gsm_tclip', mins_filt_lowf*60., mins_filt_highf*60., /lsp
;tbandpass_filter, 'ths_fgs_gsm_tclip', mins_filt_lowf*60., mins_filt_highf*60., /lsp
;split_vec, 'thr_fgs_gsm_tclip_bpfilt'
;split_vec, 'ths_fgs_gsm_tclip_bpfilt'
;split_vec, 'the_fgs_gsm_tclip'
;
;options, 'the_fgs_gsm_tclip_z', ytitle = 'P4  B!dz', ysubtitle = '!c[nT]', colors = 6, thick = l_thick
;options, 'thr_fgs_gsm_tclip_bpfilt_z', ytitle = 'RB-A  B!s!dz!r!ufilt', ysubtitle = '!c[nT]', colors = 6, thick = l_thick
;options, 'ths_fgs_gsm_tclip_bpfilt_z', ytitle = 'RB-B  B!s!dz!r!ufilt', ysubtitle = '!c[nT]', colors = 6, thick = l_thick
;ylim, 'thr_fgs_gsm_tclip_bpfilt_z', -0.3, 0.5
;
;tplot_names = ['the_fgs_gsm_tclip_z', 'thr_fgs_gsm_tclip_bpfilt_z']
;;tplot_names = ['the_fgs_gsm_tclip_z', 'ths_fgs_gsm_tclip_bpfilt_z', 'thr_fgs_gsm_tclip_bpfilt_z']
;
;case n_elements(tplot_names) of
;	2: begin
;		ysize = 3.5
;		title = 'Correlation between P3 and RB-A'
;		end
;	3: begin
;		ysize = 5.5
;		title = 'Correlations'
;		end
;endcase
;		
;
;popen, pic_folder+'/correlation'
;print_options,xsize=xsize, ysize=ysize
;tplot, tplot_names, trange = trange_show, title = title
;
;;;;;; add pointers
;case n_elements(tplot_names) of
;	2: begin
;		;; THE
;		add_pointer, 0.28, 0.64, '0', dir = 'southeast'
;		add_pointer, 0.434, 0.77, '1', dir = 'southeast'
;		add_pointer, 0.475, 0.87, '2', dir = 'southeast'
;		add_pointer, 0.5, 0.84, '3', dir = 'southwest'
;		add_pointer, 0.57, 0.79, '4', dir = 'south'
;		add_pointer, 0.64, 0.85, '5', dir = 'southeast'
;		;; RB-A
;		add_pointer, 0.315, 0.35, '0', dir = 'southeast'
;		add_pointer, 0.49, 0.45, '1', dir = 'southeast'
;		add_pointer, 0.545, 0.4, '2', dir = 'south'
;		add_pointer, 0.58, 0.39, '3', dir = 'southwest'
;		add_pointer, 0.67, 0.5, '4', dir = 'southwest'
;		add_pointer, 0.75, 0.39, '5', dir = 'southeast'
;		end
;	3: begin
;		;; THE
;		offset_top = 0.1
;		add_pointer, 0.28, 0.64+offset_top, '0', dir = 'southeast'
;		add_pointer, 0.434, 0.77+offset_top, '1', dir = 'southeast'
;		add_pointer, 0.475, 0.87+offset_top, '2', dir = 'southeast'
;		add_pointer, 0.5, 0.84+offset_top, '3', dir = 'southwest'
;		add_pointer, 0.57, 0.79+offset_top, '4', dir = 'south'
;		add_pointer, 0.64, 0.85+offset_top, '5', dir = 'southeast'
;		;; RB-A
;		offset_bottom = -0.1
;		add_pointer, 0.315, 0.35+offset_bottom, '0', dir = 'southeast'
;		add_pointer, 0.49, 0.45+offset_bottom, '1', dir = 'southeast'
;		add_pointer, 0.545, 0.4+offset_bottom, '2', dir = 'south'
;		add_pointer, 0.58, 0.39+offset_bottom, '3', dir = 'southwest'
;		add_pointer, 0.67, 0.5+offset_bottom, '4', dir = 'southwest'
;		add_pointer, 0.75, 0.39+offset_bottom, '5', dir = 'southeast'
;		end
;endcase
;
;xyouts, [0.22, 0.22], [0.48, 0.85], ['(b)', '(a)'], /normal, align = 0.5
;pclose
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

stop
end
