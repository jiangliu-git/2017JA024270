pro main_thm_plot
;;; plot themis data to check the quiet time flow
computer = 'I:'
@folders
thm_init
del_data, '*'

plot_type = 'publication' ;; for publication
;plot_type = 'view'

;;;;;;;; Define time ranges and and onset time
date = '2013 9 30 '
trange_show  = ['13 9 29 15 30', date+'1 39 30']
rtrange = '13 9 29 '+['15', '15 50']
t_onset = date+'1 20 16' ;; change by clicking
;;; original values
;mins_filt_lowf = 1. ;; minutes for high-pass_filter
;mins_filt_highf = 0.3 ;; minutes for low-pass_filter
;;; values from Saito, 1969, T = 40 - 150 s
mins_filt_lowf = 150./60. ;; minutes for high-pass_filter
mins_filt_highf = 40./60. ;; minutes for low-pass_filter
;;; trend for the electric field to see whether convection happened.
mins_filt_trend = 3.
trange_load = time_double(trange_show)+[-100, 100]
trange_load_plasma_tha  = time_double(date+['1 0 1', '1 39 30'])+[-100, 100]

;; load kyoto AL
load_bin_data, trange = trange_load, datatype = 'kyoto_al', datafolder = kyoto_al_folder
load_bin_data, trange = trange_load, datatype = 'pseudo_al', datafolder = al_folder
options, 'kyoto_al', ytitle = 'Kyoto AL', ysubtitle = '[nT]', thick = l_thick
options, 'thg_idx_al', ytitle = 'THM AL', ysubtitle = '[nT]', thick = l_thick

;;;;;;;;;;;;;;;;;;;;;;;;;;;;; THEMIS D, E, A ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; plot settings

;;; fill or not fill plasma moments
fill_type = 'fill'
;fill_type = 'nofill'

;;; size in the order of a, d, e
size_x = [6,6,6]
size_y = [8,8,8]
x_abc = 0.2

;trange_show  = date+['1 17', '1 30'] ;; change trange_show to zoom in

probes  = ['a', 'e', 'd'] ;; e must be before 

for i = 0, n_elements(probes)-1 do begin
	sc = probes[i]
	;if ~strcmp(sc, 'a') then continue ;; plot choice probes

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
	del_data, 'th*' ;; for faster load


	if strcmp(sc, 'a') then trange_load_plasma = trange_load_plasma_tha else trange_load_plasma = trange_load
	thm_load_esansst2, trange=trange_load_plasma, probe = sc, i_mask = i_mask, e_mask = e_mask, /fill

	;;;;; to keep no fill moments, use the following
	;if strcmp(fill_type, 'nofill') then begin
	;	thm_load_esansst2, trange=trange_load, probe = sc, i_mask = i_mask, e_mask = e_mask
	;	copy_data, 'th'+sc+'_ptirf_eflux_energy', 'th'+sc+'_ptix_en_eflux'
	;	copy_data, 'th'+sc+'_pterf_eflux_energy', 'th'+sc+'_ptex_en_eflux'
	;endif

	;;; load this last because thm_load_esansst2 will load efs too
	load_efi_data, trange=trange_load, probe=sc, datatype = 'efs', /tclip, e_folder = efs_folder, b_folder = fgs_folder, rtrange = rtrange, angle_good = 15.

	calc,'"th?_state_pos_gsm_RE" = "th?_state_pos_gsm"/6374.4'
	split_vec, 'th'+sc+'_Pth'
	split_vec, 'th'+sc+'_state_pos_gsm_RE'
	;;; store the thermal + magnetic variable
	;; compute Pbx
	tinterpol_mxn, 'th'+sc+'_fgs_gsm', 'th'+sc+'_Pth', newname='th'+sc+'_fgs_gsm_int'
	get_data, 'th'+sc+'_fgs_gsm_int', t, bdata
	Pbx = nTesla2_to_nPa*bdata[*,0]^2
	get_data, 'th'+sc+'_Pall', t, pdata
	store_data, 'th'+sc+'_Pbth', data = {x:t, y:[[Pbx], [pdata[*,0:1]]]} ;; Pbx, Pb, Pth

	;;; compute field-aligned B
	cotrans2fac, 'th'+sc+'_fgs_gsm', b_variable = 'th'+sc+'_fgs_gsm', pos = 'th'+sc+'_state_pos_gsm', smooth_time = 6.*60., newname = 'th'+sc+'_fgs_fac'

	;;; compute spin-plane field-alighed E
	cotrans2fac, 'th'+sc+'_efs_dot0_dsl', b_variable = 'th'+sc+'_fgs_dsl', smooth_time = 6.*60., newname = 'th'+sc+'_efs_dot0_sfac', /sc_spin

	;;; compute B strength
	get_data, 'th'+sc+'_fgs_gsm', t, bdata
	store_data, 'th'+sc+'_fgs_b', data={x:t, y:sqrt(total(bdata^2, 2))}

	;;; filter
	tbandpass_filter, 'th'+sc+'_fgs_gsm', mins_filt_lowf*60., mins_filt_highf*60.
	tbandpass_filter, 'th'+sc+'_fgs_fac', mins_filt_lowf*60., mins_filt_highf*60.
	tbandpass_filter, 'th'+sc+'_fgs_b', mins_filt_lowf*60., mins_filt_highf*60.
	tbandpass_filter, 'th'+sc+'_Pth_z', mins_filt_lowf*60., mins_filt_highf*60.
	tbandpass_filter, 'th'+sc+'_Pbth', mins_filt_lowf*60., mins_filt_highf*60.
	tbandpass_filter, 'th'+sc+'_efs_gsm', mins_filt_lowf*60., mins_filt_highf*60.
	tbandpass_filter, 'th'+sc+'_efs_dot0_dsl', mins_filt_lowf*60., mins_filt_highf*60.
	tbandpass_filter, 'th'+sc+'_efs_dot0_sfac', mins_filt_lowf*60., mins_filt_highf*60.
	;; E trend
	tlowpass_filter, 'th'+sc+'_efs_dot0_sfac', mins_filt_trend*60.
	tlowpass_filter, 'th'+sc+'_efs_gsm', mins_filt_trend*60.
	tlowpass_filter, 'th'+sc+'_efs_vperp_gsm', mins_filt_trend*60.
	tlowpass_filter, 'th'+sc+'_ptix_vperp_gsm', mins_filt_trend*60.

	;;; compute poyinting vector
	tinterpol_mxn,'th'+sc+'_efs_gsm_bpfilt','th'+sc+'_fgs_gsm_bpfilt',newname='th'+sc+'_efs_gsm_bpfilt_int'
	get_data, 'th'+sc+'_efs_gsm_bpfilt_int', t, e_wave
	get_data, 'th'+sc+'_fgs_gsm_bpfilt', t, b_wave
	poynting = transpose(crossp_long(transpose(e_wave), transpose(b_wave), n_dir = poynting_dir))
	store_data, 'th'+sc+'_poynting', data={x:t, y:poynting}, dlimits = {data_att:{coord_sys:'gsm'}}
	store_data, 'th'+sc+'_poynting_dir', data={x:t, y:transpose(poynting_dir)}, dlimits = {data_att:{coord_sys:'gsm'}}
	cotrans2fac, 'th'+sc+'_poynting', b_variable = 'th'+sc+'_fgs_gsm', pos = 'th'+sc+'_state_pos_gsm', smooth_time = 6.*60., newname = 'th'+sc+'_poynting_fac'
	cotrans2fac, 'th'+sc+'_poynting_dir', b_variable = 'th'+sc+'_fgs_gsm', pos = 'th'+sc+'_state_pos_gsm', smooth_time = 6.*60., newname = 'th'+sc+'_poynting_dir_fac'

	;;; original quantities labels
	options, 'th'+sc+'_state_pos_gsm_RE_x', ytitle = 'X [R!dE]', ysubtitle = ''
	options, 'th'+sc+'_state_pos_gsm_RE_y', ytitle = 'Y [R!dE]', ysubtitle = ''
	options, 'th'+sc+'_state_pos_gsm_RE_z', ytitle = 'Z [R!dE]', ysubtitle = ''
	options, 'th'+sc+'_fgs_gsm*',colors=[2,4,6],labels=['B!dx','B!dy','B!dz'], ytitle = 'B', ysubtitle = '!c[nT]', labflag = 1, thick = l_thick
	options, 'th'+sc+'_fgs_fac*',colors=[2,4,6],labels=["B!dr'",'B!d'+phi_letter+"'",'B!db'], ytitle = 'B!dFAC', ysubtitle = '!c[nT]', labflag = 1, thick = l_thick
	options, 'th'+sc+'_fgs_b*', ytitle = '|B|', ysubtitle = '!c[nT]', thick = l_thick
	options, 'th'+sc+'_efs_dot0_dsl*',colors=[2,4,6],labels=['E!dx','E!dy','E!dz'], ytitle = 'E!dDSL', ysubtitle = '!c[mV/m]', labflag = 1, thick = l_thick
	options, 'th'+sc+'_efs_dot0_sfac*',colors=[2,4,6],labels=['E!d-ax','E!d'+phi_letter+'"',"E!db'"], ytitle = 'E!dsFAC', ysubtitle = '!c[mV/m]', labflag = 1, thick = l_thick
	options, 'th'+sc+'_efs_gsm*',colors=[2,4,6],labels=['E!dx','E!dy','E!dz'], ytitle = 'E', ysubtitle = '!c[mV/m]', labflag = 1, thick = l_thick
	options, 'th'+sc+'_ptix_vperp_gsm*',colors=[2,4,6],labels=['V!dx','V!dy','V!dz'], ytitle = 'V!di'+perp_sign+'!n', ysubtitle = '!c[km/s]', labflag = 1, thick = l_thick
	options, 'th'+sc+'_ptix_velocity_gsm*',colors=[2,4,6],labels=['V!dx','V!dy','V!dz'], ytitle = ' V!di!n', ysubtitle = '!c[km/s]', labflag = 1, thick = l_thick
	options, 'th'+sc+'_poynting_dir_fac*',colors=[2,4,6],labels=["S!dr'",'S!d'+phi_letter+"'",'S!db'], ytitle = 'S!dFAC!n dir', ysubtitle = '', labflag = 1, thick = l_thick
	options, 'th'+sc+'_poynting_fac*',colors=[2,4,6],labels=["S!dr'",'S!d'+phi_letter+"'",'S!db'], ytitle = 'S!dFAC!n', ysubtitle = '!c[unit]', labflag = 1, thick = l_thick
	options, 'th'+sc+'_Pth_z*', ytitle = 'P!dth!n', ysubtitle = '!c[nPa]', color = 0, labels = '', thick = l_thick
	options, 'th'+sc+'_Pbth*', ytitle = 'P', ysubtitle = '!c[nPa]', colors = [2,1,4], labels = ['P!dB!Bx', 'P!dB', 'P!dth'], thick = l_thick, labflag = 1
	options, 'th'+sc+'_ptix_en_eflux', ytitle = 'p!u+!n', ysubtitle = '!c[eV]', ztitle = ''
	options, 'th'+sc+'_ptex_en_eflux', ytitle = 'e!u'+minus_sign+'!n', ysubtitle = '!c[eV]', ztitle = ''

	;;; filtered data
	options, 'th'+sc+'_fgs_gsm_bpfilt', ytitle = 'B!ufilt!n'
	options, 'th'+sc+'_fgs_fac_bpfilt', ytitle = 'B!s!ufilt!r!dFAC!n'
	options, 'th'+sc+'_fgs_b_bpfilt', ytitle = '|B|!ufilt!n'
	options, 'th'+sc+'_Pth_z_bpfilt', ytitle = 'P!dth!n', ylog=0
	options, 'th'+sc+'_Pbth_bpfilt', ytitle = 'P!ufilt', ylog=0
	options, 'th'+sc+'_efs_gsm_bpfilt', ytitle = 'E!ufilt!n'
	options, 'th'+sc+'_efs_dot0_dsl_bpfilt', ytitle = 'E!s!dDSL!r!ufilt!n'
	options, 'th'+sc+'_efs_dot0_sfac_bpfilt', ytitle = 'E!s!dsFAC!r!ufilt!n'
	options, 'th'+sc+'_efs_vperp_gsm_lpfilt', ytitle = 'V!s!dExB!r!utrend'
	options, 'th'+sc+'_ptix_vperp_gsm_lpfilt', ytitle = 'V!s!di, perpr!r!utrend'

	;;;; limits
	ylim, 'th'+sc+'_fgs_gsm', -70, 70.
	ylim, 'th'+sc+'_efs_dot0_dsl', -30, 30.
	ylim, 'th'+sc+'_efs_vperp_gsm', -150, 150.
	ylim, 'th'+sc+'_ptix_vperp_gsm', -150, 150.
	ylim, 'th'+sc+'_efs_vperp_gsm_lpfilt', -200, 200.
	ylim, 'th'+sc+'_ptix_vperp_gsm_lpfilt', -200, 200.

	;;; tplot names
	tplotnames_line = 'th'+sc+'_'+['fgs_gsm', 'efs_dot0_dsl', 'Pth_z', 'ptix_vperp_gsm_lpfilt', 'efs_vperp_gsm_lpfilt']
	tnames_spec = 'th'+sc+'_'+['ptix_en_eflux', 'ptex_en_eflux']
	tplotnames = ['kyoto_al', tplotnames_line, tnames_spec]
	;labelnames = 'th'+sc+'_'+['state_pos_gsm_RE_z', 'state_pos_gsm_RE_y', 'state_pos_gsm_RE_x'] ;; no need for this because already have location plot

	;;; plot
	popen, pic_folder+'/th'+sc
	print_options,xsize=size_x[i], ysize=size_y[i]
	tplot, tplotnames, trange = trange_show, title = thm_probe_color(sc, /number, /long);, var_label = labelnames
	timebar, t_onset, line = 1
	timebar_mass, 0, /databar, varname = tplotnames_line, line = 1
	;; seperation between ESA and SST
	timebar_mass, 23000., /databar, varname = tnames_spec, line = 2
	pclose
endfor

stop
end
