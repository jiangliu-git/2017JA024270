pro main_plot_events_rbsp
;; plot an example of multi-point observation
thm_init
computer = 'I:'
;computer = '/home/jliu'
@folders

del_data, '*'

;; events to plot
;i_check = [66,71,13,20,25,60,70]

;; what purpose to plot
;purpose = 'view'
;purpose = 'publication'
purpose = 'wave'

;; what type to plot
;plot_type = 'normal' ;; mageis only, with density and Pth
plot_type = 'wave' ;; mageis only, with density and Pth
;plot_type = 'ect_both' ;; both mageis and hope, but no moments
;plot_type = 'plasmapause' ;; examine wither the change is caused by plasmapause crossing
;plot_type = 'diagnose'

plot_both = 'yes' ;; set to plot both probes for comparison
;plot_both = 'no'

;; whether to load electric field (which is slow)
;load_e = 'yes'
load_e = 'no'

;; lower cut for computing moments
energy_min = 200 ;; eV

;; minutes to filter data (for high pass)
minutes_filter = 5.

;; some constants for loading data
case plot_type of
'wave': minutes_load = 20.
'diagnose': minutes_load = 20.
'plasmapause': minutes_load = 4*60.
else: minutes_load = 4.
endcase

;;;;;; choose list
;; the type of detrend for the list
;list_dtr = 'model' ;; using model
;list_dtr = '10' ;; using 10-min average
list_dtr = 'both' ;; using model

list_name = 'dfb_rbsp_sub'+list_dtr+'_list_good'
;list_name = 'dfb_rbsp_list_bigEnoinj' ;; events with big E but no injection
;list_name = 'dfb_rbsp_list_smallEinj' ;; events with vanishing E but have injection

;;; load events
events = load_list(list_name+'.txt', folder = list_folder)
;; diagnose
;events = events(*, [16,40]) ;; typical no injec and have injec for examples in paper
events = events(*, 54)

times = time_double(events(0,*))
probes = events(3,*)

;;; settings for plotting data
spec = 0 ;; 1 for spectrum, 2 for line plot
if strcmp(purpose, 'publication') then spec = 1

;;; for publication quality plot
size_x = 6.6
size_y = 10.3
abc_inj = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']
abc_noinj = ['i', 'j', 'k', 'l', 'm', 'n', 'o', 'p']
x_abc = replicate(0.19, n_elements(abc_inj))
y_abc = [0.93, 0.83, 0.73, 0.6, 0.48, 0.38, 0.28, 0.19]
char_abc = 1.

for i = 0, n_elements(events(0,*))-1 do begin
;for i = 1, 1 do begin
	;;; only plot for selected events
	if keyword_set(i_check) then begin
		if total(i eq i_check) eq 0 then continue
	endif

	sc = probes[i]
	case sc of
	'r': probe_rb = 'a'
	's': probe_rb = 'b'
	endcase
	if strcmp(plot_both, 'yes') then probes_load = ['r', 's'] else probes_load = sc
	trange_show = [times(i)-minutes_load*60., times(i)+minutes_load*60.] ;;; for DFB
	trange_load = [trange_show(0)-60., trange_show(1)+60.]
	timespan, trange_load[0], trange_load[1]-trange_load[0], /sec
	rbsp_load, trange=trange_load, probe=probes_load, datatype = 'fgs', /tclip, rbsp_folder = rbsp_folder
	rbsp_load, trange=trange_load, probe=probes_load, datatype = 'mageis', /tclip, rbsp_folder = rbsp_folder, level = 2
	rbsp_load, trange=trange_load, probe=probes_load, datatype = 'hope_sa', /tclip, rbsp_folder = rbsp_folder, level = 2, reduce_connect = 'algebra'
	;; official load mageis (slower)
	if keyword_set('publication') then begin
		for j = 0, n_elements(probes_load)-1 do begin
			case probes_load[j] of
			'r': probe_rb_j = 'a'
			's': probe_rb_j = 'b'
			endcase
			rbsp_load_mageis_l2, probe = probes_rb_j
			copy_data, 'rbsp'+probe_rb_j+'_ect_mageis_L2_FESA', 'th'+probes_load[j]+'_mageis_espec_official'
			copy_data, 'rbsp'+probe_rb_j+'_ect_mageis_L2_FPSA', 'th'+probes_load[j]+'_mageis_pspec_official'
		endfor
	endif

	if strcmp(plot_type, 'plasmapause') then rbsp_load, trange=trange_load, probe=probes_load, datatype = 'hfr_spectra_merged', /tclip, rbsp_folder = rbsp_folder
	rbsp_load, trange=trange_load, probe=probes_load, datatype = 'efw_density', /tclip, rbsp_folder = rbsp_folder
	if strcmp(load_e, 'yes') then begin
		catch, err
		if err eq 0 then begin
			rbsp_load_efw, trange=trange_load, probe=probes_load, /tclip
			;rbsp_load_density, trange=trange_load, probe=probes_load, /tclip
			for j = 0, n_elements(probes_load)-1 do begin
				rbsp_mgse2gse, 'th'+probes_load[j]+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit', newname='th'+probes_load[j]+'_efw_esvy_gse', /no_spice_load
				cotrans, 'th'+probes_load[j]+'_efw_esvy_gse', 'th'+probes_load[j]+'_efw_esvy_gsm', /GSE2GSM
				options, 'th'+probes_load[j]+'_efw_esvy_gsm', ytitle = 'E', ysubtitle = '!c[mV/m]', labels = ['E!dx', 'E!dy', 'E!dz'], labflag = 1, colors = [2,4,6]
				if strcmp(purpose, 'publication') then ylim, 'th'+probes_load[j]+'_efw_esvy_gsm', -9., 19.
			endfor
		endif else begin
			dprint, !error_state.msg
		endelse
		catch, /cancel
	endif

	;;;; manage tplot variables
	;; energy channels
	if spec eq 0 then begin
		if tv_exist('th'+sc+'_mageis_pspec_tclip') then begin
			get_data, 'th'+sc+'_mageis_pspec_tclip', data = pspec
			pspec_mageis_labels = strcompress(string(pspec.v(0,*),format='(i)'), /remove)
		endif
		if tv_exist('th'+sc+'_mageis_espec_tclip') then begin
			get_data, 'th'+sc+'_mageis_espec_tclip', data = espec
			espec_mageis_labels = strcompress(string(espec.v(0,*),format='(i)'), /remove)
		endif
		if tv_exist('th'+sc+'_hope_sa_pspec_tclip') then begin
			get_data, 'th'+sc+'_hope_sa_pspec_tclip', data = pspec
			pspec_hope_labels = strcompress(string(pspec.v(0,*),format='(i)'), /remove)
		endif
		if tv_exist('th'+sc+'_hope_sa_espec_tclip') then begin
			get_data, 'th'+sc+'_hope_sa_espec_tclip', data = espec
			espec_hope_labels = strcompress(string(espec.v(0,*),format='(i)'), /remove)
		endif
	endif else begin
		pspec_mageis_labels = ''
		espec_mageis_labels = ''
		pspec_hope_labels = ''
		espec_hope_labels = ''
	endelse

	for j = 0, n_elements(probes_load)-1 do begin
		;;; !!! do not change probes_load[j] to sc here or there will be errors.
		if strcmp(plot_type, 'diagnose') then begin
			;;; diagnose data (Antenna potentials)
			case probes_load[j] of
			'r': probe_rb_j = 'a'
			's': probe_rb_j = 'b'
			endcase
			catch, err
			if err eq 0 then begin
				rbsp_load_efw_waveform, trange = trange_load, probe=probe_rb_j, type='calibrated', datatype='vsvy'
				get_data, 'rbsp'+probe_rb_j+'_efw_vsvy', t, vsvy_all
				store_data, 'rbsp'+probe_rb_j+'_efw_vsvy_1234', data = {x:t, y:vsvy_all[*,0:3]}
				options, 'rbsp'+probe_rb_j+'_efw_vsvy_1234', colors = [1,2,3,4], labels = ['V1', 'V2', 'V3', 'V4'], labflag = 1
			endif else begin
				dprint, !error_state.msg
			endelse
			catch, /cancel
		endif

		;;;; band pass filter fields
		thigh_pass_filter, 'th'+probes_load[j]+'_fgs_gsm_tclip', minutes_filter*60.
		thigh_pass_filter, 'th'+probes_load[j]+'_efw_esvy_gsm', minutes_filter*60.
		;; remove spin frequency
		tsmooth2, 'th'+probes_load[j]+'_fgs_gsm_tclip_hpfilt', 8
		tsmooth2, 'th'+probes_load[j]+'_efw_esvy_gsm_hpfilt', 3

		;;;; fields
		copy_data, 'th'+probes_load[j]+'_fgs_gsm_tclip', 'th'+probes_load[j]+'_fgs_gsm_original'
		split_vec, 'th'+probes_load[j]+'_fgs_gsm_original'
		detrend_b, 'th'+probes_load[j]+'_fgs_gsm_tclip', method = 'T96', rbsp_folder = rbsp_folder, pre_hour = 2
		split_vec, 'th'+probes_load[j]+'_fgs_gsm_tclip'
		tsmooth2, 'th'+probes_load[j]+'_fgs_gsm_tclip_z', 3
		deriv_data, 'th'+probes_load[j]+'_fgs_gsm_tclip_z_sm'
		;; total field
		get_data, 'th'+probes_load[j]+'_fgs_gsm_original', t, b
		bttl = sqrt(total(b^2,2))
		store_data, 'th'+probes_load[j]+'_fgs_original_ttl', data = {x:t, y:bttl}
		time_clip, 'th'+probes_load[j]+'_fgs_original_ttl', trange_show[0], trange_show[1]
		get_data, 'th'+probes_load[j]+'_fgs_original_ttl_tclip', t, bttl_show
		bttl_max = max(bttl_show, /nan)
		bttl_min = min(bttl_show, /nan)
		bttl_span = bttl_max-bttl_min
		ylim, 'th'+probes_load[j]+'_fgs_original_ttl', bttl_min-0.1*bttl_span, bttl_max+0.1*bttl_span

		;;;; plasmapause deciding variable
		if strcmp(plot_type, 'plasmapause') then begin
			;; compute fce
			calc, "'th"+probes_load[j]+"_fce' = 28*sqrt(total('th"+probes_load[j]+"_fgs_gsm_original'^2, 2))"
			;; change density to frequency
			tinterpol_mxn,'th'+probes_load[j]+'_fce','th'+probes_load[j]+'_efw_density_tclip',newname='th'+sc+'_fce_int'
			calc, "'th"+probes_load[j]+"_efw_densfreq' = sqrt(8980.^2*'th"+probes_load[j]+"_efw_density_tclip'+'th"+probes_load[j]+"_fce_int'^2)"
			store_data, 'th'+probes_load[j]+'_hfr_spectra_efw_density', data = ['th'+probes_load[j]+'_hfr_spectra_tclip', 'th'+probes_load[j]+'_efw_densfreq']
			options, 'th'+probes_load[j]+'_efw_densfreq', colors = 1
			options,'th'+probes_load[j]+'_hfr_spectra_tclip', ylog = 1, zlog = 1
			ylim, 'th'+probes_load[j]+'_hfr_spectra_tclip', 1e4, 5e5
			ylim, 'th'+probes_load[j]+'_hfr_spectra_efw_density', 1e4, 5e5
		endif

		;;;; moments
		combine_spec, 'th'+probes_load[j]+'_mageis_pspec_tclip', 'th'+probes_load[j]+'_hope_sa_pspec_tclip', newname = 'th'+probes_load[j]+'_pspec_combined_tclip', /eV2keV_2nd
		combine_spec, 'th'+probes_load[j]+'_mageis_espec_tclip', 'th'+probes_load[j]+'_hope_sa_espec_tclip', newname = 'th'+probes_load[j]+'_espec_combined_tclip', /eV2keV_2nd
		;; all population for pressures
		compute_moments, 'th'+probes_load[j]+['_pspec_combined_tclip', '_espec_combined_tclip', '_hope_sa_hespec_tclip', '_hope_sa_ospec_tclip'], intypes = 'flux', inunits_energy = ['keV', 'keV', 'eV', 'eV'], inunits_flux = 'keV', energy_min = energy_min, particle_types = ['p', 'e', 'He+', 'O+'], /combine, newname_combine = 'all', rbsp_folder = rbsp_folder
		;; all ions
		compute_moments, 'th'+probes_load[j]+['_pspec_combined_tclip', '_hope_sa_hespec_tclip', '_hope_sa_ospec_tclip'], intypes = 'flux', inunits_energy = ['keV', 'eV', 'eV'], inunits_flux = 'keV', energy_min = energy_min, particle_types = ['p', 'He+', 'O+'], /combine, newname_combine = 'ion', rbsp_folder = rbsp_folder
		;; electron
		compute_moments, 'th'+probes_load[j]+'_espec_combined_tclip', intype = 'flux', inunits_energy = 'keV', inunits_flux = 'keV', energy_min = energy_min, particle_types = 'e', rbsp_folder = rbsp_folder

		;;;; locations
		get_data, 'th'+probes_load[j]+'_state_pos_tclip', t_pos, loc_pos
		store_data, 'th'+probes_load[j]+'_state_pos_RE_tclip', data = {x:t_pos, y:loc_pos/RE} 
		store_data, 'th'+probes_load[j]+'_state_pos_RE_tclip_rho', data = {x:t_pos, y:sqrt(total((loc_pos/RE)^2, 2))} 
		;;;; get L* values
		;; T89
		;ttrace2equator, 'th'+probes_load[j]+'_state_pos_RE_tclip', newname='th'+probes_load[j]+'_efoot',external_model='t89',par=2.0D,in_coord='gsm',out_coord='gsm'
		;; T96
		ttrace2equator96, 'th'+probes_load[j]+'_state_pos_RE_tclip'
		get_data, 'th'+probes_load[j]+'_efoot', t_efoot, efoot
		store_data, 'th'+probes_load[j]+'_Lstar', data = {x:t_efoot, y:sqrt(total(efoot[*,0:1]^2,2))}

		if strcmp(purpose, 'publication') then suf_label = '' else suf_label = probes_load[j]
		split_vec, 'th'+probes_load[j]+'_state_pos_RE_tclip'
		options, 'th'+probes_load[j]+'_state_pos_RE_tclip_x', ytitle = suf_label+'X [R!dE!n]'
		options, 'th'+probes_load[j]+'_state_pos_RE_tclip_y', ytitle = suf_label+'Y [R!dE!n]'
		options, 'th'+probes_load[j]+'_state_pos_RE_tclip_z', ytitle = suf_label+'Z [R!dE!n]'
		options, 'th'+probes_load[j]+'_state_pos_RE_tclip_rho', ytitle = suf_label+rho_letter+' [R!dE!n]'
		options, 'th'+probes_load[j]+'_Lstar', ytitle = suf_label+'L'
	endfor

	;;; fields
	options, 'th?_fgs_gsm_*',colors=[2,4,6],labels=['B!dx','B!dy','B!dz'], ytitle = 'B', ysubtitle = '!c[nT]', labflag = 1
	options, 'B_subtr',colors=[2,4,6],labels=['B!dx','B!dy','B!dz'], labflag = 1
	options, 'th?_fgs_gsm*_z*',colors=6, ytitle = 'B!dz!n', ysubtitle = '!c[nT]', labels = ''
	options, 'th?_fgs_original_ttl',colors=0, ytitle = '|B|', ysubtitle = '!c[nT]'
	options, 'th?_efw_density*', ytitle = 'n!dEFW!n'
	;; filtered
	options, 'th?_fgs_gsm_tclip_hpfilt', ytitle = 'B filtered'
	options, 'th?_efw_esvy_gsm', ytitle = 'E filtered', ysubtitle = '!c[mV/m]'

	;;; moments
	options, 'th?_ion_density*', ytitle = 'n!di!n', ysubtitle = '!c[cm!u-3!n]'
	if strcmp(purpose, 'publication') and (i eq 1) then begin
		ylim, 'th?_ion_density*', 0.6, 1.1
	endif
	options, 'th?_all_Pall*', ytitle = 'Pressures', ysubtitle = '!c[nPa]'

	;;; spectra
	if strcmp(purpose, 'publication') then ztitle = '' else ztitle = 'Flux!c[cm!u-2!ns!u-1!nsr!u-1!nkeV!u-1!n]'
	if spec then begin
		ytitle_p = 'p!u+!n energy'
		ytitle_e = 'e!u'+minus_sign+'!n energy'
	endif else begin
		ytitle_p = 'Proton flux'
		ytitle_e = 'Electron flux'
	endelse

	options,'th?_mageis_pspec_*', ylog = 1, zlog = 1, ytitle = ytitle_p, ysubtitle = '!c[keV]', ztitle = ztitle, spec = spec
	options,'th?_mageis_espec_*', ylog = 1, zlog = 1, ytitle = ytitle_e, ysubtitle = '!c[keV]', ztitle = ztitle, spec = spec
	options,'th?_hope_sa_pspec_tclip', ylog = 1, zlog = 1, ytitle = ytitle_p, ysubtitle = '!c[eV]', ztitle = ztitle, spec = spec, labels = pspec_hope_labels 
	options,'th?_hope_sa_espec_tclip', ylog = 1, zlog = 1, ytitle = ytitle_e, ysubtitle = '!c[eV]', ztitle = ztitle, spec = spec, labels = espec_hope_labels 
	options,'th?_mageis_pspec_tclip', labels = pspec_mageis_labels 
	options,'th?_mageis_espec_tclip', labels = espec_mageis_labels 
	if spec eq 1 then begin
		ylim,'th?_mageis_pspec_*', 60., 1100. 
		ylim,'th?_mageis_espec_*', 40., 4000.
		ylim,'th?_hope_sa_pspec_tclip', 24., 52000. 
		ylim,'th?_hope_sa_espec_tclip', 14., 52000.
		;if strcmp(purpose, 'publication') then begin
		;	zlim,'th?_mageis_pspec_tclip', 1.1, 1.1e5. 
		;	zlim,'th?_mageis_espec_tclip', 1.1, 1.1e5.
		;endif
	endif

	;;; decide whether this is a injection or not
	inj_p = injection_decide('th'+sc+'_mageis_pspec', t0 = times[i])
	inj_e = injection_decide('th'+sc+'_mageis_espec', t0 = times[i])
	;; injection
	if (inj_p eq 1) or (inj_e eq 1) then begin
		class_folder = 'injection'
	endif else begin
		;; no injection
		if (inj_p eq 0) and (inj_e eq 0) then begin
			class_folder = 'no_injection'
		endif else begin
			class_folder = 'unclear'
		endelse
	endelse

	;; create the save folder if not exist
	if strcmp(purpose, 'publication') then begin
		save_folder = pic_folder
		if strcmp(class_folder, 'injection') then title = 'Injection example (Probe '+strupcase(probe_rb)+')' else begin
			if strcmp(class_folder, 'no_injection') then title = 'No injection example (Probe '+strupcase(probe_rb)+')' else title = 'Unclear example (Probe '+strupcase(probe_rb)+')'
		endelse
	endif else begin
		save_folder = pic_folder+'/events_rbsp/'+class_folder
		if ~file_test(save_folder, /directory) then file_mkdir, save_folder
		title = 'DFB '+strcompress(string(i), /remove)+' RBSP-'+probe_rb
	endelse

;	popen, save_folder+'/event_'+strcompress(string(i), /remove)
;	print_options,xsize=size_x,ysize=size_y
	if strcmp(plot_both, 'yes') then begin
		case sc of
		'r': sc2 = 's'
		's': sc2 = 'r'
		endcase
		if strcmp(plot_type, 'normal') then begin
			tvnames_plot = ['th'+sc+'_fgs_gsm_original_z', 'th'+sc+'_mageis_pspec_tclip', 'th'+sc+'_mageis_espec_tclip', 'th'+sc+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit', 'th'+sc2+'_fgs_gsm_original_z', 'th'+sc2+'_mageis_pspec_tclip', 'th'+sc2+'_mageis_espec_tclip', 'th'+sc2+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit']
		endif
		if strcmp(plot_type, 'wave') then begin
			tvnames_plot = ['th'+sc+'_fgs_gsm_original', 'th'+sc+'_fgs_gsm_tclip_hpfilt_sm', 'th'+sc+'_efw_esvy_gsm_hpfilt_sm', 'th'+sc+'_mageis_pspec_tclip', 'th'+sc+'_mageis_espec_tclip', 'th'+sc2+'_fgs_gsm_original', 'th'+sc2+'_fgs_gsm_tclip_hpfilt_sm', 'th'+sc2+'_efw_esvy_gsm_hpfilt_sm', 'th'+sc2+'_mageis_pspec_tclip', 'th'+sc2+'_mageis_espec_tclip']
		endif
		varlabels = ['th'+sc2+'_state_pos_RE_tclip_z', 'th'+sc2+'_state_pos_RE_tclip_y', 'th'+sc2+'_state_pos_RE_tclip_x', 'th'+sc2+'_Lstar', 'th'+sc+'_state_pos_RE_tclip_z', 'th'+sc+'_state_pos_RE_tclip_y', 'th'+sc+'_state_pos_RE_tclip_x', 'th'+sc+'_Lstar']
	endif else begin
		varlabels = ['th'+sc+'_state_pos_RE_tclip_z', 'th'+sc+'_state_pos_RE_tclip_y', 'th'+sc+'_state_pos_RE_tclip_x', 'th'+sc+'_Lstar']
		;;; normal plot
		if strcmp(plot_type, 'normal') then begin
;			tvnames_plot = ['th'+sc+'_fgs_gsm_original', 'th'+sc+'_fgs_gsm_original_z', 'th'+sc+'_fgs_original_ttl', 'th'+sc+'_efw_density', 'th'+sc+'_ion_density', 'th'+sc+'_all_Pall', 'th'+sc+'_mageis_pspec_tclip', 'th'+sc+'_mageis_espec_tclip', 'th'+sc+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit']
			if strcmp(purpose, 'publication') then begin
				tvnames_plot = ['th'+sc+'_fgs_gsm_original', 'th'+sc+'_fgs_gsm_original_z', 'th'+sc+'_fgs_original_ttl', 'th'+sc+'_ion_density', 'th'+sc+'_all_Pall', 'th'+sc+'_mageis_pspec_official', 'th'+sc+'_mageis_espec_official', 'th'+sc+'_efw_esvy_gsm']
				options, tvnames_plot, thick = l_thick
				;;;; avoid repeating labels
				;if ~strcmp(class_folder, 'injection') then begin
				;	options, tvnames_plot, ytitle = '', ysubtitle = ''
				;	options, varlabels, ytitle = '', ysubtitle = ''
				;endif else begin
				;	options, tvnames_plot, labels = ''
				;endelse
			endif else begin
				tvnames_plot = ['th'+sc+'_fgs_gsm_original', 'th'+sc+'_fgs_gsm_original_z', 'th'+sc+'_fgs_original_ttl', 'th'+sc+'_ion_density', 'th'+sc+'_all_Pall', 'th'+sc+'_mageis_pspec_tclip', 'th'+sc+'_mageis_espec_tclip', 'th'+sc+'_efw_esvy_gsm']
			endelse
		endif
		;;; both mageis and hope
		if strcmp(plot_type, 'ect_both') then begin
			tvnames_plot = ['th'+sc+'_fgs_gsm_original', 'th'+sc+'_fgs_gsm_original_z', 'th'+sc+'_fgs_original_ttl', 'th'+sc+'_mageis_pspec_tclip', 'th'+sc+'_hope_sa_pspec_tclip', 'th'+sc+'_mageis_espec_tclip', 'th'+sc+'_hope_sa_espec_tclip', 'th'+sc+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit']
		endif
		;;; examine plasmapause crossing
		if strcmp(plot_type, 'plasmapause') then begin
			tvnames_plot = ['th'+sc+'_fgs_gsm_original', 'th'+sc+'_fgs_gsm_original_z', 'th'+sc+'_fgs_original_ttl', 'th'+sc+'_ion_density', 'th'+sc+'_hfr_spectra_efw_density', 'th'+sc+'_mageis_pspec_tclip', 'th'+sc+'_mageis_espec_tclip', 'th'+sc+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit']
		endif
		;;; diagnose plot
		if strcmp(plot_type, 'diagnose') then begin
			tvnames_plot = ['th'+sc+'_fgs_gsm_original', 'th'+sc+'_mageis_pspec_tclip', 'th'+sc+'_mageis_espec_tclip', 'th'+sc+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit', 'rbsp'+probe_rb+'_efw_vsvy_1234']
		endif
	endelse
	tplot, tvnames_plot, var_label = varlabels, trange = trange_show, title = title
	if strcmp(plot_type, 'diagnose') then timebar_mass, [-190., 190.], varname = 'rbsp'+probe_rb+'_efw_vsvy_1234', /databar
	timebar, times[i], line = 1
	if strcmp(purpose, 'publication') then begin
		if strcmp(class_folder, 'injection') then abc_this = abc_inj else abc_this = abc_noinj
		xyouts, 0.943, 0.34, 'Flux [cm!u-2!ns!u-1!nsr!u-1!nkeV!u-1!n]', /normal, orientation = 90, align = 0.5
		xyouts, x_abc, y_abc, '('+abc_this+')', charsize = char_abc, /normal
	endif
;	makepng, save_folder+'/event_'+strcompress(string(i), /remove)
;	pclose
endfor ;; for of i, element on multi_list

stop
end
