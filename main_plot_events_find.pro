pro main_plot_events_find
;; plot an example of multi-point observation
thm_init
computer = 'I:'
;computer = '/home/jliu'
@folders
pic_folder = pic_folder+'/dfbs_conj'

del_data, '*'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; settings ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;; THEMIS settings
plot_themis = 0
;plot_themis = 1 ;; plot themis

;plot_thm_type = 'one' ;; plot only the one in the list
;plot_thm_type = 'all' ;; plot all available THEMIS spacecraft.
;plot_thm_type = 'artemis' ;; plot all artemis probes. For event 1784 thb is unavailable, use 'c' instead
plot_thm_type = 'c' ;; plot these probes only, can set to be any probes

plot_thm_plasma = 1 ;; plot energy spectra (very slow)
;plot_thm_plasma = 0

;; good location of THEMIS
rho_thm_min = 2.5 ;; rho must be larger than this

;;;;;;;;;;;;;;;;;;;;;; RBSP settings
plot_rbsp = 0
;plot_rbsp = 1 ;; plot rbsp

plot_all_rb = 'yes' ;; plot all available RBSP spacecraft.
;plot_all_rb = 'no' ;; plot only one RBSP spacecraft

;rb_probes_plot = 'r' ;; comment to use all available
;rb_probes_plot = 's' ;; comment to use all available

;; whether to load RBSP electric field (which is slow)
;load_e = 'yes'
load_e = 'no'

;; good location of RBSP
rho_rb_min = 2.5 ;; rho must be larger than this

;; lower cut for computing moments
energy_min = 200 ;; eV

;;;;;;;;;;;;;;;;;;;;; GOES settings
;plot_goes = 0
plot_goes = 1 ;; whether to plot goes
probes_goes_all = '13'
;probes_goes_all = ['13', '14', '15']

;;;;;;;;;;;;;;;;;;;;; OMNI (solarwind) settings
plot_omni = 0
;plot_omni = 1 ;; plot omni solarwind data

;;;;;;;;;;;;;;;;;;;;; GBO settings
plot_gbo = 0
;plot_gbo = 1 ;; plot gruond magnetometer data (Pi2s)
;;;; required locations (between these)
;;; magnetic latitude
;mlat_suf = '' ;; all latitudes
;mlat_suf = '_low' ;; low latitudes
;mlat_suf = '_mid' ;; mid latitudes
mlat_suf = '_high' ;; high latitudes
;;; magnetic local time
mlt_suf = '' ;; all MLTs
;mlt_suf = '_night'
;mlt_suf = '_day'
;;; order in plort
;sort_suf = '_slat' ;; sort by latitude (high to low)
sort_suf = '_smlt' ;; sort by MLT (east to west)
;;;; numbers for different setting
case mlat_suf of
'_low': begin
	mlat_low = 0.
	mlat_high = 65.
	end
'_mid': begin
	mlat_low = 65. 
	mlat_high = 72.
	end
'_high': begin
	mlat_low = 72. 
	mlat_high = 999.
	end
else: begin
	mlat_low = 0. 
	mlat_high = 999.
	end
endcase

;;;;;;;;;;;;;;;;;;;;; settings for filter data
;;; whether to filter many times (for magnetic field only)
;crazy_filter_rb = [''] ;; disable
;;crazy_filter_rb = ['s'] ;; put rbsp probes here, r and/or s

;; choose themis probes to filter
;filter_thm = [''] ;; disable
;filter_thm = ['a', 'd', 'e']
filter_thm = ['b', 'c'] ;; solar wind parameters

;;;;; to reveal the 2-min period wave
;; minutes to filter data (for high pass)
minutes_filter = 3.
;; points to smooth data after filter
sm_width_b = 3
sm_width_b_goes = 8
sm_width_e = 3 ;; for rbsp only

;;;;;; to reveal the 5-min period wave
;;; minutes to filter data (for high pass)
;minutes_filter = 6.
;;; points to smooth data after filter
;sm_width_b = 30
;sm_width_b_goes = 80
;sm_width_e = 10

;;; disable smoothing
sm_width_b = 0
sm_width_b_goes = 0
sm_width_e = 0

;;;;;;;;;;;;;;;;;;;;; settings for plotting data
;; minutes before and after to plot data
minutes_show = 20.
;minutes_show = 100.
minutes_load = minutes_show+10.

;seperate_plot = 1 ;; seperate the plot
seperate_plot = 0 ;; do a whole plot

;crazy_filter = 1 ;; repeat filter for three times
crazy_filter = 0 ;; repeat filter for three times

;; whether to draw the projected locations
;project = 'yes'
project = 'no'

;;; plasma data setting
spec = 0 ;; 1 for spectrum, 0 for line plot

;;; ranges of plot
xrange = [0, -15.]
yrange = [12, -12.]

;;; parameters for plot
left_margin = 0.15
right_margin = 0.15
top_margin = 0.05
bot_margin = 0.07
space_vert = 0.01
if strcmp(project, 'yes') then n_panels = 2 else n_panels = 1
positions = panel_positions([1, n_panels], lr_margins = [left_margin, right_margin], bt_margins = [bot_margin, top_margin], space = [0, space_vert])
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; end of settings ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;; load events
list_name = 'dfb_list_lead_tail_all' ;;; themis list
;list_name = 'dfb_rbsp_subboth_list_good' ;;; rbsp list
events = load_list(list_name+'.txt', folder = list_folder)

;;;;;;;;; choose events to plot
;;;;;;; the numbers below are for THEMIS list

;;;;;; injection events
;pic_folder = pic_folder+'/injection'
;;i_plot = [1581, 1582, 1589, 1609, 1696, 1703, 1716, 1748, 1765, 1790, 1791, 1792, 1804, 1810, 1817, 1818, 1842] 
;i_plot = 1810

;;;;;;; wave events
pic_folder = pic_folder+'/wave'
;i_plot = [1710, 1754, 1759, 1772, 1773, 1774, 1775, 1776, 1777, 1783, 1784, 1787, 1789, 1813, 1836, 1837, 1838, 1839, 1869, 1876, 1877]
i_plot = 1784 ;; good wave events (1784, 1785). 1784 is for the paper.

;;;;;;;;;;;;;;;;; compute locations of GBO stations
if plot_gbo then begin
	thm_gmag_stations, gbo_stations, locations_ll, magnetic = locations_ll_mag, midnight = stations_midnight ;; locations_geo: latitude, longitude.
	gbo_stations = strlowcase(gbo_stations)
	;;;;;; settings for specific events
	if keyword_set(i_plot) then begin
		if i_plot eq 1784 then begin
			gbo_stations_this = ['atha', 'ccnv', 'fcc', 'frn', 'fsmi', 'gua', 'hon', 'larg', 'mea', 'pine', 'tpas', 'tuc', 'ukia', 'vic', 'whit', 'ykc', 'bou', 'bfe', 'dob', 'don', 'drby', 'frd', 'gbay', 'glyn', 'jck', 'kapu', 'kar', 'lrv', 'new', 'ott', 'pina', 'roe', 'rvk', 'satx', 'sjg', 'sol', 'stj', 'swno', 'tpas', 'tuc', 'wrth', 'dmh', 'kuv', 'lyr', 'svs', 'thl', 'umq', 'upn', 'abk', 'amk', 'amd', 'and','bjn', 'fcc', 'fhb', 'ghb', 'gill', 'iqa', 'kuuj', 'naq', 'nor', 'nrsq', 'sco', 'skt', 'snkq', 'sor', 'tro']
		endif
		gbo_stations_this = gbo_stations_this[uniq(gbo_stations_this, sort(gbo_stations_this))]
		locations_ll_this = fltarr(2, n_elements(gbo_stations_this))
		locations_ll_mag_this = fltarr(2, n_elements(gbo_stations_this))
		stations_midnight_this = strarr(n_elements(gbo_stations_this))
		for i = 0, n_elements(gbo_stations_this)-1 do begin
			i_match = where(strcmp(gbo_stations_this[i], gbo_stations), n_match)
			if n_match ne 1 then begin
				print, gbo_stations_this[i]
				message, 'The above station name is written wrong.'
			endif
			locations_ll_this[*,i] = locations_ll[*,i_match]
			locations_ll_mag_this[*,i] = locations_ll_mag[*,i_match]
			stations_midnight_this[i] = stations_midnight[i_match]
		endfor
		gbo_stations = gbo_stations_this
		locations_ll = locations_ll_this
		locations_ll_mag = locations_ll_mag_this
		stations_midnight = stations_midnight_this
	endif ;; if of setting i_plot
	locations_geo = latlongalt2GEO([locations_ll, replicate(1., 1, n_elements(gbo_stations))])
endif ;; if of plot GBO

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; end of all settings ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

times = time_double(events[0,*])
probes = events[3,*]

;for i = 0, n_elements(events[0,*])-1 do begin ;;; plot all
for k_plot = 0, n_elements(i_plot)-1 do begin ;;; controlled plot
	i = i_plot[k_plot]

	del_data, '*'

	sc = probes[i] ;; this is the sc of the DFB observation, can be either RBSP or THEMIS depending on the loaded list
	if times[i] lt time_double('2012 9 30') then continue

	trange_load = [times[i]-minutes_load*60., times[i]+minutes_load*60.]
	trange_show = [times[i]-minutes_show*60., times[i]+minutes_show*60.]

	;;;;;;;;;;;;;; first, decide whether to make the plot based on RBSP location (using THMIS list) or THEMIS (using RBSP list)
	if strmatch(list_name, '*rbsp*') then begin
		;;; start from RBSP list (see if THEMIS is within)
		probes_thm_test = ['a', 'd', 'e']
		load_bin_data, trange=trange_load, probe=probes_thm_test, datatype = 'pos', /tclip, datafolder = pos_folder
		if_this_event = 0
		thm_probes = ''
		for i_sc = 0, n_elements(probes_thm_test)-1 do begin
			sc_thm = probes_thm_test[i_sc]
			if tv_exist('th'+sc_thm+'_state_pos_tclip') then begin
				get_data, 'th'+sc_thm+'_state_pos_tclip', t_pos, loc_pos
				x = loc_pos[*,0]/RE
				y = loc_pos[*,1]/RE
				rho = sqrt(x^2+y^2)
				no_use = where((x lt 0) and (rho gt rho_thm_min), n_within)
				if n_within gt 0 then begin
					if_this_event = 1
					thm_probes = [thm_probes, sc_thm]
				endif
			endif
		endfor
		if n_elements(thm_probes) gt 1 then thm_probes= thm_probes[1:*]

		;;; decide RBSP probes to plot
		if strcmp(plot_all_rb, 'yes') then rb_probes_load = ['r', 's'] else rb_probes_load = sc
		rbsp_load, trange=trange_load, probe=rb_probes_load, datatype = 'fgs', /tclip, rbsp_folder = rbsp_folder
	endif else begin
		;;; start from THEMIS list
		catch, err
		if err eq 0 then begin
			rbsp_load, trange=trange_load, probe=['r', 's'], datatype = 'fgs', /tclip, rbsp_folder = rbsp_folder
			get_data, 'thr_state_pos_tclip', t_pos_r, loc_pos_r
			x_r = loc_pos_r[*,0]/RE
			y_r = loc_pos_r[*,1]/RE
			rho_r = sqrt(x_r^2+y_r^2)
		endif else begin
			dprint, !error_state.msg
		endelse
		catch, /cancel

		catch, err
		if err eq 0 then begin
			get_data, 'ths_state_pos_tclip', t_pos_s, loc_pos_s
			x_s = loc_pos_s[*,0]/RE
			y_s = loc_pos_s[*,1]/RE
			rho_s = sqrt(x_s^2+y_s^2)
		endif else begin
			dprint, !error_state.msg
		endelse
		catch, /cancel

		if_this_event = 0
		n_r_within = 0
		n_s_within = 0
		rb_probes_load = ['']
		if tv_exist('thr_state_pos_tclip') then no_use = where((x_r lt 0) and (rho_r gt rho_rb_min), n_r_within)
		if tv_exist('ths_state_pos_tclip') then no_use = where((x_s lt 0) and (rho_s gt rho_rb_min), n_s_within)
		if n_r_within gt 0 then begin
			if_this_event = 1
			rb_probes_load = [rb_probes_load, 'r']
		endif
		if n_s_within gt 0 then begin
			if_this_event = 1
			rb_probes_load = [rb_probes_load, 's']
		endif
	endelse

	if if_this_event then begin
		;;;;;;;;;;;;;;;; load THEMIS data ;;;;;;;;;;;;;;;;;;;
		if ~ strmatch(list_name, '*rbsp*') then begin
			case plot_thm_type[0] of
			'all': thm_probes = ['a', 'd', 'e']
			'one': thm_probes = sc
			'artemis': thm_probes = ['b', 'c']
			else: thm_probes = plot_thm_type
			endcase
		endif
		;;; if load themis then load
		if plot_themis and ~strcmp(thm_probes[0], '') then begin
			for j = 0, n_elements(thm_probes)-1 do begin
				sc_thm = thm_probes[j]
				if plot_thm_plasma or strcmp_or(sc_thm, ['b', 'c']) then begin
					;;; load plasma data
					thm_load_esansst2, trange=trange_load, probe = sc_thm
					time_clip, 'th'+sc_thm+'_ptix_velocity_gsm', trange_load[0], trange_load[1]
					time_clip, 'th'+sc_thm+'_ptix_density', trange_load[0], trange_load[1]
				endif
				;;; load other data
				if strcmp_or(sc_thm, ['b', 'c']) then begin
					thm_load_state, trange=trange_load, probe=sc_thm, datatype = 'pos', coord = 'gsm'
					thm_load_fit, trange=trange_load, probe=sc_thm, datatype = 'fgs', coord = 'gsm', level = 1, suffix = '_gsm'
					time_clip, 'th'+sc_thm+'_state_pos', trange_load[0], trange_load[1]
					time_clip, 'th'+sc_thm+'_fgs_gsm', trange_load[0], trange_load[1]
				endif else begin
					load_bin_data, trange=trange_load, probe=sc_thm, datatype = 'pos', /tclip, datafolder = pos_folder
					load_bin_data, trange=trange_load, probe=sc_thm, datatype = 'fgs', /tclip, datafolder = fgs_folder
					load_bin_data, trange=trange_load, probe=sc_thm, datatype = 'viperp', /tclip, datafolder = viperp_folder
				endelse

				options, 'th'+sc_thm+'_fgs_gsm*',colors=[2,4,6],labels=['B!dx','B!dy','B!dz'], ytitle = 'TH'+strupcase(sc_thm)+' B', ysubtitle = '!c[nT]', labflag = 1
				options, 'th'+sc_thm+'_ptix_vperp_gsm*',colors=[2,4,6],labels=['V!dx','V!dy','V!dz'], ytitle = 'TH'+strupcase(sc_thm)+' V!dperp!n', ysubtitle = '!c[km/s]', labflag = 1
				options, 'th'+sc_thm+'_ptix_velocity_gsm*',colors=[2,4,6],labels=['V!dx','V!dy','V!dz'], ytitle = 'TH'+strupcase(sc_thm)+' V!di!n', ysubtitle = '!c[km/s]', labflag = 1
				options, 'th'+sc_thm+'_ptix_density*', ytitle = 'TH'+strupcase(sc_thm)+' n!di!n', ysubtitle = '!c[/cc]'

				;; filter data
				if strcmp_or(sc_thm, filter_thm) then begin
					thigh_pass_filter, 'th'+sc_thm+'_fgs_gsm_tclip', minutes_filter*60., backward=0
					;; more filter
					if crazy_filter then begin ;; if contain
						thigh_pass_filter, 'th'+sc_thm+'_fgs_gsm_tclip_hpfilt', minutes_filter*60., backward=0
						thigh_pass_filter, 'th'+sc_thm+'_fgs_gsm_tclip_hpfilt_hpfilt', minutes_filter*60., backward=0, newname = 'th'+sc_thm+'_fgs_gsm_tclip_hpfilt'
					endif
					tsmooth2, 'th'+sc_thm+'_fgs_gsm_tclip_hpfilt', sm_width_b
					options, 'th'+sc_thm+'_fgs_gsm_tclip_hpfilt*', ytitle = 'TH'+strupcase(sc_thm)+' B!cfiltered'
				endif
				if plot_thm_plasma then begin
					options, 'th'+sc_thm+'_ptix_en_eflux', ytitle = 'TH'+strupcase(sc_thm)+'p!u+!n energy', ysubtitle = '!c[eV]', ztitle = 'Eflux [EFU]'
					options, 'th'+sc_thm+'_ptex_en_eflux', ytitle = 'TH'+strupcase(sc_thm)+'e!u-!n energy', ysubtitle = '!c[eV]', ztitle = 'Eflux [EFU]'
				endif
				if strcmp_or(sc_thm, ['b', 'c']) then begin
					;;; compute dynamic pressure
					split_vec, 'th'+sc_thm+'_ptix_velocity_gsm_tclip'
					tinterpol_mxn,'th'+sc_thm+'_ptix_density_tclip','th'+sc_thm+'_ptix_velocity_gsm_tclip_x',newname='th'+sc_thm+'_ptix_density_int'
					get_data, 'th'+sc_thm+'_ptix_velocity_gsm_tclip_x', t, vx
					get_data, 'th'+sc_thm+'_ptix_density_int', t, ni
					store_data, 'th'+sc_thm+'_Pdyn', data={x:t, y:ni2npa*ni*vx^2}
					options, 'th'+sc_thm+'_Pdyn', ytitle = 'TH'+strupcase(sc_thm)+' Pdyn', ysubtitle = '!c[nPa]'
				endif
			endfor
		endif

		;;; transform all locations to RE
		calc,'"th?_state_pos_RE" = "th?_state_pos_tclip"/6374.4' ;; rbsp locations are transformed to RE too

		;;;;;;;;;;;;;;;; load GOES data ;;;;;;;;;;;;;;;;;;;;;;;;
		goes_probes_load = ['']
		if if_this_event and plot_goes then begin
			goes_load_data, trange = trange_load, datatype = 'fgm', probes = probes_goes_all
			goes_load_plasma, trange = trange_load, probes = probes_goes_all
			;;; titles for the flux plot
			if spec eq 1 then begin
				ytitle_type = 'Energy'
				ysubtitle_flux = '[keV]'
			endif else begin
				ytitle_type = 'Flux'
				ysubtitle_flux = '[/(cm!U2!N-s-sr-keV)]'
			endelse
			;; position
			for j = 0, n_elements(probes_goes_all)-1 do begin
				sc_goes = probes_goes_all[j]
				;; transform coordinate system
				cotrans, 'g'+sc_goes+'_pos_gei', 'g'+sc_goes+'_pos_gse', /gei2gse
				cotrans, 'g'+sc_goes+'_pos_gse', 'g'+sc_goes+'_pos_gsm', /gse2gsm
				calc,'"g'+sc_goes+'_pos_gsm_RE" = "g'+sc_goes+'_pos_gsm"/6374.4' ;; rbsp locations are transformed to RE too
				get_data, 'g'+sc_goes+'_pos_gsm_RE', t_goes, loc_goes
				x_goes = loc_goes[*,0]
				y_goes = loc_goes[*,1]
				rho_goes = sqrt(x_goes^2+y_goes^2)
				no_use = where((x_goes lt 0) and (rho_goes gt rho_rb_min), n_goes_within)
				if (n_goes_within gt 0) and tv_exist('g'+sc_goes+'_H_enp_1') then begin
					goes_probes_load = [goes_probes_load, sc_goes]
					enp_matrix_make, 'g'+sc_goes+'_pos_gei'
					; rotate the FGM data from ENP coordinates to GEI coordinates
					tvector_rotate, 'g'+sc_goes+'_pos_gei_enp_mat', 'g'+sc_goes+'_H_enp_1', /invert
					cotrans, 'g'+sc_goes+'_H_enp_1_rot', 'g'+sc_goes+'_H_gse', /gei2gse
					cotrans, 'g'+sc_goes+'_H_gse', 'g'+sc_goes+'_H_gsm', /gse2gsm
					;;;; band pass filter fields
					thigh_pass_filter, 'g'+sc_goes+'_H_gsm', minutes_filter*60., backward=0
					;; more filter
					if crazy_filter then begin ;; if contain
						thigh_pass_filter, 'g'+sc_goes+'_H_gsm_hpfilt', minutes_filter*60., backward=0
						thigh_pass_filter, 'g'+sc_goes+'_H_gsm_hpfilt_hpfilt', minutes_filter*60., backward=0, newname = 'g'+sc_goes+'_H_gsm_hpfilt'
					endif
					tsmooth2, 'g'+sc_goes+'_H_gsm_hpfilt', sm_width_b_goes 
					;;;; copy location for later use
					copy_data, 'g'+sc_goes+'_pos_gsm_RE', 'th'+sc_goes+'_state_pos_RE'
					;;;; write labels
					options, 'g'+sc_goes+'_H_gsm*', ytitle = 'G-'+sc_goes+' H', colors=[2,4,6], labels=['B!dx','B!dy','B!dz'], ysubtitle = '[nT]', labflag = 1
					options, 'g'+sc_goes+'_H_gsm_hpfilt*', ytitle = 'G-'+sc_goes+' H!cfiltered', colors=[2,4,6], labels=['B!dx','B!dy','B!dz'], ysubtitle = '!c[nT]', labflag = 1
					options, 'g'+sc_goes+'_protflux', ytitle = 'G-'+sc_goes+' p!u+!n '+ytitle_type, ysubtitle = ysubtitle_flux, spec = spec
					options, 'g'+sc_goes+'_elecflux', ytitle = 'G-'+sc_goes+' e!u-!n '+ytitle_type, ysubtitle = ysubtitle_flux, spec = spec
					if spec then begin
						ylim, 'g'+sc_goes+'_protflux', 90., 5e5
						ylim, 'g'+sc_goes+'_elecflux', 38., 4500
					endif
				endif
			endfor
			if n_elements(goes_probes_load) gt 1 then goes_probes_load = goes_probes_load[1:*]
		endif

		;;;;;;;;;;;;;;;; load OMNI data ;;;;;;;;;;;;;;;;;;;;;;;;
		if if_this_event and plot_omni then begin
			load_bin_data, trange=trange_load, datatype = 'omni_b_gsm', /tclip, datafolder = imf_folder
			load_bin_data, trange=trange_load, datatype = 'omni_v_gsm', /tclip, datafolder = vsw_folder
			load_bin_data, trange=trange_load, datatype = 'omni_ni', /tclip, datafolder = nisw_folder
			load_bin_data, trange=trange_load, datatype = 'omni_Pdyn', /tclip, datafolder = nisw_folder
			get_data, 'omni_ni', t_omni, ni_omni
			get_data, 'omni_v_gsm', t_omni, v_omni
			store_data, 'omni_Pdyn', data={x:t_omni, y:ni2npa*ni_omni*v_omni[*,0]^2}
		endif

		;;;;;;;;;;;;;;;; load GBO data ;;;;;;;;;;;;;;;;;;;;;;;;;
		gbo_stations_load = ''
		if if_this_event and plot_gbo then begin
			mlat_load = 0.
			mlt_load = 0.
			;;; get the locations of GBO and select which to plot
			for i_stat = 0, n_elements(gbo_stations)-1 do begin
				stat = gbo_stations[i_stat]
				mlat_stat = locations_ll_mag[0, i_stat]
				midnight_stat_str = stations_midnight[i_stat]
				midnight_str_full = '2000/01/01/'+strcompress(midnight_stat_str, /remove)
				t_pos_r_str = time_string(t_pos_r)
				t_pos_r_str_full = '2000/01/01/'+strmid(t_pos_r_str, 11)
				mlt_stat = (time_double(t_pos_r_str_full)-time_double(midnight_str_full))/3600. ;; t_pos_r is the time of first loaded quantity
				i_negative = where(mlt_stat lt 0, n_negative)
				if n_negative gt 0 then mlt_stat[i_negative] = mlt_stat[i_negative]+24.
				mlt_stat_mid = mlt_stat[floor(0.5*n_elements(mlt_stat))]
				;;; tell good latitude
				case mlt_suf of
				'_night': tell_good_mlt = ((mlt_stat gt 18) and (mlt_stat lt 24)) or ((mlt_stat gt 0) and (mlt_stat lt 6))
				'_day': tell_good_mlt = (mlt_stat gt 6) and (mlt_stat lt 18)
				else: tell_good_mlt = replicate(1, n_elements(mlt_stat))
				endcase
				i_good_mlt = where(tell_good_mlt, n_good_mlt)
				;;; determine the spacecraft to load the pi2 data
				if (mlat_stat ge mlat_low) and (mlat_stat le mlat_high) and (n_good_mlt gt 0) then begin
					gbo_stations_load = [gbo_stations_load, stat]
					mlat_load = [mlat_load, mlat_stat]
					mlt_load = [mlt_load, mlt_stat_mid]
					thm_load_pi2, site = stat, trange = trange_load
					options, 'thg_mag_'+stat+'_pi2', ytitle = stat+' B!cmlat='+strmid(strcompress(string(mlat_stat), /remove), 0, 4)+'!cMLT='+strmid(strcompress(string(mlt_stat_mid), /remove), 0, 4), ysubtitle = '[nT]'
				endif
			endfor
			if n_elements(gbo_stations_load) gt 1 then begin
				gbo_stations_load = gbo_stations_load[1:*]
				mlat_load = mlat_load[1:*]
				mlt_load = mlt_load[1:*]
			endif
		endif

		;;;;;;;;;;;;;;;; load RBSP data ;;;;;;;;;;;;;;;;;;;;;;;;
		if (n_elements(rb_probes_load) gt 1) and strcmp(rb_probes_load[0], '') then rb_probes_load = rb_probes_load[1:*]
		if if_this_event and plot_rbsp then begin
			timespan, trange_load[0], trange_load[1]-trange_load[0], /sec
			rbsp_load, trange=trange_load, probe=rb_probes_load, datatype = 'mageis', /tclip, rbsp_folder = rbsp_folder, level = 2
			if ~plot_themis then begin
				rbsp_load, trange=trange_load, probe=rb_probes_load, datatype = 'hope_sa', /tclip, rbsp_folder = rbsp_folder, level = 2, reduce_connect = 'algebra'
				rbsp_load, trange=trange_load, probe=rb_probes_load, datatype = 'rept', /tclip, rbsp_folder = rbsp_folder, level = 2
;				rbsp_load, trange=trange_load, probe=rb_probes_load, datatype = 'efw_density', /tclip, rbsp_folder = rbsp_folder
			endif

			if strcmp(load_e, 'yes') then begin
				for j = 0, n_elements(rb_probes_load)-1 do begin
					case rb_probes_load[j] of
					'r': probe_rb = 'a'
					's': probe_rb = 'b'
					endcase
					catch, err
					if err eq 0 then begin
						rbsp_load_efw, trange=trange_load, probe=rb_probes_load[j], /tclip
                	
						copy_data, 'th'+rb_probes_load[j]+'_efw_esvy_mgse_vxb_removed_coro_removed_spinfit', 'th'+rb_probes_load[j]+'_efw_esvy_mgse'
						rbsp_mgse2gse, 'th'+rb_probes_load[j]+'_efw_esvy_mgse', newname='th'+rb_probes_load[j]+'_efw_esvy_gse', /no_spice_load
						cotrans, 'th'+rb_probes_load[j]+'_efw_esvy_gse', 'th'+rb_probes_load[j]+'_efw_esvy_gsm', /GSE2GSM
					endif else begin
						dprint, !error_state.msg
					endelse
					catch, /cancel
					options, 'th'+rb_probes_load[j]+'_efw_esvy_gsm', ytitle = 'RB-'+strupcase(probe_rb)+' E', ysubtitle = '!c[mV/m]', labels = ['E!dx', 'E!dy', 'E!dz'], labflag = 1, colors = [2,4,6]
				endfor
			endif

			;;;; manage tplot variables
			;; energy channels
			if spec eq 0 then begin
				for j = 0, n_elements(rb_probes_load)-1 do begin
					sc_rb = rb_probes_load[j]
					if tv_exist('th'+sc_rb+'_mageis_pspec_tclip') then begin
						get_data, 'th'+sc_rb+'_mageis_pspec_tclip', data = pspec
						pspec_mageis_labels = strcompress(string(pspec.v(0,*),format='(i)'), /remove)
					endif
					if tv_exist('th'+sc_rb+'_mageis_espec_tclip') then begin
						get_data, 'th'+sc_rb+'_mageis_espec_tclip', data = espec
						espec_mageis_labels = strcompress(string(espec.v(0,*),format='(i)'), /remove)
					endif
					if tv_exist('th'+sc_rb+'_hope_sa_pspec_tclip') then begin
						get_data, 'th'+sc_rb+'_hope_sa_pspec_tclip', data = pspec
						pspec_hope_labels = strcompress(string(pspec.v(0,*),format='(i)'), /remove)
					endif
					if tv_exist('th'+sc_rb+'_hope_sa_espec_tclip') then begin
						get_data, 'th'+sc_rb+'_hope_sa_espec_tclip', data = espec
						espec_hope_labels = strcompress(string(espec.v(0,*),format='(i)'), /remove)
					endif
					if tv_exist('th'+sc_rb+'_rept_pspec_tclip') then begin
						get_data, 'th'+sc_rb+'_rept_pspec_tclip', data = pspec
						pspec_rept_labels = strcompress(string(pspec.v(0,*),format='(i)'), /remove)
					endif
					if tv_exist('th'+sc_rb+'_rept_espec_tclip') then begin
						get_data, 'th'+sc_rb+'_rept_espec_tclip', data = espec
						espec_rept_labels = strcompress(string(espec.v(0,*),format='(i)'), /remove)
					endif
				endfor
			endif else begin
				pspec_mageis_labels = ''
				espec_mageis_labels = ''
				pspec_hope_labels = ''
				espec_hope_labels = ''
				pspec_rept_labels = ''
				espec_rept_labels = ''
			endelse

			;;; title for spectra
			ztitle = 'Flux!c[cm!u-2!ns!u-1!nsr!u-1!nkeV!u-1!n]'
			if spec then begin
				ytitle_p = 'p!u+!n energy'
				ytitle_e = 'e!u'+minus_sign+'!n energy'
				ysubtitle_flux_mageis = '[keV]'
				ysubtitle_flux_hope = '[eV]'
				ysubtitle_flux_rept = '[MeV]'
			endif else begin
				ytitle_p = 'p!u+!n flux'
				ytitle_e = 'e!u'+minus_sign+'!n flux'
				ysubtitle_flux_mageis = '[cm!u-2!ns!u-1!nsr!u-1!nkeV!u-1!n]'
				ysubtitle_flux_hope = '[cm!u-2!ns!u-1!nsr!u-1!nkeV!u-1!n]'
				ysubtitle_flux_rept = '[cm!u-2!ns!u-1!nsr!u-1!nkeV!u-1!n]'
			endelse

			for j = 0, n_elements(rb_probes_load)-1 do begin
				case rb_probes_load[j] of
				'r': probe_rb = 'a'
				's': probe_rb = 'b'
				endcase

				;;;; band pass filter fields
				thigh_pass_filter, 'th'+rb_probes_load[j]+'_fgs_gsm_tclip', minutes_filter*60., backward=0
				thigh_pass_filter, 'th'+rb_probes_load[j]+'_efw_esvy_mgse', minutes_filter*60., backward=0
				thigh_pass_filter, 'th'+rb_probes_load[j]+'_efw_esvy_gsm', minutes_filter*60., backward=0

				;;;; more filter
				if crazy_filter then begin ;; if contain
					thigh_pass_filter, 'th'+rb_probes_load[j]+'_fgs_gsm_tclip_hpfilt', minutes_filter*60., backward=0
					thigh_pass_filter, 'th'+rb_probes_load[j]+'_fgs_gsm_tclip_hpfilt_hpfilt', minutes_filter*60., backward=0, newname = 'th'+rb_probes_load[j]+'_fgs_gsm_tclip_hpfilt'
				endif

				;; remove spin frequency
				tsmooth2, 'th'+rb_probes_load[j]+'_fgs_gsm_tclip_hpfilt', sm_width_b 
				tsmooth2, 'th'+rb_probes_load[j]+'_efw_esvy_mgse_hpfilt', sm_width_e
				tsmooth2, 'th'+rb_probes_load[j]+'_efw_esvy_gsm_hpfilt', sm_width_e

				;;;;; moments
				;combine_spec, 'th'+rb_probes_load[j]+'_mageis_pspec_tclip', 'th'+rb_probes_load[j]+'_hope_sa_pspec_tclip', newname = 'th'+rb_probes_load[j]+'_pspec_combined_tclip', /eV2keV_2nd
				;combine_spec, 'th'+rb_probes_load[j]+'_mageis_espec_tclip', 'th'+rb_probes_load[j]+'_hope_sa_espec_tclip', newname = 'th'+rb_probes_load[j]+'_espec_combined_tclip', /eV2keV_2nd
				;;; all population for pressures
				;compute_moments, 'th'+rb_probes_load[j]+['_pspec_combined_tclip', '_espec_combined_tclip', '_hope_sa_hespec_tclip', '_hope_sa_ospec_tclip'], intypes = 'flux', inunits_energy = ['keV', 'keV', 'eV', 'eV'], inunits_flux = 'keV', energy_min = energy_min, particle_types = ['p', 'e', 'He+', 'O+'], /combine, newname_combine = 'all', rbsp_folder = rbsp_folder
				;;; all ions
				;compute_moments, 'th'+rb_probes_load[j]+['_pspec_combined_tclip', '_hope_sa_hespec_tclip', '_hope_sa_ospec_tclip'], intypes = 'flux', inunits_energy = ['keV', 'eV', 'eV'], inunits_flux = 'keV', energy_min = energy_min, particle_types = ['p', 'He+', 'O+'], /combine, newname_combine = 'ion', rbsp_folder = rbsp_folder
				;;; electron
				;compute_moments, 'th'+rb_probes_load[j]+'_espec_combined_tclip', intype = 'flux', inunits_energy = 'keV', inunits_flux = 'keV', energy_min = energy_min, particle_types = 'e', rbsp_folder = rbsp_folder

				;;; fields
				options, 'th'+rb_probes_load[j]+'_fgs_gsm_*',colors=[2,4,6],labels=['B!dx','B!dy','B!dz'], ytitle = 'RB-'+strupcase(probe_rb)+' B', ysubtitle = '!c[nT]', labflag = 1
				options, 'th'+rb_probes_load[j]+'_efw_density*', ytitle = 'RB-'+strupcase(probe_rb)+' n!dEFW!n'
				;; filtered
				options, 'th'+rb_probes_load[j]+'_fgs_gsm_tclip_hpfilt*', ytitle = 'RB-'+strupcase(probe_rb)+' B!cfiltered'
				options, 'th'+rb_probes_load[j]+'_efw_esvy_gsm', ytitle = 'RB-'+strupcase(probe_rb)+' E filtered', ysubtitle = '!c[mV/m]'

				;;; moments
;				options, 'th'+rb_probes_load[j]+'_ion_density*', ytitle = 'RB-'+strupcase(probe_rb)+' n!di!n', ysubtitle = '!c[cm!u-3!n]'
;				if strcmp(purpose, 'publication') and (i eq 1) then begin
;					ylim, 'th'+rb_probes_load[j]+'_ion_density*', 0.6, 1.1
;				endif
;				options, 'th'+rb_probes_load[j]+'_all_Pall*', ytitle = 'RB-'+strupcase(probe_rb)+' Pressures', ysubtitle = '!c[nPa]'

				;;; spectra
				options,'th'+rb_probes_load[j]+'_mageis_pspec_*', ylog = 1, zlog = 1, ytitle = 'RB-'+strupcase(probe_rb)+' '+ytitle_p, ysubtitle = ysubtitle_flux_mageis, ztitle = ztitle, spec = spec
				options,'th'+rb_probes_load[j]+'_mageis_espec_*', ylog = 1, zlog = 1, ytitle = 'RB-'+strupcase(probe_rb)+' '+ytitle_e, ysubtitle = ysubtitle_flux_mageis, ztitle = ztitle, spec = spec
				options,'th'+rb_probes_load[j]+'_hope_sa_pspec_tclip', ylog = 1, zlog = 1, ytitle = 'RB-'+strupcase(probe_rb)+' '+ytitle_p, ysubtitle = ysubtitle_flux_hope, ztitle = ztitle, spec = spec, labels = pspec_hope_labels 
				options,'th'+rb_probes_load[j]+'_hope_sa_espec_tclip', ylog = 1, zlog = 1, ytitle = 'RB-'+strupcase(probe_rb)+' '+ytitle_e, ysubtitle = ysubtitle_flux_hope, ztitle = ztitle, spec = spec, labels = espec_hope_labels 
				options,'th'+rb_probes_load[j]+'_rept_pspec_tclip', ylog = 1, zlog = 1, ytitle = 'RB-'+strupcase(probe_rb)+' '+ytitle_p, ysubtitle = ysubtitle_flux_rept, ztitle = ztitle, spec = spec, labels = pspec_rept_labels 
				options,'th'+rb_probes_load[j]+'_rept_espec_tclip', ylog = 1, zlog = 1, ytitle = 'RB-'+strupcase(probe_rb)+' '+ytitle_e, ysubtitle = ysubtitle_flux_rept, ztitle = ztitle, spec = spec, labels = espec_rept_labels 
			endfor

			;;; continue doing spectra
			options,'th?_mageis_pspec_tclip', labels = pspec_mageis_labels 
			options,'th?_mageis_espec_tclip', labels = espec_mageis_labels 
			if spec then begin
				ylim,'th?_mageis_pspec_*', 60., 1100. 
				ylim,'th?_mageis_espec_*', 40., 4000.
				ylim,'th?_hope_sa_pspec_tclip', 24., 52000. 
				ylim,'th?_hope_sa_espec_tclip', 14., 52000.
				ylim,'th?_rept_pspec_tclip', 20., 200. 
				ylim,'th?_rept_espec_tclip', 2., 20.
				;if strcmp(purpose, 'publication') then begin
				;	zlim,'th?_mageis_pspec_tclip', 1.1, 1.1e5. 
				;	zlim,'th?_mageis_espec_tclip', 1.1, 1.1e5.
				;endif
			endif
		endif ;; if of ploting RBSP or not

		;;;;; title of the plot
		title = 'DFB '+strcompress(string(i), /remove)

		;;;;;;;;;;;;;;;;;; make location plot ;;;;;;;;;;;;;;;;;;;;;;;;
		if ~strcmp(goes_probes_load[0], '') then probes_loc = [thm_probes, rb_probes_load, goes_probes_load] else probes_loc = [thm_probes, rb_probes_load]
		locations = [0., 0., 0.]
		if strcmp(project, 'yes') then locations_eq = [0., 0., 0.]
		for j = 0, n_elements(probes_loc)-1 do begin
			;;; store locations
			get_data, 'th'+probes_loc[j]+'_state_pos_RE', t, pos_sc
			locations = [[locations], [mean(pos_sc, dimension = 1, /nan)]]
			if strcmp(project, 'yes') then begin
				;;;; trace location to equator, out put is 'th?_efoot'
				;; T96
				ttrace2equator96, 'th'+probes_loc[j]+'_state_pos_RE'
				get_data, 'th'+probes_loc[j]+'_efoot', t, pos_sc_eq
				locations_eq = [[locations_eq], [mean(pos_sc_eq, dimension = 1, /nan)]]
			endif
		endfor
		locations = locations[*, 1:*]
		if strcmp(project, 'yes') then locations_eq = locations_eq[*, 1:*]

		probe_locations, 0.5*total(trange_load), probes = probes_loc, locations = locations, planes_plot = 'XY', title = title+' Locations', position = positions[*,0], xrange = xrange, yrange = yrange, /earth, /geo, /write_probe, size_sc = 0.5
		if strcmp(project, 'yes') then probe_locations, 0.5*total(trange_load), probes = probes_loc, locations = locations_eq, planes_plot = 'XY', title = '', position = positions[*,1], /noerase, xtitle = 'X!deq!n [R!dE!n]', ytitle = 'Y!deq!n [R!dE!n]', xrange = xrange, yrange = yrange, /earth, /geo, /write_probe, size_sc = 0.5
		makepng, pic_folder+'/event_'+strcompress(string(i), /remove)+'_location'
	
		;;;;;;;;;;;;;;;;;; make data plot ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		tvnames_plot = ['']

		;;;;;;; OMNI
		if plot_omni then begin
			tvnames_plot_omni = 'omni_'+['b_gsm', 'v_gsm', 'ni', 'Pdyn']
			tvnames_plot = [tvnames_plot, tvnames_plot_omni]
			;;; plot
			if seperate_plot then begin
				tplot, tvnames_plot_omni, trange = trange_show, title = title+' OMNI Data'
				timebar, times[i], line = 1
				makepng, pic_folder+'/event_'+strcompress(string(i), /remove)+'_omni'
			endif
		endif

		;;;;;;; THEMIS
		if plot_themis then begin
			tvnames_plot_thm = ['']
			for j = 0, n_elements(thm_probes)-1 do begin
				sc_thm = thm_probes[j]
				if strcmp_or(sc_thm, filter_thm) then b_name = ['th'+sc_thm+'_fgs_gsm_tclip', 'th'+sc_thm+'_fgs_gsm_tclip_hpfilt_sm'] else b_name = 'th'+sc_thm+'_fgs_gsm_tclip'
				if strcmp_or(sc_thm, ['b','c']) then begin
					tvnames_plot_thm = [tvnames_plot_thm, b_name, 'th'+sc_thm+'_ptix_velocity_gsm_tclip', 'th'+sc_thm+'_ptix_density_tclip', 'th'+sc_thm+'_Pdyn']
				endif else begin
					tvnames_plot_thm = [tvnames_plot_thm, b_name, 'th'+sc_thm+'_ptix_vperp_gsm_tclip', 'th'+sc_thm+'_ptix_en_eflux', 'th'+sc_thm+'_ptex_en_eflux']
				endelse
			endfor
			tvnames_plot_thm = tvnames_plot_thm[1:*]
			tvnames_plot = [tvnames_plot, tvnames_plot_thm]
			;;; plot
			if seperate_plot then begin
				tplot, tvnames_plot_thm, trange = trange_show, title = title+' THEMIS Data'
				timebar, times[i], line = 1
				makepng, pic_folder+'/event_'+strcompress(string(i), /remove)+'_thm'
			endif
		endif

		;;;;;;; RBSP
		if plot_rbsp then begin
			if ~keyword_set(rb_probes_plot) then rb_probes_plot = rb_probes_load
			tvnames_plot_rb = ['']
			for j = 0, n_elements(rb_probes_plot)-1 do begin
				sc_rb = rb_probes_plot[j]
				if plot_rbsp then begin
					tvnames_plot_rb = [tvnames_plot_rb, 'th'+sc_rb+'_fgs_gsm_tclip', 'th'+sc_rb+'_fgs_gsm_tclip_hpfilt_sm', 'th'+sc_rb+'_efw_esvy_mgse', 'th'+sc_rb+'_efw_esvy_mgse_hpfilt_sm', 'th'+sc_rb+'_rept_pspec_tclip', 'th'+sc_rb+'_mageis_pspec_tclip', 'th'+sc_rb+'_hope_sa_pspec_tclip', 'th'+sc_rb+'_rept_espec_tclip', 'th'+sc_rb+'_mageis_espec_tclip', 'th'+sc_rb+'_hope_sa_pspec_tclip']
					;; without original field
;					tvnames_plot_rb = [tvnames_plot_rb, 'th'+sc_rb+'_fgs_gsm_tclip_hpfilt_sm', 'th'+sc_rb+'_efw_esvy_mgse_hpfilt_sm', 'th'+sc_rb+'_rept_pspec_tclip', 'th'+sc_rb+'_mageis_pspec_tclip', 'th'+sc_rb+'_hope_sa_pspec_tclip', 'th'+sc_rb+'_rept_espec_tclip', 'th'+sc_rb+'_mageis_espec_tclip', 'th'+sc_rb+'_hope_sa_pspec_tclip']
				endif
			endfor
			tvnames_plot_rb = tvnames_plot_rb[1:*]
			tvnames_plot = [tvnames_plot, tvnames_plot_rb]
			;;; plot
			if seperate_plot then begin
				tplot, tvnames_plot_rb, trange = trange_show, title = title+' RBSP Data'
				timebar, times[i], line = 1
				makepng, pic_folder+'/event_'+strcompress(string(i), /remove)+'_rb'
			endif
		endif

		;;;;;;; GOES
		if plot_goes and (~strcmp(goes_probes_load[0], '')) then begin
			tvnames_plot_goes = ['']
			for j = 0, n_elements(goes_probes_load)-1 do begin
				sc_goes = goes_probes_load[j]
				tvnames_plot_goes = [tvnames_plot_goes, 'g'+sc_goes+'_H_gsm', 'g'+sc_goes+'_H_gsm_hpfilt_sm', 'g'+sc_goes+'_protflux', 'g'+sc_goes+'_elecflux']
			endfor
			tvnames_plot_goes = tvnames_plot_goes[1:*]
			tvnames_plot = [tvnames_plot, tvnames_plot_goes]
			;;; plot
			if seperate_plot then begin
				tplot, tvnames_plot_goes, trange = trange_show, title = title+' GOES Data'
				timebar, times[i], line = 1
				makepng, pic_folder+'/event_'+strcompress(string(i), /remove)+'_goes'
			endif
		endif

		;;;;;;; GBO
		if plot_gbo and (~strcmp(gbo_stations_load[0], '')) then begin
			tvnames_plot_gbo = ['']
			for j = 0, n_elements(gbo_stations_load)-1 do begin
				station = gbo_stations_load[j]
				tvnames_plot_gbo = [tvnames_plot_gbo, 'thg_mag_'+station+'_pi2']
			endfor
			tvnames_plot_gbo = tvnames_plot_gbo[1:*]
			;;; sort the tplot names
			case sort_suf of
			'_slat': i_sort = sort(mlat_load)
			'_smlt': i_sort = sort(mlt_load)
			endcase
			tvnames_plot_gbo = tvnames_plot_gbo[i_sort]
			tvnames_plot = [tvnames_plot, tvnames_plot_gbo]
			;;; plot
			if seperate_plot then begin
				tplot, tvnames_plot_gbo, trange = trange_show, title = title+' GBO Data'
				timebar, times[i], line = 1
				makepng, pic_folder+'/event_'+strcompress(string(i), /remove)+'_gbo'+mlat_suf+mlt_suf+sort_suf
			endif
		endif

		;;;; plot all
		if ~seperate_plot then begin
			tvnames_plot = tvnames_plot[1:*]
			tplot, tvnames_plot, trange = trange_show, title = title+' Data'
			timebar, times[i], line = 1
			makepng, pic_folder+'/event_'+strcompress(string(i), /remove)
		endif
	endif ;; if to plot this event
;	stop
endfor ;; for of i, element on multi_list

stop
end
