function injection_decide, tplot_var, t0=t0_in, secsrange_ahead = secsrange_ahead, secsrange_after = secsrange_after, secs_injtype = secs_injtype, max = max, n_channels = n_channels, consecutive = consecutive, crit_inj = crit_inj, crit_noinj_all = crit_noinj_all, crit_noinj_channels = crit_noinj_channels, type = type, crit_dispersion = crit_dispersion, inj_energy_all = inj_energy_all, ratios = ratio_change, inf_good = inf_good
;; to decide whether the input variable (tplot_var) is a injection
;; t0 is the zero time (e.g., the onset of a DFB) for comparing reference
;; trange_ahead, trange_after: seconds relative to t0, to compute the particle flux values.
;; return value is 1 for ture, 0 for not, -1 for not sure, -2 for no data
;; max: if set, use max value for comparison; if not using avarage value.
;; n_channels: at least this number of channels need to have injection
;; consecutive: if set, channels for injection needs to be consecutive
;; crit_inj: must exceed times of increase for injection
;; crit_noinj_all: must not exceed times of increase for noinjection (none of the channels)
;; crit_noinj_channels: must not exceed times of increase for noinjection (for n_channels); default is disabled
;; crit_dispersion: in seconds, time to tell whether is dispersed injection or dispersionless injection. if not set, data's resolution is used.
;; type: an output, 'none' for no injection or not sure, for injection has "dispersed" and "dispersionless"
;; inj_energy_all: an output, all the energies showing injection signautre, in whatever energy the input variable contains.
;; ratios: an output, the after/before ratio of all channels, has the same number of elememnts as the number of channels
;; inf_good: set this to allow infinity ratio to be also counted increased channel

if ~keyword_set(secsrange_ahead) then secsrange_ahead = [-180., -60.]
if ~keyword_set(secsrange_after) then secsrange_after = [0., 60.]
if ~keyword_set(n_channels) then n_channels = 3
if ~keyword_set(crit_inj) then crit_inj = 3.
if ~keyword_set(crit_noinj_all) then crit_noinj_all = 2. ;; all should not exceed this
if ~keyword_set(crit_noinj_channels) then crit_noinj_channels = 100. ;; n_channels of channels should not exceed this; set to a large number to disable

small_number = -1e-5 ;; used for telling if equal

t0 = time_double(t0_in)
trange_ahead = t0+secsrange_ahead
trange_after = t0+secsrange_after

if keyword_set(secs_injtype) then begin
	trange_injtype = [t0-secs_injtype, t0+secs_injtype]
endif else begin
	trange_injtype = [trange_ahead[0], trange_after[1]]
endelse

if tv_exist(tplot_var) then begin
	get_data, tplot_var, t, flux, energy_values
	;; whether all the ts are within the range of the variable
	if (t[0] gt min([t0, trange_ahead[0], trange_after[0]])) or (t[-1] lt max([t0, trange_ahead[1], trange_after[1]])) then begin
		print, 'INJECTION_DECIDE: The range of the tplot variable is not big enough.'
		type = 'nodata'
		return, -2
	endif
endif else begin
	print, 'INJECTION_DECIDE: '+tplot_var+' does not exist! Returning "no data"!'
	type = 'nodata'
	return, -2
endelse

if ~keyword_set(crit_dispersion) then begin
	crit_dispersion = 13.
endif

;;;;; calculate the before and after fluxes
;; values ahead and after
time_clip, tplot_var, trange_ahead[0], trange_ahead[1], newname = tplot_var+'_ahead'
get_data, tplot_var+'_ahead', t_ahead, flux_ahead
if keyword_set(max) then fv_ahead = max(flux_ahead, dimension = 1, /nan) else fv_ahead = mean(flux_ahead, dimension = 1, /nan)
time_clip, tplot_var, trange_after[0], trange_after[1], newname = tplot_var+'_after'
get_data, tplot_var+'_after', t_after, flux_after
if keyword_set(max) then fv_after = max(flux_after, dimension = 1, /nan) else fv_after = mean(flux_after, dimension = 1, /nan)
;; ratios of channels
ratio_change = fv_after/fv_ahead
if keyword_set(inf_good) then finite_req = 1 else finite_req = finite(ratio_change)
i_inj = where(finite_req and (ratio_change gt crit_inj), n_inj)
i_notnoinj_all = where(finite_req and (ratio_change gt crit_noinj_all), n_notnoinj_all)
i_notnoinj_channels = where(finite_req and (ratio_change gt crit_noinj_channels), n_notnoinj_channels)

;;; tell whether is injection or no injection
; injection
if n_inj ge n_channels then begin
	;; mark the energies showing injection signature
	inj_energy_all = energy_values[0, i_inj]
	;;; determine whether dispersed or dispersionless
	;; derivate the data (values strange at boundaries, so need to get ddt first)
	deriv_data, tplot_var
	time_clip, tplot_var+'_ddt', trange_injtype[0], trange_injtype[1], newname = tplot_var+'_interest_ddt'
	get_data, tplot_var+'_interest_ddt', t_ddt, fluxes_ddt, energies 
	;; determine the lowest n_channels of channels
	fluxes_ddt_inj = fluxes_ddt[*,i_inj]
	energies_inj = energies[*,i_inj]
	i_sort = sort(energies_inj[0,*])
	fluxes_ddt_compare = fluxes_ddt_inj[*, i_sort[0:n_channels-1]]
	;; the max derivative locations
	max_ddt = fltarr(n_channels)
	i_max_ddt = intarr(n_channels)
	for i = 0, n_channels-1 do begin
		max_ddt_this = max(fluxes_ddt_compare[*,i], i_max_this)
		i_max_ddt[i] = i_max_this
		max_ddt[i] = max_ddt_this
	endfor
	;;;; must contain no NaN to decide type
	if finite(total(max_ddt)) then begin
		t_max_ddt = t_ddt[i_max_ddt]
		;; the first point that exceed a ratio of max ddt
		i_first_ex = intarr(n_channels)
		;;;;;; rate_max: 
		;;;; use inj_decide range:
		;; good: 1;0.95[p15/1;e13/1]
		;; bad: 0.7[p14/3;e9/1], 0.75[p14/2;e8/1], 0.85[p15/1;e11/2], 0.9[p15/1;e13/2], 0.8[p14/1;e7/2], 0.6, 0.65, 0.5, 0.3
		rate_max = 1. ;; to be the same as max
		for i = 0, n_channels-1 do begin
			i_ex = where(fluxes_ddt_compare[*,i]-rate_max*max_ddt[i] ge small_number, n_ex)
			if n_ex gt 0 then begin
				i_first_ex[i] = i_ex[0]
			endif else begin
				print, 'INJECTION_DECIDE: this cannot happen!'
				stop
			endelse
		endfor
		t_first_ex = t_ddt[i_first_ex]
		j_last = where(i_first_ex eq n_elements(t_ddt)-1, n_last)
		j_first = where(i_first_ex eq 0, n_first)
		if (max(t_first_ex)-min(t_first_ex) lt crit_dispersion) and (n_last eq 0) and (n_first eq 0) then begin
			type = 'dispersionless'
		endif else begin
			;; if in order high first low last, dispersed; if not, unclear
			diff_tddt = t_first_ex[1:-1]-t_first_ex[0:-2] ;; should be all <0 for dispersed injection
			if (product(diff_tddt le 0.) gt 0) and (n_last le 1) and (n_first le 1) then begin
				type = 'dispersed'
			endif else begin
				;; if in order low first high last, inverse-dispersed.
				if (product(diff_tddt gt 0.) gt 0) and (n_last le 1) and (n_first le 1) then begin
					type = 'inverse-dispersed'
				endif else type = 'unclear'
			endelse
		endelse
	endif else type = 'unclear'

	;;; decide whether is consecutive
	if (n_channels gt 1) and keyword_set(consecutive) then begin
		diff_ind_inj = i_inj[1:*]-i_inj[0:-2]
		for i = 0, n_elements(diff_ind_inj)-1-(n_channels-2) do begin
			window_this = diff_ind_inj[i:i+n_channels-2]
			if product(window_this) eq 1 then begin
				return, 1
			endif else begin
				continue
			endelse
		endfor	
	endif else begin
		;;; no consequtive requirement decide injection directly
		return, 1
	endelse
endif

;;; out put energy channels
inj_energy_all = 0

;;; not injection
type = 'none'
;; no injection
if (n_notnoinj_all eq 0) and (n_notnoinj_channels lt n_channels) then return, 0
;; not sure
return, -1

end
