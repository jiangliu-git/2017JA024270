pro main_df_normal
thm_init
computer = 'I:'
@folders
;;; determine the normal directions of the DFs
probes = ['d', 'e']
trange_plot = ['13 9 30 1 15', '13 9 30 1 25']
onsets = ['13 9 30 1 21 25', '13 9 30 1 20 19']

;;; criterion
;c_mn = 10.
c_mn = 3.

;;; specify try seconds
n_sec_try = 2.
n_try = 10 ;; will try +-n_sec_try*n_try seconds around

tranges_df = '2013-09-30/'+[['01:21:24', '01:21:33'],['01:20:21', '01:20:34']] ;; for binxbout, requirement should be 25 degs. ny from this: -0.32, 0.72.

;;;;; BinxBout ;;;;;;;
;tranges_normal = '2013-09-30/'+[['01:21:19', '01:21:38'],['01:20:18', '01:20:35']] ;; for binxbout, requirement should be 25 degs. ny from this: -0.32, 0.72.
;;;;; MVA ;;;;;;;;;
;trange_n_d = ['01:18:51', '01:21:05'] ;; this gives 32.6 degs
trange_n_d = ['01:21:24', '01:21:45'] ;; this gives 83.7 degs
trange_n_e = ['01:20:03', '01:20:42'] ;; this gives 64.6 degs
tranges_normal = '2013-09-30/'+[[trange_n_d], [trange_n_e]] ;; for MVA, requirements can be 10. ny from this: -0.42, 0.81

del_data, '*'
for i = 0, n_elements(probes)-1 do begin
;for i = 1, n_elements(probes)-1 do begin
	sc = probes[i]
	trange_normal_set = time_double(tranges_normal[*,i])
	trange_df_this = time_double(tranges_df[*,i])
	event_this = [time_string(time_double(onsets[i])), time_string(time_double(onsets[i])+15.), 'm', sc]
	load_bin_data, trange = trange_plot, probe = sc, datatype = 'fgl', datafolder = fgl_folder
	;;; now manually select the range if no range is set.
	tplot, 'th'+sc+'_fgl_gsm', trange = trange_plot

	;;;;; click two points
	;ctime, trange_normal

	;;;;; use pre-set range
	trange_normal = trange_normal_set

	;;;;;; use automatic method to find the largest angle (each step is n_sec_try seconds)
	;angles_all = 0.
	;tranges_all = [0., 0.]
	;dur_df = trange_df_this[1]-trange_df_this[0]
	;t_try_start = trange_df_this[0];-0.1*dur_df
	;t_try_end = trange_df_this[1];+0.1*dur_df
	;for i_try = 0, n_try-1 do begin
	;	t_try_start_this = t_try_start-i_try*n_sec_try
	;	for j_try = 1, n_try-1 do begin
	;		t_try_end_this = t_try_end+j_try*n_sec_try
	;		trange_this = [t_try_start_this, t_try_end_this]
	;		normal_temp = df_normal(event_this, method = 1, bfolder = fgl_folder, tranges_in = trange_this, c_mn = c_mn, c_angle = 25, tranges_used = trange_used, lambda = lambda, /silent);, care = 'no')
	;		if finite(normal_temp[0]) then begin
	;			angle_this = atan2(normal_temp[1], normal_temp[0])*180.;/!pi
	;			angles_all = [angles_all, angle_this]
	;			tranges_all = [[tranges_all], [trange_this]]
	;		endif
	;	endfor
	;endfor
	;angles_all = angles_all[1:*]
	;tranges_all = tranges_all[*, 1:*]
	;;; find the suitable one
	;case sc of
	;'d': no_use = min(angles_all, i_match, /nan)
	;'e': no_use = max(angles_all, i_match, /nan)
	;endcase
	;trange_normal = tranges_all[*, i_match]
		
	normal_this = df_normal(event_this, method = 1, bfolder = fgl_folder, tranges_in = trange_normal, c_mn = c_mn, c_angle = 25, tranges_used = trange_used, lambda = lambda);, care = 'no')
	print, sc+':'
	print, 'Time range for normal direction:'
	print, time_string(trange_normal)
	print, 'Normal directions:'
	print, normal_this
	print, 'Angle from X:'
	print, atan2(normal_this[1], normal_this[0])*180./!pi

	;if keyword_set(lambda) then print, lambda
	timebar, trange_used, line = 1

	stop
endfor

end
