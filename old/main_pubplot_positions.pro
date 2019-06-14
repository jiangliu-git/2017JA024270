pro main_pubplot_positions
;;; plot the positions of the probes
thm_init
computer = 'I:'
@folders
del_data, '*'

;; Define time ranges and and onset time
kp = double(1.+2./3.) ;; looked up from kp list
time = '2013 9 30 1 21'
trange_load  = time_double(time)+[-150., 150.]


;;;; THEMIS
probes = ['a', 'd', 'e']
load_bin_data, trange=trange_load, probe=probes, datatype = 'pos', /tclip, datafolder = pos_folder

;;; RBSP
probes_rb = ['r', 's']
rbsp_load, trange=trange_load, probe=probes_rb, datatype = 'fgs', /tclip, rbsp_folder = rbsp_folder

;;; GOES
sc = '13'
trange_goes = trange_load+[-300., 300.] ;; goes resolution is low
goes_load_data, trange = trange_goes, datatype = 'fgm', probes = sc, tplotnames = tplotnames

cotrans, 'g'+sc+'_pos_gei', 'g'+sc+'_pos_gse', /gei2gse
cotrans, 'g'+sc+'_pos_gse', 'g'+sc+'_pos_gsm', /gse2gsm
copy_data, 'g'+sc+'_pos_gsm', 'th'+sc+'_state_pos'
time_clip, 'th'+sc+'_state_pos', trange_goes[0], trange_goes[1]

;;; change km to RE
calc,'"th*_state_pos_RE" = "th*_state_pos_tclip"/6374.4'

probes_loc = [probes, probes_rb, sc] ;; a, d, e, r, s, 13
locations = [0., 0., 0.]
locations_eq = [0., 0., 0.]
locations_iono = [0., 0., 0.]
for j = 0, n_elements(probes_loc)-1 do begin
	;;; store locations
	get_data, 'th'+probes_loc[j]+'_state_pos_RE', t, pos_sc
	locations = [[locations], [mean(pos_sc, dimension = 1, /nan)]]

	;;;;;;;;; trace location to equator, out put is 'th?_efoot'
	;; T89
	ttrace2equator, 'th'+probes_loc[j]+'_state_pos_RE', newname='th'+probes_loc[j]+'_efoot', external_model='t89', par=kp
	get_data, 'th'+probes_loc[j]+'_efoot', t, pos_sc_eq
	locations_eq = [[locations_eq], [mean(pos_sc_eq, dimension = 1, /nan)]]

	;;;;;;;;; trace location to ionosphere, output is 'th?_ifoot'
	;; T89
	ttrace2iono, 'th'+probes_loc[j]+'_state_pos_RE', newname='th'+probes_loc[j]+'_ifoot', external_model='t89', par=kp, out_coord = 'geo'
	get_data, 'th'+probes_loc[j]+'_ifoot', t, pos_sc_iono
	locations_iono = [[locations_iono], [mean(pos_sc_iono, dimension = 1, /nan)]]
endfor
locations = locations[*, 1:*]
locations_eq = locations_eq[*, 1:*]
locations_iono = locations_iono[*, 1:*]
locations_iono_ll = geo2latlong(locations_iono, /degrees, /full_circle)
locations_sc_ll = locations_iono_ll[1:*,*]


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; make space location plot ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;title = 'positions' ;; location plot
title = 'cartoon' ;; discussion plot

;;; ranges of plot
xrange = [0, -10.9]
yrange = [8.9, -1.4]
zrange = [-2.9, 1.9]

;;; parameters for plot
if strcmp(title, 'positions') then begin
	xsize = 6
	ysize = 4
	left_margin = 0.15
	right_margin = 0.15
	top_margin = 0.05
	bot_margin = 0.07
	space_vert = 0.01
	n_panels = 2
	positions = panel_positions([1, n_panels], lr_margins = [left_margin, right_margin], bt_margins = [bot_margin, top_margin], space = [0, space_vert])
	offset = 0.15838 ;;; change this to match up and down, bigger: upper one bigger
	positions[1,0] = positions[1,0]-offset
	positions[3,1] = positions[3,1]-offset
endif else begin
	xsize = 5.5
	ysize = 3.2
	positions = [0.1, 0.1, 0.1, 0.1]
	normal_probes = ['d', 'e']
	normals = [[0.10732344, -0.97118409, -0.21279835], [0.33533106, 0.70478590, -0.62516391]]
endelse

x_abc = 0.85
y_abc = 0.82
case title of
'positions': abc = ['d', 'e']
else: abc = ['a', 'b']
endcase

popen, pic_folder+'/'+title
print_options,xsize=xsize, ysize=ysize


;;;;;;; location plot
if strcmp(title, 'positions') then begin
	;;;; XY plane
	probe_locations, 0.5*total(trange_load), probes = probes_loc, locations = locations, planes_plot = 'XY', title = 'Probe Locations', position = positions[*,0], xrange = xrange, yrange = yrange, xtitle = ' ', /earth, /geo, /write_probe, size_sc = 0.65, reference_axis = 'X', /no_xticks
	;; add equatorial projections of the probes
	probe_locations, 0.5*total(trange_load), probes = probes_loc, locations = locations_eq, planes_plot = 'XY', size_sc = 0.3, /add, reference_axis = 'X'
	xyouts, !x.crange[0]+x_abc*(!x.crange[1]-!x.crange[0]), !y.crange[0]+y_abc*(!y.crange[1]-!y.crange[0]), '('+abc[0]+')'
	;;;; XZ plane
	probe_locations, 0.5*total(trange_load), probes = probes_loc, locations = locations, planes_plot = 'XZ', title = '', position = positions[*,1], /noerase, xrange = xrange, zrange = zrange, /earth, /write_probe, size_sc = 0.65, reference_axis = 'X'
	xyouts, !x.crange[0]+x_abc*(!x.crange[1]-!x.crange[0]), !y.crange[0]+y_abc*(!y.crange[1]-!y.crange[0]), '('+abc[1]+')'
endif



;;;;;;; cartoon plot
theta = findgen(30)/29.*!pi-0.5*!pi ;; used to draw DF
if strcmp(title, 'cartoon') then begin
	;;;;; compute the DF radius
	dfront = locations[1,1]-locations[1,2]
	alpha_d = atan2(normals[1,0], normals[0,0])
	alpha_e = atan2(normals[1,1], normals[0,1])
	r_this = abs(dfront/(sin(alpha_d)-sin(alpha_e)))
	;;;;;; DF center point
	Od = [locations[0,1]-r_this*cos(alpha_d), locations[1,1]-r_this*sin(alpha_d)]
	Oe = [locations[0,2]-r_this*cos(alpha_e), locations[1,2]-r_this*sin(alpha_e)]
	x_Od = Od(0)+r_this*cos(theta)
	y_Od = Od(1)+r_this*sin(theta)
	x_Oe = Oe(0)+r_this*cos(theta)
	y_Oe = Oe(1)+r_this*sin(theta)
	;;; average DF location
	;Oave = 0.5*(Od+Oe)
	Oave = Od-0.02
	x_ODF = Oave(0)+r_this*cos(theta)
	y_ODF = Oave(1)+r_this*sin(theta)

	positions = panel_positions([1, n_panels], lr_margins = [left_margin, right_margin], bt_margins = [top_margin, bot_margin], space = [0., vspace], ratio_vert = abs([xrange[1]-xrange[0], zrange[1]-zrange[0]]))

	;;;; XY plane only
	;probe_locations, 0.5*total(trange_load), probes = probes_loc, locations = locations, planes_plot = 'XY', title = '', xrange = xrange, yrange = yrange, /earth, /geo, /write_probe, size_sc = 0.65, reference_axis = 'X', normals = normals, normal_probes = ['d','e'], length_n_factor = 0.07, length_tan_factor = 0.05
	;;;; draw axis only
	probe_locations, 0.5*total(trange_load), probes = probes_loc, locations = locations, planes_plot = 'XY', title = '', xrange = xrange, yrange = yrange, normals = normals, normal_probes = ['d','e'], length_n_factor = 0.07, length_tan_factor = 0.05, /nodata

	;;; draw DF
	df_thick = 2.
	oplot, x_Od, y_Od, line = 1, thick = df_thick
	oplot, x_Oe, y_Oe, line = 1, thick = df_thick
	;oplot, x_ODF, y_ODF, line = 1

	;;; draw the arrow between DFs
	x_df_head_d = max(x_Od)
	y_df_head_d = 0.5*(max(y_Od)+min(y_Od))
	x_df_head_e = max(x_Oe)
	y_df_head_e = 0.5*(max(y_Oe)+min(y_Oe))
	x_distance_de = x_df_head_d-x_df_head_e
	x_arrow_start = x_df_head_e+0.2*x_distance_de
	x_arrow_end = x_df_head_d-0.2*x_distance_de-0.1
	arrow_solid, x_arrow_start, y_df_head_e-0.2, x_arrow_end, y_df_head_d+0.2, hlength = 0.2, hangle = 80., in_thick = 0.8, in_color = 255, out_color = 0, out_thick = 1.1
	xyouts, 0.5*(x_arrow_start+x_arrow_end), 0.5*(y_df_head_d+y_df_head_e)+0.1, '~1 min', align = 0.5, charsize = 0.5, orientation = 18.

	;;; draw wave cartoons
	length_period = 1.
	wave_thick = 1.7
	wave_cartoon, -6.5, 2.5, -1, 1.5, n_periods = 5, length = 5*length_period, width = 0.6, arrow_length = 0.1, thick = wave_thick ;; to RBB
	wave_cartoon, -6.5, 4.5, -0.3, 5.7, n_periods = 6, length = 6*length_period, width = 0.6, arrow_length = 0.1, thick = wave_thick ;; to G13 and RBA
	wave_cartoon, -10.5, 5.2, -7.5, 8.3, n_periods = 5, length = 5*length_period, width = 0.6, arrow_length = 0.1, thick = wave_thick ;; to P5

	;;; draw probes
	probe_locations, 0.5*total(trange_load), probes = probes_loc, locations = locations, planes_plot = 'XY', /earth, /geo, /write_probe, size_sc = 0.65, reference_axis = 'X', normals = normals, normal_probes = ['d','e'], length_n_factor = 0.07, length_tan_factor = 0.05, /add

	;;; write abc
	;xyouts, !x.crange[0]+x_abc*(!x.crange[1]-!x.crange[0]), !y.crange[0]+y_abc*(!y.crange[1]-!y.crange[0]), '(a)'
endif
pclose
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;; make ground location plot ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;xsize = 6.01
;ysize = 4.6
;
;;xsize = 12
;;ysize = 9
;
;gbos = strupcase(['kuuj', 'fhb', 'iqa', 'tpas', 'drby', 'snkq'])
;thm_gmag_stations, gbos_all, locations_all, magnetic = locations_mag_all
;locations_gbos = fltarr(2, n_elements(gbos))
;locations_mag_gbos = fltarr(2, n_elements(gbos))
;
;for i = 0, n_elements(gbos)-1 do begin
;	i_match = where(strcmp(gbos[i], gbos_all), n_match)
;	if n_match gt 0 then begin
;		locations_gbos[*,i] = locations_all[*, i_match]
;		locations_mag_gbos[*,i] = locations_mag_all[*, i_match]
;	endif else begin
;		locations_gbos[*,i] = !values.f_nan+[0.,0.]
;		locations_mag_gbos[*,i] = !values.f_nan+[0.,0.]
;	endelse
;endfor
;
;locations_scgbo = [[locations_sc_ll], [locations_gbos]]
;
;popen, pic_folder+'/map'
;print_options,xsize=xsize, ysize=ysize
;
;thm_map_set, central_lon=(max(locations_scgbo[1,*])+min(locations_scgbo[1,*]))/2, central_lat=(max(locations_scgbo[0,*])+min(locations_scgbo[0,*]))/2, /noerase, /no_color, color_continent = 255, color_background = 85, scale = 2.4e7
;
;;;;;;;; add latitude and longitude
;;;;; add latitude
;thm_map_add, invariant_lats=[50, 55, 60, 65, 70, 75], invariant_color=0, invariant_thick=0.3, invariant_linestyle = 1
;xyouts, 0.25, 0.26, '60!uo!n', /normal
;xyouts, 0.34, 0.88, '75!uo!n', /normal
;;;;; add MLT
;map_add_mlt, time, [23, 22, 21, 20, 19, 18], thick = 0.3, line = 1
;xyouts, 0.09, 0.01, '19MLT', /normal
;xyouts, 0.664, 0.01, '21MLT', /normal
;
;;;;;;;; add satellites ;;;;;;;;;;;;;;
;;; make circle symbol
;npts = 49
;angle = findgen(npts)/(npts-1)*2*!pi
;usersym, cos(angle), sin(angle), /fill ;; circle
;;; add satellites
;for i = 0, n_elements(probes_loc)-1 do begin
;	oplot, [locations_sc_ll[1,i], !values.f_nan], [locations_sc_ll[0,i], !values.f_nan], psym = 8, symsize = 2., color = thm_probe_color(probes_loc[i])
;	color_sc = thm_probe_color(probes_loc[i])
;	;; decide the increment for text
;	case probes_loc[i] of
;	'd': begin
;		inc_x = -3.2
;		inc_y = -1.1
;		end
;	'13': begin
;		inc_x = -4.2
;		inc_y = 0.5
;		color_sc = 5
;		end
;	else: begin
;		inc_x = 0.5
;		inc_y = 0.6
;		end
;	endcase
;	xyouts, locations_sc_ll[1,i]+inc_x, locations_sc_ll[0,i]+inc_y, thm_probe_color(probes_loc[i], /number), color = color_sc, charsize = 1.2
;endfor
;
;;;;;;;; add GBOs ;;;;;;;;;;
;usersym, [-1,0,1,0], [0,-1,0,1], /fill ;; solid diamond
;;; add gbos
;for i = 0, n_elements(gbos)-1 do begin
;	oplot, [locations_gbos[1,i], !values.f_nan], [locations_gbos[0,i], !values.f_nan], psym = 8, symsize = 2.2
;	case gbos[i] of
;	'TPAS': begin
;		inc_x = 0.85
;		inc_y = -1.6
;		end
;	'IQA': begin
;		inc_x = 0.
;		inc_y = 0.7
;		end
;	'FHB': begin
;		inc_x = 0.5
;		inc_y = 0.5
;		end
;	'DRBY': begin
;		inc_x = -2.7
;		inc_y = -1.
;		end
;	else: begin
;		inc_x = 3.3
;		inc_y = -1.1 ;; SNKQ and KUUJ
;		end
;	endcase
;	xyouts, locations_gbos[1,i]+inc_x, locations_gbos[0,i]+inc_y, gbos[i], align = 0.5, charsize = 1.2
;endfor
;
;xyouts, 0.08, 0.9, '(a)', charsize = 1.8, /normal
;
;;;; plot an overall box
;plots, [0., 1., 1., 0., 0.], [0., 0., 1., 1., 0.], /normal
;
;pclose
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

stop
end
