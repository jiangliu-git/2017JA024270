pro main_test_goes
goes_init
goes_lib
;
date = '2013 9 30 '
trange = date+['0 1', '1 39']

probe = '13'
suffix = ''

del_data, '*'

timespan, time_double(trange[0]), time_double(trange[1])-time_double(trange[0]), /sec

datatype = 'maged'
;datatype = 'magpd'

;load magnetic field
goes_load_data, probe=probe, datatype='fgm', trange=trange, /noeph

;load particle data
goes_load_data, probe=probe, datatype=datatype, trange=trange, /noeph
goes_part_products, probe=probe, datatype=datatype, trange=trange, output='pa gyro', g_interpolate=1, uncorrected=1;,  energy = 1000.*[250, 300]
 
tplot, 'g'+probe+'_'+datatype+'_dtc_uncor_flux_pa'
get_data, 'g'+probe+'_'+datatype+'_dtc_uncor_flux_pa', data = data
pm, data.y


;goes_lib
;goes_load_data, trange=trange, datatype='fgm',probes='13';, /avg_1m
;goes_pitch_angles, 'g13_H_enp_1', 'g13_HT_1', prefix = 'g13'
;;options, 'g15_pitch_angles', spec = 1
;tplot, 'g13_pitch_angles', trange = trange


stop
end
