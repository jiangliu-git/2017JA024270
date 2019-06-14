pro aaron_part_cmb_test

probe = 'a'
trange = '2013-09-30/' + ['1:00','1:50']

esa_datatype = 'peir'
sst_datatype = 'psif'

combined = thm_part_combine(probe=probe, trange=trange, $
                            esa_datatype=esa_datatype, sst_datatype=sst_datatype, $
                            orig_esa=esa, orig_sst=sst) 


;interpolated E spec
thm_part_products, dist_array=combined, outputs='energy', tplotnames=names

;original data
thm_part_products, dist_array=esa, outputs='energy', tplotnames=esa_names
thm_part_products, dist_array=sst, outputs='energy', tplotnames=sst_names


;create pseudo-var with original data
store_data, names[0]+'_orig', data= [esa_names, sst_names]


;set z axis for all plots
options, '*_eflux_energy*', zrange=[1e1,1e6], /zlog

;set y axis for both combined plots
options, names[0]+'*', yrange=[6,3e5]

tplot, names[0]+'*'

end 
