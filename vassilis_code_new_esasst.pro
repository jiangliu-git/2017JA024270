thm_init

tplot_options, 'xmargin', [15,9]

cwdirname='C:\My Documents\ssl\transfers\analysis\20151220_storm_time_Rx'

cwd,cwdirname

;

; Note: tha=P5, the=P4, thd=P3

mysc='e'

thsc='the'

;

; FOR SST CONTAMINATION REMOVAL USE THE VALUES BELOW

; THEY WERE DETERMINED FROM GET_SST_BINS2MASK.PRO

;

  bin_numbers_tha_e = [8,24,33,40,56,57] ; TH-A (P5)

  min_energy_sst_tha_e = 25000. ; TH-A (P5)

  bin_numbers_thd_e = [8,24,40,56] ; TH-D (P3)

  min_energy_sst_thd_e = 35000. ; TH-D (P3)

  bin_numbers_the_e = [8,24,40,56] ; TH-E (P4)

  min_energy_sst_the_e = 35000. ; TH-E (P4)

 

  bin_numbers_tha_i =  [8,24,40,55,56,57] ; TH-A (P5)

  min_energy_sst_tha_i = 25000. ; TH-A (P5)

  bin_numbers_thd_i =  [8,24,40,55,56] ; TH-D (P3)

  min_energy_sst_thd_i = 25000. ; TH-D (P3)

  bin_numbers_the_i =  [8,24,40,55,56] ; TH-E (P4)

  min_energy_sst_the_i = 25000. ; TH-E (P4)

;

bin_numbers_e = bin_numbers_the_e       ; <------------- MAKE SURE TO CHANGE THSC HERE

min_energy_sst_e = min_energy_sst_the_e ; <------------- MAKE SURE TO CHANGE THSC HERE

bin_numbers_i = bin_numbers_the_i       ; <------------- MAKE SURE TO CHANGE THSC HERE

min_energy_sst_i = min_energy_sst_the_i ; <------------- MAKE SURE TO CHANGE THSC HERE

;

timespan,¡¯mytimehere¡¯,6,/hours

trange=mytrangehere

;

Re=6378. ; km - Earth equatorial radius

mu0=4*!PI*1.e-7

pival=!PI

eVpercc_to_nPa=0.1602/1000.    ; multiply

nTesla2_to_nPa=0.01/25.132741  ; multiply

;

onearray=[1,1,1]

xprojarray=[[1,0,0],[0,0,0],[0,0,0]]

zprojarray=[[0,0,0],[0,0,0],[0,0,1]]

xunit=[1,0,0]

yunit=[0,1,0]

zunit=[0,0,1]

identarray=[[1,0,0],[0,1,0],[0,0,1]]

;

thm_load_state,probe=mysc,/get_supp,coord='gsm'

thm_load_fit,probe=mysc,coord='dsl',suff='_dsl'

thm_load_pxxm_pot4esa, probe=mysc, trange=trange

calc," 'th?_pxxm_pot_1'=1.1*'th?_pxxm_pot_0'+3. "; shift up 3eV, remember calc works with ? too !

;

; ELECTRONS

combined_e = thm_part_combine(probe=mysc, trange=trange, $

esa_datatype='peef', sst_datatype='psef', orig_esa=esa_e, orig_sst=sst_e, $

/esa_bgnd_remove, bgnd_type='anode', bgnd_npoints=3, bgnd_scale=1.5, $

sst_sun_bins=bin_numbers_e, sst_min_energy=min_energy_sst_e)

;

thm_part_products, dist_array=combined_e, outputs='moments energy', $

mag_name=thsc+'_fgs_dsl', sc_pot_name = thsc+'_pxxm_pot_1'

;

thm_part_products, dist_array=combined_e, outputs='phi', $

theta=[-45,45],energy=[30000,60000], $

mag_name=thsc+'_fgs_dsl', sc_pot_name = thsc+'_pxxm_pot_1'

;

thm_part_products, dist_array=esa_e, outputs='moments energy', $

mag_name=thsc+'_fgs_dsl', sc_pot_name = thsc+'_pxxm_pot_1'

;

; IONS

combined_i = thm_part_combine(probe=mysc, trange=trange, $

esa_datatype='peir', sst_datatype='psif', orig_esa=esa_i, orig_sst=sst_i, $

/esa_bgnd_remove, bgnd_type='anode', bgnd_npoints=3, bgnd_scale=1.5, $

sst_sun_bins=bin_numbers_i, sst_min_energy=min_energy_sst_i)

;

thm_part_products, dist_array=combined_i, outputs='moments energy', $

mag_name=thsc+'_fgs_dsl', sc_pot_name = thsc+'_pxxm_pot_1'

;

thm_part_products, dist_array=combined_i, outputs='phi', $

theta=[-45,45],energy=[30000,60000], $

mag_name=thsc+'_fgs_dsl', sc_pot_name = thsc+'_pxxm_pot_1'

;

thm_part_products, dist_array=esa_i, outputs='moments energy', $

mag_name=thsc+'_fgs_dsl', sc_pot_name = thsc+'_pxxm_pot_1'

;

thm_part_products, dist_array=sst_i, outputs='moments energy', $

mag_name=thsc+'_fgs_dsl', sc_pot_name = thsc+'_pxxm_pot_1'

;

store_data,thsc+'_pteffpot_eflux_energy',data=thsc+'_pteff_eflux_energy '+thsc+'_pxxm_pot_1'

ylim,thsc+'_pteffpot_eflux_energy',7.,1.e6,1

zlim,thsc+'_pteffpot_eflux_energy',1.e1,2.e8,1

;

store_data,thsc+'_ptxxx_density',data=thsc+'_ptirf_density '+thsc+'_pteff_density'

options,'th?_ptxxx_density','colors',[2,6]

;

ylim,'th?_ptirf_velocity',0,0,0

ylim,'th?_pt??f*eflux_energy',4.,1.e6,1

zlim,'th?_ptirf_eflux_*',1.e1,1.e7,1

zlim,'th?_pteff_eflux_*',1.e1,2.e8,1

;

tplot,[thsc+'_'+['ptxxx_density','ptirf_velocity','ptirf_ptens','ptirf_eflux_phi', $

'ptirf_eflux_energy','pteff_eflux_phi','pteffpot_eflux_energy']]

;

; HERE SAVE THE MOMENTS AND SPECTRA DATA, ONCE YOU HAVE CREATED THEM FOR ALL 3 VARIABLES

;

tplot_save,filename='get_plasma'
