
prb='13'
prb='15'
timespan,'2015-03-16',6,/days
trange=['2015-03-16', '2015-03-22']

;timespan,'2015-02-02',1,/days
;trange=['2015-02-02', '2015-02-03']
kyoto_load_dst
thm_make_AE
kyoto_load_ae,trange=trange,datatype='all'
  ; start by loading some 1-m averaged magnetometer data for GOES-15. Note that you can see the
  ; tplot variables created by the load routine in 'tplotnames'

  goes_load_data, trange=trange, datatype='fgm', probes=prb, /avg_1m, tplotnames = tplotnames

  ; we can transform the FGM data into other coordinate systems as well,
  ; first by making a transformation matrix, using the position data loaded from the load routine.
  ; the new transformation matrix should be named 'g15_pos_gei_enp_mat'

  enp_matrix_make, 'g'+prb+'_pos_gei'

  ; rotate the FGM data from ENP coordinates to GEI coordinates

  tvector_rotate, 'g'+prb+'_pos_gei_enp_mat', 'g'+prb+'_H_enp_1', /invert

  ; that rotation gives a tplot variable with the horrible name 'g15_b_enp_rot', we can copy

  ; it to a new variable with a better name, 'g15_H_gei'

  copy_data, 'g'+prb+'_H_enp_1_rot', 'g'+prb+'_H_gei'

  thm_cotrans, 'g'+prb+'_H', in_suffix='_gei', in_coord='gei', out_suffix='_gsm', out_coord='gsm'

  thm_cotrans, 'g'+prb+'_pos', in_suffix='_gei', in_coord='gei', out_suffix='_gsm', out_coord='gsm'

  ; and now we can set the labels and title appropriately

  options,'g'+prb+'_H_gsm',labels=['B!dX','B!dY','B!dZ'],ytitle='B-GSM', thick = 1.5

  ylim, 'g'+prb+'_H_gsm',-10.,150.
  tplot,'g'+prb+'_H_gsm'

  
  
  
   ; Load GOES MagED data

 ; goes_load_data, trange=trange, datatype='magpd', probes=prb, /avg_1m, tplotnames = tplotnames, /noephem

  goes_load_data, trange=trange, datatype='maged', probes=prb, /avg_1m, tplotnames = tplotnames, /noephem
  goes_load_data, trange=trange, datatype='epead', probes=prb, /avg_1m 
  
  goes_lib ; compile the GOES library routines

  goes_pitch_angles, 'g'+prb+'_H_enp_1', 'g'+prb+'_HT_1', prefix = 'g'+prb+''

  tplot, 'g'+prb+'_pitch_angles' 
  get_data, 'g'+prb+'_pitch_angles',data=pitch
  
  pitch.y[where(pitch.y gt 90.)]=180.-pitch.y[where(pitch.y gt 90.)]
  
  pitch_order=dblarr(n_elements(pitch.x),9)
  newpitch=dblarr(n_elements(pitch.x),9)
  for i=0,n_elements(pitch.x)-1 do  pitch_order[i,*]=sort(pitch.y[i,*])            ;Put pitch angles in ascending order
    for i=0,n_elements(pitch.x)-1 do newpitch[i,*]=pitch.y[i,[pitch_order[i,*]]]
    
     ;[-cos(alpha_{n+1}) + cos(alpha_{n-1})]/2  (replacing   sin(alpha_n)*[alpha_{n+1} - alpha_{n-1}]/2  )
     
     alpha_nplus1=[[newpitch[*,1:8]],[replicate(90.,n_elements(pitch.x))]]*!PI/180.
     alpha_nmins1=[[replicate(0.,n_elements(pitch.x))],[newpitch[*,0:7]]]*!PI/180.
     dalpha=(-cos(alpha_nplus1) + cos(alpha_nmins1))/2.
  
;  dalpha1=(newpitch[*,2:8]-newpitch[*,0:6])/2.
;  dalpha0=(newpitch[*,1])/2.
;  dalpha9=(90.-newpitch[*,7])/2.
;  dalpha = [[dalpha0],[dalpha1],[dalpha9]]
  
;  dalpha1=pitch.y[*,[pitch_order[*,2:8]]]-pitch.y[*,[pitch_order[*,0:6]]]
;  dalpha0=((pitch.y[*,[pitch_order[*,1]]]-pitch.y[*,[pitch_order[*,0]]]) + 360.-(pitch.y[*,[pitch_order[*,8]]]+270.+pitch.y[*,[pitch_order[*,0]]]))/2.
;  dalpha9=((pitch.y[*,[pitch_order[*,8]]]-pitch.y[*,[pitch_order[*,7]]])+((90.+pitch.y[*,[pitch_order[*,0]]])-pitch.y[*,[pitch_order[*,8]]]))/2.
;  dalpha=[[dalpha0],[dalpha1],[dalpha9]]
   
  ; MAGED

  ; 40 keV
  get_data, 'g'+prb+'_maged_40keV_dtc_cor_flux', data = d

  dp1 = dblarr(n_elements(d.x))

  dcom = {x: d.x, y: dblarr(n_elements(d.x), 5), v: [40.0, 75.0, 150.0, 275.0, 475.0]}
  dcom2 = {x: d.x, y: dblarr(n_elements(d.x), 6), v: [40.0, 75.0, 150.0, 275.0, 475.0, 2000]}

  flux=dblarr(n_elements(d.x),9)
  newy=dblarr(n_elements(d.x),9)
  for i = 0, n_elements(d.x)-1 do newy[i,*]=d.y[i,[pitch_order[i,*]]]
  for i = 0, n_elements(d.x)-1 do flux[i,*] = newy[i,*]*dalpha              
  for i = 0, n_elements(d.x)-1 do dp1[i] = total(flux[i,*])
 
  dcom.y[*,0] = dp1

  ; 75 keV

  get_data, 'g'+prb+'_maged_75keV_dtc_cor_flux', data = d

  dp2 = dblarr(n_elements(d.x))

  flux=dblarr(n_elements(d.x),9)
  newy=dblarr(n_elements(d.x),9)
  for i = 0, n_elements(d.x)-1 do newy[i,*]=d.y[i,[pitch_order[i,*]]]
  for i = 0, n_elements(d.x)-1 do flux[i,*] = newy[i,*]*dalpha; *sin(newpitch[i,*]*!pi/180.)*dalpha[i,*]*!pi/180.
  for i = 0, n_elements(d.x)-1 do dp2[i] = total(flux[i,*])
 

  dcom.y[*,1] = dp2

  ; 150 keV

  get_data, 'g'+prb+'_maged_150keV_dtc_cor_flux', data = d

  dp3 = dblarr(n_elements(d.x))

  flux=dblarr(n_elements(d.x),9)
  newy=dblarr(n_elements(d.x),9)
  for i = 0, n_elements(d.x)-1 do newy[i,*]=d.y[i,[pitch_order[i,*]]]
  for i = 0, n_elements(d.x)-1 do flux[i,*] = newy[i,*]*dalpha ;*sin(newpitch[i,*]*!pi/180.)*dalpha[i,*]*!pi/180.
  for i = 0, n_elements(d.x)-1 do dp3[i] = total(flux[i,*])
 

  dcom.y[*,2] = dp3

  ; 275 keV

  get_data, 'g'+prb+'_maged_275keV_dtc_cor_flux', data = d

  dp4 = dblarr(n_elements(d.x))

  flux=dblarr(n_elements(d.x),9)
  newy=dblarr(n_elements(d.x),9)
  for i = 0, n_elements(d.x)-1 do newy[i,*]=d.y[i,[pitch_order[i,*]]]
  for i = 0, n_elements(d.x)-1 do flux[i,*] = newy[i,*]*dalpha ;*sin(newpitch[i,*]*!pi/180.)*dalpha[i,*]*!pi/180.
  for i = 0, n_elements(d.x)-1 do dp4[i] = total(flux[i,*])
 

  dcom.y[*,3] = dp4

  ; 475 keV

  get_data, 'g'+prb+'_maged_475keV_dtc_cor_flux', data = d

  dp5 = dblarr(n_elements(d.x))

  flux=dblarr(n_elements(d.x),9)
  newy=dblarr(n_elements(d.x),9)
  for i = 0, n_elements(d.x)-1 do newy[i,*]=d.y[i,[pitch_order[i,*]]]
  for i = 0, n_elements(d.x)-1 do flux[i,*] = newy[i,*]*dalpha ;*sin(newpitch[i,*]*!pi/180.)*dalpha[i,*]*!pi/180.
  for i = 0, n_elements(d.x)-1 do dp5[i] = total(flux[i,*])
 

  dcom.y[*,4] = dp5

  store_data, 'g'+prb+'_maged_combo_flux', data = dcom

  options, 'g'+prb+'_maged_combo_flux', /ylog, ytitle = 'MAGED e!u-!n flux', ysubtitle = '[#/cm!u2!n-s-sr-keV]', thick = 1.5, labels = string([40.0, 75.0, 150.0, 275.0, 475.0], FORMAT = '(i3)')+' keV', labflag = 1

  tplot,'g'+prb+'_maged_combo_flux'

  ; GOES Positions

  if prb eq '13' then calc, "'g13_pos_gsm_re' = 'g13_pos_gsm'/6378.137"     ; Convert to [RE]
  if prb eq '14' then calc, "'g14_pos_gsm_re' = 'g14_pos_gsm'/6378.137"     ; Convert to [RE]
  if prb eq '15' then calc, "'g15_pos_gsm_re' = 'g15_pos_gsm'/6378.137"     ; Convert to [RE]

  options, 'g1?_pos_gsm_re', labels = ['X!dGSM', 'Y!dGSM', 'Z!dGSM']

  ; Ephemeris
  get_data, 'g'+prb+'_pos_gsm', data = posgsm

  xgsm = {x : posgsm.x, y : posgsm.y[*,0]/6378.137}
  ygsm = {x : posgsm.x, y : posgsm.y[*,1]/6378.137}
  zgsm = {x : posgsm.x, y : posgsm.y[*,2]/6378.137}

  store_data, 'g'+prb+'_state_pos_xgsm', data = xgsm
  store_data, 'g'+prb+'_state_pos_ygsm', data = ygsm
  store_data, 'g'+prb+'_state_pos_zgsm', data = zgsm

  options, 'g'+prb+'_state_pos_xgsm', ytitle='X!dGSM!n'
  options, 'g'+prb+'_state_pos_ygsm', ytitle='Y!dGSM!n'
  options, 'g'+prb+'_state_pos_zgsm', ytitle='Z!dGSM!n'


tplot,'g'+prb+'_H_gsm g'+prb+'_maged_combo_flux', title='GOES '+prb, $
  var_label=['g'+prb+'_state_pos_zgsm','g'+prb+'_state_pos_ygsm','g'+prb+'_state_pos_xgsm']
  

  
  
  ;Load GOES EPEAD data
  
  goes_load_data, trange=trange, datatype='epead', probes=prb, /avg_1m 
  goes_lib 
  goes_epead_center_pitch_angles, 'g'+prb+'_Bsc_1', 'g'+prb+'_BTSC_1' 
  if prb eq '13' then goes_epead_contam_cor, ['g13_elec_0.6MeV_uncor_flux','g13_elec_2MeV_uncor_flux'],['g13_prot_30.6MeV_uncor_flux','g13_prot_63.1MeV_uncor_flux','g13_prot_165MeV_uncor_flux','g13_prot_433MeV_uncor_flux']
  if prb eq '15' then goes_epead_contam_cor, ['g15_elec_0.6MeV_uncor_flux','g15_elec_2MeV_uncor_flux'],['g15_prot_30.6MeV_uncor_flux','g15_prot_63.1MeV_uncor_flux','g15_prot_165MeV_uncor_flux','g15_prot_433MeV_uncor_flux']
  
  
  if prb eq '13' then get_data,'g13_elec_2MeV_uncor_flux',data=two
  if prb eq '15' then get_data,'g15_elec_2MeV_uncor_flux',data=two
  
  get_data,'goes_epead_center_pitch_angles',data=pitch
  dp6 = dblarr(n_elements(d.x))
    
  pitch.y[where(pitch.y gt 90.)]=180.-pitch.y[where(pitch.y gt 90.)]
  pitch_order=dblarr(n_elements(pitch.x),2)
  newpitch=dblarr(n_elements(pitch.x),2)
  for i=0,n_elements(pitch.x)-1 do  pitch_order[i,*]=sort(pitch.y[i,*])            ;Put pitch angles in ascending order
    for i=0,n_elements(pitch.x)-1 do newpitch[i,*]=pitch.y[i,[pitch_order[i,*]]]
    
     ;[-cos(alpha_{n+1}) + cos(alpha_{n-1})]/2  (replacing   sin(alpha_n)*[alpha_{n+1} - alpha_{n-1}]/2  )
     
     alpha_nplus1=[[newpitch[*,1]],[replicate(90.,n_elements(pitch.x))]]*!PI/180.
     alpha_nmins1=[[replicate(0.,n_elements(pitch.x))],[newpitch[*,0]]]*!PI/180.
     dalpha=(-cos(alpha_nplus1) + cos(alpha_nmins1))/2.
  
  flux=dblarr(n_elements(two.x),2)
  newy=dblarr(n_elements(two.x),2)
  for i = 0, n_elements(two.x)-1 do newy[i,*]=d.y[i,[pitch_order[i,*]]]
  for i = 0, n_elements(two.x)-1 do flux[i,*] = newy[i,*]*dalpha              
  for i = 0, n_elements(two.x)-1 do dp6[i] = total(flux[i,*])
 
;  two_1=two.y[*,0]*sin(angles.y[*,0]*!PI/180.)
;  two_2=two.y[*,1]*sin(angles.y[*,1]*!PI/180.)
;  two_omni=(two_1+two_2)/2.
  
  store_data,'g'+prb+'_2MeV_omni',data={x:two.x,y:dp6}
  
  tplot,'g'+prb+'_2MeV_omni'
  
  
  ;dp6 = two_omni

  dcom2.y[*,0] = dp1
  dcom2.y[*,1] = dp2
  dcom2.y[*,2] = dp3
  dcom2.y[*,3] = dp4
  dcom2.y[*,4] = dp5
  dcom2.y[*,5] = dp6
  
  store_data, 'g'+prb+'_maged_combo_flux_2', data = dcom2
  options, 'g'+prb+'_maged_combo_flux_2', /ylog, ytitle = 'MAGED e!u-!n flux', ysubtitle = '[#/cm!u2!n-s-sr-keV]', thick = 1.5, $ 
    labels = string([40.0, 75.0, 150.0, 275.0, 475.0, 2000], FORMAT = '(i4)')+' keV', labflag = 1
  ylim,'g'+prb+'_maged_combo_flux_2',1.,1e6,1
  
   
   
;
;omni_hro_load
;store_data,'omni_imf',data=['OMNI_HRO_1min_BY_GSM','OMNI_HRO_1min_BZ_GSM']
;
;get_tsy_params,'kyoto_dst','omni_imf',$
;               'OMNI_HRO_1min_proton_density', $
;               'OMNI_HRO_1min_flow_speed', $
;               'T96',/speed,/imf_yz
;
;get_data,'t96_par',data=t96_par   
;tinterpol_mxn,'t96_par','g'+prb+'_pos_gsm',newname='t96_par_int'
;  
;ttrace2equator,'g'+prb+'_pos_gsm', in_coord='gsm', out_coord='gsm', newname='GOES'+prb+'_eqpts',external_model='t96',$
;  par='t96_par_int',/km
;  
;get_data,'GOES'+prb+'_eqpts',data=eqpts
;Lstar=sqrt(eqpts.y[*,0]^2+eqpts.y[*,1]^2+eqpts.y[*,2]^2)/6371.
;store_data,'GOES'+prb+'_Lstar',data={x:eqpts.x,y:Lstar}   
;tinterpol_mxn,'GOES'+prb+'_Lstar','g'+prb+'_maged_combo_flux_2',newname='GOES'+prb+'_Lstar_int'
;get_data,'GOES'+prb+'_Lstar_int',data=Lstar
;

tplot,'kyoto_dst thmAL g'+prb+'_H_gsm g'+prb+'_maged_combo_flux_2', $
  var_label=['g'+prb+'_state_pos_zgsm','g'+prb+'_state_pos_ygsm','g'+prb+'_state_pos_xgsm'],title='GOES '+prb
 
tplot,'kyoto_dst thmAL  g13_H_gsm g13_maged_combo_flux_2 g'+prb+'_H_gsm g'+prb+'_maged_combo_flux_2',var_label=['g'+prb+'_state_pos_zgsm','g'+prb+'_state_pos_ygsm','g'+prb+'_state_pos_xgsm']
   
   
   

store_data,'GOES'+prb+'_Lstar_flux_40keV',data={x:Lstar.x,y:Lstar.y/6371.,v:dcom2.y[*,0]}
store_data,'GOES'+prb+'_Lstar_flux_75keV',data={x:Lstar.x,y:Lstar.y/6371.,v:dcom2.y[*,1]}
store_data,'GOES'+prb+'_Lstar_flux_150keV',data={x:Lstar.x,y:Lstar.y/6371.,v:dcom2.y[*,2]}
store_data,'GOES'+prb+'_Lstar_flux_275keV',data={x:Lstar.x,y:Lstar.y/6371.,v:dcom2.y[*,3]}
store_data,'GOES'+prb+'_Lstar_flux_475keV',data={x:Lstar.x,y:Lstar.y/6371.,v:dcom2.y[*,4]}
store_data,'GOES'+prb+'_Lstar_flux_2000keV',data={x:Lstar.x,y:Lstar.y/6371.,v:dcom2.y[*,5]}

options,'GOES*Lstar_flux*','spec',1
zlim,'GOES*Lstar_flux*',10.,1.e6,1