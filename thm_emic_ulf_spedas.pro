pro thm_EMIC_ULF_SPEDAS
  thm_init
  !p.charsize = 1.2
  tplot_options,'xmargin',[20.,15.]
  tplot_options,'ymargin',[3,3]
  tplot_options,'xticklen',0.08
  tplot_options,'yticklen',0.02
  tplot_options,'xthick',2
  tplot_options,'ythick',2
  probes = ['a','d','e']
  prange = [time_double('2015-03-17/00:00:00'), time_double('2015-03-18/00:00:00')]
  timespan, prange
  thm_load_state, probe = probes, trange = prange,/get_support_data
  thm_load_fgm, probe = probes, trange = prange, level=2, coord='gse', data='fgl'
  thm_load_fgm, probe = probes, trange = trange, level=2, coord='dsl', data='fgs'
;  stop
;  thm_load_scm, probe=probes, trange = trange, level=1,/get_support_data, datatype = 'scf',cleanup='full'
  for ip=0., n_elements(probes)-1 do begin
    sc = probes[ip]
    tvectot, 'th'+sc+'_fgl_gse',tot='Bmag'
    get_data, 'Bmag', data = bmag
    fce_eq = 28.*abs(bmag.y)
    store_data,'th'+sc+'_fce_eq',data={x:bmag.x,y:fce_eq},dlim={colors:'w',thick:'1',linestyles:'0'}
    store_data,'th'+sc+'_fcp_eq',data={x:bmag.x,y:fce_eq/1837},dlim={colors:1,thick:'1',linestyles:'2'}
    store_data,'th'+sc+'_fche_eq',data={x:bmag.x,y:fce_eq/1837/4.0},dlim={colors:1,thick:'1',linestyles:'4'}
    store_data,'th'+sc+'_fco_eq',data={x:bmag.x,y:fce_eq/1837/16.0},dlim={colors:1,thick:'1',linestyles:'1'}
    tdegap, 'th'+sc+'_fgl_gse',/overwrite 
    tsmooth2, 'th'+sc+'_fgl_gse',20,newname='th'+sc+'_fgl_gse_sm' ;;smooth for a more stable background Bfield for the fac transformation
    nop = 1024.0 ;;number of points for fft
    ;Transform the Mag data to FAC coordinates
    thm_fac_matrix_make,'th'+sc+'_fgl_gse_sm',other_dim='Rgeo',pos_var_name='th'+sc+'_state_pos',newname='th'+sc+'_fac_mat'
    tinterpol_mxn, 'th'+sc+'_fac_mat', 'th'+sc+'_fgl_gse', newname='th'+sc+'_fac_mat_int'
    tvector_rotate,'th'+sc+'_fac_mat_int','th'+sc+'_fgl_gse',newname='th'+sc+'_fac_rgeo'
    tsmooth2, 'th'+sc+'_fac_rgeo', 400.0, newname='th'+sc+'_fac_rgeo_sm' ;;high-pass filter (keep f>0.01 Hz)
    dif_data, 'th'+sc+'_fac_rgeo', 'th'+sc+'_fac_rgeo_sm', newname='th'+sc+'_fac_rgeo_hp',copy_dlimits=1
    options, 'th'+sc+'_fac_rgeo_hp', ytitle='BFAC'
    options, ['th'+sc+'fgl_gse', 'th'+sc+'_fac_rgeo_hp', 'th'+sc+'_fac_rgeo'],colors=[2,4,6]
    get_timespan,t
    time_clip, 'th'+sc+'_fac_rgeo_hp',t(0),t(1),newname='th'+sc+'_fac_rgeo_tclip'
    twavpol, 'th'+sc+'_fac_rgeo_tclip', nopfft=nop, steplength=fix(nop*0.7)
    get_data, 'th'+sc+'_fac_rgeo_tclip_waveangle', data=wna_tmp,dlimit=dlim
    wna_tmp.y = wna_tmp.y/!dtor
    store_data, 'th'+sc+'_fac_rgeo_tclip_waveangle', data={x:wna_tmp.x, y:wna_tmp.y, v:wna_tmp.v}, dlimit=dlim
    store_data, 'th'+sc+'_Bspec', data=['th'+sc+'_fac_rgeo_tclip_powspec','th'+sc+'_fcp_eq','th'+sc+'_fche_eq','th'+sc+'_fco_eq']
    store_data, 'th'+sc+'_degpol',data=['th'+sc+'_fac_rgeo_tclip_degpol','th'+sc+'_fcp_eq','th'+sc+'_fche_eq','th'+sc+'_fco_eq']
    store_data, 'th'+sc+'_waveangle_theta',data=['th'+sc+'_fac_rgeo_tclip_waveangle','th'+sc+'_fcp_eq','th'+sc+'_fche_eq','th'+sc+'_fco_eq']
    store_data, 'th'+sc+'_elliptict',data=['th'+sc+'_fac_rgeo_tclip_elliptict','th'+sc+'_fcp_eq','th'+sc+'_fche_eq','th'+sc+'_fco_eq']
    options, 'th'+sc+['_fac_rgeo_tclip_powspec', '_fac_rgeo_tclip_degpol', '_fac_rgeo_tclip_waveangle', '_fac_rgeo_tclip_elliptict'], spec=1
    ylim, 'th'+sc+['_Bspec','_degpol','_waveangle_theta','_elliptict'],0.02,2,1
    zlim,'th'+sc+'_Bspec',1.e-5,1.e2,1
    zlim, 'th'+sc+'_waveangle_theta', 0.,90.,0
    options, 'th'+sc+'_Bspec', ztitle='(nT)!u2!n/Hz', ytitle='Frequency', ysubtitle='[Hz]',/zlog
    options, 'th'+sc+'_degpol', ztitle='degpol', ytitle='Frequency', ysubtitle='[Hz]'
    options, 'th'+sc+'_waveangle_theta', ztitle='WNA', ytitle='Frequency', ysubtitle='[Hz]',ystyle=1
    options, 'th'+sc+'_elliptict', ztitle='ellipticity', ytitle='Frequency', ysubtitle='[Hz]',zrange=[-1.0,1.0],ystyle=1
      
    time_stamp,/off
    ;;plot variables: fft-spectra from fgl data; degree of polarization; wavenormal angle in deg; ellipticity
    tplot,['th'+sc+'_Bspec', 'th'+sc+'_degpol','th'+sc+'_waveangle_theta','th'+sc+'_elliptict'],title='THEMIS-'+strupcase(sc)
    stop
    ;;for the ULF wave
    tdegap, 'th'+sc+'_fgs_dsl',/overwrite
    tdegap, 'th'+sc+'_scf',/overwrite
    tsmooth2, 'th'+sc+'_fgs_dsl',167,newname='th'+sc+'_fgs_dsl_sm' ;;smooth for a more stable background Bfield for the fac transformation
    
    split_vec,'th'+sc+'_fgs_dsl'
    get_data, 'th'+sc+'_fgs_dsl_x', data= Bx
    get_data, 'th'+sc+'_fgs_dsl_y', data= By
    get_data, 'th'+sc+'_fgs_dsl_z', data= Bz
    ;;specify the lower and upper frequency limits for ULF waves
    flow=0.002
    fhigh=0.02
    dt=3.0
    Bx_out = thm_lsp_filter_mod(Bx.y, dt, flow, fhigh)
    By_out = thm_lsp_filter_mod(By.y, dt, flow, fhigh)
    Bz_out = thm_lsp_filter_mod(Bz.y, dt, flow, fhigh)
    store_data, 'th'+sc+'_fgs_dsl_bp', data= {x:Bx.x, y:[[Bx_out], [By_out], [Bz_out]]}
    options,'th'+sc+'_fgs_dsl_bp',colors=[2,4,6]
    tsmooth2, 'th'+sc+'_fgs_dsl',167,newname='th'+sc+'_fgs_dsl_sm'
    thm_fac_matrix_make,'th'+sc+'_fgs_dsl_sm',other_dim='Rgeo',pos_var_name='th'+sc+'_state_pos',newname='th'+sc+'_fac_mat2'
    tvector_rotate,'th'+sc+'_fac_mat2','th'+sc+'_fgs_dsl_bp',newname='th'+sc+'_fgs_fac_bp'
    ;;wavelet
    wav_data,'th'+sc+'_fgs_fac_bp',trange=trange,dimennum=0
    wav_data,'th'+sc+'_fgs_fac_bp',trange=trange,dimennum=1
    wav_data,'th'+sc+'_fgs_fac_bp',trange=trange,dimennum=2
    options,'th'+sc+'_fgs_fac_bp*_wv_pow',yrnage=[0.001,0.02],/ylog,zrange=[1,1.e4],/zlog
    ;;fft
    twavpol,'th'+sc+'_fgs_fac_bp',nopfft=512.0,steplength=0.25*512.0
    ;;wave power spectra of each component
    get_data,'th'+sc+'_fgs_fac_bp_pspec3',data=pspec3
    store_data,'th'+sc+'_bfacx_spec',data={x:pspec3.x,y:reform(pspec3.y[*,*,0]),v:pspec3.v}
    store_data,'th'+sc+'_bfacy_spec',data={x:pspec3.x,y:reform(pspec3.y[*,*,1]),v:pspec3.v}
    store_data,'th'+sc+'_bfacz_spec',data={x:pspec3.x,y:reform(pspec3.y[*,*,2]),v:pspec3.v}
    options,'th'+sc+'_bfac?_spec',spec=1,yrange=[0.,0.02],zrange=[1,1.e4],/zlog
    options,'th'+sc+'_fgs_fac_bp_powspec',yrange=[0.0,0.02],ytitle='BFAC!cHz',zrange=[1,1.e4],/zlog
    stop
  endfor ;;loop for each probe
end