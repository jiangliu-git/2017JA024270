themis_folder = '/Weiyun/data/themis'
rbsp_bin_folder = '/Weiyun/data/rbsp_bin'
plasma_folder = '/plasma_fill'
ground_folder = '/Weiyun/data/ground'

if keyword_set(computer) then begin
	;;; THEMIS
	fgs_folder=computer + themis_folder + '/fgs'
	fgl_folder=computer + themis_folder + '/fgl'
	pos_folder=computer + themis_folder + '/pos'
	mode_folder=computer + themis_folder + '/mode'
	roi_folder=computer + themis_folder + '/roi'
	efs_folder=computer + themis_folder + '/efs_dsl'
	eff_folder=computer + themis_folder + '/eff_dot0_dsl'
	ni_folder=computer + themis_folder + plasma_folder + '/ni'
	Pall_folder=computer + themis_folder + plasma_folder + '/Pall'
	Pth_folder=computer + themis_folder + plasma_folder + '/Pth'
	beta_folder=computer + themis_folder + plasma_folder + '/beta'
	vi_folder = computer + themis_folder + plasma_folder + '/vi'
	viperp_folder = computer + themis_folder + plasma_folder + '/viperp'
	ve_folder = computer + themis_folder + plasma_folder + '/ve'
	vexb_folder = computer + themis_folder + plasma_folder + '/vexb'
	vexb_dsl_folder = computer + themis_folder + plasma_folder + '/vexb_dsl'
	veperp_folder = computer + themis_folder + plasma_folder + '/veperp'
	Pttl_fgs_folder=computer+themis_folder+plasma_folder+'/Pttl_fgs'
	Blobe_folder=computer+themis_folder+plasma_folder+'/Blobe'
	;;; RBSP
	rbsp_folder = computer+'/Weiyun/data/rbsp'
	pos_rbsp_folder = computer + rbsp_bin_folder + '/pos'
	;;; solar wind
	imf_folder = computer+'/Weiyun/data/sw/omni_b_gsm'
	vsw_folder = computer+'/Weiyun/data/sw/omni_v_gsm'
	nisw_folder = computer+'/Weiyun/data/sw/omni_ni'
	Pdynsw_folder = computer+'/Weiyun/data/sw/omni_Pdyn'
	dst_folder = computer+'/Weiyun/data/ground/kyoto_dst'
	;;; al
	al_folder = computer+ground_folder+'/pseudo_al'
	kyoto_al_folder = computer+ground_folder+'/kyoto_al'
endif

list_folder='../lists'
save_folder='variables'
pic_folder = '../../results/temp'

;;;;;; 
l_thick = 2.4 ;; for publication
;;; constants
Re = 6371.
mu0 = 4*!pi*1e-7
nTesla2_to_nPa = 0.01/25.132741
ni2npa = 1.67262158e-27*1e21 ;; for computation of Pdyn
mass_proton = 1.67262158e-27 ;; kg

;;; signs
perp_sign = '!9'+string("136B)+'!X'
minus_sign = '!9'+string("055B)+'!X'
cross = '!9'+string("264B)+'!X'
dot_sign = '!9'+string("327B)+'!X'
theta_letter = '!9'+string("161B)+'!X'
phi_letter = '!9'+string("152B)+'!X'
delta_letter = '!9'+string("144B)+'!X'
deltau_letter = '!9'+string("104B)+'!X' ;; upper case delta
rho_letter = '!9'+string("162B)+'!X'
l_angle = '!9'+string("341B)+'!X'
r_angle = '!9'+string("361B)+'!X'

;;; for DFB properties
seconds_check = 15.

;;; turn off the timestamp of tplot
time_stamp,/off
