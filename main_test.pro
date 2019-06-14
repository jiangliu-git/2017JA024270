pro main_test
;;; testing things
computer = 'I:'
@folders
thm_init
del_data, '*'

date = '2013 9 30 '
trange_show = date+['1 0 1', '1 39 30'] ;; regular range
trange_load = time_double(trange_show)+[-100, 100]

;;; load data
load_bin_data, trange = trange_load, datatype = 'kyoto_al', datafolder = kyoto_al_folder
load_bin_data, trange = trange_load, datatype = 'pseudo_al', datafolder = al_folder

tplot, ['kyoto_al', 'pseudo_al'], trange = trange_show

;;; test draw_in_tpanel
draw_in_tpanel, '13 9 30 '+['1 10', '1 20'], [-30, -40], varname = 'kyoto_al', line = 1, color = 1, thick = 2

end
