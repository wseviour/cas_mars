import os
import glob

ext =  'isentropic_surface_data_ls=240_to_ls=266'
day = 10
lev = 300

cons = sorted([os.path.basename(x) for x in glob.glob("input_contours_MarsWRF/%sK/pv*%s*.in" % (lev,ext))])

if not os.path.isdir('cas_output_MarsWRF/%sK' % lev):
    os.system("mkdir -p cas_output_MarsWRF/%sK" % lev)

for icon in cons:
    print "**** %s ****" % icon
    os.system("rm -r winds")
    os.system("cp -R ../MarsWRF_data/MarsWRF_winds_%sK/%s winds" % (lev,ext))
    os.system("cp input_contours_MarsWRF/%sK/%s pv.in" % (lev,icon))
    os.system("./cas < run_params_MarsWRF.in")
    os.system("cp output.cas cas_output_MarsWRF/%sK/%s.cas" % (lev,icon))
