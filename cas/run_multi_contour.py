import os
import glob

ext =  'res_test57.-70.-nu4-urlx-kt4.0.c-0020sat200.0.T85'
day = 200

cons = sorted([os.path.basename(x) for x in glob.glob("input_contours/pv*%s*.in" % ext)])

if not os.path.isdir('cas_output'):
    os.system("mkdir cas_output")

for icon in cons:
    print "**** %s ****" % icon
    os.system("rm -r winds")
    os.system("cp -R ../model_output/winds_%s winds" % ext)
    os.system("cp input_contours/%s pv.in" % icon)
    os.system("./cas < run_params_will.in")
    os.system("cp output.cas cas_output/%s.cas" % icon)
