import os

width=15
#outdir='swvac/const_width_%s' % width
outdir='swvac/instability_test_nosmooth_mars-like'

for width in range(5,26,5):
    for iclat1 in range(0,86-width,5):
        os.system('./rbvac_batch %02d %02d %s' % (iclat1, iclat1+width, outdir))
