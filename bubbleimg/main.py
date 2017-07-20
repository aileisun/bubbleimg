# import os
# import numpy as np
import astropy.table as at

import bubbleimg

reload(bubbleimg)

survey = 'hsc'
f = '/asiaa/home/ssp201718/Documents/nnarenraju/bubbleimg/bubbleimg/Database/conesearch_test.csv'
tsample = at.Table.read(f)

dir_batch = '/asiaa/home/ssp201718/Downloads/Directory/Batch/'

b = bubbleimg.batch.hscBatch(dir_batch = dir_batch, catalog = tsample, args_to_list = ['Z', 'logOIII'])
b.build(overwrite=False)
