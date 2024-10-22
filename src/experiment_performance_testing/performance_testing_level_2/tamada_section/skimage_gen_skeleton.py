import kimimaro
import numpy as np
import imageio
from skimage import measure
import os
import edt
from cloudvolume import Skeleton
from skimage.morphology import skeletonize
from skimage import data
import matplotlib.pyplot as plt
from skimage.util import invert

import os

# Run lzma -d connectomics.npy.lzma on the command line to
# obtain this 512 MB segmentation volume. Details below.
input_volume_file = "/work/boyu/Ethan_myelin_sheath/2048_myelin_sheath_ttest/t_18.tifrefined_skel_dilated_test.tif"
namex = 't_18'
skelFolder = "/work/boyu/Ethan_myelin_sheath/2048_myelin_sheath_ttest"
vol = imageio.mimread(input_volume_file, memtest=False)
labels = np.stack(vol, 2)
labels.astype(int)

# Invert the horse image
image = invert(labels)

# perform skeletonization
skeleton = skeletonize(image,method="lee")
w = imageio.get_writer(os.path.join(skelFolder, namex) + '_3Dskel_skimage.tif', mode='v')
lenx, leny, lenz = skeleton.shape
for i in range(lenz):
    w.append_data(skeleton[:, :, i])
w.close()
