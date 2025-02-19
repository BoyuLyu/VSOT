from cloudvolume import CloudVolume
from cloudvolume.lib import Bbox, Vec, yellow, mkdir
from concurrent.futures import ProcessPoolExecutor
import imageio
import numpy as np
import os
import fastremap
import json
if __name__ == '__main__':
    volSeg = CloudVolume('precomputed://https://storage.googleapis.com/iarpa_microns/minnie/minnie65/seg', progress=True,
                         mip=1,
                         fill_missing=True,
                         parallel=True)
    volCleft = CloudVolume("precomputed://https://storage.googleapis.com/iarpa_microns/minnie/minnie65/clefts",
                           mip=1,
                           fill_missing=True,
                           parallel=True)
    rootFolder = "/data_path/astro_19_minnie65"
    # chunk size varies for each astrocyte region
    # the resolution downloaded is 16x16x40
    if not os.path.isdir(rootFolder):
        os.mkdir(rootFolder)
    for m in range(5):
        if not os.path.isdir(os.path.join(rootFolder, str(m))):
            os.mkdir(os.path.join(rootFolder, str(m)))
        for n in range(5):
            if not os.path.isdir(os.path.join(rootFolder, str(m), str(n))):
                os.mkdir(os.path.join(rootFolder, str(m), str(n)))
            for k in range(5):
                if not os.path.isdir(os.path.join(rootFolder, str(m), str(n), str(k))):
                    os.mkdir(os.path.join(rootFolder, str(m), str(n), str(k)))
                print([m, n, k])
                # select the bounding box for each astrocyte region based on the coordinates of the region in Neuroglancer
                bbox = Bbox((160678 / 2 + 900 * m, 101200 / 2 + 700 * n, 20994 + 400 * (k)),
                            (160678 / 2 + 900 * (m + 1), 101200 / 2 + 700 * (n + 1), 20994 + 400 * (k + 1)))
                data = volSeg.download(bbox)
                data2 = np.squeeze(data)
                astro1 = data2 == 864691136273618701 # output the mask for the target astrocyte region
                astro1.astype(int)
                w = imageio.get_writer(os.path.join(rootFolder, str(m), str(n), str(k), "astro_segMask.tif"),
                                       mode='v')
                lenx, leny, lenz = astro1.shape
                for i in range(lenz):
                    w.append_data(astro1[:, :, i])
                w.close()
                labels, remapping = fastremap.renumber(data2, in_place=True,
                                                       preserve_zero=True)  # relabel values from 1 and refit data type
                with open(os.path.join(rootFolder, str(m), str(n), str(k), "seg_mapping.txt"), 'w') as file:
                    file.write(json.dumps(remapping))
                file.close()

                w = imageio.get_writer(os.path.join(rootFolder, str(m), str(n), str(k), "segMaskFull.tif"),
                                       mode='v')
                lenx, leny, lenz = labels.shape
                for i in range(lenz):
                    w.append_data(labels[:, :, i])
                w.close()
                cleftdata = volCleft.download(bbox)
                cleftdata2 = np.squeeze(cleftdata)
                cleftdata2, remapping = fastremap.renumber(cleftdata2, in_place=True, preserve_zero=True)
                w = imageio.get_writer(os.path.join(rootFolder, str(m), str(n), str(k), "cleft.tif"), mode='v')
                lenx, leny, lenz = cleftdata2.shape
                for i in range(lenz):
                    w.append_data(cleftdata2[:, :, i])
                w.close()
                with open(os.path.join(rootFolder, str(m), str(n), str(k), "cleft_mapping.txt"), 'w') as file:
                    file.write(json.dumps(remapping))
                file.close()