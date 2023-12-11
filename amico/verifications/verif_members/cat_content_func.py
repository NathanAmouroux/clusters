import numpy as np
from astropy.table import Table, vstack
import healpy as hp
import sys
import os
import time
import matplotlib.pyplot as plt
from astropy.io import ascii, fits
import pandas as pd
import random

curr_tile = str(sys.argv[1])

###AMICO CosmoDC2 mag_i
inpath = "/sps/lsst/groups/clusters/amico_validation_project/catalogs/AMICO/raw_amico_cats/cosmoDC2_photoz_flexzboost_v1_iband/DETECTIONS_DERIVED/"
gen_inpath = "/sps/lsst/groups/clusters/amico_validation_project/catalogs/AMICO/amico_cats/amico_map_associations_flxzb_mag/mag_i/"
gen_mb_inpath = gen_inpath + '9939_map_associations_w_mag.fits'
raw_mb_inpath = inpath + '9939_map_associations_noBuffer.fits'
cl_inpath = inpath + '9939_map_detections_refined_noBuffer.txt'

###AMICO CosmoDC2 mag_y
#inpath = "/sps/lsst/groups/clusters/amico_validation_project/catalogs/AMICO/raw_amico_cats/cosmoDC2_photoz_flexzboost_v1_yband/DETECTIONS_DERIVED/"
#gen_inpath = "/sps/lsst/groups/clusters/amico_validation_project/catalogs/AMICO/amico_cats/amico_map_associations_flxzb_mag/mag_y/"
#gen_mb_inpath = gen_inpath + '9939_map_associations_w_mag.fits'
#raw_mb_inpath = inpath + '9939_map_associations_noBuffer.fits'
#cl_inpath = inpath + '9939_map_detections_refined_noBuffer.txt'

###AMICO DC2 mag_y
#inpath = "/sps/lsst/groups/clusters/amico_validation_project/catalogs/AMICO/raw_amico_cats/DC2_v0_yband/output/DETECTIONS_DERIVED/"
#gen_inpath = "/sps/lsst/groups/clusters/amico_validation_project/catalogs/AMICO/amico_cats/DC2/mag_y/"
#gen_mb_inpath = gen_inpath + curr_tile + '_map_associations_w_mag.fits'
#raw_mb_inpath = inpath + curr_tile + '_map_associations_noBuffer.fits'
#cl_inpath = inpath + curr_tile + '_map_detections_refined_noBuffer.txt'
#oupath = "/pbs/home/n/namourou/workspace/side_codes/clusters/amico/verifications/verif_members/results/"

gen_mb = Table.read(gen_mb_inpath)
raw_mb = Table.read(raw_mb_inpath)
cl = pd.read_csv(cl_inpath, sep="\t", header = 243)
cl = Table.from_pandas(cl)


print("Variables avaible for raw clusters catalog are : ", cl.keys()[:20], "\nVariables avaible for raw members catalog are : ", raw_mb.keys(), "\nVariables avaible for generated member catalog are : ", gen_mb.keys())

x = np.mean(cl['Xphys'])
y = np.mean(cl['Yphys'])
l_x = np.std(cl['Xphys'])/2
l_y = np.std(cl['Yphys'])/2

x_range = [x-l_x, x+l_x]
y_range = [y-l_y, y+l_y]
x_corner, y_corner = np.meshgrid(x_range,y_range)

clbox = cl[(cl["Xphys"]>=x_range[0])*(cl["Xphys"]<=x_range[1])*(cl["Yphys"]>=y_range[0])*(cl["Yphys"]<=y_range[1])]

if len(clbox)>=20:
    index = random.sample(range(len(clbox)),20)
else : 
cl_test = clbox[index]

### Test function
l = 0
for id in cl_test["# ID"]:
    raw_cdt = raw_mb[raw_mb["ASSOC_ID"]==id]
    gen_cdt = gen_mb[gen_mb["ASSOC_ID"]==id]
    if len(raw_cdt) == len(gen_cdt) and len(raw_cdt) !=0:
        l += 1
        n = 0
        m = 0
        print(f"cluster {id} has the same amount of members (n={len(raw_cdt)})")
        for galid in raw_cdt["GALID"]:
            raw_mb_cdt = raw_cdt[raw_cdt["GALID"]==galid]
            gen_mb_cdt = gen_cdt[gen_cdt["GALID"]==galid]
            if raw_mb_cdt["GALID"] != gen_mb_cdt["GALID"] :
                print(f"For cluster {id} and member {galid}, ids don't match")
                n += 1
            if raw_mb_cdt["ASSOC_PROB"] != gen_mb_cdt["ASSOC_PROB"] :
                print(f"For cluster {id} and member {galid}, pmem don't match (pmem_raw = {raw_mb_cdt['ASSOC_PROB']}, pmem_gen = {gen_mb_cdt['ASSOC_PROB']})")
                m += 1
        if n == 0 and m == 0 :
            print(f"\tmembers of cluster {id} are the same")
        else :
            lines = ["There is an error with this cluster in this tile"]
            with open(outpath + f'{curr_tile}_{id}.txt', 'w') as f:
                f.writelines(lines)
print(f"Number of clusters which have members {l}")
if l != 20 :
    lines = ["There is an error somewhere with this tile"]
    with open(outpath + f'{curr_tile}.txt', 'w') as f:
        f.writelines(lines)

        