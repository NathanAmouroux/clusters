import numpy as np
from astropy.table import Table, vstack
import healpy as hp
import sys
import os
import time
import matplotlib.pyplot as plt
from astropy.io import ascii
mb_inpath = '/sps/lsst/groups/clusters/amico_validation_project/catalogs/AMICO/raw_amico_cats/test_cosmoDC2_compute_lambstar/9559_map_associations_noBuffer.fits'
outpath = '/sps/lsst/groups/clusters/amico_validation_project/catalogs/AMICO/raw_amico_cats/test_cosmoDC2_compute_lambstar/'
neighbours_path = '/sps/lsst/groups/clusters/amico_validation_project/catalogs/CosmoDC2/cosmodc2_neighbours.fits'
healpath = '/sps/lsst/users/tguillem/web/clusters/catalogs/cosmoDC2_photoz_flexzboost/v1/'
#inpath = '/sps/lsst/groups/clusters/amico_validation_project/catalogs/AMICO/amico_cats/small/amico_map_associations/'
# = '/sps/lsst/groups/clusters/amico_validation_project/catalogs/AMICO/amico_cats/small/amico_map_associations/all_maps.fits'
curr_tile = 9559
print('\nProcessing', curr_tile)   
start = time.time()
# ---Reading tables---#
print("\nReading members table")
mbn = Table.read(mb_inpath)

tile_t = Table.read(neighbours_path)

tile_l = tile_t[tile_t['tile']==int(curr_tile)]['list_of_neighbour_tiles'][0].split(',') #Because neighbours listed w/ comas in a long string --> need to convert it to a list
for i in range(len(tile_l)):
    tile_l[i] = int(tile_l[i])
print('Neighbour healpix list (including current healpix)', ':', tile_l, '\nNow merging all input healpix files')

large_heal = Table.read(healpath + str(tile_l[0]) + '/galaxies.fits')
for tile in tile_l[1:]:
    heal = Table.read(healpath + str(tile) + '/galaxies.fits')
    large_heal = vstack([large_heal, heal])
print('Done. \nIncorporation of magnitudes')


# ---Add Magnitudes---#

heal_dict = {}

#CosmoDC2
for i, ID in enumerate(large_heal['galaxy_id']):
    heal_dict[ID] = heal_dict.get(ID, [])+[i]

#DC2
#for i, ID in enumerate(large_heal['id']):
#    heal_dict[ID] = heal_dict.get(ID, [])+[i]
mbn['matched'] = np.array([i in heal_dict for i in mbn['GALID']])
mbn['mt_ids'] = None
for i, ID in enumerate(mbn['GALID']):
    mbn['mt_ids'][i] = heal_dict.get(ID, [])

mbn['mag_g'] = 0.0
mbn['mag_i'] = 0.0
mbn['mag_r'] = 0.0
mbn['mag_z'] = 0.0
mbn['mag_y'] = 0.0
mbn['ra'] = 0.0
mbn['dec'] = 0.0
mbn['redshift'] = 0.0
mags = ['mag_g','mag_i','mag_r','mag_z','mag_y','ra','dec','redshift']
per_list = [.1,.2,.4,.6,.8,1.0]
k = 0
mbn['mt_id_final'] = None

for i, amg in enumerate(mbn):
    per = round(i/len(mbn),2)
    if per in per_list and k != per:
        print(i, per)
        k = per
    if len(amg['mt_ids'])==1:
        mbn['mt_id_final'][i] = amg['mt_ids'][0]
        healsh = large_heal[amg['mt_ids']]
        for mag in mags:
            mbn[mag][i] =  healsh[mag]
    elif len(amg['mt_ids'])>1:
        healsh = large_heal[amg['mt_ids']]
        j = np.arange(healsh.size, dtype=float)[0]
        #print(dc2h['mt_ids'], mtmask, j)
        mbn['mt_id_final'][i] = amg['mt_ids'][j]
        for mag in mags:
            mbn[mag][i] =  healsh[mag][j]


print(mbn[:5])
print('There are :', len(mbn[mbn['mag_g']==0]['mag_g']), 'galaxies which are not found inside input catalog')


mbn_mag = mbn['GALID','FIELD_PROB','ASSOC_ID','ASSOC_PROB', 'ra', 'dec', 'redshift', 'mag_g', 'mag_i', 'mag_r', 'mag_z', 'mag_y']
mbn_mag.write(outpath + str(curr_tile) + '_map_associations_w_mag'  + '.fits', overwrite = True)