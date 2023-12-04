from astropy.table import Table, vstack
import numpy as np
import matplotlib.pyplot as plt

healpath = "/sps/lsst/users/tguillem/web/clusters/catalogs/DC2_photoz_flexzboost/v0/"
heal_raw = Table.read(healpath + '3631/galaxies.fits')
heal_raw2 = Table.read(healpath + '3632/galaxies.fits')
heal_raw3 = Table.read(healpath + '3633/galaxies.fits')

c2_mb = Table.read('/sps/lsst/groups/clusters/amico_validation_project/catalogs/AMICO/amico_cats/DC2/mag_y/3631_map_associations_w_mag.fits')
l = [3632,3633]
for i in l :
    cur_t = Table.read(f'/sps/lsst/groups/clusters/amico_validation_project/catalogs/AMICO/amico_cats/DC2/mag_y/{str(i)}_map_associations_w_mag.fits')
    c2_mb = vstack([c2_mb,cur_t])

all_t = Table.read('/sps/lsst/groups/clusters/amico_validation_project/catalogs/AMICO/amico_cats/DC2/mag_y/all_maps.fits')

plt.scatter(all_t['ra'], all_t['dec'], s= 1, alpha = .3, label = 'large_cat')
#plt.scatter(heal_raw['ra'], heal_raw['dec'], s= 1, alpha = .5, label = '3631')
#plt.scatter(heal_raw2['ra'], heal_raw2['dec'], s= 1, alpha = .5, label = '3632')
#plt.scatter(heal_raw3['ra'], heal_raw3['dec'], s= 1, alpha = .5, label = '3633')
plt.scatter(c2_mb['ra'], c2_mb['dec'], s= 1, alpha = .1, label = 'AMICO')
plt.legend()
plt.xlim([47,76])
plt.ylim([-46,-26])
plt.savefig('verif_1_vs_all.png')