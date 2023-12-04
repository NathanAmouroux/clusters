import healsparse as hsp
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table

filename_r = "/sps/lsst/groups/desc/shared/DC2-prod/Run2.2i/addons/supreme/dr6/supreme_dc2_dr6d_v3_u_maglim_psf_wmean.hs"
filename_i = "/sps/lsst/groups/desc/shared/DC2-prod/Run2.2i/addons/supreme/dr6/supreme_dc2_dr6d_v3_y_maglim_psf_wmean.hs"

mask_r = hsp.HealSparseMap.read(filename_r)
mask_i = hsp.HealSparseMap.read(filename_i)
valid_pixels_indices_r = mask_r.valid_pixels
valid_pixels_indices_i = mask_i.valid_pixels

nside_sparse_r = mask_r._nside_sparse
nside_sparse_i = mask_i._nside_sparse
print("processing vpix")
vpix_r, ra_r, dec_r = mask_r.valid_pixels_pos(return_pixels=True)
vpix_i, ra_i, dec_i = mask_i.valid_pixels_pos(return_pixels=True)

val_min = 30000000
val_max = 40000000

gal = Table.read("/sps/lsst/users/namourou/web/desc/clusters/DC2_mask/galaxies.fits")
ra_min, ra_max = min(gal['ra']), max(gal['ra'])
dec_min, dec_max = min(gal['dec']), max(gal['dec'])
#ra_min, ra_max = 60,60.5
#dec_min, dec_max = -36.4, -36

# Définir la taille des carrés (en degrés)
cote_carré = 0.01  # Vous pouvez ajuster cette valeur

# Créer une grille de carrés
ra_bins = np.arange(ra_min, ra_max, cote_carré)
dec_bins = np.arange(dec_min, dec_max, cote_carré)

# Créer des tableaux pour stocker les densités de galaxies
densités = np.zeros((len(ra_bins) - 1, len(dec_bins) - 1), dtype=float)
dens = np.histogram2d(gal['ra'], gal['dec'], bins = (ra_bins,dec_bins))[0]
ra_pt = []
for i in range(len(ra_bins)-1):
    ra_pt.append((ra_bins[i]+ra_bins[i+1])/2) 
dec_pt = []
for i in range(len(dec_bins)-1):
    dec_pt.append((dec_bins[i]+dec_bins[i+1])/2) 
ra_pt_grid, dec_pt_grid = np.meshgrid(ra_pt, dec_pt)
mu = np.mean(dens)
sigma = np.std(dens)
print(mu,sigma)
indices_inférieurs_à_la_moyenne = np.where(dens.T < mu-2*sigma)
# Maintenant, extrayez les valeurs de RA et Dec correspondant aux indices inférieurs à la moyenne
ra_mask = ra_pt_grid[indices_inférieurs_à_la_moyenne]
dec_mask = dec_pt_grid[indices_inférieurs_à_la_moyenne]
print("first plot")
fig, ax = plt.subplots(figsize=(10,6))
#plt.hexbin(ra[val_min:val_max], dec[val_min:val_max], C=mask[vpix[val_min:val_max]], gridsize = 2500, vmin=24.9, vmax=25.1)
plt.scatter(ra_r[val_min:val_max], dec_r[val_min:val_max], c=mask_r[vpix_r[val_min:val_max]], vmin=26, vmax=27)
plt.colorbar()
#plt.scatter(gal_nextd['ra'],gal_nextd['dec'], s=1, alpha = 1, color = 'red')
#plt.scatter(ra_mask, dec_mask, s=5, color = 'white')
#plt.scatter(max(ra_bins),max(dec_bins), color = 'black' )
plt.xlabel("ra")
plt.ylabel("dec")
plt.title('galaxies in tract')
#plt.xlim([60,60.5])
#plt.ylim([-36.4,-36])
#for i in range(len(ra_mask)):
#    rectangle = plt.Rectangle((ra_mask[i]-0.5*cote_carré, dec_mask[i]-0.5*cote_carré), cote_carré, cote_carré, fill=True, color='black', alpha = .1)
#    ax.add_patch(rectangle)
plt.xlim([60.45,60.1])
plt.ylim([-36.35,-36.10])
#plt.legend()
plt.savefig("u_band_mag_fig.png")
plt.close()
print('sec plot')
fig, ax = plt.subplots(figsize=(10,6))
#plt.hexbin(ra[val_min:val_max], dec[val_min:val_max], C=mask[vpix[val_min:val_max]], gridsize = 2500, vmin=24.9, vmax=25.1)
plt.scatter(ra_i[val_min:val_max], dec_i[val_min:val_max], c=mask_i[vpix_i[val_min:val_max]])#, vmin=24.9, vmax=25.1)
plt.colorbar()
#plt.scatter(gal_nextd['ra'],gal_nextd['dec'], s=1, alpha = 1, color = 'red')
#plt.scatter(ra_mask, dec_mask, s=5, color = 'white')
#plt.scatter(max(ra_bins),max(dec_bins), color = 'black' )
plt.xlabel("ra")
plt.ylabel("dec")
plt.title('galaxies in tract')
#plt.xlim([60,60.5])
#plt.ylim([-36.4,-36])
#for i in range(len(ra_mask)):
#    rectangle = plt.Rectangle((ra_mask[i]-0.5*cote_carré, dec_mask[i]-0.5*cote_carré), cote_carré, cote_carré, fill=True, color='black', alpha = .1)
#    ax.add_patch(rectangle)
plt.xlim([60.45,60.1])
plt.ylim([-36.35,-36.10])
plt.savefig("y_band_mag_fig.png")
#plt.legend()
plt.close()