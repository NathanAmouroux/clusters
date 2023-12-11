# !/usr/bin/env python
# coding: utf-8
from astropy.table import Table, vstack, hstack
import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
import numpy.ma as ma
import sys
curr_tile = str(sys.argv[1])
#healpath_cdc = '/sps/lsst/users/tguillem/web/clusters/catalogs/cosmoDC2_photoz_flexzboost/v1/'
healpath_dc =  "/sps/lsst/users/tguillem/web/clusters/catalogs/DC2_photoz_flexzboost/v0/" + curr_tile + "/galaxies.fits"
halospath_dc = "/sps/lsst/groups/clusters/cluster_comparison_project/before_matching/halos/cosmoDC2/DC2.masked/Catalog.fits"
clusterspath_dc = "/sps/lsst/groups/clusters/amico_validation_project/catalogs/AMICO/raw_amico_cats/DC2_v0_yband/map_detections_refined_noBuffer_all_noDoubles.fits"
outpath = "/sps/lsst/groups/clusters/amico_validation_project/catalogs/DC2/masks/"

halos_dc = Table.read(halospath_dc)
clusters_dc = Table.read(clusterspath_dc)
gals_dc = Table.read(healpath_dc)
ra_min, ra_max = min(gals_dc["ra"]), max(gals_dc["ra"])
dec_min, dec_max = min(gals_dc["dec"]), max(gals_dc["dec"])

cote_carré = 0.01  # Vous pouvez ajuster cette valeur

# Créer une grille de carrés
ra_bins = np.arange(ra_min, ra_max + cote_carré, cote_carré)
dec_bins = np.arange(dec_min, dec_max + cote_carré, cote_carré)

# Créer des tableaux pour stocker les densités de galaxies
dens_dc = np.histogram2d(gals_dc['ra'], gals_dc['dec'], bins = (ra_bins,dec_bins))[0]
ra_grid, dec_grid = np.meshgrid(ra_bins, dec_bins)
mu_dc = np.mean(dens_dc)
sigma_dc = np.std(dens_dc)
print(mu_dc,sigma_dc)

d_flat_dc=dens_dc.ravel()

ra_pt = []
for i in range(len(ra_bins)-1):
    ra_pt.append((ra_bins[i]+ra_bins[i+1])/2) 

dec_pt = []
for i in range(len(dec_bins)-1):
    dec_pt.append((dec_bins[i]+dec_bins[i+1])/2) 
    
ra_pt_grid, dec_pt_grid = np.meshgrid(ra_pt, dec_pt)
indices_inférieurs_à_la_moyenne_dc = np.where(dens_dc.T < mu_dc-3*sigma_dc)

# Maintenant, extrayez les valeurs de RA et Dec correspondant aux indices inférieurs à la moyenne
ra_mask_dc = ra_pt_grid[indices_inférieurs_à_la_moyenne_dc]
dec_mask_dc = dec_pt_grid[indices_inférieurs_à_la_moyenne_dc]

# Exemple de coordonnées (RA, Dec) en degrés
ra1, dec1 = min(gals_dc["ra"]), min(gals_dc["dec"])
ra2, dec2 = max(gals_dc["ra"]), max(gals_dc["dec"])
area_tot = abs(ra2-ra1)*abs(dec2-dec1)
print(f"Surface sur le ciel du tract : {area_tot:.2f} degrés carrés")

surface_mask_dc = len(ra_mask_dc)*(cote_carré**2)
s_frac_dc = surface_mask_dc/area_tot
    
masked_gals_dc = []
masked_halos_dc = []
masked_clusters_dc = []
for i in range(len(ra_mask_dc)):
    ra_min = ra_mask_dc[i]-0.5*cote_carré
    ra_max = ra_mask_dc[i]+0.5*cote_carré
    dec_min = dec_mask_dc[i]-0.5*cote_carré
    dec_max = dec_mask_dc[i]+0.5*cote_carré
    masked_gals = gals_dc[(gals_dc['ra'] >= ra_min) & (gals_dc['ra'] <= ra_max) & \
                    (gals_dc['dec'] >= dec_min) & (gals_dc['dec'] <= dec_max)]
    masked_halos = halos_dc[(halos_dc['ra_cl'] >= ra_min) & (halos_dc['ra_cl'] <= ra_max) & \
                    (halos_dc['dec_cl'] >= dec_min) & (halos_dc['dec_cl'] <= dec_max)]
    masked_clusters = clusters_dc[(clusters_dc['Xphys'] >= ra_min) & (clusters_dc['Xphys'] <= ra_max) & \
                    (clusters_dc['Yphys'] >= dec_min) & (clusters_dc['Yphys'] <= dec_max)]
    masked_gals_dc.append(len(masked_gals))
    masked_halos_dc.append(len(masked_halos))
    masked_clusters_dc.append(len(masked_clusters))
    

n_gal_frac_dc = sum(masked_gals_dc)/len(gals_dc)
n_halo_frac_dc = sum(masked_halos_dc)/len(halos_dc)
n_clusters_frac_dc = sum(masked_clusters_dc)/len(clusters_dc)
sumary = {"cat_name" : ["DC2"], "masked surface (deg**2)" : [surface_mask_dc], "fraction of surface masked" : [round(s_frac_dc,4)]
          ,"n_masked_gals" : [sum(masked_gals_dc)], "frac_masked_gals" : [round(n_gal_frac_dc,4)],
         "n_masked_halos" : [sum(masked_halos_dc)], "frac_masked_halos" : [round(n_halo_frac_dc,4)],
         "n_masked_clusters" : [sum(masked_clusters_dc)], "frac_masked_clusters" : [round(n_clusters_frac_dc,4)],
         "box_size" : [0.01]}
sumary = Table(sumary)
sumary.write(outpath + curr_tile +"_infos.fits")
mask_coords = Table({'ra':ra_mask_dc, 'dec':dec_mask_dc})
mask_coords.write(outpath + curr_tile+"_mask_coords.fits")