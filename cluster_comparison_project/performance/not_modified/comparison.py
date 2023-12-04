#!/usr/bin/env python
# coding: utf-8

###import
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import sys
import os

###clevar
import clevar
from clevar.catalog import ClCatalog

inpath = '/sps/lsst/groups/clusters/amico_validation_project/catalogs/matching_cats/'
cat_path = inpath + str(sys.argv[1])

print('Using catalogs located at :', cat_path)

# #---Select catalogs to match---##

outpath = '/pbs/home/n/namourou/test_jupyter/cluster_challenge/plots/'
outpath += str(sys.argv[1])

print('Saving plots at :', outpath)

if os.path.exists(outpath):
     shutil.rmtree(outpath)
os.makedirs(outpath)

##########select case
if str(sys.argv[1]) == 'p_matching':
    matching = 'p'
elif str(sys.argv[1]) == 'mb_matching':
    matching = 'mb'
catalog1 = 'c1_' + matching + '.fits'
catalog2 = 'c2_' + matching + '.fits'
##########

#load c1, c2, c3 and c4
c1 = ClCatalog.read_full(cat_path + catalog1)
c2 = ClCatalog.read_full(cat_path + catalog2)

#create a merged catalog for the cross-matched pairs
output_matched_catalog(matching_folder_1+catalog1, matching_folder_1+catalog2,matching_folder_1+'output_catalog_12.fits', c1, c2, matching_type='cross', overwrite=True)
c_merged_12 = ClCatalog.read(matching_folder_1+'output_catalog_12.fits', 'merged',  z_cl='cat1_z', richness = 'cat1_mass', z_halo='cat2_z', mass = 'cat2_mass', log_mass = 'cat2_log_mass', m200c = 'cat2_m200c', log_m200c = 'cat2_log_m200c',)
c_merged_12 = c_merged_12[c_merged_12.data['z_cl']<1.15]

output_matched_catalog(matching_folder_2+catalog1, matching_folder_2+catalog2,matching_folder_2+'output_catalog_34.fits', c3, c4, matching_type='cross', overwrite=True)
c_merged_34 = ClCatalog.read(matching_folder_2+'output_catalog_34.fits', 'merged',  z_cl='cat1_z', richness = 'cat1_mass', z_halo='cat2_z', mass = 'cat2_mass', log_mass = 'cat2_log_mass', m200c = 'cat2_m200c', log_m200c = 'cat2_m200c',)

#mass plot
bin_range = [0,60]
nbins = 60
plot_name = ""
for i in range(0,3):
     if(i==0):
          plot_name = 'mass_all.png'
          c_cut_12 = c_merged_12
          c_cut_34 = c_merged_34
     if(i==1):
          plot_name = 'mass_lt1014.png'
          c_cut_12 = c_merged_12[c_merged_12.data['m200c']<10**14]
          c_cut_34 = c_merged_34[c_merged_34.data['m200c']<10**14]
     if(i==2):
          plot_name = 'mass_gt1014.png'
          c_cut_12 = c_merged_12[c_merged_12.data['m200c']>10**14]
          c_cut_34 = c_merged_34[c_merged_34.data['m200c']>10**14]
     plt.figure()
     plt.hist(c_cut_12.data['richness'], range=bin_range, bins=nbins, label='RM', histtype='step', color = 'red')
     plt.hist(c_cut_34.data['richness'], range=bin_range, bins=nbins, label='WaZP', histtype='step', color = 'black')
     plt.xlabel("alg. richness");
     plt.ylabel("clusters")
     plt.grid(which='major', axis='both', linestyle='-', linewidth='0.1', color='grey')
     plt.grid(which='minor', axis='both', linestyle=':', linewidth='0.1', color='grey')
     plt.legend(title = '', loc='upper right')
     plt.title('cosmoDC2')
     plt.savefig(outpath+plot_name, bbox_inches='tight')
     plt.close()

#redshift plot
bin_range = [0,1.6]
nbins = 32 
plot_name = ""
for i in range(0,3):
     if(i==0):
          plot_name = 'redshift_all.png'
          c_cut_12 = c_merged_12
          c_cut_34 = c_merged_34
     if(i==1):
          plot_name = 'redshift_lt1014.png'
          c_cut_12 = c_merged_12[c_merged_12.data['m200c']<10**14]
          c_cut_34 = c_merged_34[c_merged_34.data['m200c']<10**14]
     if(i==2):
          plot_name = 'redshift_gt1014.png'
          c_cut_12 = c_merged_12[c_merged_12.data['m200c']>10**14]
          c_cut_34 = c_merged_34[c_merged_34.data['m200c']>10**14]
     plt.figure()
     plt.hist(c_cut_12.data['z_cl'], range=bin_range, bins=nbins, label='RM', histtype='step', color = 'red')
     plt.hist(c_cut_34.data['z_cl'], range=bin_range, bins=nbins, label='WaZP', histtype='step', color = 'black')
     plt.xlabel("redshift");
     plt.ylabel("clusters / 0.05 dz")
     plt.grid(which='major', axis='both', linestyle='-', linewidth='0.1', color='grey')
     plt.grid(which='minor', axis='both', linestyle=':', linewidth='0.1', color='grey')
     plt.legend(title = '', loc='upper right')
     plt.title('cosmoDC2')
     plt.savefig(outpath+plot_name, bbox_inches='tight')
     plt.close()
