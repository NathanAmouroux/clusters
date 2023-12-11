from matplotlib import pylab as plt
from astropy.io import fits
import numpy as np
import os, glob, sys
from astropy.table import Table
from astropy.io import fits

tract = sys.argv[1]

inpath = "/sps/lsst/groups/clusters/amico_validation_project/catalogs/AMICO/raw_amico_cats/DC2_v0_yband/masks/output/MASK/"
raw = fits.open(inpath + f"{tract}_mask_binary.fits")
buff = fits.open(inpath + f"{tract}_mask_binary_buffer.fits")

gals = Table.read(f"/sps/lsst/users/tguillem/web/clusters/catalogs/DC2_photoz_flexzboost/v0/{tract}/galaxies.fits")


ra_min, ra_max, step1 = raw[0].header["START_1"], raw[0].header["END_1"], raw[0].header["STEP_1"]
n_bin1 = int((ra_max-ra_min)/step1)
dec_min, dec_max, step2 = raw[0].header["START_2"], raw[0].header["END_2"], raw[0].header["STEP_2"]
n_bin2 = int((dec_max-dec_min)/step2)
ra = np.linspace(ra_min, ra_max, n_bin1)
dec =  np.linspace(dec_min, dec_max, n_bin2)
ra_m, dec_m = np.meshgrid(ra, dec)

plt.figure(figsize=(10,10))
plt.imshow(raw[0].data, extent = [ra_min, ra_max, dec_max, dec_min])
plt.scatter(gals['ra'], gals["dec"],s=.5, alpha = .1, label = "galaxies")
plt.xlim([ra_min, ra_max])
plt.ylim([dec_min, dec_max])
plt.legend()
plt.savefig(f"output/{tract}_masks.png", bbox_inches="tight")
plt.close()

ra_min, ra_max, step1 = buff[0].header["START_1"], buff[0].header["END_1"], buff[0].header["STEP_1"]
n_bin1 = int((ra_max-ra_min)/step1)
dec_min, dec_max, step2 = buff[0].header["START_2"], buff[0].header["END_2"], buff[0].header["STEP_2"]
n_bin2 = int((dec_max-dec_min)/step2)
ra = np.linspace(ra_min, ra_max, n_bin1)
dec =  np.linspace(dec_min, dec_max, n_bin2)
ra_m, dec_m = np.meshgrid(ra, dec)

plt.figure(figsize=(10,10))
plt.imshow(buff[0].data, extent = [ra_min, ra_max, dec_max, dec_min])
plt.scatter(gals['ra'], gals["dec"],s=.5, alpha = .1, label = "galaxies")
plt.xlim([ra_min, ra_max])
plt.ylim([dec_min, dec_max])
plt.legend()
plt.savefig(f"output/{tract}_masks_Buffer.png", bbox_inches="tight")
plt.close()
