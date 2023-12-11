import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
gal = Table.read("/sps/lsst/users/namourou/web/desc/clusters/DC2_mask/galaxies.fits")
plt.figure(figsize=(8,6))
plt.scatter(gal['ra'],gal['dec'], s=0.01, alpha = .5, color = 'red')
plt.xlabel("ra")
plt.ylabel("dec")
plt.title('galaxies in tract')
plt.savefig('copy_hehe.png', bbox_inches='tight')
#plt.legend()