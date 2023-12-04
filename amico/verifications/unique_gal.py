from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
import pandas as pd

inpath = "/sps/lsst/groups/clusters/cluster_comparison_project/before_matching/"
outpath = "/sps/lsst/groups/clusters/amico_validation_project/catalogs/"

c2_mb = Table.read(inpath + "amico/DC2.fzb.magy/v0/Catalog_members.fits")
# Convertissez la table Astropy en un DataFrame Pandas.
df = c2_mb.to_pandas()

# Utilisez Pandas pour sélectionner les entrées avec des IDs uniques.
df_unique = df.drop_duplicates(subset='id_mb')

# Convertissez le DataFrame Pandas résultant en une table Astropy.
nouvelle_table = Table.from_pandas(df_unique)
nouvelle_table.write(outpath + 'testsubset.fits')