"""Map the 204 GBD countries and territories to their GBD super-region and region.
Inputs: the location hierarchy and English location metadata of the GBD results tool (JSON), and one GBD 2023 results file
that lists the 204 analysed countries and territories (IHME-GBD_2023_DATA-aaf6002d-1.zip).
Usage: python3 -I build_location_map.py <hierarchy_json> <metadata_json> <gbd_2019_2023_dir> <out_csv>
"""
import sys, io, json, zipfile
from pathlib import Path
import pandas as pd

hier, meta, gdir, out = sys.argv[1:5]
tree = json.load(open(hier, encoding='utf-8'))['data']['locations']
names = json.load(open(meta, encoding='utf-8'))['data']['location']
z = zipfile.ZipFile(Path(gdir) / 'IHME-GBD_2023_DATA-aaf6002d-1.zip')
g = pd.read_csv(io.BytesIO(z.read([n for n in z.namelist() if n.endswith('.csv')][0])), usecols=['location_id'])
countries = set(g.location_id) - {1}
assert len(countries) == 204
SUPER = {31, 64, 103, 137, 158, 4, 166}
rows = {}
def walk(node, sr=None, region=None):
    i = node['id']
    if i in SUPER: sr = i
    elif sr is not None and region is None: region = i
    if i in countries and i not in rows:
        rows[i] = dict(location_id=i, location_name=names[str(i)]['name'], super_region=names[str(sr)]['name'], super_region_id=sr, region=names[str(region)]['name'])
    for k in node.get('children', []): walk(k, sr, region)
for top in tree:
    if top['id'] == 1:
        for c in top['children']: walk(c)
m = pd.DataFrame(rows.values()).sort_values(['super_region_id', 'region', 'location_name'])
assert len(m) == 204, len(m)
m.to_csv(out, index=False)
print('locations', len(m), m.super_region.value_counts().to_dict())
