import pandas as pd
import requests
import sys
import json

ids = pd.read_csv('/home/pelmo/sscrofa_variant_remapping/data/starting_ids.csv', header = None)
ids.columns = ['ensembl_id_10']
ids

server = "http://rest.ensembl.org"
id_dict = {}
for id in ids['ensembl_id_10']:
    ext = "/archive/id/" + id + "?"
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    if not r.ok:
      r.raise_for_status()
      sys.exit()
    decoded = r.json()
    id_dict[id] = decoded
id_df = pd.DataFrame.from_dict(id_dict).transpose()
id_df.to_csv('/home/pelmo/sscrofa_variant_remapping/data/id_investigation.csv')
