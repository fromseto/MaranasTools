import pandas as pd
import json

df = pd.read_excel("metabolite_details_to_dict_input.xlsx",header=None)

met_ids = df[0].tolist()
attribute_names = df[1].tolist()
attribute_values = df[2].tolist()

met_detail_dict = {}

for idx,met in enumerate(met_ids):
	if met not in met_detail_dict.keys():
		met_detail_dict[met] = {}
	met_detail_dict[met][attribute_names[idx]] = attribute_values[idx]
# dictionary = {k: {k2: g2 for k2,g2 in df.groupby(1)} for k,g in df.groupby(0)}

elements = ['C','H','N','O','P','S','F','Cl','Mg','Fe','Se','Co','As','Br',
                'I','R']
# nonsense_detail = []
for met in met_detail_dict.keys():
	val = 0
	for elm in elements:
		val = val + met_detail_dict[met][elm]
	if val ==0:
		# nonsense_detail.append(met)
		met_detail_dict.pop(met, None)

with open('met_detials_dict_v2.json', 'w') as fp:
    json.dump(met_detail_dict, fp)
