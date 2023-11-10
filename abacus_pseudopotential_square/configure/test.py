import psedupotential_parse as psp
import json
dict1, dict2 = psp.StandardHtmlParser(
    FileName='../pseudopotentials/resources/sg15_10/In_ONCV_PBE-1.0.upf',
    header='<UPF version="2.0.1">',
    footer='</UPF>'
)

with open('test1.json', 'w') as f:
    json.dump(dict1, f, indent=4)
with open('test2.json', 'w') as f:
    json.dump(dict2, f, indent=4)