import rdkit
from rdkit import Chem
import os
import pandas as pd
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import DataStructs
from rdkit.Chem import CanonSmiles
from rdkit.Chem import MACCSkeys
from rdkit.Chem import AllChem

os.chdir('C:\\Users\mradaeva\Desktop\SOFT\Similarity search')

database = Chem.SmilesMolSupplier("merged_all.smi")
query = Chem.SmilesMolSupplier("query.smi")

query_smile=Chem.MolToSmiles(query[0])
query = Chem.MolFromSmiles(query_smile)

query_fp = AllChem.GetMorganFingerprint(query,1)


print(query_fp)

na, ta, sim = [],[], []
count = 0
for n in database:
    try:
        smile = Chem.MolToSmiles(n)
        mol=Chem.MolFromSmiles(smile)
        mol_fp =AllChem.GetMorganFingerprint(mol,1)
        similarity = DataStructs.TanimotoSimilarity(mol_fp, query_fp)
        count = count + 1
        if similarity>0.55:
            smile= Chem.MolToSmiles(n)
            name = n.GetProp("_Name")
            ta.append(smile)
            na.append(name)
            sim.append(similarity)
            print(name, smile, similarity, count)
    except:
        continue

d = {'Name':na,'target':ta, 'Similarity':sim}
df_final = pd.DataFrame(data=d)
df_final = df_final.sort_values('Similarity', ascending=False)

df_final.to_csv('Morgan.csv', index=False, sep=',')
