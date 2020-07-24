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

name =[x.GetProp("_Name") for x in query]
print(name[0])

query_fp = MACCSkeys.GenMACCSKeys(query[0])

print(query_fp)

ta, sim, na = [], [], []
count = 0
for n in database:
    try:
        smile = Chem.MolToSmiles(n)
        mol=Chem.MolFromSmiles(smile)
        mol_fp =MACCSkeys.GenMACCSKeys(mol)
        similarity = DataStructs.TanimotoSimilarity(mol_fp, query_fp)
        count = count + 1
        if similarity>0.75:
            smile= Chem.MolToSmiles(n)
            name = n.GetProp("_Name")
            na.append(name)
            ta.append(smile)
            sim.append(similarity)
            print(name,smile, similarity, count)
    except:
        continue

d = {'Name':na,'Smile':ta, 'Similarity':sim}
df_final = pd.DataFrame(data=d)
df_final = df_final.sort_values('Similarity', ascending=False)

df_final.to_csv('MACCS.csv', index=False, sep=',')