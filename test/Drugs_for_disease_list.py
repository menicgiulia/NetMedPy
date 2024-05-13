# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 19:14:22 2024

@author: aalda
"""

import numpy as np
import random
import pandas as pd
import matplotlib.pyplot as plt

dbFile = "C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/data/Drug_Bank/DB_Drug_Targets_2023.csv"

synom_file = "C:/Users/aalda/Dropbox (Personal)/Barabasi Lab/data/Drug_Bank/DB_Drug_Dictionary_2023.csv"
synoms = pd.read_csv(synom_file)
synoms = synoms[['DB_id','Synonym']]

db = pd.read_csv(dbFile)
db = db.query("organism=='Humans'")
db = db[db['Status'].str.contains('approved',case=False)]

db = pd.merge(db,synoms, how='left',on='DB_id')


dis = 'lung cancer'


dis_dg = db[db['Indication'].str.contains(dis,case=False,na=False)]
dis_dg = dis_dg[['Gene_Target','Name']]
dis_dg = dis_dg.drop_duplicates()


targets_per_drug = dis_dg.groupby('Name').size().reset_index(name='target_count')
targets_per_drug = targets_per_drug.sort_values(by='target_count', ascending=False)


# for d in targets_per_drug.Name:
#     q = db.query("Name==@d")
#     q = q.iloc[0,]
    
#     print(f"******************{d}******************")
#     print("Description ")
#     print(q.Description)
#     print("Indication")
#     print(q.Indication)

drug = 'Methylphenidate'
dinfo = db.query('Name==@drug or Synonym==@drug')

print(f"******************{drug}******************")
print("Description ")
print(dinfo.iloc[0].Description)
print("Indication")
print(dinfo.iloc[0].Indication)

targets = dinfo[['Gene_Target','Name']]
targets = targets.drop_duplicates()

diseases = ['inflammation','coronary artery disease',
            'chronic obstructive pulmonary disease','lung neoplasm',
            'factor ix', 'fragile x']

drugs = ['Ibuprofen','Dexamethasone',
         'Atorvastatin','Amlodipine',
         'Salmeterol','Fluticasone furoate',
         'Pralsetinib','Nintedanib',
         'Coagulation Factor IX Human','Nonacog beta pegol',
         'Risperidone','Fluoxetine']


prescribed_drugs = {'inflammation':['Ibuprofen','Dexamethasone'],
                    'coronary artery disease':['Atorvastatin','Amlodipine'],
                    'chronic obstructive pulmonary disease':['Salmeterol','Fluticasone furoate'],
                    'lung neoplasm':['Pralsetinib','Nintedanib'],
                    'factor ix':['Coagulation Factor IX Human','Nonacog beta pegol'],
                    'fragile x':['Risperidone','Fluoxetine']
                    }


df = pd.DataFrame(columns = ['Disease','Prescription','Drug_Targets'])

for d in prescribed_drugs.keys():
    for p in prescribed_drugs[d]:
        
        targets = db.query('Name==@p')
        targets = targets[['Name','Gene_Target']]
        targets = targets.drop_duplicates()
        l = len(targets.Gene_Target)
        
        df.loc[len(df.index)] = [d,p,l]
        
df.to_csv('data/drug_disease_list/drug_list.csv',index=False)
df.to_excel('data/drug_disease_list/drug_list.xlsx', index=False, engine='openpyxl')
