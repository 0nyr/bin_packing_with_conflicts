
#%%

import os
import pandas as pd


folder ='Elhedhli'


elhedhli_dictionary = {1: 'u120', 2: 'u250', 3: 'u500', 5:'t60', 6:'t120', 7:'t249'}

def elhedhli_get_set_from_name(s):
    data = s.split('.')[0].split('_')
    
    inst_set = elhedhli_dictionary[int(data[1])]
    return inst_set

def elhedhli_get_density_from_name(s):
    data = s.split('.')[0].split('_')
    
    d = int(data[2])*10
    return d


#%%

df = pd.read_csv('time.csv', header=None)
df.columns = ['file_name', 'optimal', 'time']
df['density'] = df['file_name'].apply(elhedhli_get_density_from_name)
df['set'] = df['file_name'].apply(elhedhli_get_set_from_name)

df

# %%

instance_sets = [elhedhli_dictionary[k] for k in elhedhli_dictionary]
set_dfs = {k: df.loc[df['set'] == k] for k in instance_sets}

set_dfs

for k in set_dfs:
    if len(set_dfs[k]) != 0:
        print(f'set {k}:')
        print(set_dfs[k])
        print('\n')

# %%

set_dfs['u120']

#%%

set_dfs['u120'].loc[set_dfs['u120']['optimal'] == 1][['time']].describe()

#%% 

set_dfs['u250']

# %%
