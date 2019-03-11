import pandas as pd
import re  # find digits in a string
import matplotlib.pyplot as plt
import aflowlib as aflow
import asyncio
from concurrent.futures import ThreadPoolExecutor
import nest_asyncio
nest_asyncio.apply()
      

# specify your path where to save your data
path = 'E:\\work\\Results_aflow\\GitHub\\'
# and desierd anion and cations
cation = 'Mn'
anion = 'F'
#if you already downloaded data after first launch, make it True
data_saved = False  

if not data_saved:
  db_links = aflow.get_links(cation, anion)
  #better to save links
  db_links.to_csv(path+cation+anion+'_links'+'.csv', index=False)
  #uncomment if saved and comment upper one
#        db_links = pd.read_csv(path+cation+anion+'_links'+'.csv')

  #get full data for each compound
  db = asyncio.run(aflow.get_data_asynchronous(db_links)) 

  #assign oxidation states based on the number of anions
  db['oxidation_state'] = db.apply(lambda row: aflow.oxidation_state(row, cation, anion), axis=1)

  #finds all digits in 'aurl' and takes only last one
  db['ICSD'] = db.aurl.apply(lambda row: re.findall(r'\d+', row)[-1])

  # next string splits Bader charges for anion's and cation's charges
  db[['cation_charges', 'anion_charges']]  = db.apply(lambda row: aflow.split_bader(row), axis=1).apply(pd.Series)  # bader returns tuple. "apply" makes it data series
  #drop records without Bader charge
  db = db.dropna(subset=['cation_charges']) 
  #save data   
  db.to_csv(path+cation+anion+'_full'+'.csv', index=False)

else:
  db = pd.read_csv(path+cation+anion+'_full'+'.csv')
  #converts string to list of floats
  db.cation_charges  = db.cation_charges.apply(lambda row: aflow.bader_to_list_of_floats(row))

#converting database of compounds to database of atoms and drop columns
db_for_each = aflow.to_db_with_bader_for_each_atom(db, atom_type='cation')

#plot    
aflow.plot_os_vs_bader(db_for_each, cation, anion)
    
    
