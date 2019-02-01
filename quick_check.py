import pandas as pd
import requests  # read url
import re  # find digits in a string
import matplotlib.pyplot as plt
from mendeleev import element


def get_links(cation, anion):
    '''get_links()
    return dataframe of links to compounds in AFLOW database
    for all given available compounds'''
    # we need to download first page, to get total number of compounds
    print('Download all links to available binaries in the AFLOW...')
    link = "http://aflowlib.duke.edu/search/API/?species({},{}),$nspecies(2),paging(1)".format(cation, anion)
    db_links = pd.read_json(link, orient='index')
    total_number_of_compounds = db_links.index[0].split()[2]
    number_of_pages = 1 + int(int(total_number_of_compounds) / 64)  # website gives only 64 compounds per page
    if number_of_pages > 1:  # for fancy looking
        for page in range(2, number_of_pages+1):
            print('page number', page, 'from', number_of_pages)
            link_temp = link[:-2] + str(page) + ')'
            db_links = db_links.append(pd.read_json(link_temp, orient='index'))
    print('Done.')
    print(" ".join(['In total', total_number_of_compounds,cation,anion,'binaries']))
    return db_links


def bader_for_quick(row):
    '''download Bader charges'''
    cation_bader = []
    anion_bader = []  
    
    link_to_get_charges = 'http://aflowlib.duke.edu/' + row['aurl'][18:] + '/?bader_net_charges'
    charges = requests.get(link_to_get_charges)
    
    # skip compounds without computed Bader charges
    if charges.headers['Content-Length']=='0':
        return None
    
    for charge in charges.text.strip().split(sep=','):  # list of string of bader charges
        if float(charge) > 0:
            cation_bader.append(float(charge))
        else:
            anion_bader.append(float(charge))

    return cation_bader, anion_bader


def oxidation_state(row, cation, anion):
    '''calculates mean oxidation state in a compound, 
    based on number of anion and cations'''
    
    electrons_taken_by_groupid = {15:3, 16:2, 17:1}
    number_of_electron_taken = electrons_taken_by_groupid[element(anion).group_id] #oxidation state of anion.
    
    # in 'compound' elements in alphabetic order
    spicies = sorted([cation, anion])
    o = 0
    me = 0
    #map returns object, we make it list
    consequiteve_digits = list(map(int, re.findall(r'\d+', row['compound'])))
    if anion == spicies[0]:
        o = consequiteve_digits[0]
        me = consequiteve_digits[1]
    elif cation == spicies[0]:
        o = consequiteve_digits[1]
        me = consequiteve_digits[0]
    return number_of_electron_taken * o/me


def bader_to_list_of_floats(bader):
    charges = []
    # if problems try str(bader) instead bader in for loop
    for charge in bader.split(sep=','):
        try:
            charges.append(float(charge.replace('\'', '').replace('[', '').replace(']', '').strip()))
        except ValueError:
            print("ERROR")
    return charges


def to_db_with_bader_for_each_atom(db, atom_type='cation'):
    '''Returns Dataframe whrere each row correspons to 
    an atom with bader charge and estimated oxidation state'''
    rows_list = []
    for db_i in db.index:
        if atom_type == 'cation':
            number_of_atoms = len(db['cation_charges'][db_i])
        elif atom_type == 'anion':
            number_of_atoms = len(db['anion_charges'][db_i])
            
        for chrg_j in range(number_of_atoms):
            dict1 = {}
            dict1.update({'charge': db[atom_type+'_charges'][db_i][chrg_j],
                          'oxidation_state': db['oxidation_state'][db_i],
                          'auid': db['auid'][db_i],
                          'aurl': db['aurl'][db_i],
                          'species': db['species'][db_i],
                          'compound': db['compound'][db_i]})
            rows_list.append(dict1)
    return pd.DataFrame(rows_list)

def plot_os_vs_bader(db):

    size=100
    #dict for drawing scatter plot for each oxidation state. 9 for mixed
    dict_of_df =  {1:{'marker':'o', 'linewidth':1, 'color':'b', 'edgecolors':'None', 'size':size*0.7},
                   2:{'marker':'+', 'linewidth':2, 'color':'r', 'edgecolors':'r', 'size':size},
                   3:{'marker':'2', 'linewidth':2, 'color':'g', 'edgecolors':'r', 'size':size*1.3},
                   4:{'marker':'d', 'linewidth':2, 'color':'y', 'edgecolors':'blue', 'size':size*0.7},
                   5:{'marker':'p', 'linewidth':1, 'color':'k', 'edgecolors':'r', 'size':size*0.7},
                   6:{'marker':'H', 'linewidth':1, 'color':'c', 'edgecolors':'r', 'size':size},
                   7:{'marker':'*', 'linewidth':1, 'color':'m', 'edgecolors':'r', 'size':size},
                   8:{'marker':'$8$', 'linewidth':1, 'color':'purple', 'edgecolors':'black', 'size':size},
                   9:{'marker':'s', 'linewidth':1, 'color':'b', 'edgecolors':'r', 'size':size}}

    db_single_valence = db[db.oxidation_state % 1 == 0]
    db_mixed_valence  = db[db.oxidation_state % 1 != 0]
    
    for os in db_single_valence.oxidation_state.unique():
        os = int(os)
        number_of_atoms = db_single_valence[db_single_valence.oxidation_state==os].shape[0]
        plt.scatter([0]*number_of_atoms,
                    db_single_valence[db_single_valence.oxidation_state==os].charge,
                    marker = dict_of_df[os]['marker'],
                    s=dict_of_df[os]['size'], alpha=0.5, c=dict_of_df[os]['color'], 
                    label = str(os) + ' oxidation state, '+ str(number_of_atoms) +' atoms ')


if __name__ == '__main__':
    # specify your path where to save your data
    path = 'E:\\work\\Results_aflow\\GitHub\\'
    # and desierd anion and cations
    cation = 'Fe'
    anion = 'F'
    
    
    db = get_links(cation, anion)
    #better to save links
    db.to_csv(path+cation+anion+'links'+'.csv')
    #uncomment if saved and comment upper one
    db = pd.read_csv(path+cation+anion+'links'+'.csv')
    db['oxidation_state'] = db.apply(lambda row: oxidation_state(row, cation, anion), axis=1)
    
    #    # next string gets bader charges for each compound
    db[['cation_charges', 'anion_charges']]  = db.apply(lambda row: bader_for_quick(row), axis=1).apply(pd.Series)  # bader returns tuple. "apply" makes it data series
    
    db = db.dropna(subset=['cation_charges'])   
    db.to_csv(path+cation+anion+'.csv')
    
    db_for_each = to_db_with_bader_for_each_atom(db, atom_type='cation')
    
    fig, ax = plt.subplots()
    plot_os_vs_bader(db_for_each)
    plt.legend()
    plt.show()
