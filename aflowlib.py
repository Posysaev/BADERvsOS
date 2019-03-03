import pandas as pd
import requests  # read url
import re  # find digits in a string
import matplotlib.pyplot as plt
from mendeleev import element
import asyncio
from concurrent.futures import ThreadPoolExecutor


def get_links(cation, anion):
    '''get_links()
    return dataframe of links to compounds in AFLOW database
    for all given available compounds'''
    # we need to download first page, to get total number of compounds
    print('Downloading all links to available binaries in the AFLOW...')
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
    return db_links.reset_index(drop=True)


def bader_for_quick(row):
    '''download Bader charges'''
    
    cation_bader = []
    anion_bader = []  
    
    link_to_get_charges = 'http://aflowlib.duke.edu/' + row['aurl'][18:] + '/?bader_net_charges'
    charges = requests.get(link_to_get_charges)
    
    
    # skip compounds without computed Bader charges
    if charges.headers['Content-Length']=='0':
        print('{:11d} ICSD: {:7s} No Bader charges for this record'.format(row.name, row['ICSD']))
        return None
    
    for charge in charges.text.strip().split(sep=','):  # list of string of bader charges
        if float(charge) > 0:
            cation_bader.append(float(charge))
        else:
            anion_bader.append(float(charge))
    print('{:11d} ICSD: {:7s} Bader charges are recieved'.format(row.name, row['ICSD']))
    return cation_bader, anion_bader


def split_bader(row):
    '''split Bader charges to anion's and cation's charges'''
    
    cation_bader = []
    anion_bader = []  
    
    charges = row.bader_net_charges  
    
    # !!!!!ADD HOW TO SKIP RECORDS WITHOUT BADER CHARGES skip compounds without computed Bader charges
#    if charges.headers['Content-Length']=='0':
#        print('{:11d} ICSD: {:7s} No Bader charges for this record'.format(row.name, row['ICSD']))
#        return None
    
    for charge in charges.strip().split(sep=','):  # list of string of bader charges
        if float(charge) > 0:
            cation_bader.append(float(charge))
        else:
            anion_bader.append(float(charge))
    print('{:11d} ICSD: {:7s} Bader charges are recieved'.format(row.name, row['ICSD']))
    return cation_bader, anion_bader


def get_data(db_list, link):
    '''Get all data available for a compound. Recieves json dictionary'''
#    print(link[18:])
    db_list.append(pd.read_json('http://aflowlib.duke.edu/' + link[18:] + '/?format=json', orient='columns')[:1])  
#    try:
#        return data[data['bader_net_charges'].notnull()]
#    except:
#        print('All records are without bader')
    return None
        

async def get_data_asynchronous(db_links):
    '''Aflow can send you data only for one compound per respond. In order to speed up the process we use 
    asynchronous programming and send up to 50 requests. Each respons is appended to db_list. Returns dataframe'''
    db_list = []
    with ThreadPoolExecutor(max_workers=50) as executor:    
        loop = asyncio.get_event_loop()
        print('Downloading all data to available binaries in the AFLOW...')
        tasks  = [loop.run_in_executor(executor, get_data, *(db_list, link)) for link in db_links['aurl'] ]
        for response in await asyncio.gather(*tasks):
            pass
        print('Done.')
    return pd.concat(db_list, ignore_index=True)


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
    an atom with bader charge and estimated oxidation state. The dataframe is sorted by OS'''
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
#                          'auid': db['auid'][db_i],
#                          'aurl': db['aurl'][db_i],
                          'species': db['species'][db_i],
                          'compound': db['compound'][db_i],
                          'ICSD': db['ICSD'][db_i]})
            rows_list.append(dict1)
    
    return pd.DataFrame(rows_list).sort_values(by=['oxidation_state', 'ICSD'])

def plot_os_vs_bader(db, cation, anion):

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
    number_of_compounds = db_single_valence.ICSD.unique().shape[0]
    fig, ax = plt.subplots()
#    plt.style.use('bmh')
#    plt.rcParams["font.family"] = "serif"
#    plt.rcParams["font.serif"] = "Times new roman"
#    plt.rcParams['axes.facecolor']='white'

    #plot each compound on x axis
    for j, icsd in enumerate(db_single_valence.ICSD.unique()):
        number_of_atoms = db_single_valence[db_single_valence.ICSD==icsd].shape[0]
        os = int(db_single_valence[db_single_valence.ICSD==icsd].oxidation_state.iloc[0])
#        print(int(db_single_valence[db_single_valence.ICSD==icsd].oxidation_state.iloc[0]))
        plt.scatter([j]*number_of_atoms,
                    db_single_valence[db_single_valence.ICSD==icsd].charge,
                    marker = dict_of_df[os]['marker'],
                    s=dict_of_df[os]['size'], alpha=0.5, c=dict_of_df[os]['color'], 
                    label = '')
        plt.axvline(j, color='k', linewidth = 0.1, linestyle='dashed')
        
    # plot all compunds
    for os in db_single_valence.oxidation_state.unique():
        os = int(os)
        number_of_all_atoms = db_single_valence[db_single_valence.oxidation_state==os].shape[0]
        plt.scatter([number_of_compounds]*number_of_all_atoms,
                    db_single_valence[db_single_valence.oxidation_state==os].charge,
                    marker = dict_of_df[os]['marker'],
                    s=dict_of_df[os]['size'], alpha=0.5, c=dict_of_df[os]['color'], 
                    label = str(os) + ' oxidation state, '+ str(number_of_all_atoms) +' atoms ')
#    plt.axvline(0, color='k', linewidth = 0.1, linestyle='dashed') 
            
    plt.xlabel('ICSD', fontsize=12, color='k')
    plt.ylabel('Bader charge ($e$)', fontsize=12)
    plt.title(cation+'$_x$' + anion+ '$_y$'+' binaries' + '\n'+
              'Total number of cations is {} in ({}) compounds'.format(
              db_single_valence.shape[0],
              db_single_valence.ICSD.unique().shape[0]))
    
    x1, x2, y1, y2 = plt.axis()
    x1 = -1
#    x2 = 0.4
    y1 = 0
    y2 = y2+0.5
    plt.axis((x1, x2, y1, y2))
    plt.xticks([i for i in range(db_single_valence.ICSD.unique().shape[0]+1)], 
                list(db_single_valence.ICSD.unique())+['All'], rotation=80)
    plt.legend(title='Legend', bbox_to_anchor=(1, 1))
    fig.set_size_inches(12,6)
    #plot intervals
    text=''
    for os in db_single_valence.oxidation_state.unique():
        text+='\n'+ str(os)+'  '+str(db_single_valence[db_single_valence.oxidation_state==os].charge.min())+'-'\
        +str(db_single_valence[db_single_valence.oxidation_state==os].charge.max())
        
    # these are matplotlib.patch.Patch properties
    props = dict(boxstyle='round', facecolor='green', alpha=0.15)
    # place a text box 
    plt.text(1.02, 0.5, 'Intervals'+text, transform=ax.transAxes, fontsize=14, bbox=props)

    plt.tight_layout()
    plt.show()
