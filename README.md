# BADERvsOS
Explore trend between Bader charges and oxidation states for binaries with AFLOW database and Python.
File 'Full_for_SV_and_MV.csv' includes final data set used in the article Comput. Mater. Sci., 2019, Accepted for 
publication. Currently available on chemrxiv.org/s/c9564778b1aa7523f855 

Note: asyncio library has been added. Asynchronous programming allows making multiple requests to AFLOWLIB server which allows receiving needed data much faster.

quick_results.py can be used to analyze trends between OS and Bader charges.

aflowlib.py contains all needed functions, make sure it imported properly.

Needed libraries: pandas, requests, re, matplotlib, mendeleev, asyncio, nest_asyncio

How to:
In quick_results.py specify 
1) path to save data 
2) desired anion and cation. 
In the end, the graph will appear.

TODO:
1. Analysis of mixed valence compounds
2. And analysis of structures by calculating coordination numbers, cure from duplicates, make sure that each atom has anion-cation bond
3. Calculation of oxidation states in mixed-valence compounds. 
![mv_mn_new](https://user-images.githubusercontent.com/43289846/52565604-f87f7300-2e0f-11e9-87cd-889a66b39084.png)


If you have any questions, please write to posysaev.sergey@gmail.com
