# BADERvsOS
Explore trend between Bader charges and oxidation states for binaries with AFLOW database and Python.

File 'Full_for_SV_and_MV.csv' includes final data set used in the article Posysaev, S., Miroshnichenko, O., Alatalo, M., Le, D. & Rahman, T. S. Oxidation states of binary oxides from data analytics of the electronic structure. Comput. Mater. Sci. 161, 403–414 (2019) https://doi.org/10.1016/j.commatsci.2019.01.046. You are welcome to read or download the article by April 23 2019 free of charge here https://authors.elsevier.com/c/1Yf~k3In-ul718.

quick_results.py can be used to analyze trends between OS and Bader charges.

aflowlib.py contains all needed functions, make sure it imported properly.

Needed libraries: pandas, requests, re, matplotlib, mendeleev, asyncio, nest_asyncio
Note: asyncio library has been added. Asynchronous programming allows making multiple requests to AFLOWLIB server which allows receiving needed data much faster.

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
