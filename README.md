# BADERvsOS
Explore trend between Bader charges and oxidation states for binaries with AFLOW database and Python.
File 'Full_for_SV_and_MV.csv' includes final data set used in article ... 

quick_check.py can be used to analyses trends between OS and Bader charges. It is very first version and it needs a lot of work to be usable for wider community.

Needed libraries: pandas, requests, re, matplotlib, mendeleev

How to:
In quick_check specify 
1) path to save data 
2) desired anion and cation. 
In the end, graph will appear.

TODO:
1. quick_check analyses only single valence compounds
2. The script which analysis structures by calculating coordination numbers, cures from duplicates, make sure that each atom has anion-cation bond


If you have any question, please write to posysaev.sergey@gmail.com
