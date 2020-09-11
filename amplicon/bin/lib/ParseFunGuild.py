# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 09:07:31 2017

@author: songz
"""

#%% funguild module in FAST
def parse_funguild_database(database_type):
    # Get the url for the database (fungi or nematode)
    if database_type == 'fungi':
    url = 'http://www.stbates.org/funguild_db.php'
    elif database_type == 'nematode':
    url = 'http://www.stbates.org/nemaguild_db.php'
    
    