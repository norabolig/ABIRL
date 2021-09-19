#!/usr/bin/env python

## [Includes]
import sys
import os
import time
import re
import json

import requests
from urllib.parse import quote as urlencode            
## [Includes]

## [Mast Query]
def mast_query(request):
    """Perform a MAST query.
    
        Parameters
        ----------
        request (dictionary): The MAST request json object
        
        Returns head,content where head is the response HTTP headers, and content is the returned data"""
    
    # Base API url
    request_url='https://mast.stsci.edu/api/v0/invoke'    
    
    # Grab Python Version 
    version = ".".join(map(str, sys.version_info[:3]))

    # Create Http Header Variables
    headers = {"Content-type": "application/x-www-form-urlencoded",
               "Accept": "text/plain",
               "User-agent":"python-requests/"+version}

    # Encoding the request as a json string
    req_string = json.dumps(request)
    req_string = urlencode(req_string)
    
    # Perform the HTTP request
    resp = requests.post(request_url, data="request="+req_string, headers=headers)
    
    # Pull out the headers and response content
    head = resp.headers
    content = resp.content.decode('utf-8')

    return head, content
## [Mast Query]


## [Download Request]
def download_request(payload, filename, download_type="file"):
    request_url='https://mast.stsci.edu/api/v0.1/Download/' + download_type
    resp = requests.post(request_url, data=payload)

    with open(filename,'wb') as FLE:
        FLE.write(resp.content)

    return filename
## [Download Request]    


## [Json to csv]
def mast_json2csv(json):    
    csv_str =  ",".join([x['name'] for x in json['fields']])
    csv_str += "\n"
    csv_str += ",".join([x['type'] for x in json['fields']])
    csv_str += "\n"

    col_names = [x['name'] for x in json['fields']]  
    for row in json['data']:
        csv_str += ",".join([str(row.get(col,"nul")) for col in col_names]) + "\n"
        
    return csv_str
## [Json to csv]


## [Json to astropy]
from astropy.table import Table
import numpy as np

def mast_json2table(json_obj):

    data_table = Table()

    for col,atype in [(x['name'],x['type']) for x in json_obj['fields']]:
            if atype=="string":
                atype="str"
            if atype=="boolean":
                atype="bool"
            
            try:data_table[col] = np.array([x.get(col,None) for x in json_obj['data']],dtype=atype)
            except:pass
        
    return data_table
## [Json to astropy]


## [Name Resolver]
def resolve_name():

    request = {'service':'Mast.Name.Lookup',
               'params':{'input':'M101',
                         'format':'json'},
    }

    headers, out_string = mast_query(request)

    out_data = json.loads(out_string)

    return out_data
## [Name Resolver]

## [List Missions]
def list_caom_missions():

    request = {
        'service':'Mast.Missions.List',
        'params':{},
        'format':'json'
    }
    
    headers, out_string = mast_query(request)

    out_data = [x['distinctValue'] for x in json.loads(out_string)['data']]

    return out_data
## [List Missions]

## [CAOM Cone Search]
def caom_cone_search():
    cachebreaker = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
    request = {'service':'Mast.Caom.Cone',
               'params':{'ra':254.28746,
                         'dec':-4.09933,
                         'radius': 0.2,
                         'cachebreaker': cachebreaker},
               'format':'json',
               'pagesize':5000,
               'page':1,
               'removenullcolumns':True,
               'timeout':3,
               'cachebreaker': cachebreaker}

    status = 'EXECUTING'
    while status == 'EXECUTING':
        headers, out_string = mast_query(request)
        out_data = json.loads(out_string)
        status = out_data['status']
        print(status)

    headers, out_string = mast_query(request)
    out_data = json.loads(out_string)

    return out_data
## [CAOM Cone Search]

## [VO Cone Search]
def vo_cone_search():

    request = {'service':'Vo.Hesarc.DatascopeListable',
               'params':{'ra':254.28746,
                         'dec':-4.09933,
                         'radius':0.2},
               'format':'json',
               'removenullcolumns':True}   

    all_data = []
    start_time = time.time()
    
    while True:
        headers, out_string = mast_query(request)
        out_data = json.loads(out_string)
        all_data.append(out_data)
        if out_data['status'] != "EXECUTING":
            break
        if time.time() - start_time > 30:
            print("Working...")
            start_time = time.time()
        time.sleep(10)
   
    return all_data
## [VO Cone Search]

## [HSC V2 Cone Search]
def hscV2_cone_search():

    request = {'service':'Mast.Hsc.Db.v2',
               'params':{'ra':254.287,
                         'dec':-4.09933,
                         'radius':0.2,
                         'nr':5000,
                         'ni':1,
                         'magtype':1},
               'format':'json',
               'pagesize':1000,
               'page':1,
               'removenullcolumns':True}   

    headers, out_string = mast_query(request)

    out_data = json.loads(out_string)

    return out_data
## [HSC V2 Cone Search]

## [HSC V3 Cone Search]
def hscV3_cone_search():

    request = {'service':'Mast.Hsc.Db.v3',
               'params':{'ra':254.287,
                         'dec':-4.09933,
                         'radius':0.2,
                         'nr':5000,
                         'ni':1,
                         'magtype':1},
               'format':'json',
               'pagesize':1000,
               'page':1,
               'removenullcolumns':True}   

    headers, out_string = mast_query(request)

    out_data = json.loads(out_string)

    return out_data
## [HSC V3 Cone Search]

### [GAIA DR2 Cone Search]
def gaiaDR2_cone_search(ra,dec,rad,pgsize=1000,page=1):
    request = {'service':'Mast.Catalogs.GaiaDR2.Cone',
               'params':{'ra':ra,
                         'dec':dec,
                         'radius':rad},
               'format':'json',
               'pagesize':pgsize,
               'page':page}   

    headers, out_string = mast_query(request)

    out_data = json.loads(out_string)

    return out_data
## [GAIA DR2 Cone Search]

### [GAIA DR3 Cone Search]
def gaiaDR3_cone_search(ra,dec,rad,pgsize=1000,page=1):
    request = {'service':'Mast.Catalogs.GaiaDR3.Cone',
               'params':{'ra':ra,
                         'dec':dec,
                         'radius':rad},
               'format':'json',
               'pagesize':pgsize,
               'page':page}   

    headers, out_string = mast_query(request)

    out_data = json.loads(out_string)

    return out_data
## [GAIA DR3 Cone Search]


## [TGAS Cone Search]
def tgas_cone_search():
    request = { "service":"Mast.Catalogs.Tgas.Cone",
                "params":{
                    "ra":254.28746,
                    "dec":-4.09933,
                    "radius":0.2},
                "format":"json",
                "timeout":10}   

    headers, out_string = mast_query(request)

    out_data = json.loads(out_string)

    return out_data
## [TGAS Cone Search]

## [GALEX Cone Search]
def galex_cone_search():
    request = { "service":"Mast.Galex.Catalog",
                "params":{
                    "ra":254.28746,
                    "dec":-4.09933,
                    "radius":0.2},
                "format":"json",
                "timeout":10}   

    headers, out_string = mast_query(request)

    out_data = json.loads(out_string)

    return out_data
## [GALEX Cone Search]

## [TIC Cone Search]
def tic_cone_search():
    request = { "service":"Mast.Catalogs.Tic.Cone",
                "params":{
                    "ra":254.28746,
                    "dec":-4.09933,
                    "radius":0.2},
                "format":"json",
                "timeout":10}   

    headers, out_string = mast_query(request)

    out_data = json.loads(out_string)

    return out_data
## [TIC Cone Search]

## [TIC Advanced Search]
def tic_advanced_search():
    request = {"service":"Mast.Catalogs.Filtered.Tic",
               "format":"json",
               "params":{
                   "columns":"c.*",
                   "filters":[
                       {"paramName":"dec",
                        "values":[{"min":-90.,"max":-30.}]},
                       {"paramName":"Teff",
                        "values":[{"min":4250.,"max":4500.}]},
                       {"paramName":"logg",
                        "values":[{"min":4.4,"max":5.0}]},
                       {"paramName":"Tmag",
                        "values":[{"min":8.,"max":10.}]}]
               }}

    headers, out_string = mast_query(request)
    out_data = json.loads(out_string)

    return out_data
## [TIC Advanced Search]

## [TIC Advanced Search (Data Only)]
def tic_advanced_search_rows():
    request = {"service":"Mast.Catalogs.Filtered.Tic.Rows",
               "format":"json",
               "params":{
                   "filters":[
                       {"paramName":"dec",
                        "values":[{"min":-90.,"max":-30.}]},
                       {"paramName":"Teff",
                        "values":[{"min":4250.,"max":4500.}]},
                       {"paramName":"logg",
                        "values":[{"min":4.4,"max":5.0}]},
                       {"paramName":"Tmag",
                        "values":[{"min":8.,"max":10.}]}]
               }}

    headers, out_string = mast_query(request)
    out_data = json.loads(out_string)

    return out_data
## [TIC Advanced Search (Data Only)]

## [TIC Advanced Search Position]
def tic_advanced_search_position():
    request = {"service":"Mast.Catalogs.Filtered.Tic.Position",
               "format":"json",
               "params":{
                   "columns":"c.*",
                   "filters":[
                       {"paramName":"Teff",
                        "values":[{"min":4250.,"max":4500.}]},
                       {"paramName":"logg",
                        "values":[{"min":4.4,"max":5.0}]},
                       {"paramName":"Tmag",
                        "values":[{"min":8.,"max":10.}]}],
                   "ra": 210.8023,
                   "dec": 54.349,
                   "radius": .2
               }}

    headers, out_string = mast_query(request)
    out_data = json.loads(out_string)

    return out_data
## [TIC Advanced Search Position]

## [TIC Advanced Search Position (Data Only)]
def tic_advanced_search_position_rows():
    request = {"service":"Mast.Catalogs.Filtered.Tic.Position.Rows",
               "format":"json",
               "params":{
                   "filters":[
                       {"paramName":"Teff",
                        "values":[{"min":4250.,"max":4500.}]},
                       {"paramName":"logg",
                        "values":[{"min":4.4,"max":5.0}]},
                       {"paramName":"Tmag",
                        "values":[{"min":8.,"max":10.}]}],
                   "ra": 210.8023,
                   "dec": 54.349,
                   "radius": .2
               }}

    headers, out_string = mast_query(request)
    out_data = json.loads(out_string)

    return out_data
## [TIC Advanced Search Position (Data Only)]

## [CTL Cone Search]
def ctl_cone_search():
    request = { "service":"Mast.Catalogs.Ctl.Cone",
                "params":{
                    "ra":254.28746,
                    "dec":-4.09933,
                    "radius":0.2},
                "format":"json",
                "timeout":10}

    headers, out_string = mast_query(request)

    out_data = json.loads(out_string)

    return out_data
## [CTL Cone Search]

## [CTL Advanced Search]
def ctl_advanced_search():
    request = {"service":"Mast.Catalogs.Filtered.Ctl",
               "format":"json",
               "params":{
                   "columns":"c.*",
                   "filters":[
                       {"paramName":"dec",
                        "values":[{"min":-90.,"max":-30.}]},
                       {"paramName":"Teff",
                        "values":[{"min":4250.,"max":4500.}]},
                       {"paramName":"logg",
                        "values":[{"min":4.4,"max":5.0}]},
                       {"paramName":"Tmag",
                        "values":[{"min":8.,"max":10.}]}]
               }}

    headers, out_string = mast_query(request)
    out_data = json.loads(out_string)

    return out_data
## [CTL Advanced Search]

## [CTL Advanced Search (Data Only)]
def ctl_advanced_search_rows():
    request = {"service":"Mast.Catalogs.Filtered.Ctl.Rows",
               "format":"json",
               "params":{
                   "filters":[
                       {"paramName":"dec",
                        "values":[{"min":-90.,"max":-30.}]},
                       {"paramName":"Teff",
                        "values":[{"min":4250.,"max":4500.}]},
                       {"paramName":"logg",
                        "values":[{"min":4.4,"max":5.0}]},
                       {"paramName":"Tmag",
                        "values":[{"min":8.,"max":10.}]}]
               }}

    headers, out_string = mast_query(request)
    out_data = json.loads(out_string)

    return out_data
## [CTL Advanced Search (Data Only)]

## [CTL Advanced Search Position]
def ctl_advanced_search_position():
    request = {"service":"Mast.Catalogs.Filtered.Ctl.Position",
               "format":"json",
               "params":{
                   "columns":"c.*",
                   "filters":[
                       {"paramName":"Teff",
                        "values":[{"min":4250.,"max":4500.}]},
                       {"paramName":"logg",
                        "values":[{"min":4.4,"max":5.0}]},
                       {"paramName":"Tmag",
                        "values":[{"min":8.,"max":10.}]}],
                   "ra": 210.8023,
                   "dec": 54.349,
                   "radius": .2
               }}

    headers, out_string = mast_query(request)
    out_data = json.loads(out_string)

    return out_data
## [CTL Advanced Search Position]

## [CTL Advanced Search Position (Data Only)]
def ctl_advanced_search_position_rows():
    request = {"service":"Mast.Catalogs.Filtered.Ctl.Position.Rows",
               "format":"json",
               "params":{
                   "filters":[
                       {"paramName":"Teff",
                        "values":[{"min":4250.,"max":4500.}]},
                       {"paramName":"logg",
                        "values":[{"min":4.4,"max":5.0}]},
                       {"paramName":"Tmag",
                        "values":[{"min":8.,"max":10.}]}],
                   "ra": 210.8023,
                   "dec": 54.349,
                   "radius": .2
               }}

    headers, out_string = mast_query(request)
    out_data = json.loads(out_string)

    return out_data
## [CTL Advanced Search Position (Data Only)]

## [DD Cone Search]
def dd_cone_search():
    request = { "service":"Mast.Catalogs.DiskDetective.Cone",
                "params":{
                    "ra":254.28746,
                    "dec":-4.09933,
                    "radius":0.2},
                "format":"json",
                "timeout":10}   

    headers, out_string = mast_query(request)

    out_data = json.loads(out_string)

    return out_data
## [DD Cone Search]

## [DD Advanced Search]
def dd_advanced_search():
    request = {"service":"Mast.Catalogs.Filtered.DiskDetective",
               "format":"json",
               "params":{
                   "filters":[
                       {"paramName":"classifiers",
                        "values":[{"min":10,"max":18}]},
                       {"paramName":"oval",
                        "values":[{"min":15,"max":76}]}]
               }}

    headers, out_string = mast_query(request)
    out_data = json.loads(out_string)

    return out_data
## [DD Advanced Search]

## [DD Advanced Search Position]
def dd_advanced_search_position():
    request = {"service":"Mast.Catalogs.Filtered.DiskDetective.Position",
               "format":"json",
               "params":{
                   "filters":[
                       {"paramName":"classifiers",
                        "values":[{"min":10,"max":18}]}],
                   "ra": 86.6909,
                   "dec": 0.079,
                   "radius": .2
               }}

    headers, out_string = mast_query(request)
    out_data = json.loads(out_string)

    return out_data
## [DD Advanced Search Position]

## [DD Advanced Search Counts]
def dd_advanced_search_counts():
    request = {"service":"Mast.Catalogs.Filtered.DiskDetective.Count",
               "format":"json",
               "params":{
                   "filters":[
                       {"paramName":"classifiers",
                        "values":[{"min":10,"max":18}]},
                       {"paramName":"oval",
                        "values":[{"min":15,"max":76}]}]
               }}

    headers, out_string = mast_query(request)
    out_data = json.loads(out_string)

    return out_data
## [DD Advanced Search Counts]

## [DD Advanced Search Position Counts]
def dd_advanced_search_position_counts():
    request = {"service":"Mast.Catalogs.Filtered.DiskDetective.Position.Count",
               "format":"json",
               "params":{
                   "filters":[
                       {"paramName":"classifiers",
                        "values":[{"min":10,"max":18}]}],
                   "ra": 86.6909,
                   "dec": 0.079,
                   "radius": .2
               }}

    headers, out_string = mast_query(request)
    out_data = json.loads(out_string)

    return out_data
## [DD Advanced Search Position Counts]

## [WFC3 PSF UVIS Filtered]
def get_wfc3uvis_matches():
    request = {"service":"Mast.Catalogs.Filtered.Wfc3.Uvis",
               "format":"json",
               "params":{
                   "filters":[
                       {"paramName":"filter",
                        "values": ["F814W", "F390W"]},
                       {"paramName": "psf_x_center",
                        "values":[{"min":6,"max":18}]}]
               }}

    headers, out_string = mast_query(request)
    out_data = json.loads(out_string)

    return out_data
## [WFC3 PSF UVIS Filtered]

## [WFC3 PSF IR Filtered]
def get_wfc3ir_matches():
    request = {"service":"Mast.Catalogs.Filtered.Wfc3.Ir",
               "format":"json",
               "params":{
                   "filters":[
                       {"paramName":"filter",
                        "values": ["F814W", "F390W"]},
                       {"paramName": "psf_x_center",
                        "values":[{"min":6,"max":18}]}]
               }}

    headers, out_string = mast_query(request)
    out_data = json.loads(out_string)

    return out_data
## [WFC3 PSF IR Filtered]

## [Advanced Search]
def advanced_search_counts():
    request = {"service":"Mast.Caom.Filtered",
               "format":"json",
               "params":{
                   "columns":"COUNT_BIG(*)",
                   "filters":[
                       {"paramName":"filters",
                        "values":["NUV","FUV"],
                        "separator":";"
                       },
                       {"paramName":"t_max",
                        "values":[{"min":52264.4586,"max":54452.8914}], #MJD
                       },
                       {"paramName":"obsid",
                        "values":[],
                        "freeText":"%200%"}
                   ]}}

    headers, out_string = mast_query(request)
    out_data = json.loads(out_string)

    return out_data


def advanced_search():
    request = {"service":"Mast.Caom.Filtered",
               "format":"json",
               "params":{
                   "columns":"*",
                   "filters":[
                       {"paramName":"dataproduct_type",
                        "values":["image"],
                       },
                       {"paramName":"proposal_pi",
                        "values":["Osten"]
                       }
                   ],
                   "obstype":"all"
               }}

    headers, out_string = mast_query(request)
    out_data = json.loads(out_string)

    return out_data
## [Advanced Search]


## [Advanced Search Position]
def advanced_search_with_position_counts():
    request = { "service":"Mast.Caom.Filtered.Position",
                "format":"json",
                "params":{
                    "columns":"COUNT_BIG(*)",
                    "filters":[
                        {"paramName":"dataproduct_type",
                         "values":["cube","image"]
                        }],
                    "position":"210.8023, 54.349, 5"
                }}

    headers, out_string = mast_query(request)
    out_data = json.loads(out_string)

    return out_data


def advanced_search_with_position():
    request = {"service":"Mast.Caom.Filtered.Position",
               "format":"json",
               "params":{
                   "columns":"*",
                   "filters":[
                       {"paramName":"dataproduct_type",
                        "values":["cube"]
                       }],
                   "position":"210.8023, 54.349, 0.24"
               }}

    headers, out_string = mast_query(request)
    out_data = json.loads(out_string)

    return out_data
## [Advanced Search Position]


## [HSC Spectra]
def hsc_spectra_search():

    request = {'service':'Mast.HscSpectra.Db.All',
               'format':'votable'}   

    headers, out_string = mast_query(request)

    return out_string
## [HSC Spectra]

## [HSC Spectra Download]
def download_hsc_spectra():

    # grab all the hsc spectra
    request = {'service':'Mast.HscSpectra.Db.All',
               'format':'json'}   
    headers, out_string = mast_query(request)
    out_data = json.loads(out_string)

    # download the first 3 spects
    for spec in out_data['data'][:3]:
        # build the url
        if spec['SpectrumType'] < 2:
            dataurl = 'https://hla.stsci.edu/cgi-bin/getdata.cgi?config=ops&dataset=' \
                      + spec['DatasetName']
            filename = spec['DatasetName']
        else:
            dataurl = 'https://hla.stsci.edu/cgi-bin/ecfproxy?file_id=' \
                      + spec['DatasetName'] + '.fits'
            filename = spec['DatasetName'] + '.fits'

        # download
        out_file = download_request(dataurl, filename)

        print(out_file)   
## [HSC Spectra Download]


## [HSC V2 Matches]
def get_hscv2_matches():
    # perform the HSC search
    result = hscV2_cone_search()
    data = result['data']

    # get the match id
    matchId = data[0]['MatchID']

    # get detailed results for chosen match id
    request = {'service':'Mast.HscMatches.Db.v2',
               'params':{'input':matchId},
               'format':'json',
               'page':1,
               'pagesize':4}   

    headers, out_string = mast_query(request)

    out_data = json.loads(out_string)

    return out_data
## [HSC V2 Matches]

## [HSC V3 Matches]
def get_hscv3_matches():
    # perform the HSC search
    result = hscV3_cone_search()
    data = result['data']

    # get the match id
    matchId = data[0]['MatchID']

    # get detailed results for chosen match id
    request = {'service':'Mast.HscMatches.Db.v3',
               'params':{'input':matchId},
               'format':'json',
               'page':1,
               'pagesize':4}   

    headers, out_string = mast_query(request)

    out_data = json.loads(out_string)

    return out_data
## [HSC V3 Matches]


## [Get VO Data]
def get_vo_data():

    # perform vo cone search
    vo_data = vo_cone_search()
    vo_json = vo_data[0]
    
    row = json['data'][2]  
    url = row['tableURL']
        
    request = {'service':'Vo.Generic.Table',
               'params':{'url':url},
               'format':'json',
               'page':1,
               'pagesize':1000}   

    headers, out_string = mast_query(request)

    out_data = json.loads(out_string)

    return out_data
## [Get VO Data]


## [CAOM Crossmatch]
def crossmatch_from_cone_search():
    
    # This is a json object
    crossmatch_input = caom_cone_search()
    
    request =  {"service":"Mast.Caom.Crossmatch",
                "data":crossmatch_input,
                "params":{
                    "raColumn":"s_ra",
                    "decColumn":"s_dec",
                    "radius":0.001
                },
                "pagesize":1000,
                "page":1,
                "format":"json",
                "removecache":True}

    headers, out_string = mast_query(request)

    out_data = json.loads(out_string)

    return out_data

def crossmatch_from_minimal_json():
    
    crossmatch_input = {"fields":[{"name":"ra","type":"float"},
                                  {"name":"dec","type":"float"}],
                        "data":[{"ra":210.8,"dec":54.3}]}
    
    request =  {"service":"Mast.Caom.Crossmatch",
                "data":crossmatch_input,
                "params":{
                    "raColumn":"ra",
                    "decColumn":"dec",
                    "radius":0.001
                },
                "pagesize":1000,
                "page":1,
                "format":"json",
                "clearcache":True}

    headers, out_string = mast_query(request)

    out_data = json.loads(out_string)

    return out_data
## [CAOM Crossmatch]

## [Galex Crossmatch]
def galex_crossmatch():

    # This is a json object
    crossmatch_input = {"fields":[{"name":"s_ra","type":"float"},
                                  {"name":"s_dec","type":"float"}],
                        "data":[{"s_ra":210.8,"s_dec":54.3}]}
    
    request =  {"service":"Mast.Galex.Crossmatch",
                "data":crossmatch_input,
                "params":{
                    "raColumn":"s_ra",
                    "decColumn":"s_dec",
                    "radius":0.01
                },
                "pagesize":1000,
                "page":1,
                "format":"json"}

    headers, out_string = mast_query(request)

    out_data = json.loads(out_string)

    return out_data
## [Galex Crossmatch]    

## [sdss Crossmatch]
def sdss_crossmatch():

    # This is a json object
    crossmatch_input = {"fields":[{"name":"ra","type":"float"},
                                  {"name":"dec","type":"float"}],
                        "data":[{"ra":337.10977,"dec":30.30261}]} 

    request ={"service":"Mast.Sdss.Crossmatch",
              "data":crossmatch_input,
              "params": {
                  "raColumn":"ra",
                  "decColumn":"dec",
                  "radius":0.01
              },
              "format":"json",
              "pagesize":1000,
              "page":1}

    headers, out_string = mast_query(request)

    out_data = json.loads(out_string)

    return out_data
## [sdss Crossmatch]


## [2mass Crossmatch]
def twoMass_crossmatch():

    # This is a json object
    crossmatch_input = {"fields":[{"name":"ra","type":"float"},
                                  {"name":"dec","type":"float"}],
                        "data":[{"ra":210.88447,"dec":54.332}]}
    
    request =  {"service":"Mast.2Mass.Crossmatch",
                "data":crossmatch_input,
                "params":{
                    "raColumn":"ra",
                    "decColumn":"dec",
                    "radius":0.04
                },
                "pagesize":1000,
                "page":1,
                "format":"json"}

    headers, out_string = mast_query(request)

    out_data = json.loads(out_string)

    return out_data
## [2mass Crossmatch]


## [hsc2 Crossmatch]
def hscMagAper2_crossmatch():

    # This is a json object
    crossmatch_input = {"fields":[{"name":"ra","type":"float"},
                                  {"name":"dec","type":"float"}],
                        "data":[{"ra":210.8,"dec":54.3}]}
    
    request =  {"service":"Mast.Hsc.Crossmatch.MagAper2v3",
                "data":crossmatch_input,
                "params":{
                    "raColumn":"ra",
                    "decColumn":"dec",
                    "radius":0.001
                },
                "pagesize":1000,
                "page":1,
                "format":"json"}

    headers, out_string = mast_query(request)

    out_data = json.loads(out_string)

    return out_data
## [hsc2 Crossmatch]

## [hscauto Crossmatch]
def hscMagAuto_crossmatch():

    # This is a json object
    crossmatch_input = {"fields":[{"name":"ra","type":"float"},
                                  {"name":"dec","type":"float"}],
                        "data":[{"ra":210.8,"dec":54.3}]}
    
    request =  {"service":"Mast.Hsc.Crossmatch.MagAutov3",
                "data":crossmatch_input,
                "params":{
                    "raColumn":"ra",
                    "decColumn":"dec",
                    "radius":0.001
                },
                "pagesize":1000,
                "page":1,
                "format":"json"}

    headers, out_string = mast_query(request)

    out_data = json.loads(out_string)

    return out_data
## [hscauto Crossmatch]

## [gaia DR1 Crossmatch]
def gaiaDR1_crossmatch():

    # This is a json object
    crossmatch_input = {"fields":[{"name":"ra","type":"float"},
                                  {"name":"dec","type":"float"}],
                        "data":[{"ra":210.8,"dec":54.3}]}
    
    request =  {"service":"Mast.GaiaDR1.Crossmatch",
                "data":crossmatch_input,
                "params":{
                    "raColumn":"ra",
                    "decColumn":"dec",
                    "radius":0.1
                },
                "format":"json"}

    headers, out_string = mast_query(request)

    out_data = json.loads(out_string)

    return out_data
## [gaia DR1 Crossmatch]

## [gaia DR2 Crossmatch]
def gaiaDR2_crossmatch():

    # This is a json object
    crossmatch_input = {"fields":[{"name":"ra","type":"float"},
                                 {"name":"dec","type":"float"}],
                       "data":[{"ra":210.8,"dec":54.3}]}
    
    request =  {"service":"Mast.GaiaDR2.Crossmatch",
                "data":crossmatch_input,
                "params":{
                    "raColumn":"ra",
                    "decColumn":"dec",
                    "radius":0.1
                },
                "format":"json"}

    headers, out_string = mast_query(request)

    out_data = json.loads(out_string)

    return out_data
## [gaia DR2 Crossmatch]

## [tgas Crossmatch]
def tgas_crossmatch():

    # This is a json object
    crossmatch_input = {"fields":[{"name":"ra","type":"float"},
                                 {"name":"dec","type":"float"}],
                       "data":[{"ra":211.09,"dec":54.3228}]}
    
    request =  {"service":"Mast.Tgas.Crossmatch",
                "data":crossmatch_input,
                "params":{
                    "raColumn":"ra",
                    "decColumn":"dec",
                    "radius":0.2
                },
                "format":"json"}

    headers, out_string = mast_query(request)

    out_data = json.loads(out_string)

    return out_data
## [tgas Crossmatch]

## [tic Crossmatch]
def tic_crossmatch():

    # This is a json object
    crossmatch_input = {"fields":[{"name":"ra","type":"float"},
                                 {"name":"dec","type":"float"}],
                       "data":[{"ra":211.09,"dec":54.3228}]}
    
    request =  {"service":"Mast.Tic.Crossmatch",
                "data":crossmatch_input,
                "params":{
                    "raColumn":"ra",
                    "decColumn":"dec",
                    "radius":0.2
                },
                "format":"json"}

    headers,out_string = mast_query(request)

    out_data = json.loads(out_string)

    return out_data
## [tic Crossmatch]

## [ctl Crossmatch]
def ctl_crossmatch():

    # This is a json object
    crossmatch_input = {"fields":[{"name":"ra","type":"float"},
                                 {"name":"dec","type":"float"}],
                       "data":[{"ra":211.09,"dec":54.3228}]}

    request =  {"service":"Mast.Ctl.Crossmatch",
                "data":crossmatch_input,
                "params":{
                    "raColumn":"ra",
                    "decColumn":"dec",
                    "radius":0.2
                },
                "format":"json"}

    headers, out_string = mast_query(request)

    out_data = json.loads(out_string)

    return out_data
## [ctl Crossmatch]

## [Product Query]
def get_caom_products():

    # perform the CAOM search
    result = caom_cone_search()
    data = result['data']

    # get the product group id (obsid)
    obsid = data[1]['obsid']

    # get detailed results for chosen match id
    request = {'service':'Mast.Caom.Products',
               'params':{'obsid':obsid},
               'format':'json',
               'pagesize':4,
               'page':1}   

    headers, out_string = mast_query(request)

    out_data = json.loads(out_string)

    return out_data
## [Product Query]


## [Download Product]
def download_multiple_products():

    # get data products
    result = get_caom_products()
    product_list = result['data']

    # Collecting the parameters you need
    url_list = [("uri", url) for url in product_list['dataURI'][:2]]
    extension = ".tar.gz"
    filename = "mastData"

    # download the file
    bundle = download_request(url_list, filename=local_path, download_type="bundle")

    return bundle
## [Download Product]


## [Direct Download]
def download_individual_products():

    # get data products
    result = get_caom_products()
    data = result['data']
    product = data[0]

    # Make the local file name and make sure the directory is there
    local_path = os.path.join("mastFiles", product['obs_collection'], product['obs_id'])
    if not os.path.exists(local_path):
        os.makedirs(local_path)
    print(local_path)
    local_path = os.path.join(local_path,  os.path.basename(product['productFilename']))
    print(local_path)
    
    # download the file
    download_filename = download_request(product['dataURI'], filename=local_path)
    
    return download_filename
## [Direct Download]

