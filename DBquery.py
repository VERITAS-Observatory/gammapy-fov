import argparse
import time
from multiprocessing import Pool
import os
import json
import argparse
import pymysql
import datetime
import math
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import ascii
from astropy.table import Table

# https://veritas.sao.arizona.edu/wiki/VOFFLINE_Database_Tables
# https://veritas.sao.arizona.edu/wiki/Ryan_Dickherber%27s_Wiki
# Log Gen script: http://veritash.sao.arizona.edu:8081/OfflineAnalysis-WG/230319_221625/query_night
# Partially borrowed from Ruo's veritas_db_query.py

datadir = os.environ['VTS_GP']+'/data'
targetlist = f'{datadir}/VERITAS_targets.json'
# epoch date ranges stored in json can be found at https://veritas.sao.arizona.edu/wiki/Epoch_boundary_run_number
epochs = f'{datadir}/epochs.json'
with open(epochs, 'r') as fp: epochs = json.load(fp)

db_host = os.environ['DB_HOST']
db_user = os.environ['DB_USER']


def get_epoch(run_id):
    # setup database connection
    dbcnx=pymysql.connect(host=db_host, db='VERITAS', user=db_user, cursorclass=pymysql.cursors.DictCursor, charset="utf8")
    # connect to database
    crs=dbcnx.cursor()
    query = f'SELECT run_id,data_start_time,data_end_time FROM tblRun_Info WHERE run_id={run_id}'
    crs.execute(query)
    res = crs.fetchall()
    dbcnx.close()

    start_time = res[0]['data_start_time']
    date_format = '%Y-%m-%d %H:%M:%S'

    epoch = None
    for ep in epochs.keys():
        start, end = epochs[ep]
        after = datetime.datetime.strptime(end, date_format) > start_time
        before = start_time > datetime.datetime.strptime(start, date_format)
        if after and before:
            epoch = ep
        else: continue
    if epoch is None:
        print("Run with run_id", run_id, "is missing epoch info. Is epochs.json up to date?")
    return epoch


def ConvertRadRaDecToGalactic(ra, dec):
    # Convert RA, Dec in radian to Glat and Glon in deg
    
    ra *= 180./math.pi
    dec *= 180./math.pi
    my_sky = SkyCoord(ra*u.deg, dec*u.deg, frame='icrs')
    return my_sky.galactic.l.deg, my_sky.galactic.b.deg


def get_runs_in_epoch(ep):
    # setup database connection
    dbcnx=pymysql.connect(host=db_host, db='VERITAS', user=db_user, cursorclass=pymysql.cursors.DictCursor, charset="utf8")
    # connect to database
    crs=dbcnx.cursor()
    start_time, end_time = epochs[ep]
    start_time, end_time = [datetime.datetime.strptime(t, '%Y-%m-%d %H:%M:%S') for t in [start_time, end_time]]
    query = f'SELECT run_id,data_start_time,data_end_time FROM tblRun_Info WHERE data_start_time>\'{start_time}\' AND data_end_time<\'{end_time}\''
    crs.execute(query)
    res = crs.fetchall()
    dbcnx.close()
    runs = [run['run_id'] for run in res]
    return runs


def get_src_ra_dec():  
    dbcnx=pymysql.connect(host=db_host, db='VERITAS', user=db_user, cursorclass=pymysql.cursors.DictCursor, charset="utf8")
    crs=dbcnx.cursor()
    query = 'SELECT source_id,ra,decl FROM tblObserving_Sources'
    crs.execute(query)
    sources = crs.fetchall()
    dbcnx.close()
    for src in sources:
        glon, glat = ConvertRadRaDecToGalactic(src['ra'], src['decl'])
        src.update({'glon': glon, 'glat': glat})
    with open(targetlist, 'w+') as fp:
        json.dump(sources, fp)
            

def get_run_dur(run_id):
    # setup database connection
    dbcnx=pymysql.connect(host=db_host, db='VOFFLINE', user=db_user, cursorclass=pymysql.cursors.DictCursor, charset="utf8")
    # connect to database
    crs=dbcnx.cursor()
    query = f'SELECT run_id,usable_duration FROM tblRun_Analysis_Comments WHERE run_id={run_id}'
    crs.execute(query)
    res = crs.fetchall()
    dbcnx.close()
    if len(res) == 0 or res[0]['usable_duration'] is None:
        dbcnx=pymysql.connect(host=db_host, db='VERITAS', user=db_user, cursorclass=pymysql.cursors.DictCursor, charset="utf8")
        crs=dbcnx.cursor()
        query = f'SELECT run_id,data_start_time,data_end_time FROM tblRun_Info WHERE run_id={run_id}'
        crs.execute(query)
        res = crs.fetchall()[0]
        return res['data_end_time']-res['data_start_time']
    else: return res[0]['usable_duration']
     

def get_wobble(run_id):
    dbcnx=pymysql.connect(host=db_host, db='VERITAS', user=db_user, cursorclass=pymysql.cursors.DictCursor, charset = "utf8")
    crs=dbcnx.cursor()
    query = f'SELECT run_id,offset_distance FROM tblRun_Info WHERE run_id={run_id}'
    crs.execute(query)
    res = crs.fetchall()
    dbcnx.close()
    return res[0]['offset_distance']


def get_wobble_angle_and_distance(run_id):
    dbcnx=pymysql.connect(host=db_host, db='VERITAS', user=db_user, cursorclass=pymysql.cursors.DictCursor, charset="utf8")
    crs=dbcnx.cursor()
    query = f'SELECT run_id,offset_distance,offset_angle FROM tblRun_Info WHERE run_id={run_id}'
    crs.execute(query)
    res = crs.fetchall()
    dbcnx.close()
    angle = {0: 'N', 90: 'E', 180: 'S', 270: 'W'}
    return angle[res[0]['offset_angle']], res[0]['offset_distance']

    
def get_run_el_az_current(run_id):

    # setup database connection
    dbcnx=pymysql.connect(host=db_host, db='VERITAS', user=db_user, cursorclass=pymysql.cursors.DictCursor, charset="utf8")
    # connect to database
    crs=dbcnx.cursor()
    query = f'SELECT run_id,data_start_time,data_end_time FROM tblRun_Info WHERE run_id={run_id}'
    crs.execute(query)
    res = crs.fetchall()
    
    tstart = str(res[0]['data_start_time'])
    tend = str(res[0]['data_end_time'])
    if tstart=='None': return 0, 0
    if tend=='None': return 0, 0
    timestamp_start = tstart.replace(' ','').replace('-','').replace(':','')
    timestamp_end = tend.replace(' ','').replace('-','').replace(':','')
    timestamp_start += '000'
    timestamp_end += '000'
    el_avg = 0.
    az_avg = 0.
    total_entries = 0.
    
    query = f'SELECT elevation_target,azimuth_target FROM tblPositioner_Telescope1_Status WHERE timestamp>{timestamp_start} AND timestamp<{timestamp_end}'
    crs.execute(query)
    res = crs.fetchall()
    for x in res:
        el_avg += x['elevation_target']
        az_avg += x['azimuth_target']
        total_entries += 1.
    if total_entries==0.:
        return 0., 0.
    el_avg = el_avg/total_entries*180./math.pi
    az_avg = az_avg/total_entries*180./math.pi
    
    current_avg = 0.
    total_channels = 0.
    for ch in range(1,500):
        current_avg_ch = 0.
        total_entries = 0.
        query = f'SELECT current_meas FROM tblHV_Telescope1_Status WHERE db_start_time>\'{tstart}\' AND db_start_time<\'{tend}\' AND channel={ch}'
        crs.execute(query)
        res = crs.fetchall()
        for x in res:
            current_avg_ch += x['current_meas']
            total_entries += 1.
        current_avg_ch = current_avg_ch/total_entries
        current_avg += current_avg_ch
        total_channels += 1.
    current_avg = current_avg/total_channels
    dbcnx.close()
    
    return el_avg, az_avg, current_avg


def get_run_srcid_el_az_current(run_id):

    # setup database connection
    dbcnx=pymysql.connect(host=db_host, db='VERITAS', user=db_user, cursorclass=pymysql.cursors.DictCursor, charset="utf8")
    # connect to database
    crs=dbcnx.cursor()
    query = f'SELECT run_id,source_id,data_start_time,data_end_time FROM tblRun_Info WHERE run_id={run_id}'
    crs.execute(query)
    res = crs.fetchall()
    
    srcid = res[0]['source_id']
    tstart = str(res[0]['data_start_time'])
    tend = str(res[0]['data_end_time'])
    if tstart=='None': return 0, 0
    if tend=='None': return 0, 0
    timestamp_start = tstart.replace(' ','').replace('-','').replace(':','')
    timestamp_end = tend.replace(' ','').replace('-','').replace(':','')
    timestamp_start += '000'
    timestamp_end += '000'
    el_avg = 0.
    az_avg = 0.
    total_entries = 0.
    
    query = f'SELECT elevation_target,azimuth_target FROM tblPositioner_Telescope1_Status WHERE timestamp>{timestamp_start} AND timestamp<{timestamp_end}'
    crs.execute(query)
    res = crs.fetchall()
    for x in res:
        el_avg += x['elevation_target']
        az_avg += x['azimuth_target']
        total_entries += 1.
    if total_entries==0.:
        return 0., 0.
    el_avg = el_avg/total_entries*180./math.pi
    az_avg = az_avg/total_entries*180./math.pi
    
    current_avg = 0.
    total_channels = 0.
    for ch in range(1,500):
        current_avg_ch = 0.
        total_entries = 0.
        query = f'SELECT current_meas FROM tblHV_Telescope1_Status WHERE db_start_time>\'{tstart}\' AND db_start_time<\'{tend}\' AND channel={ch}'
        crs.execute(query)
        res = crs.fetchall()
        for x in res:
            current_avg_ch += x['current_meas']
            total_entries += 1.
        current_avg_ch = current_avg_ch/total_entries
        current_avg += current_avg_ch
        total_channels += 1.
    current_avg = current_avg/total_channels
    dbcnx.close()
    
    return srcid, el_avg, az_avg, current_avg


def get_run_timecut(runs):

    # setup database connection
    dbcnx=pymysql.connect(host=db_host, db='VOFFLINE', user=db_user, cursorclass=pymysql.cursors.DictCursor, charset="utf8")
    # connect to database
    crs=dbcnx.cursor()
    time_cut_mask = []
    for run_id in runs:
        query = 'SELECT run_id,time_cut_mask FROM tblRun_Analysis_Comments WHERE run_id=%s'%(run_id)
        crs.execute(query)
        # fetch from cursor
        res = crs.fetchall()
        timecuts = res[0]['time_cut_mask']
        if timecuts is not None:
            print(f'Run {run_id} time mask {timecuts}')
            time_cut_mask.append([run_id, timecuts])
    dbcnx.close()
    return time_cut_mask


def write_anasum_timemask(timemaskfile, runs):
    
    time_cut_mask = get_run_timecut(runs)
    if len(time_cut_mask) > 0:
        with open(timemaskfile, 'a') as f:
            for mask in time_cut_mask:
                run, cuts = mask; cuts = cuts.split(',');
                for cut in cuts:
                    cutstart, cutend = cut.split('/')
                    cutdur = int(cutend) - int(cutstart)
                    f.write('\n* %s %s %s 0'%(run, cutstart, cutdur))
        print('%s time cuts written in %s'%(len(time_cut_mask), timemaskfile))
    else: print('No time cuts to apply')