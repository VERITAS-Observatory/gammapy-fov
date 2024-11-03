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

condition={'elevation': 5, 'azimuth': 10, 'current': 0.2} # el, az, current
# mindur = 3 # times on run duration

glat_threshold = 10
dur_threshold = datetime.timedelta(seconds=600)

targetlist = '/a/data/tehanu/jw3855/ed_analysis/VERITAS_target_list_coord.dat'
good_weather = ['A+', 'A', 'A-', 'B+', 'B', 'B-']
epochs = [["oa", "2004-01-01 00:00:00", "2009-07-31 23:59:59"],
          ["na", "2009-08-01 00:00:00", "2012-08-31 23:59:59"],
          ["12s", "2012-09-01 00:00:00", "2012-10-28 23:59:59"],
          ["1213w", "2012-10-29 00:00:00", "2013-05-22 23:59:59"],
          ["13s", "2013-05-23 00:00:00", "2013-11-16 23:59:59"],
          ["1314w", "2013-11-17 00:00:00", "2014-05-12 23:59:59"],
          ["14s", "2014-05-13 00:00:00", "2014-11-07 23:59:59"],
          ["1415w", "2014-11-08 00:00:00", "2015-05-31 23:59:59"],
          ["15s", "2015-06-01 00:00:00", "2015-11-25 23:59:59"],
          ["1516w", "2015-11-26 00:00:00", "2016-05-20 23:59:59"],
          ["16s", "2016-05-21 00:00:00", "2016-11-14 23:59:59"],
          ["1617w", "2016-11-15 00:00:00", "2017-05-10 23:59:59"],
          ["17s", "2017-05-11 00:00:00", "2017-12-05 23:59:59"],
          ["1718w", "2017-12-06 00:00:00", "2018-04-30 23:59:59"],
          ["18s", "2018-05-01 00:00:00", "2018-10-24 23:59:59"],
          ["1819w", "2018-10-25 00:00:00", "2019-05-18 23:59:59"],
          ["19s", "2019-05-19 00:00:00", "2019-11-12 23:59:59"],
          ["1920w", "2019-11-13 00:00:00", "2020-05-07 23:59:59"],
          ["20s", "2020-05-08 00:00:00", "2020-11-04 23:59:59"],
          ["2021w", "2020-11-05 00:00:00", "2021-04-27 23:59:59"],
          ["21s", "2021-04-28 00:00:00", "2021-11-16 23:59:59"],
          ["2122w", "2021-11-17 00:00:00", "2022-04-15 23:59:59"],
          ["22s", "2022-04-16 00:00:00", "2022-11-08 23:59:59"],
          ["2223w", "2022-11-09 00:00:00", "2023-05-01 23:59:59"]]

timemaskfile = '/a/data/tehanu/jw3855/software/ED490.1/Eventdisplay_AnalysisFiles_VTS/ParameterFiles/ANASUM.timemask.dat'


def get_epoch(run_id):
    # setup database connection
    dbcnx=pymysql.connect(host='romulus.ucsc.edu', db='VERITAS', user='readonly', cursorclass=pymysql.cursors.DictCursor)
    # connect to database
    crs=dbcnx.cursor()
    query = f'SELECT run_id,data_start_time,data_end_time FROM tblRun_Info WHERE run_id={run_id}'
    crs.execute(query)
    res = crs.fetchall()

    start_time = res[0]['data_start_time']
    date_format = '%Y-%m-%d %H:%M:%S'
    for ep, start, end in epochs:
        after = datetime.datetime.strptime(end, date_format) > start_time
        before = start_time > datetime.datetime.strptime(start, date_format)
        if after and before:
            epoch = ep
        else: continue
    return epoch


def ConvertRadRaDecToGalactic(ra, dec):
    # Convert RA, Dec in radian to Glat and Glon in deg
    
    ra *= 180./math.pi
    dec *= 180./math.pi
    my_sky = SkyCoord(ra*u.deg, dec*u.deg, frame='icrs')
    return my_sky.galactic.l.deg, my_sky.galactic.b.deg


def get_runs_in_epoch(ep):
    # setup database connection
    dbcnx=pymysql.connect(host='romulus.ucsc.edu', db='VERITAS', user='readonly', cursorclass=pymysql.cursors.DictCursor)
    # connect to database
    crs=dbcnx.cursor()
    start_time, end_time = [epoch[1:] for epoch in epochs if epoch[0] == ep][0]
    start_time, end_time = [datetime.datetime.strptime(t, '%Y-%m-%d %H:%M:%S') for t in [start_time, end_time]]
    query = f'SELECT run_id,data_start_time,data_end_time FROM tblRun_Info WHERE data_start_time>\'{start_time}\' AND data_end_time<\'{end_time}\''
    crs.execute(query)
    res = crs.fetchall()
    runs = [run['run_id'] for run in res]
    return runs


def get_src_ra_dec():  
    # setup database connection
    dbcnx=pymysql.connect(host='romulus.ucsc.edu', db='VERITAS', user='readonly', cursorclass=pymysql.cursors.DictCursor)
    # connect to database
    crs=dbcnx.cursor()
    query = 'SELECT source_id,ra,decl FROM tblObserving_Sources'
    crs.execute(query)
    sources = crs.fetchall()
    for src in sources:
        glon, glat = ConvertRadRaDecToGalactic(src['ra'], src['decl'])
        src.update({'glon': glon, 'glat': glat})
    Table(sources).write(targetlist, format='ascii', overwrite=True)
    

def get_epoch_good_egruns(ep):
    
    # science, good run, weather => B, usable_durateion => 10 min, abs(glat) => glat_threshold
    
    start_time, end_time = [epoch[1:] for epoch in epochs if epoch[0] == ep][0]
    start_time, end_time = [datetime.datetime.strptime(t, '%Y-%m-%d %H:%M:%S') for t in [start_time, end_time]]
    
    # setup database connection
    dbcnx=pymysql.connect(host='romulus.ucsc.edu', db='VERITAS', user='readonly', cursorclass=pymysql.cursors.DictCursor)
    # connect to database
    crs=dbcnx.cursor()
    query = f'SELECT run_id,data_start_time,data_end_time,weather,config_mask,source_id FROM tblRun_Info WHERE data_start_time>\'{start_time}\' AND data_end_time<\'{end_time}\' AND config_mask=15'
    crs.execute(query)
    good_runs = crs.fetchall()
    good_runs = [run for run in good_runs if run['weather'] in good_weather]
    
    eg_targets = ascii.read(targetlist)
    eg_targets = list(eg_targets[abs(eg_targets['glat'])>glat_threshold]['source_id'])
    good_eg_runs = []
    # setup database connection
    dbcnx=pymysql.connect(host='romulus.ucsc.edu', db='VOFFLINE', user='readonly', cursorclass=pymysql.cursors.DictCursor)
    # connect to database
    crs=dbcnx.cursor()
    
    print(f'Usable extragalactic (glat > +- {glat_threshold} deg) runs in epoch {ep}')
    for run in good_runs:
        run_id = run['run_id']
        query = f'SELECT run_id,data_category,status,usable_duration FROM tblRun_Analysis_Comments WHERE run_id={run_id}'
        crs.execute(query)
        res = crs.fetchall()
        if len(res) == 0: continue;
        if res[0]['data_category'] != 'science': continue;
        if res[0]['status'] != 'good_run': continue;
        if res[0]['usable_duration'] is None: continue;
        if res[0]['usable_duration'] < dur_threshold: continue;
        if run['source_id'] not in eg_targets: continue;
        else:
            good_eg_runs.append(run_id)
            print(run_id)
            
    return good_eg_runs


def multi_get_run_el_az_current(egrun):
    try:
        tinit = datetime.datetime.now()
        off_el, off_az, off_current = get_run_el_az_current(egrun)
        egrun_dict = {'run_id': egrun, 'elevation': off_el, 'azimuth': off_az, 'current': off_current}
        tfin = datetime.datetime.now()
        print(f'{tfin-tinit} {egrun_dict}')
        return egrun_dict
    except ZeroDivisionError as e:
        print(e)
        return


def multi_get_offruns(on_off_dict):
    onrun, offrun = on_off_dict
    if 'elevation' in condition.keys():
        if abs(onrun['elevation'] - offrun['elevation']) > condition['elevation']: return
    if 'azimuth' in condition.keys():
        az = onrun['azimuth']
        az_min = onrun['azimuth'] - condition['azimuth']
        az_max = onrun['azimuth'] + condition['azimuth']
        if az_min < 0: az_min += 360
        if az_max >= 360: az_max -= 360
        if az_min > az_max:
            if offrun['azimuth'] < az_min and offrun['azimuth'] > az_max: return
        else:
            if offrun['azimuth'] < az_min or offrun['azimuth'] > az_max: return
    if 'current' in condition.keys():
        if abs(onrun['current'] - offrun['current'])/onrun['current'] > condition['current']: return
    print(f'{onrun} {offrun}')
    return offrun['run_id']


def get_offruns(onruns):
    
    previous_epoch = ''
    previous_epoch_runs = []
    
    with open(onrun_list, 'r') as f:
        onruns = f.read().splitlines()
    onruns.sort()
    
    on_off_pairs = []
    for onrun in onruns:
        epoch = get_epoch(onrun)
        el, az, current = get_run_el_az_current(onrun)
        text = f'Extragalactic runs in {epoch} with '
        if 'elevation' in condition.keys(): text += 'el {} +- {} '.format(el, condition['elevation'])
        if 'azimuth' in condition.keys(): text += 'az {} +- {} '.format(az, condition['azimuth'])
        if 'current' in condition.keys(): text += 'current {} +- {}%'.format(current, condition['current']*100)
        print(text)

        if epoch == previous_epoch:
            good_egruns_dict = previous_epoch_runs
        else:
            previous_epoch = epoch
            egrun_dict_file = f'glat10_good_runs/{epoch}_good_egruns_el_az_current.json'
            if os.path.isfile(egrun_dict_file):
                with open(egrun_dict_file, 'r') as fp:
                    good_egruns_dict = json.load(fp)
            else:
                good_egruns = get_epoch_good_egruns(epoch)
                good_egruns_dict = []
                with Pool() as pool:
                    for result in pool.imap(multi_get_run_el_az_current, good_egruns):
                        good_egruns_dict.append(result)
                good_egruns_dict = [egruns_dict for egruns_dict in good_egruns_dict if egruns_dict is not None]
                with open(egrun_dict_file, 'w+') as fp:
                    json.dump(good_egruns_dict, fp)
            previous_epoch_runs = good_egruns_dict
        
        on_off_pair = [({'run_id': onrun, 'elevation': el, 'azimuth': az, 'current': current},
                        egrun_dict) for egrun_dict in good_egruns_dict]
        offruns = []
        with Pool() as pool:
            for result in pool.imap_unordered(multi_get_offruns, on_off_pair):
                offruns.append(result)
        offruns = [offrun for offrun in offruns if offrun is not None]
        if not len(offruns) > 0: offruns = [None]
        with open(offrun_list, 'a+')as f:
            for offrun in offruns: f.write(f'{onrun} {offrun}\n')
            

def get_run_dur(run_id):
    # setup database connection
    dbcnx=pymysql.connect(host='romulus.ucsc.edu', db='VOFFLINE', user='readonly', cursorclass=pymysql.cursors.DictCursor)
    # connect to database
    crs=dbcnx.cursor()
    query = f'SELECT run_id,usable_duration FROM tblRun_Analysis_Comments WHERE run_id={run_id}'
    crs.execute(query)
    res = crs.fetchall()
    return res[0]['usable_duration']
     
        
def get_run_el_az_current(run_id):

    # setup database connection
    dbcnx=pymysql.connect(host='romulus.ucsc.edu', db='VERITAS', user='readonly', cursorclass=pymysql.cursors.DictCursor)
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
    
    return el_avg, az_avg, current_avg
                       

def get_run_el_az(run_id):

    # setup database connection
    dbcnx=pymysql.connect(host='romulus.ucsc.edu', db='VERITAS', user='readonly', cursorclass=pymysql.cursors.DictCursor)
    # connect to database
    crs=dbcnx.cursor()
    
    query = f'SELECT run_id,data_start_time,data_end_time FROM tblRun_Info WHERE run_id={run_id}'
    crs.execute(query)
    res = crs.fetchall()
    
    timestamp_start = '%s'%(res[0]['data_start_time'])
    timestamp_end = '%s'%(res[0]['data_end_time'])
    if timestamp_start=='None': return 0, 0
    if timestamp_end=='None': return 0, 0
    timestamp_start = timestamp_start.replace(' ','').replace('-','').replace(':','')
    timestamp_end = timestamp_end.replace(' ','').replace('-','').replace(':','')
    timestamp_start += '000'
    timestamp_end += '000'
    #print ('timestamp_start = %s'%(timestamp_start))
    el_avg_run = 0.
    az_avg_run = 0.
    total_entries = 0.
    query = "SELECT elevation_target,azimuth_target FROM tblPositioner_Telescope1_Status WHERE timestamp>%s AND timestamp<%s"%(timestamp_start,timestamp_end)
    crs.execute(query)
    res = crs.fetchall()
    for x in res:
        el_avg_run += x['elevation_target']
        az_avg_run += x['azimuth_target']
        total_entries += 1.
    if total_entries==0.:
        return 0., 0.
    el_avg_run = el_avg_run/total_entries*180./math.pi
    az_avg_run = az_avg_run/total_entries*180./math.pi
    #print ('run_id = %s, el_avg_run = %s, az_avg_run = %s'%(run_id,el_avg_run,az_avg_run))
    return epoch, el_avg_run, az_avg_run


def get_run_timecut(runs):

    # setup database connection
    dbcnx=pymysql.connect(host='romulus.ucsc.edu', db='VOFFLINE', user='readonly', cursorclass=pymysql.cursors.DictCursor)
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
            time_cut_mask.append([run_id, timecuts])
    return time_cut_mask


def write_anasum_timemask(runs):
    
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
        
        
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Query the DB and perform various tasks.")
    parser.add_argument('mode', help="timemask / offrun")
    parser.add_argument('runlist', help="path to input runlist")
    
    args = parser.parse_args()
    
    if args.mode == 'timemask':
        # Query and write the time maskes for the input runs
        with open(args.runlist, 'r') as f:
            runs = f.read().splitlines()
        write_anasum_timemask(runs)
    
    if args.mode == 'offrun':
        # Get off runs for runwise background
        offrun_list_name = args.runlist.split('.')[0]
        offrun_list = f'{offrun_list_name}_off.txt'
        i = 0
        while i < 5:
            try:
                with open(args.runlist, 'r') as f:
                    onruns = f.read().splitlines()
                onruns = [int(run) for run in onruns]
                if os.path.isfile(offrun_list):
                    with open(offrun_list, 'r') as f:
                        finished = f.read().splitlines()
                    processed = [int(run.split(' ')[0]) for run in finished]
                    processed = list(set(processed))
                    processed.sort()
                    runs_to_process = [onrun for onrun in onruns if onrun not in processed]
                else: runs_to_process = onruns
                get_offruns(runs_to_process)
            except pymysql.err.OperationalError as e:
                i += 1
                print(e)
                time.sleep(60)
            else:
                i = 5
                with open(offrun_list, 'r') as f:
                    pairs = f.read().splitlines()
                offruns = [pair.split(' ')[-1] for pair in pairs]
                offruns = [pair for pair in offruns if pair != 'None']
                offruns = list(set(offruns))
                offruns.sort()
                with open(f'{offrun_list_name}_off_download.txt', 'w+') as f:
                    for run in offruns: f.write(f'{run}\n')