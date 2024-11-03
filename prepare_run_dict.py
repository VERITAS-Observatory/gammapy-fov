import argparse
import time
from multiprocessing import Pool
import os
import json
import pymysql
import datetime
import traceback
from DBquery import *

glat_threshold = 10
dur_threshold = datetime.timedelta(seconds=600)
good_weather = ['A+', 'A', 'A-', 'B+', 'B', 'B-']

datadir = os.environ['VTS_GP']+'/data'
targetlist = f'{datadir}/VERITAS_targets.json'
# timemaskfile = os.environ['VERITAS_EVNDISP_AUX_DIR']+'/ParameterFiles/ANASUM.timemask.dat'

epochs = f'{datadir}/epochs.json'
with open(epochs, 'r') as fp: epochs = json.load(fp)

db_host = os.environ['DB_HOST']
db_user = os.environ['DB_USER']

srcdict_suffix = '_ep_wob_dur_el_az_cur.json'
egdict_suffix = '_good_egruns_dur_el_az_cur_srcid.json'


def multi_good_runs(run):
    dbcnx=pymysql.connect(host=db_host, db='VOFFLINE', user=db_user, cursorclass=pymysql.cursors.DictCursor, charset="utf8")
    crs=dbcnx.cursor()
    run_id = run['run_id']
    query = f'SELECT run_id,data_category,status,usable_duration FROM tblRun_Analysis_Comments WHERE run_id={run_id}'
    crs.execute(query)
    res = crs.fetchall()
    dbcnx.close()
    if len(res) == 0: return
    elif res[0]['data_category'] != 'science': return
    elif res[0]['status'] != 'good_run': return
    elif res[0]['usable_duration'] is None: return
    elif res[0]['usable_duration'] < dur_threshold: return
    else:
        print(run)
        return run

def multi_get_run_info(run):
    try:
        tinit = datetime.datetime.now()
        el, az, current = get_run_el_az_current(run['run_id'])
        dur = get_run_dur(run['run_id'])
        run.update({'duration': dur.total_seconds(), 'elevation': el, 'azimuth': az, 'current': current})
        tfin = datetime.datetime.now()
        print(f'{tfin-tinit} {run}')
        return run
    except ZeroDivisionError as e:
        traceback.print_exc()
        return
    
def get_epoch_good_srcruns_info(ep, source_id=None):
    
    # science, good run, weather => B, usable_durateion => 10 min, abs(glat) => glat_threshold
    # 4 tel, run_type = observing (not obsLowHV or obsFilter)
    start_time, end_time = epochs[ep]
    start_time, end_time = [datetime.datetime.strptime(t, '%Y-%m-%d %H:%M:%S') for t in [start_time, end_time]]
    dbcnx=pymysql.connect(host=db_host, db='VERITAS', user=db_user, cursorclass=pymysql.cursors.DictCursor, charset="utf8")
    crs=dbcnx.cursor()
    if source_id is None:
        query = f'SELECT run_id,run_type,data_start_time,data_end_time,weather,config_mask,source_id FROM tblRun_Info WHERE run_type=\'observing\' AND data_start_time>\'{start_time}\' AND data_end_time<\'{end_time}\' AND config_mask=15'
    else:
        query = f'SELECT run_id,run_type,data_start_time,data_end_time,weather,config_mask,offset_distance,source_id FROM tblRun_Info WHERE run_type=\'observing\' AND data_start_time>\'{start_time}\' AND data_end_time<\'{end_time}\' AND config_mask=15 AND source_id=\'{source_id}\''
    crs.execute(query)
    epoch_runs = crs.fetchall()
    dbcnx.close()
    epoch_runs = [run for run in epoch_runs if run['weather'] in good_weather]
    if source_id is None:
        with open(targetlist, 'r') as fp: eg_targets = json.load(fp)
        eg_targets = [target for target in eg_targets if target['glat'] > glat_threshold]
        eg_targets = [target['source_id'] for target in eg_targets]
        epoch_runs = [run for run in epoch_runs if run['source_id'] in eg_targets]
        epoch_runs = [{'run_id': run['run_id'], 'source_id': run['source_id']} for run in epoch_runs]
    else:
        epoch_runs = [{'run_id': run['run_id'], 'epoch': ep, 'wobble': run['offset_distance']} for run in epoch_runs]
    good_runs = []
    with Pool() as pool:
        for result in pool.imap(multi_good_runs, epoch_runs):
            good_runs.append(result)
    good_runs = [run for run in good_runs if run != None]
    run_info = []
    with Pool() as pool:
        for result in pool.imap(multi_get_run_info, good_runs):
            run_info.append(result)
    run_info = [run for run in run_info if run != None]
    return run_info


def get_any_srcrun_info(run_id):
    tinit = datetime.datetime.now()
    try:
        epoch = get_epoch(run_id)
        wobble = get_wobble(run_id)
        dur = get_run_dur(run_id)
        if dur is not None: dur = dur.total_seconds()
        srcid, el, az, current = get_run_srcid_el_az_current(run_id)
        tfin = datetime.datetime.now()
        print(f'{tfin-tinit} {run_id}')
        return {'run_id': run_id, 'source_id': srcid, 'epoch': epoch, 'wobble': wobble,
             'duration': dur, 'elevation': el, 'azimuth': az, 'current': current}
    except IndexError as e:
        print(run_id)
        traceback.print_exc()
        return
    except ZeroDivisionError as e:
        print(run_id)
        traceback.print_exc()
        return

    
            
            
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create run info dictionary for a source / runs (and for epochs).")
    parser.add_argument('-s', '--srcid', help="get all runs for particular source")
    parser.add_argument('-r', '--runid', nargs='+', help="get particular runs (run numbers or text file including runs")
    parser.add_argument('-p', '--prefix', help="prefix of the output json file")
    parser.add_argument('-e', '--epoch', choices=epochs.keys(), nargs='+', help="get extragalactic runs for particular epochs")
    args = parser.parse_args()
    
    if args.srcid is None and args.runid is None and args.epoch is None:
        raise Exception('Provide either srcid, runid, or epoch.')
        
    elif args.srcid is not None:
        if ' ' in args.srcid: srcid = args.srcid.replace(' ', '')
        if '+' in args.srcid: srcid = args.srcid.replace('+', 'p')
        if '-' in args.srcid: srcid = args.srcid.replace('-', 'n')
        dictname = f'{datadir}/{srcid}{srcdict_suffix}'
        if args.epoch is not None: epoch2run = args.epoch
        else: epoch2run = epochs.keys()
        for epoch in epoch2run:
            run_info = get_epoch_good_srcruns_info(epoch, args.srcid)
            if len(run_info) > 0:
                with open(dictname, 'a+') as fp: json.dump(run_info, fp)
            print('{} runs in {} written in {}'.format(args.srcid, epoch, dictname))
#             write_anasum_timemask([run['run_id'] for run in run_info])
        with open(dictname, 'r') as f: srcdict = f.read()
        srcdict = srcdict.replace('][', ', ')
        with open(dictname, 'w+') as f: f.write(srcdict)            
    
    elif args.runid is not None:
        if len(args.runid)==1 and os.path.exists(args.runid[0]):
            with open(args.runid[0], 'r') as f: args.runid = f.read().splitlines()
        args.runid = [int(run) for run in args.runid]
        if args.prefix is None: args.prefix = 'runs'
        dictname = f'{datadir}/{args.prefix}{srcdict_suffix}'
        run_info = []
        with Pool() as pool:
            for result in pool.imap(get_any_srcrun_info, args.runid):
                run_info.append(result)
        run_info = [run for run in run_info if run != None]
        with open(dictname, 'a+') as fp: json.dump(run_info, fp)
        print('{} runs written in {}'.format(len(args.runid), dictname))
        with open(dictname, 'r') as f: srcdict = f.read()
        srcdict = srcdict.replace('][', ', ')
        with open(dictname, 'w+') as f: f.write(srcdict)            
#         write_anasum_timemask([run['run_id'] for run in run_info])  
    
    else:
        dictname = [f'{datadir}/{epoch}{egdict_suffix}' for epoch in args.epoch]
        for epoch, file in zip(args.epoch, dictname):
            run_info = get_epoch_good_srcruns_info(epoch)
            with open(file, 'w+') as fp: json.dump(run_info, fp)
            print('Extragalactic runs in {} written in {}'.format(epoch, file))
#             write_anasum_timemask([run['run_id'] for run in run_info])