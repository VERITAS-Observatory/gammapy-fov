import json
import argparse
from multiprocessing import Pool
from gammapy.data import Observations
from gammapy_bkg2 import *


dl3dir = os.environ['DL3DIR']
datadir = os.environ['VTS_GP']+'/data'
egdict_suffix = '_good_egruns_dur_el_az_cur_srcid.json'

epochs = f'{datadir}/epochs.json'
with open(epochs, 'r') as fp: epochs = json.load(fp)


def check_blank_sky(obs):
    return blank_sky(obs, check_only=True)

def select_empty_runs(observations):
    copy_obs = Observations()
    with Pool() as pool:
        for result in tqdm(pool.imap(check_blank_sky, observations)):
            if result is not None: copy_obs.append(result)
    return copy_obs


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Find extragalactic runs that can be a blank sky by excising and filling holes.')
    parser.add_argument('-e', '--epoch', nargs='+', choices=epochs.keys(), help='Epochs to process.')
    args = parser.parse_args()
    
    egdict = [f'{datadir}/{epoch}{egdict_suffix}' for epoch in args.epoch]
    ds = DataStore.from_dir(dl3dir)
    
    for idx, epoch in enumerate(args.epoch):
        print(f'{epoch} started.')
        with open(egdict[idx], 'r') as fp: eginfo = json.load(fp)
        egobs = get_available_obs(ds, [run['run_id'] for run in eginfo])
        egobs = select_empty_runs(egobs)
        eginfo = [run for run in eginfo if str(run['run_id']) in egobs.ids]
        with open(egdict[idx], 'w+') as fp: json.dump(eginfo, fp)
        print(f'{epoch} finished.')