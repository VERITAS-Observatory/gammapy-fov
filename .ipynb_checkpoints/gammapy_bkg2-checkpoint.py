import math
from regions import CircleSkyRegion, PixCoord
import astropy.units as u
from astropy.io import fits
from astropy.table import Table, unique, vstack
from astropy.coordinates import SkyCoord
from gammapy.irf import Background3D, Background2D
from gammapy.maps import WcsGeom, WcsNDMap, RegionGeom, Map
from gammapy.data import DataStore, Observations, EventList
from gammapy.analysis import Analysis, AnalysisConfig
from gammapy.datasets import MapDatasetEventSampler, Datasets, MapDatasetOnOff, SpectrumDataset
from gammapy.makers import ReflectedRegionsFinder, RingBackgroundMaker, SpectrumDatasetMaker, ReflectedRegionsBackgroundMaker
from gammapy.visualization import plot_spectrum_datasets_off_regions
from gammapy.stats import CashCountsStatistic
from copy import deepcopy
import numpy as np
from scipy.stats import norm
import json
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from tqdm import tqdm
from astropy.convolution import convolve, Gaussian2DKernel
import warnings
from multiprocessing import Pool
import os

wdir = os.environ['VTS_GP']
starcat = f'{wdir}/data/Hipparcos_MAG8_1997.dat'
vtscat_fits = f'{wdir}/data/vtscat.fits'


def star_mask(center_pos, minmag=6, r_star=0.3):
    geom = WcsGeom.create(skydir=center_pos, width="3.5 deg")
    with open(starcat, 'r') as f:
        stars = f.read().splitlines()
    stars = [star for star in stars if '#' not in star]
    stars = stars[2:]
    stars = [star.split('\t') for star in stars]
    for i in range(len(stars)):
        stars[i] = [star.strip() for star in stars[i]]
    stars = [star for star in stars if '' not in star]
    for i in range(len(stars)):
        stars[i] = [float(star) for star in stars[i]]
    stars = [star for star in stars if star[-1] + star[-2] <= 6]
    stars = [SkyCoord(star[0], star[1], unit='deg') for star in stars]
    mask_star = [geom.contains(star) for star in stars]
    mask_star = [star_coord[0] for star_coord in zip(stars, mask_star) if star_coord[1][0] == True]
    return [[pos.ra.value, pos.dec.value, r_star] for pos in mask_star]

def src_mask(center_pos, r_src=0.4):
    vtscat = Table.read(vtscat_fits, hdu=1)
    vtscat = unique(vtscat, keys='NAME')
    vtscat = [SkyCoord(vts['RA'], vts['DEC'], unit='deg') for vts in vtscat]
    geom = WcsGeom.create(skydir=center_pos, width="3.5 deg")
    mask_src = [geom.contains(vtscat_coord) for vtscat_coord in vtscat]
    mask_src = [vtscat_coord[0] for vtscat_coord in zip(vtscat, mask_src) if vtscat_coord[1][0] == True]
    regions = [[pos.ra.value, pos.dec.value, r_src] for pos in mask_src]
    return regions

def region_from_mask(mask):
    return [CircleSkyRegion(center=SkyCoord(reg[0], reg[1], unit='deg'), radius=reg[2]*u.deg) for reg in mask]

def prepare_mask(geom, center_pos, star=True, minmag=6, r_star=0.3, src=True, r_src=0.4, r_exclude=None):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")    
        mask = []
        if star: mask += star_mask(center_pos, minmag, r_star)
        if src: mask += src_mask(center_pos, r_src)
        if r_exclude is not None:
            mask += [[center_pos.ra.value, center_pos.dec.value, r_exclude]]
        mask = region_from_mask(mask)
        mask = ~geom.region_mask(mask)
        return mask

def prepare_on_region(center_pos, theta):
    return CircleSkyRegion(center=center_pos, radius=theta)



def weight_by_duration(runinfo, name, mode):
    data = np.array([run[name] for run in runinfo])
    weights = np.array([run['duration'] for run in runinfo])
    if mode == 'average': return np.average(data, weights=weights)
    elif mode == 'data': return data * weights / np.average(weights)

def plot_runs_info(runinfo, thresh_az=None, prefix=None, elmin=50, elmax=90, ):
    src_id = np.array([run['run_id'] for run in runinfo])
    src_el = np.array([run['elevation'] for run in runinfo])
    src_nsb = np.array([run['current'] for run in runinfo])
    src_az = np.array([run['azimuth'] for run in runinfo])
    
    f = plt.figure(figsize=(24, 4))
    
    ax1 = plt.subplot(141)
    src_el_weighted = weight_by_duration(runinfo, 'elevation', mode='data')
    x = np.linspace(np.min(src_el),np.max(src_el),30)
    ax1.hist(src_el,density=True,bins=x)
    mu, std = norm.fit(src_el_weighted)
    p = norm.pdf(x, mu, std)
    ax1.plot(x, p, lw=2, color="black")
    ax1.set_title(f'elevation (mu: {mu:.4}, std: {std:.4})')
    ax1.set_xlim(elmin,elmax)
    ax1.annotate(prefix, (0.1,0.9), xycoords='axes fraction', fontsize=16)
    
    ax2 = plt.subplot(142)
    ax2.scatter(src_id, src_el, s=10)
    ax2.hlines(mu,np.min(src_id),np.max(src_id),lw=3,color='r')
    ax2.fill_between([np.min(src_id),np.max(src_id)], mu-std, mu+std, color='r', alpha=0.2)
    ax2.set_title(f'elevation ({mu:.4} +- {std:.4})')
    ax2.set_ylim(elmin,elmax)
    
    ax3 = plt.subplot(143)
    src_az_avg = weight_by_duration(runinfo, 'azimuth', mode='average')
    ax3.scatter(src_id, src_az, s=10)
    ax3.hlines(src_az_avg,np.min(src_id),np.max(src_id),lw=3,color='r')
    title = f'azimuth ({src_az_avg:.4}'
    if thresh_az is not None:
        ax3.fill_between([np.min(src_id),np.max(src_id)], src_az_avg-thresh_az, src_az_avg+thresh_az, color='r', alpha=0.2)
        title += f' +- {thresh_az}'
    ax3.set_title(f'{title})')
    ax3.set_ylim(0,360)
    
    ax4 = plt.subplot(144)
    x = np.linspace(np.min(src_nsb),np.max(src_nsb),30)
    ax4.hist(src_nsb,density=True,bins=x)
    mu, std = norm.fit(src_nsb)
    p = norm.pdf(x, mu, std)
    ax4.plot(x, p, lw=2, color="black")
    ax4.set_title(f'current (mu: {mu:.3}, std: {std:.2})')
    ax4.set_xlim(0,16)
    
    if prefix is not None: f.savefig(f'{prefix}_runinfo.png', bbox_inches='tight')
#     plt.close(f)


def write_idx(bkg_model, run_id, src_name, hdu_name, dl3_dir, flush=True):
    hdu = bkg_model.to_table_hdu()
    hdu.name = hdu_name
    file_name = f'{src_name}.bkg'
    with fits.open(f'{dl3_dir}/{file_name}', mode='append') as hdu_list:
        hdu_list.append(hdu)
    
    ds = DataStore.from_dir(dl3_dir)
    hdu_table = deepcopy(ds.hdu_table)
    if flush:
        if 'bkg_3d' in hdu_table['HDU_CLASS']:
                hdu_table.remove_rows(hdu_table['HDU_CLASS'] == 'bkg_3d')
    rows = []
    for obs in run_id:
        idfilter = hdu_table['OBS_ID'] == int(obs)
        hdufilter = hdu_table['HDU_CLASS'] == 'bkg_3d'
        hdu_table.remove_rows(idfilter*hdufilter)
        row = {
            "OBS_ID": int(obs),
            "HDU_TYPE": "bkg",
            "HDU_CLASS": "bkg_3d",
            "FILE_DIR": "",
            "FILE_NAME": file_name,
            "HDU_NAME": hdu_name
        }
        rows.append(row)

    hdu_table_bkg = Table(rows=rows)
    hdu_table = vstack([hdu_table, hdu_table_bkg])
    hdu_table.sort("OBS_ID")

    hdu_list = fits.HDUList()
    hdu = fits.BinTableHDU(hdu_table)
    hdu.name = "HDU_INDEX"
    hdu_list.append(hdu)
    hdu_list.writeto(dl3_dir+'/hdu-index.fits.gz', overwrite=True)
    Table.read(dl3_dir+'/hdu-index.fits.gz')
    
def get_available_obs(data_store, obs_ids):
    obs_ids = deepcopy(obs_ids)
    i = True
    while i == True:
        try:
            if len(obs_ids) > 0: obs = data_store.get_observations(obs_ids, required_irf='all-optional')
            else: return None
        except ValueError as e:
            print(e)
            obs_ids.remove(int(str(e).split(' ')[2]))
        else: i = False
    return obs

def select_hadrons(obs, energy=None):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        hadrons = obs.events
        skydir = obs.pointing_radec
        hadrons = hadrons.select_region(f'fk5;circle({skydir.ra.value},{skydir.dec.value},1.75)')
        hadrons = hadrons.select_row_subset(hadrons.table['IS_GAMMA']==False)
        if energy is not None: hadrons = hadrons.select_energy(energy.bounds)
    return hadrons

def select_gammas(obs, energy=None):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        gammas = obs.events
        skydir = obs.pointing_radec
        gammas = gammas.select_region(f'fk5;circle({skydir.ra.value},{skydir.dec.value},1.75)')
        gammas = gammas.select_row_subset(gammas.table['IS_GAMMA']==True)
        if energy is not None: gammas = gammas.select_energy(energy.bounds)
    return gammas

def hadron_rate(observations, energy):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        hadron_counts = 0
        observation_time = 0
        for obs in tqdm(observations):
            hadrons = select_hadrons(obs, energy)
            hadron_counts += len(hadrons.energy)
            observation_time += obs.observation_time_duration
    return hadron_counts / observation_time

def convolve_nonzero(data2d, kernelsz=1):
    if kernelsz == 0: convolved_data = data2d
    else:
        kernel = Gaussian2DKernel(kernelsz)
        convolved_data = convolve(data2d, kernel)
    if np.any(convolved_data <= 0):
        if np.max(convolved_data) <= 0: return None
        else:
            convolved_data[np.where(convolved_data <= 0)] = np.sort(list(set(convolved_data.flatten())))[1]/100
    return convolved_data



def blank_sky(obs, check_only=False, r_src=0.4, binsz=0.05, width=3.5):
        
    """
    Excise sources and stars, then fill the hole with events from a nearby region with the same offset.
    1. Find holes (star, source).
    2. Find available patches to fill the holes.
    3. Rotate the events in each patch in the pixel coordinates to fill the corresponding hole.
    """

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        
        pointing = obs.pointing_radec
        regions = star_mask(pointing)
        regions += src_mask(pointing, r_src)
        
        if len(regions) == 0:
            copy_obs = obs.copy()
            events = select_gammas(copy_obs)
            events.stack(select_hadrons(copy_obs))
            copy_obs = copy_obs.copy(in_memory=True, events=events)
            return copy_obs
        
        else:
            copy_obs = obs.copy()
            geom = WcsGeom.create(binsz=binsz, frame='icrs', skydir=pointing, width=width*u.deg)
            regions_hole = region_from_mask(regions)
            mask = geom.region_mask(regions_hole)
            rrfinder = ReflectedRegionsFinder(binsz=binsz*u.deg)
            regions_patch = []
            for reg in regions_hole:
                rr = rrfinder.run(reg, pointing, ~mask)[0]
                if len(rr) == 0:
                    print(f'{obs.obs_id}: No region available to fill the hole.')
                    return None
                else:
                    regions_patch.append(rr[-1])
#                     mask = geom.region_mask(regions_hole+regions_patch)            
            if check_only: return copy_obs
            
            events_excise = select_gammas(copy_obs)
            for reg in regions:
                ra, dec, radius = reg
                events_excise = events_excise.select_region(f'fk5;-circle({ra},{dec},{radius})')
                
            events_patch = events_excise.select_region(regions_patch)
            center_pix = PixCoord(geom.center_pix[0],geom.center_pix[1])
            for idx, (hole, patch) in enumerate(zip(regions_hole, regions_patch)):
                hole_pix = geom.coord_to_pix(hole.center)
                patch_pix = geom.coord_to_pix(patch.center)
                hole_ang = math.atan2(center_pix.xy[1]-hole_pix[1],center_pix.xy[0]-hole_pix[0])
                patch_ang = math.atan2(center_pix.xy[1]-patch_pix[1],center_pix.xy[0]-patch_pix[0])
                diff_ang = hole_ang - patch_ang

                events_fill = events_patch.select_region(patch)
                fill_pix = [geom.coord_to_pix(radec) for radec in events_fill.radec]
                fill_pix = [PixCoord(pix[0],pix[1]).rotate(center_pix,diff_ang) for pix in fill_pix]
                fill_pix = [geom.pix_to_coord(pix.xy) for pix in fill_pix]
                events_fill_table = events_fill.table
                for i in range(len(events_fill_table)):
                    events_fill_table[i]['RA'] = fill_pix[i][0].value
                    events_fill_table[i]['DEC'] = fill_pix[i][1].value
                events_fill = EventList(events_fill_table)
                if idx > 0:
                    events_fill = events_fill.select_mask(~geom.region_mask(regions_hole[:idx]))
                events_excise.stack(events_fill)
#                 events_patch = events_patch.select_mask(~geom.region_mask(patch))

            events_excise.stack(select_hadrons(copy_obs))    
            copy_obs = copy_obs.copy(in_memory=True, events=events_excise)
            return copy_obs

def shifted_sky(obspair, r_src=0.4, binsz=0.05, width=3.5):
    """
    Rotate the events in the pixel coordinates.
    Change the pointing information to match the target pointing coordinates.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        obs, refobs = obspair
        pointing = obs.pointing.radec
        refpointing = refobs.pointing.radec
        copy_obs = obs.copy()
        geom_old = WcsGeom.create(binsz=binsz, frame='icrs', skydir=pointing, width=width)
        geom_new = WcsGeom.create(binsz=binsz, frame='icrs', skydir=refpointing, width=width)
        evlist = select_gammas(copy_obs)
        evlist = evlist.table

        for idx, evt in enumerate(evlist):
            old_pix = geom_old.coord_to_pix(SkyCoord(evt['RA'], evt['DEC'], unit='deg'))
            new_sky = geom_new.pix_to_coord(old_pix)
            evlist[idx]['RA'] = new_sky[0].value
            evlist[idx]['DEC'] = new_sky[1].value

        evlist.meta['RA_PNT'] = refpointing.ra.value
        evlist.meta['DEC_PNT'] = refpointing.dec.value
        copy_obs = copy_obs.copy(in_memory=True, events=EventList(evlist), pointing=refobs.pointing)
        return copy_obs



def plot_obs(observations, gamma_only=True, mask=None, binsz=0.05, width=3.5, **kwargs):
    nrows = int(np.ceil(len(observations)/8))
    f = plt.figure(figsize=(8*2,nrows*2))
    for idx, observation in enumerate(observations):
        pointing = observation.pointing.radec
        geom = WcsGeom.create(binsz=binsz, frame='icrs', skydir=pointing, width=width)
        skymap = WcsNDMap(geom)
        if gamma_only: events = select_gammas(observation)
        else: events = observation.events
        skymap.fill_events(events)
        if mask is not None: skymap *= mask
        ax = plt.subplot(nrows, 8, idx+1, projection=skymap.geom.wcs)
        skymap.plot(ax=ax, **kwargs)
        ax.annotate(observation.obs_id,(0.03,0.9),xycoords='axes fraction',color='w')
    f.tight_layout()

    
def lima17(Non, Noff, alpha):
    # Eq. 17 in Li & Ma https://ui.adsabs.harvard.edu/abs/1983ApJ...272..317L/abstract
    minusMLLR = Non*np.log((1+alpha)/alpha*Non/(Non+Noff)) + Noff*np.log((1+alpha)*Noff/(Non+Noff))
    minusMLLR = np.sqrt(minusMLLR)
    return np.sign(Non-alpha*Noff)*np.sqrt(2)*minusMLLR

    
def sigmap(dataset, kernel=None, stat='cash'):
    mask = dataset.mask_safe
    n_on = (dataset.counts * mask).sum_over_axes()  
    if isinstance(dataset, MapDatasetOnOff):
        raise Exception("This function is meant for the FOV method. Otherwise use ExcessMapEstimator!")
    else:
        background = (dataset.npred() * mask).sum_over_axes()
        if stat=='cash':
            sqrt_ts = CashCountsStatistic(n_on.data, background.data).sqrt_ts
        elif stat=='lima':
            with np.errstate(invalid="ignore", divide="ignore"):
                sqrt_ts = np.nan_to_num(lima17(n_on.data, background.data, np.ones(n_on.data.shape)))
        sqrt_ts = Map.from_geom(dataset.counts.geom.squash('energy'), data=sqrt_ts)
        if kernel is not None:
            pixel_size = np.mean(np.abs(dataset.counts.geom.wcs.wcs.cdelt))
            size = kernel.to(u.deg).value / pixel_size
            sqrt_ts_conv = sqrt_ts.convolve(Gaussian2DKernel(size).array)
        else: sqrt_ts_conv = sqrt_ts.copy()
        sqrt_ts.data[~mask.sum_over_axes()] = np.nan
        sqrt_ts_conv.data[~mask.sum_over_axes()] = np.nan
        return {'sqrt_ts': sqrt_ts, 'sqrt_ts_conv': sqrt_ts_conv}
    

def plot_sig(significance_map, exclusion_mask, exclusion_radius, significance_map_conv=None,
             vmin=-5, vmax=5, prefix=None):
    
    mask_copy = exclusion_mask.copy()*1.0
    mask_copy.data[mask_copy.data==0.0] = np.inf
    significance_map_off = significance_map * mask_copy
    significance_all = significance_map.data[np.isfinite(significance_map.data)]
    significance_off = significance_map_off.data[np.isfinite(significance_map_off.data)]

    f = plt.figure(figsize=(13, 10))

    ax1 = plt.subplot(221, projection=significance_map.geom.wcs)
    if significance_map_conv is not None:
        significance_map_conv.plot(ax=ax1, add_cbar=True, vmin=vmin, vmax=vmax, cmap='coolwarm')
    else: significance_map.plot(ax=ax1, add_cbar=True, vmin=vmin, vmax=vmax, cmap='coolwarm')
    CircleSkyRegion(
        significance_map.geom.center_skydir,exclusion_radius).to_pixel(
        significance_map.geom.wcs).plot(
        ax=ax1, lw=2)


    x = np.linspace(-6, 6, 60)
    mu, std = norm.fit(significance_off)
    p = norm.pdf(x, 0, 1)
    significance_all[significance_all <= -6] = -6
    significance_all[significance_all >= 6] = 6
    significance_off[significance_off <= -6] = -6
    significance_off[significance_off >= 6] = 6   

    ax2 = plt.subplot(222)
    ax2.hist(
        significance_all,
        density=True,
        alpha=0.5,
        label="all bins",
        bins=x,
    )

    ax2.hist(
        significance_off,
        density=True,
        facecolor="None",
        alpha=0.5,
        edgecolor="k",
        label="off bins",
        bins=x,
    )

    ax2.plot(x, p, lw=2, color="black")
    ax2.legend()
    ax2.set_xlabel("Significance")
    ax2.set_yscale("log")
    ax2.set_ylim(1e-5, 1)
    ax2.set_xlim(-6, 6)
    ax2.annotate(f"mu = {mu:.2f}\nstd = {std:.2f}", xy=(0.05, 0.9), xycoords="axes fraction")
    
    if prefix is not None: f.savefig(f'{prefix}_sig_dist.png', bbox_inches='tight')  
    
    
def plot_bkg(bkgrate, prefix=None):
    nrows = int(np.rint(len(bkgrate.data)/2))
    extent_lon = bkgrate.axes['fov_lon'].edges
    extent_lon = np.linspace(extent_lon[0],extent_lon[-1],6)
    extent_lat = bkgrate.axes['fov_lat'].edges
    extent_lat = np.linspace(extent_lat[0],extent_lat[-1],6)
    extent = [extent_lon[0].value, extent_lon[-1].value, extent_lat[0].value, extent_lat[-1].value]
    energy_edges = np.round(bkgrate.axes['energy'].edges,1)
    f, axes = plt.subplots(ncols=2, nrows=nrows, figsize=(10,nrows*4))
    for idx, ax in enumerate(axes.flat[:len(bkgrate.data)]):
        im = ax.imshow(bkgrate.data[idx], origin='lower', extent=extent)
        ax.set_xticks(extent_lon.value)
        ax.set_yticks(extent_lat.value)
        ax.set_xlabel(f'DETX ({extent_lon.unit})')
        ax.set_ylabel(f'DETY ({extent_lat.unit})')
        divider = make_axes_locatable(ax)
        colorbar_axes = divider.append_axes("right",size="5%",pad=0.2)
        cbar = f.colorbar(im, cax=colorbar_axes)
        cbar.set_label(f'{bkgrate.unit}')
        ax.set_title(energy_edges[idx:idx+2])
    f.tight_layout()
    

def bias_correction(stacked_mimics, kernelsz=1, prefix=None):
    stacked_counts = WcsNDMap(stacked_mimics[0].counts.geom)
    stacked_background = WcsNDMap(stacked_mimics[0].background.geom)
    for stacked in stacked_mimics:
        stacked_counts += stacked.counts
        stacked_background += stacked.npred()
    stacked_counts = stacked_counts.convolve(Gaussian2DKernel(kernelsz).array)
    stacked_background = stacked_background.convolve(Gaussian2DKernel(kernelsz).array)
    with np.errstate(invalid="ignore", divide="ignore"):
        correction = stacked_counts / stacked_background
    energy_edges = np.round(bkgrate.axes['energy'].edges,1)
    f, axes = plt.subplots(ncols=3, nrows=2, figsize=(15,8), subplot_kw={'projection': correction.geom.wcs})
    for idx, ax in enumerate(axes.flat):
        correction.slice_by_idx({'energy': idx}).plot(ax=ax, add_cbar=True, vmin=0, vmax=2, cmap='coolwarm')
        ax.set_title(energy_edges[idx:idx+2])
    correction.data = np.nan_to_num(correction.data)
    if prefix is not None:
        f.savefig(f'{prefix}_bias.png', bbox_inches='tight')
        np.save(f'{prefix}_bias.npy', correction)
    return correction    
    
    
def plot_fov(fov_norms, dataname, prefix=None):
    fov_norms = np.array(fov_norms)
    f = plt.figure(figsize=(6, 5))
    ax = plt.gca()
    x_lo, x_hi = [0.6, 1.8]
    x = np.linspace(x_lo,x_hi,21)
    mu, std = norm.fit(fov_norms)
    p = norm.pdf(x, mu, std)
    fov_norms[fov_norms < x_lo] = x_lo
    fov_norms[fov_norms > x_hi] = x_hi
    ax.hist(fov_norms,density=True,bins=x)
    ax.plot(x, p, lw=2, color="black")
    ax.annotate(f'{dataname} FoV norm\n{mu:.2} +- {std:.2}',(0.7,0.9),xycoords='axes fraction')
    if prefix is not None: f.savefig(f'{prefix}_fov_norms.png', bbox_inches='tight')    
    
    
def load_src_info(srcname):
    with open(f'{wdir}/data/{srcname}_ep_wob_dur_el_az_cur.json', 'r') as fp:
        srcinfo = json.load(fp)
    return srcinfo

def load_eg_info(oa=False, na=False):
    eg_suffix = '_good_egruns_dur_el_az_cur_srcid.json'
    if oa and na:
        raise Exception('Do not analyze OA and NA runs together.')
    if oa:
        with open(f'{wdir}/data/oa{eg_suffix}', 'r') as fp:
            eginfo = json.load(fp)
    elif na:
        with open(f'{wdir}/data/na{eg_suffix}', 'r') as fp:
            eginfo = json.load(fp)
    else:
        eginfo = []
        files = os.listdir(f'{wdir}/data')
        files = [file for file in files if eg_suffix in file]
        files = [file for file in files if 'oa' not in file and 'na' not in file]
        for file in files:
            with open(f'{wdir}/data/{file}', 'r') as fp:
                eginfo += json.load(fp)
    return eginfo

    
class MimicDatasetMaker:
    
    def __init__(self, data_store, srcinfo, eginfo, az_thresh):
        
        self.data_store = data_store
        self.srcinfo = srcinfo
        self.eginfo = eginfo
        self.az_thresh = az_thresh
        self.mimic_datasets = []
                
    def find_mimic_runs(self, update=True, times=1):
        eginfo_copy = deepcopy(self.eginfo)
        eg_dur = np.sum([eg['duration'] for eg in eginfo_copy])
        if eg_dur < np.sum([run['duration'] for run in self.srcinfo]):
            raise Exception('Not enough candidate runs.')

        mimic_info = []
        for src in tqdm(self.srcinfo):
            for idx, run in enumerate(eginfo_copy):
                eginfo_copy[idx]['elevation_abdiff'] = np.absolute(src['elevation'] - run['elevation'])
            eginfo_copy = sorted(eginfo_copy, key=lambda x: x['elevation_abdiff'])
            goal_dur = src['duration']*times
            while goal_dur > 0:
                for idx, run in enumerate(eginfo_copy):
                    if src['azimuth'] - self.az_thresh < 0:
                        if run['azimuth'] >= src['azimuth']+self.az_thresh:
                            if run['azimuth'] <= src['azimuth']-self.az_thresh+360: continue
                    elif src['azimuth'] + self.az_thresh > 360:
                        if run['azimuth'] <= src['azimuth']-self.az_thresh:
                            if run['azimuth'] >= src['azimuth']+self.az_thresh-360: continue
                    else:
                        if run['azimuth'] >= src['azimuth']+self.az_thresh: continue
                        if run['azimuth'] <= src['azimuth']-self.az_thresh: continue
                    selected = eginfo_copy.pop(idx)
                    selected.update({'onrun_id': src['run_id']})
                    mimic_info.append(selected)
                    goal_dur -= selected['duration']
                    break
                                
        if update:
            self.eginfo = [run for run in eginfo_copy if run not in mimic_info]
            self.mimic_datasets.append({'information': mimic_info})
        else: return (mimic_info, eginfo_copy)
   
    def run_blank_sky(self, idx, update=True):
        obsids = [run['run_id'] for run in self.mimic_datasets[idx]['information']]
        obs = get_available_obs(self.data_store, obsids)
        obs = self._run_blank_sky(obs)
        if update: self.mimic_datasets[idx]['observations'] = obs
        else: return obs
    
    def run_shifted_sky(self, idx, update=True):
        src_observations = get_available_obs(self.data_store, [run['run_id'] for run in self.srcinfo])
        obs = self.mimic_datasets[idx]['observations']
        obs = self._run_shifted_sky(src_observations, obs,
                                    self.mimic_datasets[idx]['information'])
        if update: self.mimic_datasets[idx]['observations'] = obs
        else: return obs

    @staticmethod
    def _run_blank_sky(observations):
        copy_obs = Observations()
        with Pool() as pool:
            for result in tqdm(pool.imap(blank_sky, observations)):
                if result is not None: copy_obs.append(result)
        return copy_obs
    
    @staticmethod
    def _run_shifted_sky(src_observations, shift_observations, obsinfo):
        shift_ids = shift_observations.ids
        refobs = Observations()
        for obs in shift_ids:
            refobs_id = [run['onrun_id'] for run in obsinfo if str(run['run_id']) == obs][0]
            refobs.append(src_observations[src_observations.index(str(refobs_id))])
        obspair = zip(shift_observations, refobs)

        copy_obs = Observations()
        with Pool() as pool:
            for result in tqdm(pool.imap(shifted_sky, obspair)):
                copy_obs.append(result)
        return copy_obs


        
class Background3DModelEstimator:
    
    def __init__(self, energy, fov_lon, fov_lat):
        self.counts = self._make_bkg3d(energy, fov_lon, fov_lat, unit="")
        self.exposure = self._make_bkg3d(energy, fov_lon, fov_lat, unit="s TeV sr")

    @staticmethod
    def _make_bkg3d(energy, fov_lon, fov_lat, unit):
        return Background3D(axes=[energy, fov_lon, fov_lat], unit=unit)

    def run(self, observations):
        for obs in tqdm(observations):
            self.fill_counts(obs)
            self.fill_exposure(obs)

    def fill_counts(self, obs):
        energy, fov_lon, fov_lat = self.counts.axes
        events = select_gammas(obs, energy)
        events = MapDatasetEventSampler.event_det_coords(obs, events)
        counts = np.histogramdd(
            (events.energy.to('TeV'), events.table['DETX'], events.table['DETY']),
            (energy.edges, fov_lon.edges, fov_lat.edges))[0]        
        self.counts.data += counts
            
    def fill_exposure(self, obs):
        time = obs.observation_time_duration
        bin_volume = self.exposure.axes.bin_volume()
        self.exposure.quantity += time * bin_volume

    @property
    def background_rate(self):
        rate = deepcopy(self.counts)
        rate.quantity /= self.exposure.quantity
        return rate
        
        
        
class Background2DModelEstimator:
    """"""

    def __init__(self, energy, offset):
        self.counts = self._make_bkg2d(energy, offset, unit="")
        self.exposure = self._make_bkg2d(energy, offset, unit="s TeV sr")

    @staticmethod
    def _make_bkg2d(energy, offset, unit):
        return Background2D(axes=[energy, offset], unit=unit)

    def run(self, observations):
        for obs in tqdm(observations):
            self.fill_counts(obs)
            self.fill_exposure(obs)

    def fill_counts(self, obs):
        energy, offset = self.counts.axes
        events = select_gammas(obs, energy)
        counts = np.histogram2d(
            x=events.energy.to("TeV"),
            y=events.offset.to("deg"),
            bins=(energy.edges, offset.edges))[0]
        self.counts.data += counts

    def fill_exposure(self, obs):
        axes = self.exposure.axes
        offset = axes["offset"].center
        time = obs.observation_time_duration
        exposure = 2 * np.pi * offset * time * axes.bin_volume()
        self.exposure.quantity += exposure
        
    @property
    def background_rate(self):
        rate = deepcopy(self.counts)
        rate.quantity /= self.exposure.quantity
        return rate        
        
        
        
def run_rbm(geom, geom_image, energy_axis_true, radius_in, radius_width,
            mask_1d, mask_2d, on_region, datasets, containment_correction=True):
    
    spec_on_off = Datasets()
    ring_maker_1D = RingBackgroundMaker(
        r_in=radius_in, width=radius_width, exclusion_mask=mask_1d)
    for dataset in datasets:
        dataset_on_off = ring_maker_1D.run(dataset)
        dataset_on_off = dataset_on_off.to_spectrum_dataset(on_region, containment_correction=containment_correction)
        spec_on_off.append(dataset_on_off)
    spec_on_off = spec_on_off.stack_reduce()
    
    image_on_off = MapDatasetOnOff.create(
        geom=geom_image, energy_axis_true=energy_axis_true, name='stacked')
    ring_maker = RingBackgroundMaker(
        r_in=radius_in, width=radius_width, exclusion_mask=mask_2d)
    for dataset in datasets:
        dataset = dataset.to_image()
        dataset_on_off = ring_maker.run(dataset)
        image_on_off.stack(dataset_on_off)
        
    return (spec_on_off, image_on_off)

def plot_onoff(image_on_off, filename=None):
    names = ['counts', 'counts_off', 'alpha', 'background']
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10,7))
    for name, ax in zip(names, axes.flat):
        if name == 'counts':
            image_on_off.counts.plot(ax=ax, add_cbar=True)
        if name == 'counts_off':
            image_on_off.counts_off.plot(ax=ax, add_cbar=True)
        if name == 'alpha':
            image_on_off.alpha.plot(ax=ax, add_cbar=True)
        if name == 'background':
            image_on_off.background.plot(ax=ax, add_cbar=True)
        ax.set_title(name)
        if filename is not None: f.savefig(filename)

def run_rr(srcobs, on_region, exclusion_mask, energy_axis, energy_axis_true,
           safe_mask_maker, reduce=True, containment_correction=True):
    geom_region = RegionGeom.create(region=on_region, axes=[energy_axis])
    empty = SpectrumDataset.create(geom=geom_region, energy_axis_true=energy_axis_true)
    spec_dataset_maker = SpectrumDatasetMaker(containment_correction=containment_correction, use_region_center=False,
                                              selection=['counts', 'exposure', 'edisp'])
    rr_maker = ReflectedRegionsBackgroundMaker(exclusion_mask=exclusion_mask)
    spec_datasets = Datasets()
    for observation in srcobs:
        dataset = spec_dataset_maker.run(empty.copy(name=str(observation.obs_id)), observation)
        dataset = safe_mask_maker.run(dataset, observation)
        dataset = rr_maker.run(dataset, observation)
        spec_datasets.append(dataset)
    f = plt.figure()
    ax = exclusion_mask.plot()
    on_region.to_pixel(ax.wcs).plot(ax=ax, edgecolor="k")
    plot_spectrum_datasets_off_regions(ax=ax, datasets=spec_datasets)
    plt.show()
    plt.close(f)
    if reduce: spec_datasets = spec_datasets.stack_reduce()
    return spec_datasets