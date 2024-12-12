import pickle
from pathlib import Path

import numpy as np

data_path = Path(__file__).resolve().parents[1]
data_path = data_path.joinpath('data', 'dp0')


class BandSampler:

    full_labels = ['median', 'rms', 'seeing', 'exp_time', 'zero_point']

    def __init__(self, gmm_sampler, gmm_labels, constants=None):
        self.gmm_sampler = gmm_sampler
        self.gmm_labels = gmm_labels
        if constants is None:
            constants = {}
        self.constants = constants
        assert set(gmm_labels) | set(constants.keys()) == set(self.full_labels)

    def sample(self):
        samples = self.gmm_sampler.sample()[0][0]
        samples = {label: samples[i] for i, label in enumerate(self.gmm_labels)}
        samples.update({k: v for k, v in self.constants.items() if k not in samples})
        samples['rms'] = np.maximum(samples['rms'], 1e-6)
        return samples


def dp0_gmm_sampler(dp0_coadd="1"):
    if dp0_coadd == "single_visit":
        gmm_labels = ['median', 'rms', 'seeing', 'zero_point']
        gmm_sampler = pickle.load(open(data_path.joinpath('dp0_gmm_sampler_single_visit.pkl'), 'rb'))
        return {band: BandSampler(gmm_sampler[band], gmm_labels, {'exp_time': 30.0})
                for band in 'grizy'}
    elif dp0_coadd == "1":
        gmm_labels = ['median', 'rms', 'seeing', 'exp_time']
        gmm_sampler = pickle.load(open(data_path.joinpath('dp0_gmm_sampler_1yr_coadd.pkl'), 'rb'))
        return {band: BandSampler(gmm_sampler[band], gmm_labels, {'zero_point': 28.75})
                for band in 'grizy'}
    else:
        raise ValueError("dp0_coadd must be '1' or 'single_visit'")
