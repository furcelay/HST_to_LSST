import pickle
from pathlib import Path


data_path = Path(__file__).resolve().parents[1]
data_path = data_path.joinpath('data', 'dp0')
dp0_gmm_sampler = pickle.load(open(data_path.joinpath('dp0_gmm_sampler.pkl'), 'rb'))
