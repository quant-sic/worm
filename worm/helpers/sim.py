import subprocess
from dataclasses import dataclass
from pathlib import Path
import numpy as np
import h5py
from typing import Optional
from collections import defaultdict

@dataclass
class WormInputParameters:

    mu: np.ndarray
    t_hop:float = 1.0
    U_on:float = 4.0
    V_nn:float = 0.0
    model: str = "BoseHubbard"
    runtimelimit: int = 10000
    sweeps: int = 25000
    thermalization: int = 100
    Lx: int = 4
    Ly: int = 4
    Lz: int = 1
    pbcx: int = 1
    pbcy: int = 1
    pbcz: int = 1
    beta: float = 20.0
    nmax: int = 3
    E_off: float = 1.0
    canonical: int = -1
    seed: int = 30
    Ntest: int = 10000000
    Nsave: int = 100000000
    Nmeasure: int = 1
    Nmeasure2: int = 10
    C_worm: float = 2.0
    p_insertworm: float = 1.0
    p_moveworm: float = 0.3
    p_insertkink: float = 0.2
    p_deletekink: float = 0.2
    p_glueworm: float = 0.3

    h5_path: Optional[Path] = None

    @classmethod
    def from_ini(cls, ini_path: Path):

        # Read ini file
        with open(ini_path, "r") as f:
            lines = f.readlines()

        # Create dictionary
        params = {}
        for line in lines:
            if not line.startswith("#"):
                key, value = line.split("=")
                params[key.strip()] = value.strip()

        # Create input parameters
        return cls(**params)
    
    def save_h5(self):
    
        if self.h5_path is None:
            raise RuntimeError("h5_path must be set")        

        # create parent directory if it does not exist
        self.h5_path.parent.mkdir(parents=True, exist_ok=True)

        # Create h5 file
        with h5py.File(self.h5_path, 'w') as file:
            file["/mu"] = self.mu.flatten()
            file["/t_hop"] = self.t_hop
            file["/U_on"] = self.U_on
            file["/V_nn"] = self.V_nn

    def to_ini(self,checkpoint,outputfile,save_path:Path):

        # create parent directory if it does not exist
        save_path.parent.mkdir(parents=True, exist_ok=True)

        # Create ini file
        with open(save_path, "w") as f:

            for key in self.__dataclass_fields__.keys():
                if not key=="mu":
                    f.write(f"{key} = {self.__getattribute__(key)}\n")
            
            if self.h5_path is None:
                raise RuntimeError("h5_path must be set")        
            else:        
                f.write(f"site_arrays = {self.h5_path}\n")

            f.write(f"outputfile = {outputfile}\n")
            f.write(f"checkpoint = {checkpoint}\n")


@dataclass
class WormOutput:

    out_file_path: Path

    @property
    def observables(self):
        
        h5_file = h5py.File(self.out_file_path,'r')

        observables_dict = {}

        observables_dict = defaultdict(dict)

        for obs, obs_dataset in h5_file["simulation/results"].items():
            for measure,value in obs_dataset.items():
                if isinstance(value,h5py.Dataset):
                    observables_dict[obs][measure] = value[()]
                
                elif isinstance(value,h5py.Group):
                    observables_dict[obs][measure] = {}
                    for sub_measure, sub_value in value.items():
                        observables_dict[obs][measure][sub_measure] = sub_value[()]

        return observables_dict


class WormSimulation(object):

    def __init__(self, input_parameters: WormInputParameters, worm_executable: Path, save_dir, njobs=1):
        
        self.input_parameters = input_parameters
        self.worm_executable = worm_executable
        self.njobs = njobs
        self.save_dir = save_dir

        self.ini_path = self.save_dir/"parameters.ini"

    def save_parameters(self):

        h5_path = self.save_dir/"parameters.h5"
        self.input_parameters.h5_path = h5_path

        outputfile = self.save_dir/"output.h5"
        checkpoint = self.save_dir/"checkpoint.h5"

        self.input_parameters.to_ini(save_path=self.ini_path,checkpoint=checkpoint,outputfile=outputfile)
        self.input_parameters.save_h5()

    
    def run(self):

        subprocess.run(
        f"export TMPDIR=/tmp && mpirun -np {self.njobs} {self.worm_executable} {self.ini_path}", shell=True)


    def get_results(self):
        
        outputfile = self.save_dir/"output.h5"
        output = WormOutput(out_file_path=outputfile)

        return output