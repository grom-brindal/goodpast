import h5py
import numpy as np
import os
import torch
import torchani
import ase
from ase import Atom
from ase import Atoms
from ase import units
from ase.io import read, write
from ase.optimize import BFGS
from ase.md import langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.io import read

device = torch.device('cpu')
path = os.getcwd()

const_file = './rHCNO-5.2R_16-3.5A_a4-8.params'
sae_file = './sae_linfit_CAS.dat'

consts = torchani.neurochem.Constants(const_file)
aev_computer = torchani.AEVComputer(**consts)
energy_shifter = torchani.neurochem.load_sae(sae_file)

try:
    os.mkdir('./output/')
except OSError:
    pass

network = '../best.pt'
network_prefix = '/home/goodpast/johan406/sorted_ethylene/2round_start_scratch/10round/training/1round_original_AEVs/gs/1e-1_displacement/'
network_prefix = '/home/goodpast/johan406/sorted_ethylene/3round_cohesive_molecule/1round_paring_down/training/1round/gs/1e-1_displacement/'
network_prefix = '/home/goodpast/johan406/sorted_ethylene/5round_nms_exploding/training/2round_1700datasets/1round_original_AEVs/gs/1e-1_displacement/'
network = [network_prefix + str(i) + 'round/best.pt' for i in range(10)]

print(network)

def atomic(linear=[(384,128), (128,128), (128,64), (64,1)]):
    model = torch.nn.Sequential(
        torch.nn.Linear(*linear[0]),
        torch.nn.CELU(0.1),
        torch.nn.Linear(*linear[1]),
        torch.nn.CELU(0.1),
        torch.nn.Linear(*linear[2]),
        torch.nn.CELU(0.1),
        torch.nn.Linear(*linear[3])
    )
    return model
nn = torchani.ANIModel([atomic(v) for v in (
[(384,160), (160,128), (128,96), (96,1)],
[(384,144), (144,112), (112,96), (96,1)],
[(384,128), (128,112), (112,96), (96,1)],
[(384,128), (128,112), (112,96), (96,1)]
)])


#if os.path.isfile(network):
#    nn.load_state_dict(torch.load(network,map_location=torch.device('cpu')))
#else:
#    print('Error: No Potentials Detected')
#    exit()

#model = torch.nn.Sequential(aev_computer, nn, energy_shifter)


mol = [ Atom('C', (0.795463, 0.0, 0.0)),
        Atom('C', (-0.795463, 0.0, 0.0)),
        Atom('H', (1.138234, 1.022448, 0.0)),
        Atom('H', (1.138234, -1.022448, 0.0)),
        Atom('H', (-1.138234, 1.022448, 0.0)),
        Atom('H', (-1.138234, -1.022448, 0.0)) ]
#mol = read('25000_velocities.traj')
#mol = [ Atom('C', (-0.704225, -0.030802, -0.021961)),
#        Atom('C', (0.713126, -0.019103, -0.038462)),
#        Atom('H', (-1.201186, 0.957176, -0.031633)),
#        Atom('H', (-1.184440, -1.026849, -0.000917)),
#        Atom('H', (1.193350, 0.976946, -0.059477)),
#        Atom('H', (1.210095, -1.007073, -0.028778))]
       
mol = [ Atom('C', (  -0.000000,    0.663585,   -0.000000)),
        Atom('C', (  -0.000000,   -0.663585,    0.000000)),
        Atom('H', (  -0.000000,    1.235740,    0.924684)),
        Atom('H', (  -0.000000,    1.235740,   -0.924684)),
        Atom('H', (  -0.000000,   -1.235740,    0.924684)),
        Atom('H', (  -0.000000,   -1.235740,   -0.924684)) ]

mol = Atoms(mol)
def printenergy(a=mol):
    """Function to print the potential, kinetic and total energy."""
    epot = a.get_potential_energy() / len(a)
    ekin = a.get_kinetic_energy() / len(a)
    print('Energy per atom: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK)  '
          'Etot = %.3feV' % (epot, ekin, ekin / (1.5 * units.kB), epot + ekin))



def printxyz(a=mol):
    """Function to print out xyz coordinates"""
    check_out=os.getcwd() + '/output'
    os.chdir(check_out)
    i=0
    file_name = str(i) + '_checkpoint.pdb'
    file_name = str(file_name)
    while (os.path.isfile(file_name)):
        i=i+1
        file_name = str(i)+'_checkpoint.pdb'
        file_name = str(file_name)
    checkpoint_out = str(i) +  '_checkpoint.pdb'
    write(str(checkpoint_out), mol)
    os.chdir('..')



def printtraj(a=mol):
    """Function to print out trajectory file"""
    check_out=os.getcwd() + '/output'
    os.chdir(check_out)
    i=0
    file_name = str(i) + '_velocities.traj'
    file_name = str(file_name)
    while (os.path.isfile(file_name)):
        i=i+1000
        file_name = str(i)+'_velocities.traj'
        file_name = str(file_name)
    checktraj_out = str(i) +  '_velocities.traj'
    write(str(checktraj_out), mol)
    os.chdir('..')

device = torch.device('cpu')

null = []
personal_directory = os.getcwd()
with open(personal_directory + '/outfile.txt', 'w') as f:
    f.write('')
    f.close()

#calculator = torch.nn.Sequential(aev_computer, nn, energy_shifter)
calculator = torchani.ase.Calculator(personal_directory, consts.species, aev_computer, network, nn, energy_shifter, device)
mol.set_calculator(calculator)
mol.set_cell((20.0*np.identity(3)))
mol.set_pbc=(True)
opt = BFGS(mol)
#opt.run(fmax=0.001)
MaxwellBoltzmannDistribution(mol, units.kB * 5.0, force_temp=False, rng=np.random)

langevinsim = langevin.Langevin(mol, 0.5 * units.fs, temperature = 5.0 * units.kB, friction=0.01)
langevinsim.attach(printxyz, interval=1)
langevinsim.attach(printenergy, interval=100)
langevinsim.attach(printtraj, interval=1000)
langevinsim.run(500)

exit()
