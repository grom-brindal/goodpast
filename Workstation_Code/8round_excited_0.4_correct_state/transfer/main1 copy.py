import os
import h5py
import numpy as np
import matplotlib
from matplotlib import pyplot as plt



folders = ['1.2k_0.1', '12k_0.1', '1.2k_0.2', '12k_0.2']
folders1 = ['1.2k_0.1', '1.2k_0.2']
folders2 = ['12k_0.1', '12k_0.2']
coordinates = []
energy1s, energy2s, energy3s = [], [], []
allcoords = []
atype = np.array((b'C',b'C',b'H',b'H',b'H',b'H'))
natoms = 6

for i in range(1, 1201):
    #print(i)
    try:
        with open(str(i) + '.out', 'r') as f:
            lines = f.readlines()
            for y in range(len(lines)):
                if 'Results for state 1.1' in lines[y]:
                    energy1 = lines[y+3].split()
                    energy1 = energy1[4]
                    energy1s.append(energy1)
                if 'Results for state 2.1' in lines[y]:
                    energy2 = lines[y+3].split()
                    energy2 = energy2[4]
                    energy2s.append(energy2)
                if 'Results for state 3.1' in lines[y]:
                    energy3 = lines[y+3].split()
                    energy3 = energy3[4]
                    energy3s.append(energy3)
                if 'geometry={' in lines[y]:
                    for u in range(6):
                        coordinates.append(lines[y+1+u].split())
                if len(coordinates) != 0:
                    for e in range(len(coordinates)):
                        if len(coordinates[e]) != 0:
                            if len(coordinates[e]) == 4:
                                coordinates[e] = [coordinates[e][1], coordinates[e][2], coordinates[e][3]]
                            else:
                                pass
            allcoords.append(coordinates)
    except FileNotFoundError:
        pass
ener1 = ["ener%d" % x for x in range(len(energy1s))]
ener2 = ["ener%d" % x for x in range(len(energy2s))]
ener3 = ["ener%d" % x for x in range(len(energy3s))]
coor0 = ["coor%d" % x for x in range(len(allcoords))]

for p in range(len(allcoords)):
    ener1[p] = float(energy1s[p])
    ener2[p] = energy2s[p]
    ener3[p] = energy3s[p]
    coor0[p] = allcoords[p]
x_axis = []
for i in range(1, 1201):
    x_axis.append(round(float(i * 0.002),3))

line1 = plt.plot(x_axis, ener1)
plt.savefig('energy1.png')
exit()
ener1 = np.array(ener1).astype(float)
ener2 = np.array(ener2).astype(float)
ener3 = np.array(ener3).astype(float)
coor0 = np.array(coor0).astype(float)
print(ener1)
print(coor0)
newfile1 = 'energy1.h5'
newfile2 = 'energy2.h5'
newfile3 = 'energy3.h5'

newfile1 = h5py.File(newfile1, 'w')
newfile2 = h5py.File(newfile2, 'w')
newfile3 = h5py.File(newfile3, 'w')

layer1 = newfile1.create_group('firstlayer')
group1 = layer1.create_group('secondlayer')
group1.create_dataset('coordinates', data=coor0)

group1.create_dataset('coordinatesHE', (0,natoms,3), dtype=np.float32)
group1.create_dataset('energies', data=ener1)                           
group1.create_dataset('energiesHE', (0,), dtype=np.float64)
group1.create_dataset('smiles', (3,), dtype='S1', data=(b'C',b'(',b'='))
group1.create_dataset('species',(natoms,),dtype='S1',data=atype)

layer2 = newfile2.create_group('firstlayer')
group2 = layer2.create_group('secondlayer')
group2.create_dataset('coordinates', data=coor0)
group2.create_dataset('coordinatesHE', (0,natoms,3), dtype=np.float32)
group2.create_dataset('energies', data=ener2)                               
group2.create_dataset('energiesHE', (0,), dtype=np.float64)
group2.create_dataset('smiles', (3,), dtype='S1', data=(b'C',b'(',b'='))
group2.create_dataset('species',(natoms,),dtype='S1',data=atype)

layer3 = newfile3.create_group('firstlayer')
group3 = layer3.create_group('secondlayer')
group3.create_dataset('coordinates', data=coor0)
group3.create_dataset('coordinatesHE', (0,natoms,3), dtype=np.float32)
group3.create_dataset('energies', data=ener3)                           
group3.create_dataset('energiesHE', (0,), dtype=np.float64)
group3.create_dataset('smiles', (3,), dtype='S1', data=(b'C',b'(',b'='))
group3.create_dataset('species',(natoms,),dtype='S1',data=atype)
 
    

exit()


