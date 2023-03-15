import os

os.chdir('datagen')
for i in range(0, 500):
    os.system('molpro ' + str(i) + '.in')
