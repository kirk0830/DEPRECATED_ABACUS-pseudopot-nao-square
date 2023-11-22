"""
in each folder whose name starts with "t_", sed the calculation = 'scf' to calculation = 'vc-relax'
"""
import os

for folder in os.listdir('.'):
    if folder.startswith('t_') and not folder.endswith('.zip'):
        os.chdir(folder)
        print(os.getcwd())
        os.system("sed -i 's/scf/vc-relax/g' scf.in")
        os.chdir('..')
