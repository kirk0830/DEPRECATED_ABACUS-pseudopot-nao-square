import os
import json

def generate_siab_input(
        exe_mpi: str, exe_pw: str, element: str, ecut: float, rcut: float,
        pseudo_dir: str, pseudo_name: str, sigma: float,
        reference_structure: dict,
        max_steps: int, 
        orbital_configurations: dict,
        save_orbitals: dict
):

    return_str = '''#--------------------------------------------------------------------------------
#1. CMD & ENV
 EXE_mpi      mpiexec.hydra -np 20 
 EXE_pw       /home/nic/wszhang/abacus/abacus222_intel-2018u4/ABACUS.mpi

#-------------------------------------------------------------------------------- 
#2. Electronic calculatation
 element     As  # element name 
 Ecut        100  # cutoff energy (in Ry)
 Rcut        7  # cutoff radius (in a.u.)
 Pseudo_dir  /home/nic/wszhang/abacus/delta_dft/CIF_POT/SG15_ONCV_PBE-1.0/
 Pseudo_name As_ONCV_PBE-1.0.upf
 sigma       0.01 # energy range for gauss smearing (in Ry)

#--------------------------------------------------------------------------------
#3. Reference structure related parameters for PW calculation
#For the built-in structure types (including 'dimer', 'trimer' and 'tetramer'):
#STRU Name   #STRU Type  #nbands #MaxL   #nspin  #Bond Length list 
 STRU1       dimer       9       2       1      1.8 2.1 2.5 3.0 4.0
 STRU2       trimer      12      2       1      2.0 2.3 2.7

#-------------------------------------------------------------------------------- 
#4. SIAB calculatation
 max_steps    9000
#Orbital configure and reference target for each level
#LevelIndex  #Ref STRU name  #Ref Bands  #InputOrb    #OrbitalConf 
 Level1      STRU1           auto        none        1s1p      
 Level2      STRU1           auto        fix         2s2p1d    
 Level3      STRU2           auto        fix         3s3p2d    

#--------------------------------------------------------------------------------
#5. Save Orbitals
#Index    #LevelNum   #OrbitalType 
 Save1    Level1      SZ
 Save2    Level2      DZP
 Save3    Level3      TZDP
'''
    return return_str

def automatical_ptg_dpsi(work_status: dict, element: str, pseudopotential: str):

    """ Generate numerical atomic orbitals in batch, use as-prepared SG15-type pseudopotential 
    corresponding orbital generation input script as reference """
    os.chdir(work_status["numerical_orbitals"]["resources"]) # check-in
    if not os.path.isdir("sg15_10"):
        print("warning: sg15_10 not found, will download from web")
        if not os.path.exists("resources.json"):
            print("error: resources.json not found, please check")
        else:
            with open(work_status["work_folder"] + "\\resources.json", "r") as f:
                resources = json.load(f)
                website = resources["numerical_orbitals"]["resources"]["sg15_10"]
                os.system("wget " + website)
                os.system("unzip SG15-Version1p0__AllOrbitals-Version2p0.zip")
                os.system("rm SG15-Version1p0__AllOrbitals-Version2p0.zip")
                os.system("mv SG15-Version1p0__AllOrbitals-Version2p0 sg15_10")
        # after download, check again
        if not os.path.isdir("sg15_10"):
            print("error: sg15_10 not found, please check")
            return
    # check-in
    os.chdir(work_status["work_folder"])