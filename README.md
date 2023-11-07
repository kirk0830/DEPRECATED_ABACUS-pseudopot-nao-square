# ABACUS psesudopotentials and numerical orbitals test workflow
*Author: kirk0830*
## Environments and dependencies
### Executables
* Quantum ESPRESSO `I1d.x` executable is needed to generate PSlibrary 0.3.1 and 1.0.0 version PAW and US pseudopotential files. The executable is not included in the repository.
* PTG_dpsi executable is needed to generate numerical orbitals for specific kinds of pseudopotential.
* ABACUS `abacus` executable is needed to generate numerical orbitals that provides reference wavefunction.
### Python packages
* `numpy` is needed for handling numerical data.
* `scipy` is needed for some special functions.
* `matplotlib` is needed for plotting.
* `json` and `time` are needed for storing work logs.
### Python program list
* `pseudopot_test_generation.py`: generate pseudopotential files for specific kinds of pseudopotential, use Quantum ESPRESSO HSE functional.
* `nao_test_generation.py`: generate numerical orbitals for specific kinds of pseudopotential. Use ABACUS.
* `pseudopot_test_analysis.py`: band structure analysis, calculates band gap and plot
* `pseudopot_test_fold_band.py`: band folding, calculates DOS
### File structure
```
├── readme.md
├── resources.json
├── pseudopot_test_generation.py
├── nao_test_generation.py
├── pseudopot_test_analysis.py
├── pseudopot_test_fold_band.py
├── numerical_orbitals/
│  ├── resources/
│  │   ├── sg15_10 ([pseudopot_kind]_[pseudopot_version]_[appendix])
│  │   │   ├── 33_As_DZP ([atomic_number]_[element]_[basis_set])
│  │   │   │   ├── As_gga_7au_100Ry_2s2p1d.orb
│  │   │   │   └── ...
│  │   │   └── ...
│  │   ├── ...
│  │   └── tools/
│  │       ├── atom_in.py
│  │       ├── orb_io.py
│  │       ├── cal_ovlp_jjtilde.py
│  │       ├── ABACUS_orb_parser.py
│  │       ├── nao_oscillation_remote.py
│  │       └── ...
│  ├── results/
│  │   ├── ...
│  │   └── ...
│  ├── templates/
│  │   ├── INPUT
│  │   ├── KPT
│  │   └── STRU_* (for each system to test)
│  ├── config.json (for ABACUS Test setting)
│  └── ...
├── pseudopotentials/
│  ├── resources/
│  │   ├── sg15_10 ([pseudopot_kind]_[pseudopot_version]_[appendix])
│  │   ├── pd_03
│  │   │   ├── As_*.upf
│  │   │   └── ...
│  │   ├── pd_04
│  │   ├── pd_04_d
│  │   └── ...
│  ├── results/
│  │   ├── t_*
│  │   └── ...
│  ├── templates/
│  │   └── template_*.in (for each system to test)
│  ├── config.json (for ABACUS Test setting)
│  └── ...
├── analysis/
│  └── [developing]
└── ...
```
## Usage
### Configuration
#### Files auto-download
In `./resources.json`, add/change `resources` section to download files automatically. `jsons` section will be useful when preparing input scripts.
#### Pseudopotential generation
For cases such as `kjpaw`, `rrkjus` or `pslnc` those pseudopotential kinds belonging to PSlibrary, pseudopotential may need to generate with Quantum ESPRESSO `Id1.x` code. Therefore this executable must be ready. If so, change the variable `qe_path` in `resources.json` to the path of folder, where there is `Id1.x` executable.
#### Numerical orbital generation
For cases such as non-SG15, numerical atomic orbital is needed to generate. Please follow the handson of PTG_dpsi to install executables of ABACUS, PTG_dpsi and other dependencies. Change the variable `abacus_path` in `resources.json` to the path of folder, where there is `abacus` executable.  
Numerical atomic orbitals of other kinds of pseudopotentials will be generated based on parameters set in `SIAB_INPUT` of SG15-kind pseudopotential corresponding numerical atomic orbital generation.
### Input scripts generation
* execute `pseudopot_test_generation.py` to generate test input scripts.
* execute `nao_test_generation.py` to generate numerical orbitals.
### Test run
visit ABACUS Test website developed by ABACUS team to submit test jobs.  
website: https://labs.dp.tech/projects/abacustest/
#### Pseudopotential test
* use `normal| rundft, postdft` module, upload a zip file. For more test configurations, please refer to `./pseudpotentials/config.json`.
#### Numerical orbital test
* use `expert| predft, rundft, postdft` module, upload a zip file. For more test configurations, please refer to `./numerical_orbitals/config.json`.
### Analysis
* execute `pseudopot_test_analysis.py` to analyze band structure and calculate band gap. For more functions of it, please read in-file comments.
## Phylosophy of organizing information in jsons
### resources.json
TBD...
### config.json
TBD...
### work_status
`work_status` is the input variable of `pseudopot_test_generation.py`, `pseudopot_test_analysis.py` and `nao_test_generation.py`. It always has such following keywords:
```
{
    "work_folder": "path/to/work/folder",
    "pseudopot_kinds": ["pseudopot_kind_1", "pseudopot_kind_2", ...],
    "versions": ["version_1", "version_2", ...],
    "appendices": ["appendix_1", "appendix_2", ...],
    ...
    "systems": {
        "system1": {
            "experimental_value": 1.0,
        }
        "system2": {
            "experimental_value": 2.0,
        }
        ...
    }
}
```
In this json, user should give enough information for the whole test, included but not limited to systems description, sometimes the postprocessing such as plot is also needed.
### test_status
test_status is a on-the-fly or in-build json, which will pass from one function to another inside program. If there is no special case, it will not be a real output. test_status describes with more details not only for systems, but may also carry data.
The design of test_status is shown in following:
```
{
    "work_folder": "./numerical_orbitals",
    "InAs": {
        "sg15_10_7_DZP_pbe": {
            "elements": [
                "In",
                "As"
            ],
            "pseudopotentials": [
                "In_ONCV_PBE-1.0.upf",
                "As_ONCV_PBE-1.0.upf"
            ],
            "numerical_orbitals": [
                "In_gga_7au_100Ry_2s2p2d1f.orb",
                "As_gga_7au_100Ry_2s2p1d.orb"
            ]
        }
}
```
Therefore except `work_folder` and `root`, other keys are all system names. Under each system, there are test names, and under each pseudopotential, there are elements and numerical orbitals.  
- For pseudopotential tests, the name of test is organized as [pseudopot_kind]_ [pseudopot_version]_ [appendix]_ [functional].  
- For numerical atomic orbital tests, the name of test is organized as [pseudopot_kind]_ [pseudopot_version]_ [appendix]_ [cutoff radius]_ [zeta information]_ [functional].
## Other tools
### Numerical atomic orbital analysis
Due to present method to generate numerical atomic orbitals has some problems, it is necessary to check the quality of generated numerical atomic orbitals. `./numerical_orbitals/resources/tools/ABACUS_orb_parser.py` is a tool to check the oscillation of numerical atomic orbitals.  
If indeed there is high frequency, unphysical oscillation in orbitals generated, use `./numerical_orbitals/resources/tools/nao_oscillation_remote.py` to remove high frequency oscillation and generate new numerical atomic orbitals.
