# ABACUS psesudopotentials and numerical orbitals test workflow
*Author: kirk0830*
## To use
```bash
python -m main
```
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
├── analysis/
│  ├── __init__.py
│  ├── analysis.py
│  ├── astr.py
│  ├── abg.py
│  ├── fold_band.py
│  └── ...
├── base/
│  ├── __init__.py
│  ├── atom_in.py
│  ├── file_io.py
│  ├── scans.py
│  └── ...
├── configure/
│  ├── __init__.py
│  ├── configure.py
│  ├── nao_generation.py
│  ├── resources.py
│  ├── ...
│  └── support/
│      ├── resources.json
│      └── ...
├── examples/
│  ├── test_status_template.json
│  ├── work_status_template.json
│  └── ...
├── generation/
│  ├── __init__.py
│  ├── generation.py
│  ├── glpt.py
│  ├── gppt.py
│  └── ...
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
│  │   ├── Quantum_ESPRESSO/
│  │   │   ├── template_*.in (for each system to test)
│  │   │   └── ...
│  │   ├── ABACUS/
│  │   │   ├── INPUT
│  │   │   ├── KPT
│  │   │   └── STRU_* (for each system to test)
│  └── ...
├── main.py
├── __init__.py
└── ...
```
## Usage
### Configuration
#### Files auto-download (under development)
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
`work_status` is the placeholder of input parameters. For test generation and analysis, different input files are used. For test generation, it always has such following keywords:
```json
{
    "work_folder": ".",
    "basis_type": "lcao",
    "test_mode": "pseudopotential",
    "software": "ABACUS",
    "functionals": ["PBE", "HSE06"],
    "systems": {
        "InAs": {
            "experimental_values": {
                "band_gap": 1.35,
                "stress": -0.0001
            }
        },
        "InSb": {
            "experimental_values": {
                "band_gap": 0.7,
                "stress": -0.0001
            }
        }
    },
    "pseudopotentials": {
        "pseudopot_paths": {
            "work": "./pseudopotentials",
            "resources": "./pseudopotentials/resources"
        },
        "kinds": {
            "In": ["sg15", "pd"],
            "As": ["sg15", "pd"],
            "Sb": ["sg15", "pd"]
        },
        "versions": {
            "In": ["10", "04"],
            "As": ["10", "04"],
            "Sb": ["10", "04"]
        },
        "appendices": {
            "In": [""],
            "As": [""],
            "Sb": [""]
        }
    },
    "numerical_orbitals": {
        "nao_paths": {
            "work": "./numerical_orbitals",
            "resources": "./numerical_orbitals/resources"
        },
        "types": {
            "In": ["DZP", "TZDP"],
            "As": ["DZP"],
            "Sb": ["DZP", "TZDP"]
        },
        "rcuts": {
            "In": [7, 8, 9, 10],
            "As": [10],
            "Sb": [10]
        },
        "appendices": {
            "In": ["", "osci_rm"],
            "As": [""],
            "Sb": [""]
        }
    }
}
```
In this json, user should give enough information for the whole test, included but not limited to systems description, sometimes the postprocessing such as plot is also needed.
For test analysis tasks, the input json always has the following information:
```json
{
    "test_status": "test_status.json",
    "work_folder": ".",
    "basis_type": "lcao",
    "test_mode": "pseudopotential",
    "software": "ABACUS",
    "systems": {
        "InAs": {
            "experimental_values": {
                "band_gap": 1.35,
                "stress": -0.0001
            }
        },
        "InSb": {
            "experimental_values": {
                "band_gap": 0.7,
                "stress": -0.0001
            }
        }
    },
    "analysis": {
        "band_gap": {
            "statistics": {
                "draw": false,
                "draw_2d": true,
                "2d": {
                    "mode": "imshow",
                    "cmap": "jet",
                    "ncol": 2
                }
            }
        },
        "dos": {
            "draw": false,
            "window_width": 4,
            "delta_e": 0.01,
            "smear": 0.02
        },
        "stress": {

        }
    }
}
```
### test_status
test_status is a on-the-fly or in-build json, which will pass from one function to another inside program. If there is no special case, it will not be a real output. test_status describes with more details not only for systems, but may also carry data.
The design of test_status is shown in following:
```json
    "systems": {
        "InAs": {
            "sg1510sg1510_D7D10": {
                "elements": [
                    "In",
                    "As"
                ],
                "pseudopotentials": {
                    "files": {
                        "In": "In_ONCV_PBE-1.0.upf",
                        "As": "As_ONCV_PBE-1.0.upf"
                    },
                    "info": {
                        "In": {
                            "kind": "sg15",
                            "version": "10",
                            "appendix": ""
                        },
                        "As": {
                            "kind": "sg15",
                            "version": "10",
                            "appendix": ""
                        }
                    }
                },
                "numerical_orbitals": {
                    "files": {
                        "In": "In_gga_7au_100Ry_2s2p2d1f.orb",
                        "As": "As_gga_10au_100Ry_2s2p1d.orb"
                    },
                    "info": {
                        "In": {
                            "type": "DZP",
                            "rcut": 7,
                            "appendix": ""
                        },
                        "As": {
                            "type": "DZP",
                            "rcut": 10,
                            "appendix": ""
                        }
                    }
                }
            },
        ...
    }
    ...
}
```
Under each system, there are test names, and under each pseudopotential, there are elements and numerical orbitals.  
- For pseudopotential tests, the name of test is organized as [pseudopot_kind]_ [pseudopot_version]_ [appendix]_ [functional].  
- For numerical atomic orbital tests, the name of test is organized as [pseudopot_kind]_ [pseudopot_version]_ [appendix]_ [cutoff radius]_ [zeta information]_ [functional].
## Other tools
### Numerical atomic orbital analysis
Due to present method to generate numerical atomic orbitals has some problems, it is necessary to check the quality of generated numerical atomic orbitals. `./numerical_orbitals/resources/tools/ABACUS_orb_parser.py` is a tool to check the oscillation of numerical atomic orbitals.  
If indeed there is high frequency, unphysical oscillation in orbitals generated, use `./numerical_orbitals/resources/tools/nao_oscillation_remote.py` to remove high frequency oscillation and generate new numerical atomic orbitals.