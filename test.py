import abacus_pseudopotential_square.base.quantum_espresso as qe

pseudopotentials_in = {
    "files": {
        "La": "La_pseudopot",
        "O": "O_pseudopot",
    },
}
qe.cif_to_quantum_espresso("mp-26.cif", pseudopotentials_in, {})