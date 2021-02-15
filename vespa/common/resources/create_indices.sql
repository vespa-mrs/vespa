CREATE INDEX
    simulations_metabolite_id_index
ON
    simulations (metabolite_id);

CREATE INDEX
    experiment_dims_experiment_id_index
ON
    experiment_dims (experiment_id);

CREATE INDEX
    simulations_dims_id_index
ON
    simulations (dims_id);

CREATE INDEX
    metabolite_spins_metabolite_id_index
ON
    metabolite_spins (metabolite_id);

CREATE INDEX
    j_couplings_spin1_id_index
ON
    j_couplings (spin1_id);

CREATE INDEX
    j_couplings_spin2_id_index
ON
    j_couplings (spin2_id);


