CREATE TABLE 'experiments' (
'id' TEXT NOT NULL  PRIMARY KEY REFERENCES 'experiment_metabolites' ('experiment_id'),
'created' TIMESTAMP NOT NULL  DEFAULT CURRENT_TIMESTAMP,
'name' TEXT DEFAULT NULL,
'is_public' BOOLEAN NOT NULL  DEFAULT '0',
'comment' TEXT DEFAULT NULL,
'investigator' TEXT DEFAULT NULL,
'b0' REAL NOT NULL ,
'isotope' TEXT NOT NULL ,
'peak_search_ppm_low' REAL NOT NULL  DEFAULT 0.00,
'peak_search_ppm_high' REAL NOT NULL  DEFAULT 10.0,
'blend_tolerance_ppm' REAL NOT NULL  DEFAULT 0.00150,
'blend_tolerance_phase' REAL NOT NULL  DEFAULT 50.0,
'pulse_sequence_id' TEXT DEFAULT NULL REFERENCES 'pulse_sequences' ('id'),
UNIQUE (name)
);

CREATE TABLE 'pulse_sequences' (
'id' TEXT NOT NULL  PRIMARY KEY,
'name' TEXT DEFAULT NULL,
'is_public' BOOLEAN NOT NULL  DEFAULT '0',
'created' TIMESTAMP NOT NULL  DEFAULT CURRENT_TIMESTAMP,
'creator' TEXT DEFAULT NULL,
'comment' TEXT DEFAULT NULL,
'sequence_code' TEXT DEFAULT NULL,
'binning_code' TEXT DEFAULT NULL
);

CREATE TABLE 'metabolites' (
'id' TEXT NOT NULL  PRIMARY KEY REFERENCES 'experiment_metabolites' ('metabolite_id'),
'name' TEXT NOT NULL ,
'is_public' BOOLEAN NOT NULL  DEFAULT '0',
'created' TIMESTAMP NOT NULL  DEFAULT CURRENT_TIMESTAMP,
'creator' TEXT DEFAULT NULL,
'comment' TEXT DEFAULT NULL,
'deactivated' TIMESTAMP DEFAULT NULL
);

CREATE TABLE 'metabolite_spins' (
'id' INTEGER NOT NULL  PRIMARY KEY AUTOINCREMENT,
'metabolite_id' TEXT NOT NULL  REFERENCES 'metabolites' ('id'),
'isotope' TEXT NOT NULL ,
'chemical_shift' REAL NOT NULL ,
'display_order' INTEGER NOT NULL ,
UNIQUE (metabolite_id, display_order)
);

CREATE TABLE 'j_couplings' (
'id' INTEGER NOT NULL  PRIMARY KEY AUTOINCREMENT,
'value' REAL NOT NULL ,
'spin1_id' INTEGER DEFAULT NULL REFERENCES 'metabolite_spins' ('id'),
'spin2_id' INTEGER DEFAULT NULL REFERENCES 'metabolite_spins' ('id')
);

CREATE TABLE 'simulations' (
'metabolite_id' INTEGER NOT NULL  REFERENCES 'metabolites' ('id'),
'dims_id' INTEGER NOT NULL  REFERENCES 'experiment_dims' ('id'),
'started' TIMESTAMP DEFAULT NULL,
'completed' TIMESTAMP DEFAULT NULL,
'ppms' BLOB DEFAULT NULL,
'areas' BLOB DEFAULT NULL,
'phases' BLOB DEFAULT NULL,
PRIMARY KEY (metabolite_id, dims_id)
);

CREATE TABLE 'isotopes' (
'name' TEXT NOT NULL  PRIMARY KEY,
'display_order' INTEGER NOT NULL 
);

CREATE TABLE 'pulse_sequence_user_static_parameters' (
'id' INTEGER NOT NULL  PRIMARY KEY AUTOINCREMENT,
'pulse_sequence_id' TEXT NOT NULL  REFERENCES 'pulse_sequences' ('id'),
'type' TEXT NOT NULL ,
'name' TEXT NOT NULL ,
'default_value' TEXT NOT NULL ,
'display_order' INTEGER NOT NULL 
);

CREATE TABLE 'pulse_sequence_loops' (
'pulse_sequence_id' TEXT NOT NULL  REFERENCES 'pulse_sequences' ('id'),
'label' TEXT DEFAULT NULL,
'display_order' INTEGER NOT NULL ,
PRIMARY KEY (pulse_sequence_id, display_order)
);

CREATE TABLE 'b0_bins' (
'left' REAL NOT NULL ,
'center' REAL NOT NULL ,
'right' REAL NOT NULL 
);

CREATE TABLE 'experiment_user_static_parameters' (
'experiment_id' TEXT NOT NULL  REFERENCES 'experiments' ('id'),
'value' TEXT NOT NULL ,
'display_order' INTEGER NOT NULL ,
PRIMARY KEY (experiment_id, display_order)
);

CREATE TABLE 'vespa' (
'database_version' INTEGER NOT NULL 
);

CREATE TABLE 'pulse_projects' (
'id' TEXT NOT NULL  PRIMARY KEY,
'is_public' BOOLEAN NOT NULL ,
'name' TEXT NOT NULL ,
'creator' TEXT DEFAULT NULL,
'created' TIMESTAMP NOT NULL ,
'comment' TEXT DEFAULT NULL,
'machine_settings_id' INTEGER NOT NULL  REFERENCES 'machine_settings' ('id'),
'master_parameters_id' INTEGER NOT NULL  REFERENCES 'master_parameters' ('id')
);

CREATE TABLE 'machine_settings' (
'id' INTEGER NOT NULL  PRIMARY KEY AUTOINCREMENT,
'name' TEXT DEFAULT NULL,
'is_template' BOOLEAN NOT NULL ,
'is_default' BOOLEAN DEFAULT NULL,
'machine_type' TEXT NOT NULL  DEFAULT '''',
'max_b1_field' REAL NOT NULL  DEFAULT 0.0,
'field_strength' REAL NOT NULL  DEFAULT 0.0,
'zero_padding' INTEGER NOT NULL  DEFAULT 0,
'min_dwell_time' REAL NOT NULL  DEFAULT 0.0,
'dwell_time_increment' REAL NOT NULL  DEFAULT 0.0,
'gradient_raster_time' REAL NOT NULL  DEFAULT 0.0,
'gradient_slew_rate' REAL NOT NULL  DEFAULT 0.0,
'gradient_maximum' REAL NOT NULL  DEFAULT 0.0
);

CREATE TABLE 'master_parameters' (
'id' INTEGER NOT NULL  PRIMARY KEY AUTOINCREMENT,
'calc_resolution' INTEGER NOT NULL  DEFAULT 0,
'pulse_bandwidth_type' TEXT NOT NULL  DEFAULT ''''
);

CREATE TABLE 'transformations' (
'id' INTEGER NOT NULL  PRIMARY KEY AUTOINCREMENT,
'pulse_project_id' TEXT NOT NULL  REFERENCES 'pulse_projects' ('id'),
'progression' INTEGER NOT NULL  DEFAULT 0,
'transformation_type' TEXT NOT NULL  DEFAULT '''',
'parameters_id' INTEGER NOT NULL  REFERENCES 'interpolate_rescale_parameters' ('id') REFERENCES 'hs_pulse_parameters' ('id') REFERENCES 'slr_pulse_parameters' ('id') REFERENCES 'root_reflect_parameters' ('id') REFERENCES 'gaussian_pulse_parameters' ('id') REFERENCES 'randomized_pulse_parameters' ('id') REFERENCES 'import_pulse_parameters' ('id') REFERENCES 'ocn_parameters' ('id'),
'result_id' INTEGER NOT NULL  REFERENCES 'results' ('id'),
UNIQUE (pulse_project_id, progression)
);

CREATE TABLE 'slr_pulse_parameters' (
'id' INTEGER NOT NULL  PRIMARY KEY AUTOINCREMENT,
'tip_angle' REAL NOT NULL  DEFAULT 0.0,
'time_steps' INTEGER NOT NULL  DEFAULT 0,
'duration' REAL NOT NULL  DEFAULT 0.0,
'bandwidth' REAL NOT NULL  DEFAULT 0.0,
'separation' REAL NOT NULL  DEFAULT 0.0,
'is_single_band' BOOLEAN NOT NULL ,
'nc_phase_subtype' TEXT NOT NULL  DEFAULT '''',
'slr_filter_type' TEXT NOT NULL  DEFAULT '''',
'pass_ripple' REAL NOT NULL  DEFAULT 0.0,
'reject_ripple' REAL NOT NULL  DEFAULT 0.0
);

CREATE TABLE 'interpolate_rescale_parameters' (
'id' INTEGER NOT NULL  PRIMARY KEY AUTOINCREMENT,
'do_interpolate' BOOLEAN NOT NULL ,
'interpolation_factor' INTEGER NOT NULL  DEFAULT 0,
'new_dwell_time' REAL NOT NULL  DEFAULT 0.0,
'do_rescaling' BOOLEAN NOT NULL ,
'angle' REAL NOT NULL  DEFAULT 0.0
);

CREATE TABLE 'results' (
'id' INTEGER NOT NULL  PRIMARY KEY AUTOINCREMENT,
'created' TIMESTAMP NOT NULL ,
'gradient_id' INTEGER NOT NULL  DEFAULT 0 REFERENCES 'gradients' ('id'),
'ocn_state_id' INTEGER NOT NULL  DEFAULT 0 REFERENCES 'ocn_states' ('id')
);

CREATE TABLE 'hs_pulse_parameters' (
'id' INTEGER NOT NULL  PRIMARY KEY AUTOINCREMENT,
'total_rotation' REAL NOT NULL  DEFAULT 0.0,
'time_steps' INTEGER NOT NULL  DEFAULT 0,
'dwell_time' REAL DEFAULT 0.0,
'was_bandwidth_specified' BOOLEAN NOT NULL ,
'quality_cycles' REAL NOT NULL  DEFAULT 0.0,
'power_n' INTEGER NOT NULL  DEFAULT 0,
'sharpness_mu' REAL NOT NULL  DEFAULT 0.0,
'filter_type' TEXT NOT NULL  DEFAULT '''',
'filter_application' REAL NOT NULL  DEFAULT 0.0
);

CREATE TABLE 'rf_waveforms' (
'id' INTEGER NOT NULL  PRIMARY KEY AUTOINCREMENT,
'result_id' INTEGER NOT NULL  REFERENCES 'results' ('id'),
'time_point' REAL NOT NULL  DEFAULT 0,
'real_amplitude' REAL NOT NULL  DEFAULT 0.0,
'imaginary_amplitude' REAL NOT NULL  DEFAULT 0.0
);

CREATE TABLE 'gradients' (
'id' INTEGER NOT NULL  PRIMARY KEY AUTOINCREMENT REFERENCES 'gradient_waveforms' ('gradient_id'),
'linear_gradient_value' REAL NOT NULL  DEFAULT 0.0,
'refocused_gradient' REAL NOT NULL  DEFAULT 0.0,
'frequency_offset' REAL NOT NULL  DEFAULT 0.0
);

CREATE TABLE 'gradient_waveforms' (
'id' INTEGER NOT NULL  PRIMARY KEY AUTOINCREMENT,
'gradient_id' INTEGER NOT NULL ,
'time_point' REAL NOT NULL  DEFAULT 0.0,
'gradient_value' REAL NOT NULL  DEFAULT 0.0,
'f2_value' REAL NOT NULL  DEFAULT 0.0
);

CREATE TABLE 'root_reflect_parameters' (
'id' INTEGER NOT NULL  PRIMARY KEY AUTOINCREMENT,
'a_roots_only' BOOLEAN NOT NULL ,
'graph_angle' REAL NOT NULL  DEFAULT 0.0,
'x_axis_start' REAL NOT NULL  DEFAULT 0.0,
'anorm_real' REAL NOT NULL ,
'anorm_imaginary' REAL NOT NULL ,
'bnorm_real' REAL NOT NULL ,
'bnorm_imaginary' REAL NOT NULL ,
'leading_zeros' INTEGER NOT NULL ,
'trailing_zeros' INTEGER NOT NULL 
);

CREATE TABLE 'a_roots' (
'id' INTEGER NOT NULL  PRIMARY KEY AUTOINCREMENT,
'root_reflect_id' INTEGER NOT NULL  REFERENCES 'root_reflect_parameters' ('id'),
'aroot_real' REAL NOT NULL ,
'aroot_imaginary' REAL NOT NULL ,
'was_flipped' BOOLEAN NOT NULL 
);

CREATE TABLE 'b_roots' (
'id' INTEGER NOT NULL  PRIMARY KEY AUTOINCREMENT,
'root_reflect_id' INTEGER NOT NULL  REFERENCES 'root_reflect_parameters' ('id'),
'broot_real' REAL NOT NULL ,
'broot_imaginary' REAL NOT NULL ,
'was_flipped' BOOLEAN NOT NULL 
);

CREATE TABLE 'experiment_dims' (
'id' INTEGER NOT NULL  PRIMARY KEY AUTOINCREMENT,
'experiment_id' TEXT NOT NULL  REFERENCES 'experiments' ('id'),
'dim1' REAL NOT NULL ,
'dim2' REAL NOT NULL ,
'dim3' REAL NOT NULL ,
UNIQUE (experiment_id, dim1, dim2, dim3)
);

CREATE TABLE 'experiment_metabolites' (
'experiment_id' TEXT NOT NULL ,
'metabolite_id' TEXT NOT NULL ,
PRIMARY KEY (experiment_id, metabolite_id)
);

CREATE TABLE 'gaussian_pulse_parameters' (
'id' INTEGER NOT NULL  PRIMARY KEY AUTOINCREMENT,
'tip_angle' REAL NOT NULL  DEFAULT 0.0,
'time_steps' INTEGER NOT NULL  DEFAULT 0,
'duration' REAL NOT NULL  DEFAULT 0.0,
'bandwidth' REAL NOT NULL  DEFAULT 0.0,
'filter_type' TEXT NOT NULL  DEFAULT '''',
'filter_application' REAL NOT NULL  DEFAULT 0.0
);

CREATE TABLE 'randomized_pulse_parameters' (
'id' INTEGER NOT NULL  PRIMARY KEY AUTOINCREMENT,
'time_steps' INTEGER NOT NULL  DEFAULT 0,
'duration' REAL NOT NULL  DEFAULT 0.0
);

CREATE TABLE 'import_pulse_parameters' (
'id' INTEGER NOT NULL  PRIMARY KEY AUTOINCREMENT,
'file_path' TEXT NOT NULL  DEFAULT '''',
'comment' TEXT DEFAULT NULL,
'file_format' TEXT NOT NULL  DEFAULT '''',
'dwell_time' REAL NOT NULL  DEFAULT 0.0,
'use_max_intensity' BOOLEAN NOT NULL  DEFAULT '0',
'max_intensity' REAL NOT NULL  DEFAULT 0.0,
'scale_factor' REAL NOT NULL  DEFAULT 0.0
);

CREATE TABLE 'pulse_sequence_pulse_projects' (
'pulse_sequence_id' TEXT NOT NULL  REFERENCES 'pulse_sequences' ('id'),
'pulse_project_id' TEXT NOT NULL  REFERENCES 'pulse_projects' ('id'),
'progression' INTEGER NOT NULL 
);

CREATE TABLE 'ocn_states' (
'id' INTEGER NOT NULL  PRIMARY KEY AUTOINCREMENT,
'multiplier' REAL NOT NULL ,
'met_max_iterations' BOOLEAN NOT NULL ,
'met_residual_error' BOOLEAN NOT NULL ,
'met_differential_error' BOOLEAN NOT NULL ,
'met_max_time' BOOLEAN NOT NULL ,
'met_increasing_error' BOOLEAN NOT NULL ,
'run_time' REAL NOT NULL ,
'iterations' INTEGER NOT NULL ,
'decreases' INTEGER NOT NULL 
);

CREATE TABLE 'deltab1_points' (
'ocn_state_id' INTEGER NOT NULL  REFERENCES 'ocn_states' ('id'),
'progression' INTEGER NOT NULL ,
'real_amplitude' REAL NOT NULL ,
'imaginary_amplitude' REAL NOT NULL 
);

CREATE TABLE 'ocn_parameters' (
'id' INTEGER NOT NULL  PRIMARY KEY AUTOINCREMENT,
'pulse_type' TEXT NOT NULL  DEFAULT '''',
'phase_type' TEXT NOT NULL  DEFAULT '''',
'gradient_refocusing_value' REAL NOT NULL  DEFAULT 0.0,
'tip_angle' REAL NOT NULL  DEFAULT 0.0,
'bandwidth' REAL NOT NULL  DEFAULT 0.0,
'step_size_multiplier' REAL NOT NULL ,
'step_size_modification' TEXT NOT NULL  DEFAULT '''',
'excite_band_points' INTEGER NOT NULL  DEFAULT 0,
'b1_immunity_range' REAL NOT NULL  DEFAULT 0.0,
'steps' INTEGER NOT NULL  DEFAULT 0,
'b1_maximum' REAL NOT NULL ,
'limit_sar' BOOLEAN NOT NULL  DEFAULT '0',
'sar_factor' REAL NOT NULL  DEFAULT 0.0,
'error_increase_tolerance' REAL NOT NULL  DEFAULT 0.0,
'max_iteration_check' BOOLEAN NOT NULL  DEFAULT '0',
'max_iterations' INTEGER NOT NULL  DEFAULT 0,
'residual_error_check' BOOLEAN NOT NULL  DEFAULT '0',
'residual_error_tolerance' REAL NOT NULL  DEFAULT 0.0,
'differential_error_check' BOOLEAN NOT NULL  DEFAULT '0',
'differential_error_tolerance' REAL NOT NULL  DEFAULT 0.0,
'halt_if_error_increasing' BOOLEAN NOT NULL  DEFAULT '0',
'halt_on_max_time' BOOLEAN NOT NULL  DEFAULT '0',
'max_time' REAL NOT NULL  DEFAULT 0.0,
'enforce_symmetry' BOOLEAN NOT NULL  DEFAULT '0'
);

CREATE TABLE 'ocn_residual_errors' (
'ocn_state_id' INTEGER NOT NULL  REFERENCES 'ocn_states' ('id'),
'value' REAL NOT NULL ,
'progression' INTEGER NOT NULL 
);

CREATE TABLE 'pulse_designs' (
'id' TEXT NOT NULL  PRIMARY KEY,
'is_public' BOOLEAN NOT NULL ,
'name' TEXT NOT NULL ,
'creator' TEXT DEFAULT NULL,
'created' TIMESTAMP NOT NULL ,
'comment' TEXT DEFAULT NULL,
'calc_resolution' INTEGER NOT NULL  DEFAULT 0,
'pulse_bandwidth_type' TEXT NOT NULL  DEFAULT '''',
'machine_specs_id' INTEGER NOT NULL  REFERENCES 'machine_specs' ('id'),
'gyromagnetic_nuclei' TEXT NOT NULL  DEFAULT '1H'
);

CREATE TABLE 'machine_specs' (
'id' INTEGER NOT NULL  PRIMARY KEY AUTOINCREMENT,
'name' TEXT DEFAULT NULL,
'is_template' BOOLEAN NOT NULL ,
'is_default' BOOLEAN DEFAULT NULL,
'machine_type' TEXT NOT NULL  DEFAULT '''',
'max_b1_field' REAL NOT NULL  DEFAULT 0.0,
'field_strength' REAL NOT NULL  DEFAULT 0.0,
'zero_padding' INTEGER NOT NULL  DEFAULT 0,
'min_dwell_time' REAL NOT NULL  DEFAULT 0.0,
'dwell_time_increment' REAL NOT NULL  DEFAULT 0.0,
'gradient_raster_time' REAL NOT NULL  DEFAULT 0.0,
'gradient_slew_rate' REAL NOT NULL  DEFAULT 0.0,
'gradient_maximum' REAL NOT NULL  DEFAULT 0.0
);

CREATE TABLE 'transforms' (
'id' INTEGER NOT NULL  PRIMARY KEY AUTOINCREMENT,
'pulse_design_id' TEXT NOT NULL  REFERENCES 'pulse_designs' ('id'),
'progression' INTEGER NOT NULL  DEFAULT 0,
'transform_kernel_id' TEXT NOT NULL  REFERENCES 'transform_kernels' ('id'),
'rf_result_id' INTEGER NOT NULL  REFERENCES 'rf_results' ('id'),
UNIQUE (pulse_design_id, progression)
);

CREATE TABLE 'transform_kernels' (
'id' TEXT NOT NULL  PRIMARY KEY,
'type' TEXT DEFAULT NULL,
'name' TEXT DEFAULT NULL,
'menu_label' TEXT DEFAULT NULL,
'is_public' BOOLEAN NOT NULL ,
'created' TIMESTAMP NOT NULL ,
'creator' TEXT DEFAULT NULL,
'comment' TEXT DEFAULT NULL,
'algorithm_code' TEXT DEFAULT NULL,
'time_steps' TEXT DEFAULT NULL,
'duration' TEXT DEFAULT NULL,
'hide_file1' BOOLEAN NOT NULL ,
'hide_file2' BOOLEAN NOT NULL ,
'file1_label' TEXT DEFAULT NULL,
'file2_label' TEXT DEFAULT NULL,
'hide_time_steps' BOOLEAN NOT NULL ,
'hide_duration' BOOLEAN NOT NULL ,
'hide_tip_angle' BOOLEAN NOT NULL ,
'hide_bandwidth' BOOLEAN NOT NULL ,
'tip_angle' TEXT DEFAULT NULL,
'bandwidth' TEXT DEFAULT NULL
);

CREATE TABLE 'transform_parameters' (
'transform_id' INTEGER NOT NULL  REFERENCES 'transforms' ('id'),
'variable' TEXT NOT NULL ,
'type' TEXT NOT NULL ,
'value' TEXT NOT NULL ,
'sort_order' INTEGER NOT NULL 
);

CREATE TABLE 'rf_results' (
'id' INTEGER NOT NULL  PRIMARY KEY AUTOINCREMENT,
'created' TIMESTAMP NOT NULL ,
'rf_waveform' BLOB DEFAULT NULL,
'rf_xaxis' BLOB DEFAULT NULL,
'gradient' BLOB DEFAULT NULL,
'grad_xaxis' BLOB DEFAULT NULL,
'opcon_state_id' INTEGER NOT NULL  DEFAULT 0
);

CREATE TABLE 'transform_kernel_controls' (
'id' INTEGER NOT NULL  PRIMARY KEY AUTOINCREMENT,
'transform_kernel_id' TEXT NOT NULL  REFERENCES 'transform_kernels' ('id'),
'name_' TEXT NOT NULL ,
'type_' TEXT NOT NULL ,
'default_' TEXT NOT NULL ,
'variable_' TEXT NOT NULL ,
'display_order' INTEGER DEFAULT NULL
);

CREATE TABLE 'pulse_sequence_pulse_designs' (
'pulse_sequence_id' TEXT NOT NULL  REFERENCES 'pulse_sequences' ('id'),
'pulse_design_id' TEXT NOT NULL  REFERENCES 'pulse_designs' ('id'),
'progression' INTEGER NOT NULL 
);