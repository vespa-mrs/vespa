CREATE TABLE 'pulse_designs' (
'id' TEXT NOT NULL  PRIMARY KEY,
'is_public' BOOLEAN NOT NULL ,
'name' TEXT NOT NULL ,
'creator' TEXT DEFAULT NULL,
'created' TIMESTAMP NOT NULL ,
'comment' TEXT DEFAULT NULL,
'calc_resolution' INTEGER NOT NULL  DEFAULT 0,
'pulse_bandwidth_type' TEXT NOT NULL  DEFAULT '''',
'machine_specs_id' INTEGER NOT NULL  REFERENCES 'machine_specs' ('id')
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
'opcon_state_id' INTEGER NOT NULL  DEFAULT 0 REFERENCES 'opcon_states' ('id')
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



