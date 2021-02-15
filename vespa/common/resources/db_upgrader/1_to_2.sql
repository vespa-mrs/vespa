CREATE TABLE pulse_projects (
id TEXT NOT NULL  PRIMARY KEY,
is_public BOOLEAN NOT NULL ,
name TEXT NOT NULL ,
creator TEXT DEFAULT NULL,
created TIMESTAMP NOT NULL ,
comment TEXT DEFAULT NULL,
machine_settings_id INTEGER NOT NULL  REFERENCES machine_settings (id),
master_parameters_id INTEGER NOT NULL  REFERENCES master_parameters (id)
);

CREATE TABLE machine_settings (
id INTEGER NOT NULL  PRIMARY KEY AUTOINCREMENT,
name TEXT DEFAULT NULL,
is_template BOOLEAN NOT NULL ,
is_default BOOLEAN DEFAULT NULL,
machine_type TEXT NOT NULL  DEFAULT '''',
max_b1_field REAL NOT NULL  DEFAULT 0.0,
field_strength REAL NOT NULL  DEFAULT 0.0,
zero_padding INTEGER NOT NULL  DEFAULT 0,
min_dwell_time REAL NOT NULL  DEFAULT 0.0,
dwell_time_increment REAL NOT NULL  DEFAULT 0.0,
gradient_raster_time REAL NOT NULL  DEFAULT 0.0,
gradient_slew_rate REAL NOT NULL  DEFAULT 0.0,
gradient_maximum REAL NOT NULL  DEFAULT 0.0
);

CREATE TABLE master_parameters (
id INTEGER NOT NULL  PRIMARY KEY AUTOINCREMENT,
calc_resolution INTEGER NOT NULL  DEFAULT 0,
pulse_bandwidth_type TEXT NOT NULL  DEFAULT ''''
);

CREATE TABLE transformations (
id INTEGER NOT NULL  PRIMARY KEY AUTOINCREMENT,
pulse_project_id TEXT NOT NULL  REFERENCES pulse_projects (id),
progression INTEGER NOT NULL  DEFAULT 0,
transformation_type TEXT NOT NULL  DEFAULT '''',
parameters_id INTEGER NOT NULL  REFERENCES interpolate_rescale_parameters (id) REFERENCES hs_pulse_parameters (id) REFERENCES slr_pulse_parameters (id) REFERENCES root_reflect_parameters (id),
result_id INTEGER NOT NULL  REFERENCES results (id)
);

CREATE TABLE slr_pulse_parameters (
id INTEGER NOT NULL  PRIMARY KEY AUTOINCREMENT,
tip_angle REAL NOT NULL  DEFAULT 0.0,
time_steps INTEGER NOT NULL  DEFAULT 0,
duration REAL NOT NULL  DEFAULT 0.0,
bandwidth REAL NOT NULL  DEFAULT 0.0,
separation REAL NOT NULL  DEFAULT 0.0,
is_single_band BOOLEAN NOT NULL ,
nc_phase_subtype TEXT NOT NULL  DEFAULT '''',
slr_filter_type TEXT NOT NULL  DEFAULT '''',
pass_ripple REAL NOT NULL  DEFAULT 0.0,
reject_ripple REAL NOT NULL  DEFAULT 0.0
);

CREATE TABLE interpolate_rescale_parameters (
id INTEGER NOT NULL  PRIMARY KEY AUTOINCREMENT,
do_interpolate BOOLEAN NOT NULL ,
interpolation_factor INTEGER NOT NULL  DEFAULT 0,
new_dwell_time REAL NOT NULL  DEFAULT 0.0,
do_rescaling BOOLEAN NOT NULL ,
angle REAL NOT NULL  DEFAULT 0.0
);

CREATE TABLE results (
id INTEGER NOT NULL  PRIMARY KEY AUTOINCREMENT,
created TIMESTAMP NOT NULL ,
gradient_id INTEGER NOT NULL  REFERENCES gradients (id)
);

CREATE TABLE hs_pulse_parameters (
id INTEGER NOT NULL  PRIMARY KEY AUTOINCREMENT,
total_rotation REAL NOT NULL  DEFAULT 0.0,
time_steps INTEGER NOT NULL  DEFAULT 0,
dwell_time REAL DEFAULT 0.0,
was_bandwidth_specified BOOLEAN NOT NULL ,
quality_cycles REAL NOT NULL  DEFAULT 0.0,
power_n INTEGER NOT NULL  DEFAULT 0,
sharpness_mu REAL NOT NULL  DEFAULT 0.0,
filter_type TEXT NOT NULL  DEFAULT '''',
filter_application REAL NOT NULL  DEFAULT 0.0
);

CREATE TABLE rf_waveforms (
id INTEGER NOT NULL  PRIMARY KEY AUTOINCREMENT,
result_id INTEGER NOT NULL  REFERENCES results (id),
time_point REAL NOT NULL  DEFAULT 0,
real_amplitude REAL NOT NULL  DEFAULT 0.0,
imaginary_amplitude REAL NOT NULL  DEFAULT 0.0
);

CREATE TABLE gradients (
id INTEGER NOT NULL  PRIMARY KEY AUTOINCREMENT REFERENCES gradient_waveforms (gradient_id),
linear_gradient_value REAL NOT NULL  DEFAULT 0.0,
refocused_gradient REAL NOT NULL  DEFAULT 0.0,
frequency_offset REAL NOT NULL  DEFAULT 0.0
);

CREATE TABLE gradient_waveforms (
id INTEGER NOT NULL  PRIMARY KEY AUTOINCREMENT,
gradient_id INTEGER NOT NULL ,
time_point REAL NOT NULL  DEFAULT 0.0,
gradient_value REAL NOT NULL  DEFAULT 0.0,
f2_value REAL NOT NULL  DEFAULT 0.0
);

CREATE TABLE root_reflect_parameters (
id INTEGER NOT NULL  PRIMARY KEY AUTOINCREMENT,
a_roots_only BOOLEAN NOT NULL ,
graph_angle REAL NOT NULL  DEFAULT 0.0,
x_axis_start REAL NOT NULL  DEFAULT 0.0,
anorm_real REAL NOT NULL ,
anorm_imaginary REAL NOT NULL ,
bnorm_real REAL NOT NULL ,
bnorm_imaginary REAL NOT NULL ,
leading_zeros INTEGER NOT NULL ,
trailing_zeros INTEGER NOT NULL 
);

CREATE TABLE a_roots (
id INTEGER NOT NULL  PRIMARY KEY AUTOINCREMENT,
root_reflect_id INTEGER NOT NULL  REFERENCES root_reflect_parameters (id),
aroot_real REAL NOT NULL ,
aroot_imaginary REAL NOT NULL ,
was_flipped BOOLEAN NOT NULL 
);

CREATE TABLE b_roots (
id INTEGER NOT NULL  PRIMARY KEY AUTOINCREMENT,
root_reflect_id INTEGER NOT NULL  REFERENCES root_reflect_parameters (id),
broot_real REAL NOT NULL ,
broot_imaginary REAL NOT NULL ,
was_flipped BOOLEAN NOT NULL 
);


