# Add "$@" to the end of each generate_data.py call, and then test the script
# before running by calling runsimschedule.bat forced_size=1 to generate a
# a single value for each export of parameters before doing the full run

# NOTE: you need to make a copy of data without acoustic noise into the
# directories with acoustic noise, run make_noise_copy.py after this script
# completes


export SAMPLERATE=44100

python generate_data.py "modelname='joris_cat'" "stimuli_type='white'" "acoustic_noise_type='none'" "$@"
python generate_data.py "modelname='joris_cat'" "stimuli_type='coloured'" "acoustic_noise_type='none'" "$@"
python generate_data.py "modelname='joris_cat'" "stimuli_type='white'" "acoustic_noise_type='independent'" "$@"
python generate_data.py "modelname='joris_cat'" "stimuli_type='tone'" "acoustic_noise_type='independent'" "sn_level_range=(0*dB,0*dB)" "$@"
python generate_data.py "modelname='joris_cat'" "stimuli_type='bandpassnoise'" "acoustic_noise_type='independent'" "sn_level_range=(0*dB,0*dB)" "$@" 

python generate_data.py "modelname='mcalpine_guinea_pig'" "stimuli_type='white'" "acoustic_noise_type='none'" "$@"
python generate_data.py "modelname='mcalpine_guinea_pig'" "stimuli_type='coloured'" "acoustic_noise_type='none'" "$@"
python generate_data.py "modelname='mcalpine_guinea_pig'" "stimuli_type='white'" "acoustic_noise_type='independent'" "$@"
python generate_data.py "modelname='mcalpine_guinea_pig'" "stimuli_type='tone'" "acoustic_noise_type='independent'" "sn_level_range=(0*dB,0*dB)" "$@"
python generate_data.py "modelname='mcalpine_guinea_pig'" "stimuli_type='bandpassnoise'" "acoustic_noise_type='independent'" "sn_level_range=(0*dB,0*dB)" "$@" 

python generate_data.py "modelname='ircam_human_uniform_bipd'" "stimuli_type='white'" "acoustic_noise_type='none'" "$@"
python generate_data.py "modelname='ircam_human_uniform_bipd'" "stimuli_type='coloured'" "acoustic_noise_type='none'" "$@"
python generate_data.py "modelname='ircam_human_uniform_bipd'" "stimuli_type='white'" "acoustic_noise_type='independent'" "$@"
python generate_data.py "modelname='ircam_human_uniform_bipd'" "stimuli_type='tone'" "acoustic_noise_type='independent'" "sn_level_range=(0*dB,0*dB)" "$@"
python generate_data.py "modelname='ircam_human_uniform_bipd'" "stimuli_type='bandpassnoise'" "acoustic_noise_type='independent'" "sn_level_range=(0*dB,0*dB)" "$@"

python generate_data.py "modelname='ircam_human_mcalpine_bd'" "stimuli_type='white'" "acoustic_noise_type='none'" "$@"
python generate_data.py "modelname='ircam_human_mcalpine_bd'" "stimuli_type='coloured'" "acoustic_noise_type='none'" "$@"
python generate_data.py "modelname='ircam_human_mcalpine_bd'" "stimuli_type='white'" "acoustic_noise_type='independent'" "$@"
python generate_data.py "modelname='ircam_human_mcalpine_bd'" "stimuli_type='tone'" "acoustic_noise_type='independent'" "sn_level_range=(0*dB,0*dB)" "$@"
python generate_data.py "modelname='ircam_human_mcalpine_bd'" "stimuli_type='bandpassnoise'" "acoustic_noise_type='independent'" "sn_level_range=(0*dB,0*dB)" "$@"

python generate_data.py "modelname='ircam_human_tollin_bd'" "stimuli_type='white'" "acoustic_noise_type='none'" "$@"
python generate_data.py "modelname='ircam_human_tollin_bd'" "stimuli_type='coloured'" "acoustic_noise_type='none'" "$@"
python generate_data.py "modelname='ircam_human_tollin_bd'" "stimuli_type='white'" "acoustic_noise_type='independent'" "$@"
python generate_data.py "modelname='ircam_human_tollin_bd'" "stimuli_type='tone'" "acoustic_noise_type='independent'" "sn_level_range=(0*dB,0*dB)" "$@"
python generate_data.py "modelname='ircam_human_tollin_bd'" "stimuli_type='bandpassnoise'" "acoustic_noise_type='independent'" "sn_level_range=(0*dB,0*dB)" "$@"

export SAMPLERATE=192000

python generate_data.py "modelname='ircam_mcalpine_guinea_pig'" "stimuli_type='white'" "acoustic_noise_type='none'" "$@"
python generate_data.py "modelname='ircam_mcalpine_guinea_pig'" "stimuli_type='coloured'" "acoustic_noise_type='none'" "$@"
python generate_data.py "modelname='ircam_mcalpine_guinea_pig'" "stimuli_type='white'" "acoustic_noise_type='independent'" "$@" 
python generate_data.py "modelname='ircam_mcalpine_guinea_pig'" "stimuli_type='tone'" "acoustic_noise_type='independent'" "sn_level_range=(0*dB,0*dB)" "$@" 
python generate_data.py "modelname='ircam_mcalpine_guinea_pig'" "stimuli_type='bandpassnoise'" "acoustic_noise_type='independent'" "sn_level_range=(0*dB,0*dB)" "$@" 

export SAMPLERATE=97656

python generate_data.py "modelname='joris_tollin_cat'" "stimuli_type='white'" "acoustic_noise_type='none'" "$@"
python generate_data.py "modelname='joris_tollin_cat'" "stimuli_type='coloured'" "acoustic_noise_type='none'" "$@"
python generate_data.py "modelname='joris_tollin_cat'" "stimuli_type='white'" "acoustic_noise_type='independent'" "$@" 
python generate_data.py "modelname='joris_tollin_cat'" "stimuli_type='tone'" "acoustic_noise_type='independent'" "sn_level_range=(0*dB,0*dB)" "$@" 
python generate_data.py "modelname='joris_tollin_cat'" "stimuli_type='bandpassnoise'" "acoustic_noise_type='independent'" "sn_level_range=(0*dB,0*dB)" "$@" 

export SAMPLERATE=96000

python generate_data.py "modelname='wagner_owl'" "stimuli_type='white'" "acoustic_noise_type='none'" "$@"
python generate_data.py "modelname='wagner_owl'" "stimuli_type='coloured'" "acoustic_noise_type='none'" "$@"
python generate_data.py "modelname='wagner_owl'" "stimuli_type='white'" "acoustic_noise_type='independent'" "$@" 
python generate_data.py "modelname='wagner_owl'" "stimuli_type='tone'" "acoustic_noise_type='independent'" "sn_level_range=(0*dB,0*dB)" "$@" 
python generate_data.py "modelname='wagner_owl'" "stimuli_type='bandpassnoise'" "acoustic_noise_type='independent'" "sn_level_range=(0*dB,0*dB)" "$@" 

export SAMPLERATE=44100

# *** Heterogeneity testing, no acoustic noise ***

export STD_FACTOR=0.0
python generate_data.py "modelname='mcalpine_guinea_pig'" "stimuli_type='white'" "acoustic_noise_type='none'" "$@" 
export STD_FACTOR=0.05
python generate_data.py "modelname='mcalpine_guinea_pig'" "stimuli_type='white'" "acoustic_noise_type='none'" "$@"
export STD_FACTOR=0.15
python generate_data.py "modelname='mcalpine_guinea_pig'" "stimuli_type='white'" "acoustic_noise_type='none'" "$@"
export STD_FACTOR=0.25
python generate_data.py "modelname='mcalpine_guinea_pig'" "stimuli_type='white'" "acoustic_noise_type='none'" "$@"
export STD_FACTOR=0.5
python generate_data.py "modelname='mcalpine_guinea_pig'" "stimuli_type='white'" "acoustic_noise_type='none'" "$@"
export STD_FACTOR=0.75
python generate_data.py "modelname='mcalpine_guinea_pig'" "stimuli_type='white'" "acoustic_noise_type='none'" "$@"
export STD_FACTOR=1.5
python generate_data.py "modelname='mcalpine_guinea_pig'" "stimuli_type='white'" "acoustic_noise_type='none'" "$@"

# *** Heterogeneity testing, with acoustic noise ***

export STD_FACTOR=0.0
python generate_data.py "modelname='mcalpine_guinea_pig'" "stimuli_type='white'" "acoustic_noise_type='independent'" "$@" 
export STD_FACTOR=0.05
python generate_data.py "modelname='mcalpine_guinea_pig'" "stimuli_type='white'" "acoustic_noise_type='independent'" "$@"
export STD_FACTOR=0.15
python generate_data.py "modelname='mcalpine_guinea_pig'" "stimuli_type='white'" "acoustic_noise_type='independent'" "$@"
export STD_FACTOR=0.25
python generate_data.py "modelname='mcalpine_guinea_pig'" "stimuli_type='white'" "acoustic_noise_type='independent'" "$@"
export STD_FACTOR=0.5
python generate_data.py "modelname='mcalpine_guinea_pig'" "stimuli_type='white'" "acoustic_noise_type='independent'" "$@"
export STD_FACTOR=0.75
python generate_data.py "modelname='mcalpine_guinea_pig'" "stimuli_type='white'" "acoustic_noise_type='independent'" "$@"
export STD_FACTOR=1.5
python generate_data.py "modelname='mcalpine_guinea_pig'" "stimuli_type='white'" "acoustic_noise_type='independent'" "$@"
export STD_FACTOR=0.8
python generate_data.py "modelname='mcalpine_guinea_pig'" "stimuli_type='white'" "acoustic_noise_type='independent'" "$@"
export STD_FACTOR=0.9
python generate_data.py "modelname='mcalpine_guinea_pig'" "stimuli_type='white'" "acoustic_noise_type='independent'" "$@"
export STD_FACTOR=1.1
python generate_data.py "modelname='mcalpine_guinea_pig'" "stimuli_type='white'" "acoustic_noise_type='independent'" "$@"
export STD_FACTOR=1.2
python generate_data.py "modelname='mcalpine_guinea_pig'" "stimuli_type='white'" "acoustic_noise_type='independent'" "$@"
export STD_FACTOR=1.3
python generate_data.py "modelname='mcalpine_guinea_pig'" "stimuli_type='white'" "acoustic_noise_type='independent'" "$@"

export STD_FACTOR=

python make_noise_copy.py

echo Finished!
