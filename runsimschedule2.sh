export SAMPLERATE=96000

python generate_data.py "modelname='wagner_owl'" "stimuli_type='white'" "acoustic_noise_type='none'" "$@"
python generate_data.py "modelname='wagner_owl'" "stimuli_type='coloured'" "acoustic_noise_type='none'" "$@"
python generate_data.py "modelname='wagner_owl'" "stimuli_type='white'" "acoustic_noise_type='independent'" "$@" 
python generate_data.py "modelname='wagner_owl'" "stimuli_type='tone'" "acoustic_noise_type='independent'" "sn_level_range=(0*dB,0*dB)" "$@" 
python generate_data.py "modelname='wagner_owl'" "stimuli_type='bandpassnoise'" "acoustic_noise_type='independent'" "sn_level_range=(0*dB,0*dB)" "$@" 

python make_noise_copy.py

echo Finished!
