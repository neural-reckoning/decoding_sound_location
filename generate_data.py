from analyse_model import *

if __name__=='__main__':
    # defaults
    cpu = -2
    report = 'stderr'
    forced_size = 6400
    stimuli_type = 'white'
    stimuli_duration = 100*ms
    acoustic_noise_type = 'none'
    # overrides
    for arg in sys.argv[1:]:
        exec arg
    # load the model
    _cpu = cpu
    exec 'from models.'+modelname+' import *'
    cpu = _cpu

    ### STIMULI
    if stimuli_type=='white':
        stimuli_filters = ['type=="whitenoise"']
        stimuli = WhiteNoiseStimuli(space, duration=stimuli_duration,
                                    level=level)
    elif stimuli_type=='coloured':
        stimuli_filters = ['type=="powerlawnoise"']
        stimuli = PowerLawNoiseStimuli(space, duration=stimuli_duration,
                                       level=level)
    elif stimuli_type=='natural':
        stimuli_filters = ['type=="natural"']
        raise ValueError('Natural sounds not implemented yet (cutoff initial segment)')
        stimuli = NatureSoundsStimuli(space, duration=stimuli_duration,
                                      level=level, noiselevel=0*dB)
    elif stimuli_type=='tone':
        stimuli_filters = ['type=="tone"']
        stimuli = ToneStimuli(space, sort(cf), duration=stimuli_duration,
                             level=level)
    elif stimuli_type=='bandpassnoise':
        stimuli_filters = ['type=="bandpassnoise"']
        stimuli = BandPassNoiseStimuli(space, sort(cf),
                                       duration=stimuli_duration,
                                       level=level)
    else:
        raise ValueError('Stimuli type not known')

    ### ACOUSTIC NOISE
    if acoustic_noise_type=='none':
        acousticnoisemodel = NoAcousticNoise()
    elif acoustic_noise_type=='independent':
        acousticnoisemodel = IndependentWhiteAcousticNoise(tuple(level-sn for sn in sn_level_range[::-1]))
    else:
        raise ValueError('acoustic_noise_type unknown')
    extraname['acousticnoisemodel'] = acousticnoisemodel
    
    # Stupid stuff to make Analysis() work even though we don't use it
    training_filters = testing_filters = sample_filters = []
    use_ideal_responses = False
    num_shuffles = 1
    training_size = testing_size = extra_size = 0
    
    extended_name = generate_extended_datamanager_name(basename, extraname)    
    analysis = get_analysis_from_namespace()

    logstr = str(datetime.datetime.now())+':\n'
    logstr += '  Starting '+' '.join(sys.argv[1:])+'\n\n'
    open('generate_data.log', 'a').write(logstr)
    
    # Check how many we have already done
    numdone = len(analysis.response_subset(analysis.responses(), stimuli_filters)[0])
    print 'Already have', numdone, 'matching responses.'
    forced_size -= numdone
    if forced_size<0:
        status = 'Skipped (already completed)'
    else:        
        if analysis.generate_data(stimuli, cpu, extra_size, forced_size, doexit=False):
            status = 'Completed'
        else:
            status = 'Aborted'
    
    logstr = str(datetime.datetime.now())+':\n'
    logstr += '  '+status+'\n\n'
    open('generate_data.log', 'a').write(logstr)
