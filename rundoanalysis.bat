@echo off
start /min do_analysis.py analysis_type='cells' %*
start /min do_analysis.py analysis_type='frequency' %*
start /min do_analysis.py analysis_type='noise' %*
start /min do_analysis.py analysis_type='colour' %*
