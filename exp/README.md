# Experimental control of somatosensory detection task

Control of DS5 Bipolar Constant Current Stimulator via MATLAB Data Acquisition Toolbox and Psychtoolbox

## General

- Stimulation intensity relates to adjusted voltage and current (e.g., 10V:10mA -> signal of 5 relates to 5 mA)
- DS5 current and voltage output have to multiplied by 10 (1V signal equals 10mA resp. 10V)
- If input positive, then red output is positive, so current flows from black to red
- IMPORTANT: Wait until output is done "talking" (e.g., trigger(ao); wait(ao,1000);)

## Hardware Requirements

- Digitimer DS5
- National Instruments NI USB-6343
  - Manual: https://www.ni.com/en-us/support/model.usb-6343.html
- Parallel port for button box

## Software Requirements

- Windows 7
- MATLAB 32-bit [8.5.1 (R2015a) Service Pack 1]
- MATLAB Toolboxes:
  - Data Acquisition Toolbox 32-bit 3.7 (http://de.mathworks.com/products/daq/)
  - Psychtoolbox 3.0.11 (http://psychtoolbox.org)
  - Palamedes Toolbox 1.9.0 (http://palamedestoolbox.org) [included in the repository]
  - Mex-file plug-in for fast 32-bit MATLAB port I/O access on 64-bit Windows (http://apps.usd.edu/coglab/psyc770/IO32on64.html)

## Repository Structure

- respirationCA.m - main script to control experiment by executing sections step-by-step
- exp_init_NI.m - helper function to initialize experiment
- respirationCA/ - directory with code for experimental block
- thr1F/ - directory with code for threshold assessment
- assets/func - general functions used in multiple scripts
- assets/plot - functions to plot threshold assessment and analog input data
- assets/vendor - third-party code, e.g. Palamedes & io32
