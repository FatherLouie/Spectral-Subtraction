# Spectral Subtraction
The aim of this project was to explore the possibilty of speaker separation (from a recording with multiple) using noise cancellation techniques that were presented in the paper here:
[1979, SF Boll](https://ieeexplore.ieee.org/document/1163209). The outputs obtained showed excellent clarity of the separated speakers. 

## Pitch Detection
One of the files uploaded tested 2 custom functions against 5 in-built functions for detection of pitch of a speech signal. This was crucial to perform spectral subtraction, and hence a thorough review of all
available methods had to be conducted. Plots on running the code for pitch detection compare the frame-wise pitches for the input signal.

## Dual Source
Once the frame-wise pitch of the speaker is known, the act of spectral subtraction is done in the file labelled [dual_source]()

## Analysis and Synthesis
Another approach using filter banks was also tried out, where not adding certain filter banks in order to recover certain speakers was implemented. The results for these too can be found in the file labellled
[analysis_and_synthesis]()

## Transfer Functions
In order to aid the filter bank creation, the theory of transfer functions was learnt and relevant properties of transfer functions were coded with plots confirming the same.
