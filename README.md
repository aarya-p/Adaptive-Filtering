# Wind Noise Cancellation using Hybrid Adaptive Filtering 

This repository contains the MATLAB implementation and documentation for our mini-project titled **"Wind Noise Cancellation using Hybrid Adaptive Filtering"**, developed as part of the 6th semester Electronics and Communication Engineering curriculum at **Ramaiah Institute of Technology**, Bangalore.

## Project Overview

Wind noise poses a major challenge in real-time audio systems, reducing signal quality and communication clarity. Traditional fixed filters are not effective for dynamic noise scenarios.

To address this, we developed a **hybrid adaptive filtering system** that combines:
- **Least Mean Squares (LMS)**: low-complexity, efficient under stable noise
- **Recursive Least Squares (RLS)**: fast convergence, better under rapidly changing noise

The hybrid model dynamically switches between LMS and RLS based on real-time SNR estimation, leveraging the best of both approaches.

## Methodology

1. **Data Preparation**:
   - Clean signals: Piano tones (octave 2–4)
   - Noise signals: Simulated and recorded wind noise (Weibull-distributed)

2. **Algorithm Design**:
   - LMS and RLS filters run in parallel
   - Real-time SNR estimated over sliding window
   - Dynamic switching based on SNR difference threshold

3. **Evaluation Metrics**:
   - SNR Gain
   - Mean Squared Error (MSE)
   - Convergence rate
   - Visual waveform and frequency analysis

## Technologies Used

- MATLAB with DSP Toolbox
- .WAV audio samples
- Signal processing techniques
- Adaptive filtering algorithms

## Results Summary

|   Filter   |  SNR Gain  |
|------------|------------|
|    LMS     |   ~24 dB   |
|    RLS     |   ~49 dB   |
| **Hybrid** | **~52 dB** |

> The hybrid filter consistently outperformed LMS and RLS individually, with better SNR, faster convergence, and efficient computation.

## 🧩 Future Work

- Adaptive tuning of LMS/RLS parameters
- Hardware implementation on FPGA/DSP
- Application to real-world audio (hearing aids, teleconferencing, etc.)
