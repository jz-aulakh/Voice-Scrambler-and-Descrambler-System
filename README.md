# Voice-Scrambler-and-Descrambler-System


Voice scrambling system using frequency inversion technique in MATLAB.

## Description

Scrambles audio by inverting its frequency spectrum using carrier modulation and Butterworth filtering. Applying the process twice recovers the original signal.

## Requirements

- MATLAB (R2016b or later)
- Signal Processing Toolbox

## Usage

```matlab
voice_scrambler_and_descrambler
```

Select a .wav file when prompted. The system will:
1. Play original audio
2. Scramble and play unintelligible version
3. Descramble and play recovered audio
4. Generate visualization plots

## System Parameters

- **Carrier Frequency:** 3700 Hz
- **Filter:** 8th-order Butterworth Low-Pass
- **Cutoff Frequency:** 3000 Hz

## How It Works

```
Original Audio → Modulation → Low-Pass Filter → Scrambled Audio
Scrambled Audio → Modulation → Low-Pass Filter → Recovered Audio
```

The same process is applied twice to recover the original signal.

## Output

- `scrambled_audio.wav` - Unintelligible audio
- `descrambled_audio.wav` - Recovered audio
- Time-domain, frequency-domain, and spectrogram plots

## Team

Rana Muhammad Bilal, Muhammad Jahan Zaib, Umer Nisar  
EE-232: Signals and Systems | NUST Islamabad
