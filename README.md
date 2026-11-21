# Decoding Two TDMA Messages in MATLAB

How does noise affect whether a receiver can recover a message? This project decodes two French text messages from a bundled time-division multiple access (TDMA) recording. It then adds white Gaussian noise at four signal-to-noise ratios (SNRs) and decodes the same messages again.

## Run the experiment

You’ll need MATLAB with **Signal Processing Toolbox**. Download or clone this repository, open its root folder in MATLAB, and run the following in the Command Window.

```matlab
rng(0, 'twister');
projet_TDMA
```

The script prints both decoded messages before adding noise and at 20, 10, 0, and −5 dB. The seed makes the noise reproducible for the example below.

## Example result

In MATLAB R2025b, both messages matched the baseline at 20 and 10 dB. At 0 dB, each had two changed characters. At −5 dB, both had more errors.

| Added noise (SNR) | User 1 changed characters | User 2 changed characters |
| --- | ---: | ---: |
| 20 dB | 0 | 0 |
| 10 dB | 0 | 0 |
| 0 dB | 2 | 2 |
| −5 dB | 9 | 7 |

Each message contains 60 characters. The table counts positions that differ from the same message decoded before adding noise. These results are from one seeded run.

## How the receiver works

1. Extract slot 2 for user 1 (15 kHz carrier) and slot 5 for user 2 (45 kHz carrier).
2. Mix each slot down to baseband and apply a 10 kHz low-pass filter with zero-phase filtering.
3. Use matched filtering, symbol sampling, and descrambling to recover bits, then convert them to text.

## Files

| File | Role |
| --- | --- |
| [`projet_TDMA.m`](projet_TDMA.m) | Receiver and noise experiment |
| [`Signal_MF_TDMA.mat`](Signal_MF_TDMA.mat) | Bundled received signal |
| [`Demod_BdB.m`](Demod_BdB.m) | Baseband demodulation and decoding, credited in the source to Nathalie Thomas |
| [`bin2str.m`](bin2str.m) | Bit-to-text conversion, credited in the source to Rick Rosson |
