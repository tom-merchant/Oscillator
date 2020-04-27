
#pragma once

#ifndef OSCILLATOR_OSCILLATOR_H
#define OSCILLATOR_OSCILLATOR_H

#include "wavetable.h"

typedef double (signal)( double, double, double );


typedef enum oscillator_type
{
    FUNCTION,
    WAVETABLE,
    PULSEWIDTH,
    MORPH_WAVETABLE
} oscillator_type;


typedef struct oscillator
{
    double init_phase;
    double phase;
    double frequency;
    double pulse_width;
    double sample_radians_per_hertz; /* e.g M_PI * 2 / 44100 */
    oscillator_type type;
    signal *source;
    wavetable_mipmap *wavetable;
    morph_wavetable_mipmap *morph_table;
    double morph_pos;
    unsigned int oversample;
} oscillator;


double next_sample ( oscillator* osc );

double next_n_samples ( oscillator* osc, long long nsamples );

double reset_oscillator ( oscillator* osc );

signal sine;
signal unlimited_square;
signal unlimited_ramp;
signal unlimited_inverse_ramp;
signal unlimited_triangle;

oscillator* new_oscillator ( double freq, double phase, double pulsewidth, long samplerate, signal source );

oscillator* new_sine_oscillator               ( double freq, double phase, long samplerate );

oscillator* new_unlimited_square_oscillator   ( double freq, double phase, double pw, long samplerate );

oscillator* new_unlimited_sawtooth_oscillator ( double freq, double phase, long samplerate );

oscillator* new_unlimited_triangle_oscillator ( double freq, double phase, long samplerate );

oscillator* new_wavetable_oscillator          ( double freq, double phase, long samplerate, wavetable_mipmap* wt );

oscillator* new_morph_wavetable_oscillator    ( double freq, double phase, long samplerate, morph_wavetable_mipmap* mwm );

oscillator* new_pulsewidth_oscillator         ( double freq, double phase, wavetable_mipmap* saw_mipmap );

#endif //OSCILLATOR_OSCILLATOR_H