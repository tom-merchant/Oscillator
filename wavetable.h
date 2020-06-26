#pragma once

#ifndef OSCILLATOR_WAVETABLE_H
#define OSCILLATOR_WAVETABLE_H

#define OVERSAMPLE 32
#define TABLES_PER_OCTAVE 8
#define NUM_HARMONICS 512

// 1024 harmonics / 3 tables per octave also works well but the way harmonics drop out 8 tables per octave actually
// ends up sounding better over the midrange frequencies

#define MORPH_TABLE_HARMONICS 256
#define MORPH_TABLE_OVERSAMPLE 32
#define MORPH_TABLES_PER_OCTAVE 6
#define TABLES_PER_MORPH 32
#define MAX_WAVEFORMS_PER_MORPH 16


#include <stdint.h>
#include "envelope.h"


extern const char MWT_MAGIC_NUMBER[7];
extern const char MWH_MAGIC_NUMBER[7];

typedef struct wavetable
{
    uint32_t num_samples;
    uint32_t num_harmonics;
    double boundary_harmonic_cosine_amplitude;
    double boundary_harmonic_sine_amplitude;
    double* table;
} wavetable;

typedef struct wavetable_mipmap
{
    uint32_t num_wavetables;
    int *frequency_map;
    wavetable* tables;
} wavetable_mipmap;

typedef struct wavetable_harmonics
{
    uint32_t num_harmonics;
    double* cos_table;
    double* sin_table;
} wavetable_harmonics;


typedef enum morph_type
{
    MORPH_LINEAR,         // What is says on the tin
    MORPH_CONSTANT_POWER, // Angular morph using the identity cos(theta)^2 + sin(theta)^2 = 1 to preserve amplitude
    MORPH_POLYNOMIAL,     // For 3 or more wavetable morphs use an interpolating polynomial
    MORPH_RAMP_UP,        // parabolic curve between each bin with the first as the vertex
    MORPH_RAMP_DOWN,      // parabolic curve between each bin with the second as the vertex
    /* NOTE MORPH ENVELOPES MUST START AT 0 AND END AT 1 AND HAVE A LENGTH OF ONE SECOND */
    MORPH_ENVELOPE,       // shape of morph defined by a user-provided envelope
} morph_type;


typedef struct morph_wavetable_harmonics
{
    uint32_t num_harmonics;
    uint32_t num_wavetables;
    uint32_t tables_per_morph;
    uint32_t oversample;
    uint32_t tables_per_octave;
    morph_type morph_style;
    envelope* morph_env;
    wavetable_harmonics* spectrums;
} morph_wavetable_harmonics;


typedef struct morph_wavetable_mipmap
{
    uint32_t num_wavetables;
    uint32_t tables_per_morph;
    uint32_t total_tables;
    wavetable_mipmap* tables;
} morph_wavetable_mipmap;


void generate_wavetable_freqmap(wavetable_mipmap* wt, long samplerate);

void generate_morph_wavetable_map ( morph_wavetable_mipmap* mwm );

wavetable_mipmap* generate_wavetables(wavetable_harmonics* harmonics, int oversample, int tables_per_octave);
wavetable_mipmap* generate_static_wavetables(wavetable_harmonics* harmonics, int oversample, int tables_per_octave);

wavetable_mipmap* generate_wavetables_from_signal ( const double* signal, int num_samples, int num_harmonics, int oversample, int tables_per_octave );

morph_wavetable_mipmap* generate_morph_wavetables( morph_wavetable_harmonics* harmonics, long samplerate );


void save_wavetable(wavetable_mipmap* wt, const char* path);

wavetable_mipmap* load_wavetable(const char* path, int samplerate);

morph_wavetable_mipmap* load_serum_morph_wavetable ( char* path, int tables_per_octave, int oversample,
                                                     int tables_per_morph, long samplerate );

wavetable_mipmap* load_serum_wavetable (char* path, int tables_per_octave, int oversample, long samplerate);

void save_morph_wavetable_mipmap                          ( morph_wavetable_mipmap *mwm,          const char* path );
void save_morph_wavetable_harmonics                       ( const morph_wavetable_harmonics *mwh, const char* path );
morph_wavetable_mipmap* load_morph_wavetable_mipmap       (                                       const char* path );
morph_wavetable_harmonics* load_morph_wavetable_harmonics (                                       const char* path );

void free_wavetable_mipmap(wavetable_mipmap* wt);
void free_morph_wavetable_mipmap(morph_wavetable_mipmap* mwm);

#endif //OSCILLATOR_WAVETABLE_H
