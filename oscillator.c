#include "oscillator.h"

#include <math.h>
#include <stdlib.h>

#define M_2PI ( 2 * M_PI )

long _samplerate;

double sine ( double phase, double frequency, double x )
{
    return sin ( phase );
}

double unlimited_square ( double phase, double frequency, double pulsewidth )
{
    return copysign (1, phase - (M_PI + (0.5 - pulsewidth) * M_2PI) );
}

double unlimited_ramp ( double phase, double frequency, double x )
{
    return ( phase - M_PI ) / M_PI;
}

double unlimited_inverse_ramp ( double phase, double frequency, double x )
{
    return ( M_PI - phase ) / M_PI;
}

double unlimited_triangle ( double phase, double frequency, double x )
{
    return 2 * fabs ( unlimited_ramp ( phase, frequency, x ) ) - 1;
}


oscillator* new_oscillator ( double freq, double phase, double pulsewidth, long samplerate, signal source )
{
    oscillator* new_osc = malloc ( sizeof ( oscillator ) );

    new_osc->sample_radians_per_hertz = M_2PI / (double)samplerate;

    new_osc->init_phase  = phase;
    new_osc->phase       = phase;
    new_osc->frequency   = freq;
    new_osc->pulse_width = pulsewidth;
    new_osc->source      = source;
    new_osc->type        = FUNCTION;

    return new_osc;
}


oscillator* new_sine_oscillator ( double freq, double phase, long samplerate )
{
    return new_oscillator ( freq, phase, 0, samplerate, sine );
}

oscillator* new_unlimited_square_oscillator ( double freq, double phase, double pulsewidth, long samplerate )
{
    return new_oscillator ( freq, phase, pulsewidth, samplerate, unlimited_square );
}

oscillator* new_unlimited_sawtooth_oscillator ( double freq, double phase, long samplerate )
{
    return new_oscillator ( freq, phase, 0, samplerate, unlimited_ramp );
}

oscillator* new_unlimited_triangle_oscillator ( double freq, double phase, long samplerate )
{
    return new_oscillator ( freq, phase, 0, samplerate, unlimited_triangle );
}

oscillator* new_morph_wavetable_oscillator ( double freq, double phase, long samplerate, morph_wavetable_mipmap* mwm )
{
    oscillator* out = new_oscillator ( freq, phase, 50, samplerate, sine );

    out->type = MORPH_WAVETABLE;
    out->morph_table = mwm;

    return out;
}

oscillator* new_wavetable_oscillator ( double freq, double phase, long samplerate, wavetable_mipmap* wt )
{
    oscillator* out = new_oscillator ( freq, phase, 50, samplerate, sine );

    out->type = WAVETABLE;
    out->wavetable = wt;

    return out;
}

oscillator* new_wavetable_oscillator_from_harmonics( double freq, double phase, long samplerate, wavetable_harmonics* harmonics)
{
    oscillator* out = new_oscillator ( freq, phase, 50, samplerate, sine );

    out->type = WAVETABLE;
    out->wavetable = generate_wavetables( harmonics, OVERSAMPLE, TABLES_PER_OCTAVE );

    return out;
}

oscillator* new_wavetable_oscillator_from_signal( double freq, double phase, long samplerate, double* signal, int num_samples )
{
    oscillator* out = new_oscillator ( freq, phase, 50, samplerate, sine );

    out->type = WAVETABLE;
    out->wavetable = generate_wavetables_from_signal ( signal, num_samples, NUM_HARMONICS, OVERSAMPLE, TABLES_PER_OCTAVE );

    return out;
}

double wavetable_signal (oscillator* osc)
{
    int table_index;
    wavetable *table;
    double sample_index, x;
    double interp_samples[2];

    table_index = osc->wavetable->frequency_map[(int)osc->frequency];

    if(table_index == -1)
    {
        return sin(osc->phase);
    }

    table = &osc->wavetable->tables[table_index];

    sample_index = (osc->phase / M_2PI) * table->num_samples;

    interp_samples[0] = table->table[(int)sample_index % table->num_samples];
    interp_samples[1] = table->table[(int)ceil(sample_index) % table->num_samples];

    x = (sample_index - floor(sample_index));

    return interp_samples[0] + (interp_samples[1] - interp_samples[0]) * x;
}

#ifndef CLAMP
#define CLAMP(x, low, high)  (((x) > (high)) ? (high) : (((x) < (low)) ? (low) : (x)))
#endif

double morph_wavetable_signal ( oscillator* osc )
{
    int gradations = osc->morph_table->total_tables - 1;
    osc->morph_pos = CLAMP(osc->morph_pos, 0, 1);
    osc->wavetable = &osc->morph_table->tables[(int)(osc->morph_pos * gradations)];

    return wavetable_signal (osc);
}

double next_sample ( oscillator* osc )
{
    double out;

    switch(osc->type)
    {
        case FUNCTION:
            out = osc->source ( osc->phase, osc->frequency, osc->pulse_width );
            break;
        case WAVETABLE:
            out = wavetable_signal (osc);
            break;
        case MORPH_WAVETABLE:
            out = morph_wavetable_signal (osc);
            break;
        default:
            out = sin(osc->phase);
    }

    osc->phase += osc->sample_radians_per_hertz * osc->frequency;
    osc->phase  = fmod ( osc->phase, M_2PI );

    return out;
}

void skip_n_samples ( oscillator* osc, long nsamples )
{
    osc->phase += (double)nsamples * osc->frequency * osc->sample_radians_per_hertz;
    osc->phase   = fmod ( osc->phase, M_2PI );
}

double reset_oscillator ( oscillator* osc )
{
    osc->phase = osc->init_phase;
}
