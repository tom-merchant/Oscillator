//
// Created by tmerchant on 28/03/2020.
//

#include "wavetable.h"

#include <stdlib.h>
#include <string.h>
#ifdef _GNUC_
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#endif
#include <math.h>
#include <sys/param.h>
#include <fftw3.h>
#include <assert.h>
#include <complex.h>
#include <sndfile.h>
#include "AudioUtils.h"


const char MWT_MAGIC_NUMBER[7] = { 'M', 'W', 'T', (char)0xFA, (char)0x5D, (char)0x58, (char)0xA3 };
const char MWH_MAGIC_NUMBER[7] = { 'M', 'W', 'H', (char)0xB2, (char)0xA7, (char)0xC9, (char)0x41 };

wavetable_mipmap* generate_wavetables( wavetable_harmonics* harmonics, int oversample, int tables_per_octave)
{
    double reduction_factor = pow (2, 1.0f / (double)tables_per_octave);
    uint32_t i, table_size;
    uint32_t k;

    fftw_complex *freqs;
    fftw_complex *signal;
    fftw_plan plan;

    wavetable_mipmap *output = calloc (1, sizeof(wavetable_mipmap));
    wavetable_harmonics *harmonics_copy = malloc(sizeof(wavetable_harmonics));

    harmonics_copy->num_harmonics = harmonics->num_harmonics;
    harmonics_copy->cos_table = calloc (harmonics_copy->num_harmonics * 2 * oversample + 1, sizeof(double));
    harmonics_copy->sin_table = calloc (harmonics_copy->num_harmonics * 2 * oversample + 1, sizeof(double));

    output->num_wavetables = tables_per_octave * log2(harmonics->num_harmonics);
    output->tables = calloc (output->num_wavetables, sizeof(wavetable));

    memcpy ( harmonics_copy->cos_table + 1, harmonics->cos_table, harmonics->num_harmonics * sizeof (double) );
    memcpy ( harmonics_copy->sin_table + 1, harmonics->sin_table, harmonics->num_harmonics * sizeof (double));
    harmonics_copy->cos_table[0] = 0;
    harmonics_copy->sin_table[0] = 0;

    for(k = 0; k < output->num_wavetables; k++)
    {
        if(harmonics_copy->num_harmonics < 2)
        {
            output->num_wavetables = k + 1;
            break;
        }

        table_size = harmonics_copy->num_harmonics * 2 * oversample + 1;

        output->tables[k].num_harmonics = harmonics_copy->num_harmonics;
        output->tables[k].num_samples = table_size;
        output->tables[k].table = calloc(table_size, sizeof(double));
        output->tables[k].boundary_harmonic_cosine_amplitude = harmonics_copy->cos_table[table_size - 1];
        output->tables[k].boundary_harmonic_sine_amplitude = harmonics_copy->sin_table[table_size - 1];

        freqs = fftw_alloc_complex (table_size);
        signal = fftw_alloc_complex (table_size);

        for(i = 0; i < table_size; i++)
        {
            freqs[i][0] = harmonics_copy->cos_table[i];
            freqs[i][1] = harmonics_copy->sin_table[i];
        }

        plan = fftw_plan_dft_1d (table_size, freqs, signal, FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute (plan);

        for(i = 0; i < table_size; i++)
        {
            output->tables[k].table[i] = signal[i][0];
        }

        fftw_destroy_plan (plan);
        fftw_free (freqs);
        fftw_free (signal);

        remove_dc_offset ( output->tables[k].table, output->tables[k].num_samples );
        normalise ( output->tables[k].table, output->tables[k].num_samples );

        harmonics_copy->num_harmonics /= reduction_factor;

        for(i = harmonics_copy->num_harmonics + 1; i < (harmonics->num_harmonics * 2 * oversample + 1); i++)
        {
            harmonics_copy->cos_table[i] = 0;
            harmonics_copy->sin_table[i] = 0;
        }
    }

    free(harmonics_copy->cos_table);
    free(harmonics_copy->sin_table);
    free(harmonics_copy);

    return output;
}

//TODO: Does essentially the same thing as generate_wavetables except each mipmap is a constant size
wavetable_mipmap* generate_static_wavetables(wavetable_harmonics* harmonics, int oversample, int tables_per_octave)
{

}

morph_wavetable_mipmap* generate_morph_wavetables ( morph_wavetable_harmonics* harmonics, long samplerate )
{
    wavetable_harmonics* generated_harmonics;
    morph_wavetable_mipmap* out;
    uint32_t num_spectrums, i, k, waveform_index;
    wavetable_harmonics* harm1;
    wavetable_harmonics* harm2;
    wavetable_mipmap* tmp;
    double morph_time, theta, psin, pcos;
    double vertex1[2], point1[2], vertex2[2], point2[2];
    double a1, h1, c1, a2, h2, c2;
    double* cos_polynomial_points = NULL;
    double* sin_polynomial_points = NULL;
    double* polynomial_abscissae = NULL;
    double** cos_polynomial_coefficients = NULL;
    double** sin_polynomial_coefficients = NULL;
    uint32_t polynomial_order = 0;

    num_spectrums = ( harmonics->num_wavetables - 1 ) * harmonics->tables_per_morph;

    out = malloc ( sizeof ( morph_wavetable_mipmap ) );
    out->num_wavetables = harmonics->num_wavetables;
    out->total_tables = num_spectrums;
    out->tables_per_morph = harmonics->tables_per_morph;
    out->tables = malloc ( sizeof( wavetable_mipmap ) * num_spectrums );

    generated_harmonics = malloc ( num_spectrums * sizeof(wavetable_harmonics) );

    if(harmonics->morph_style == MORPH_POLYNOMIAL)
    {
        cos_polynomial_coefficients = calloc (harmonics->num_harmonics, sizeof(double*));
        sin_polynomial_coefficients = calloc (harmonics->num_harmonics, sizeof(double*));
        cos_polynomial_points = malloc ( harmonics->num_wavetables * sizeof ( double ) );
        sin_polynomial_points = malloc ( harmonics->num_wavetables * sizeof ( double ) );
        polynomial_abscissae     = malloc ( harmonics->num_wavetables * sizeof ( double ) );
        polynomial_order = harmonics->num_wavetables - 1;

        for ( i = 0; i < harmonics->num_wavetables; i++ )
        {
            polynomial_abscissae [ i ] = (double)i / (double)harmonics->num_wavetables;
        }

        for( k = 0; k < harmonics->num_harmonics; k++ )
        {

            for ( i = 0; i < harmonics->num_wavetables; i++ )
            {
                cos_polynomial_points [ i ] = harmonics->spectrums [ i ].cos_table [ k ];
                sin_polynomial_points [ i ] = harmonics->spectrums [ i ].sin_table [ k ];
            }

            cos_polynomial_coefficients[k] = malloc ( harmonics->num_wavetables * sizeof ( double ) );
            sin_polynomial_coefficients[k] = malloc ( harmonics->num_wavetables * sizeof ( double ) );

            calculate_newton_polynomial ( polynomial_order, cos_polynomial_points, polynomial_abscissae,
                    cos_polynomial_coefficients [ k ] );
            calculate_newton_polynomial ( polynomial_order, sin_polynomial_points, polynomial_abscissae,
                    sin_polynomial_coefficients [ k ] );
        }

        free(cos_polynomial_points);
        free(sin_polynomial_points);
    }

    for( i = 0; i < num_spectrums; i++ )
    {
        generated_harmonics[i].cos_table = calloc ( harmonics->num_harmonics, sizeof( double ) );
        generated_harmonics[i].sin_table = calloc ( harmonics->num_harmonics, sizeof( double ) );
        generated_harmonics[i].num_harmonics = harmonics->num_harmonics;

        waveform_index = i / harmonics->tables_per_morph;

        if( i % harmonics->tables_per_morph != 0)
        {
            harm1 = &harmonics->spectrums[waveform_index];
            harm2 = &harmonics->spectrums[waveform_index + 1];

            morph_time = (double)(i % harmonics->tables_per_morph) / (double)harmonics->tables_per_morph;

            switch(harmonics->morph_style)
            {
                default:
                case MORPH_LINEAR:
                    for ( k = 0; k < MIN(harm1->num_harmonics, harmonics->num_harmonics); k++ )
                    {
                        generated_harmonics[i].cos_table[k] = harm1->cos_table[k] + morph_time * ( harm2->cos_table[k] - harm1->cos_table[k] );
                        generated_harmonics[i].sin_table[k] = harm1->sin_table[k] + morph_time * ( harm2->sin_table[k] - harm1->sin_table[k] );
                    }
                    break;
                case MORPH_CONSTANT_POWER:
                    //I'm not actually sure if it make sense to do this, but we can, so why not?
                    theta = (M_PI_2 * morph_time) - M_PI_4;
#ifdef __GNUC__
                    sincos ( theta, &psin, &pcos );
#else
                    psin = sin ( theta );
                    pcos = cos ( theta );
#endif
                    for( k = 0; k < MIN(harm1->num_harmonics, harmonics->num_harmonics); k++ )
                    {
                        generated_harmonics[ i ].cos_table[ k ] =
                                harm1->cos_table[ k ] * ( ( M_SQRT2 / 2 ) * ( pcos - psin )) +
                                harm2->cos_table[ k ] * ( ( M_SQRT2 / 2 ) * ( pcos + psin ));
                        generated_harmonics[ i ].sin_table[ k ] =
                                harm1->sin_table[ k ] * ( ( M_SQRT2 / 2 ) * ( pcos - psin )) +
                                harm2->sin_table[ k ] * ( ( M_SQRT2 / 2 ) * ( pcos + psin ));
                    }
                    break;
                case MORPH_POLYNOMIAL:
                    morph_time = (double) i / (double)num_spectrums;

                    for( k = 0; k < MIN(harm1->num_harmonics, harmonics->num_harmonics); k++ )
                    {
                        generated_harmonics [ i ].cos_table [ k ] =
                                newton_horner_eval ( polynomial_order, cos_polynomial_coefficients [ k ],
                                        polynomial_abscissae, morph_time );
                        generated_harmonics [ i ].sin_table [ k ] =
                                newton_horner_eval ( polynomial_order, sin_polynomial_coefficients [ k ],
                                        polynomial_abscissae, morph_time );
                    }
                    break;
                case MORPH_RAMP_UP:
                    for ( k = 0; k < MIN(harm1->num_harmonics, harmonics->num_harmonics); k++ )
                    {
                        vertex1[0] = 0;
                        vertex1[1] = harm1->cos_table[k];
                        vertex2[0] = 0;
                        vertex2[1] = harm1->sin_table[k];
                        point1[0] = 1;
                        point1[1] = harm2->cos_table[k];
                        point2[0] = 1;
                        point2[1] = harm2->sin_table[k];

                        h1 = vertex1[0];
                        c1 = vertex1[1];
                        a1 = (point1[1] - c1) / pow(point1[0] - h1, 2);
                        h2 = vertex2[0];
                        c2 = vertex2[1];
                        a2 = (point2[1] - c2) / pow(point2[0] - h2, 2);

                        generated_harmonics[i].cos_table[k] = a1 * pow ( morph_time - h1, 2  ) + c1;
                        generated_harmonics[i].sin_table[k] = a2 * pow ( morph_time - h2, 2  ) + c2;
                    }
                    break;
                case MORPH_RAMP_DOWN:
                    for ( k = 0; k < MIN(harm1->num_harmonics, harmonics->num_harmonics); k++ )
                    {
                        point1[0] = 0;
                        point1[1] = harm1->cos_table[k];
                        point2[0] = 0;
                        point2[1] = harm1->sin_table[k];
                        vertex1[0] = 1;
                        vertex1[1] = harm2->cos_table[k];
                        vertex2[0] = 1;
                        vertex2[1] = harm2->sin_table[k];

                        h1 = vertex1[0];
                        c1 = vertex1[1];
                        a1 = (point1[1] - c1) / pow(point1[0] - h1, 2);
                        h2 = vertex2[0];
                        c2 = vertex2[1];
                        a2 = (point2[1] - c2) / pow(point2[0] - h2, 2);

                        generated_harmonics[i].cos_table[k] = a1 * pow ( morph_time - h1, 2  ) + c1;
                        generated_harmonics[i].sin_table[k] = a2 * pow ( morph_time - h2, 2  ) + c2;
                    }
                    break;
                case MORPH_ENVELOPE:
                    for ( k = 0; k < MIN(harm1->num_harmonics, harmonics->num_harmonics); k++ )
                    {
                        generated_harmonics[i].cos_table[k] =
                                ( 1 - value_at ( harmonics->morph_env, morph_time ) ) * harm1->cos_table[k] +
                                      value_at ( harmonics->morph_env, morph_time )   * harm2->cos_table[k];
                        generated_harmonics[i].sin_table[k] =
                                ( 1 - value_at ( harmonics->morph_env, morph_time ) ) * harm1->sin_table[k] +
                                      value_at ( harmonics->morph_env, morph_time )   * harm2->sin_table[k];
                    }
                    break;
            }
        }
        else
        {
            memcpy ( generated_harmonics[i].cos_table, harmonics->spectrums[waveform_index].cos_table,
                    MIN(harmonics->spectrums[waveform_index].num_harmonics, harmonics->num_harmonics) * sizeof( double ));
            memcpy ( generated_harmonics[i].sin_table, harmonics->spectrums[waveform_index].sin_table,
                     MIN(harmonics->spectrums[waveform_index].num_harmonics, harmonics->num_harmonics) * sizeof( double ));
        }

        tmp = generate_wavetables ( &generated_harmonics[i], harmonics->oversample, harmonics->tables_per_octave );
        generate_wavetable_freqmap ( tmp, samplerate );
        memcpy ( &out->tables[i], tmp, sizeof(wavetable_mipmap) );
        free ( tmp );

        free ( generated_harmonics[i].cos_table );
        free ( generated_harmonics[i].sin_table );
    }

    if(harmonics->morph_style == MORPH_POLYNOMIAL)
    {
        assert(cos_polynomial_coefficients);
        assert(sin_polynomial_coefficients);

        for( k = 0; k < harmonics->num_harmonics; k++ )
        {
            free ( cos_polynomial_coefficients [ k ] );
            free ( sin_polynomial_coefficients [ k ] );
        }
        free(polynomial_abscissae);
        free(cos_polynomial_coefficients);
        free(sin_polynomial_coefficients);
    }

    free ( generated_harmonics );

    return out;
}

wavetable_mipmap* generate_wavetables_from_signal ( const double* signal, int num_samples, int num_harmonics, int oversample, int tables_per_octave )
{
    int i;
    fftw_complex *freqs;
    fftw_complex *signal_array;
    fftw_plan plan;
    wavetable_harmonics *signal_spectrum;
    wavetable_mipmap *out;

    freqs = fftw_alloc_complex (num_samples);
    signal_array = fftw_alloc_complex (num_samples);

    for(i = 0; i < num_samples; i++)
    {
        signal_array[i][0] = signal[i];
        signal_array[i][1] = 0;
    }
    plan = fftw_plan_dft_1d (num_samples, signal_array, freqs, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute (plan);

    signal_spectrum = malloc (sizeof(signal_spectrum));
    signal_spectrum->cos_table = malloc ( num_harmonics );
    signal_spectrum->sin_table = malloc ( num_harmonics );

    for(i = 1; i < MIN(num_samples, num_harmonics + 1); i++)
    {
        signal_spectrum->cos_table[i - 1] = freqs[i][0];
        signal_spectrum->sin_table[i - 1] = freqs[i][1];
    }

    fftw_free (freqs);
    fftw_free (signal_array);
    fftw_destroy_plan (plan);

    out = generate_wavetables ( signal_spectrum, oversample, tables_per_octave );

    free(signal_spectrum->cos_table);
    free(signal_spectrum->sin_table);
    free(signal_spectrum);

    return out;
}

void save_wavetable(wavetable_mipmap* wt, const char* path)
{
    uint32_t i, j;
    int sample_buf_ptr;
    uint32_t* wavetable_samples;
    uint32_t* wavetable_harmonics;
    unsigned char* sample_buf;
    FILE* out;

    uint32_t out_sample;

    sample_buf = malloc (3 * wt->tables[0].num_samples);
    wavetable_samples = malloc (wt->num_wavetables * 4);
    wavetable_harmonics = malloc (wt->num_wavetables * 4);

    for(i = 0; i < wt->num_wavetables; i++)
    {
        wavetable_samples[i] = wt->tables[i].num_samples;
        wavetable_harmonics[i] = wt->tables[i].num_harmonics;
    }

    out = fopen (path, "w");

    fwrite (MWT_MAGIC_NUMBER, 1, 7, out);

    /*TODO: Endianness*/
    fwrite(&wt->num_wavetables, 4, 1, out);
    fwrite (wavetable_samples, 4, wt->num_wavetables, out);
    fwrite (wavetable_harmonics, 4, wt->num_wavetables, out);

    for(i = 0; i < wt->num_wavetables; i++)
    {
        sample_buf_ptr = 0;

        for(j = 0; j < wt->tables[i].num_samples; j++)
        {
            out_sample = (uint32_t) ((( wt->tables[ i ].table[ j ] + 1 ) / 2.0f ) * 16777215.0f ); //16777216 = 2^16 - 1

            sample_buf[ sample_buf_ptr++ ] = out_sample >> 16u;
            sample_buf[ sample_buf_ptr++ ] = out_sample >> 8u;
            sample_buf[ sample_buf_ptr++ ] = out_sample & 0xFFu;

        }

        fwrite (sample_buf, 1, sample_buf_ptr, out);
    }

    fclose (out);
    free(sample_buf);
    free(wavetable_harmonics);
    free(wavetable_samples);
}

//TODO: save and load wavetables in other wavetable formats

long get_serum_wavetable_framesize ( char* path )
{
    FILE *wtfile = fopen(path, "rb");
    int found = 0;
    long sz;

    int morph_table = 0;

    fseek(wtfile, 0L, SEEK_END);
    sz = ftell(wtfile);
    rewind (wtfile);

    char c;
    char c1[5];

    char clmdata[48];

    while (!found && ftell (wtfile) < sz)
    {
        fread (&c, 1, 1, wtfile);

        if ( c == 0x63 )
        {
            fread ( &c1, 1, 2,wtfile );

            if ( c1[0] == 0x6C && c1[1] == 0x6D )
            {
                found = 1;
            }
        }
    }

    fread (c1, 1, 5, wtfile);
    fread (clmdata, 1, 48, wtfile);

    fclose ( wtfile );

    return strtol (&clmdata[3], NULL, 10);
}

long serum_wavetable_get_num_tables ( char* path )
{
    long framesize = get_serum_wavetable_framesize ( path );
    long ntables = 0;

    SF_INFO info = {};

    SNDFILE* wtfile = sf_open (path, SFM_READ, &info);

    ntables = info.frames / framesize;

    sf_close (wtfile);

    return ntables;
}

morph_wavetable_mipmap* load_serum_morph_wavetable ( char* path, int tables_per_octave, int oversample,
        int tables_per_morph, long samplerate )
{
    long framesize = get_serum_wavetable_framesize ( path );
    long ntables = serum_wavetable_get_num_tables ( path );
    long i, j;

    morph_wavetable_harmonics mwh;

    mwh.num_wavetables = ntables;
    mwh.morph_style = MORPH_LINEAR;
    mwh.morph_env = 0;
    mwh.spectrums = calloc (ntables, sizeof (wavetable_harmonics));
    mwh.num_harmonics = framesize;
    mwh.tables_per_octave = tables_per_octave;
    mwh.oversample = oversample;
    mwh.tables_per_morph = tables_per_morph;

    SF_INFO info = {};
    SNDFILE* wtfile = sf_open (path, SFM_READ, &info);

    fftw_complex *freqs = fftw_alloc_complex ( framesize + 1 );
    double* signal = malloc ( sizeof (double) * framesize );
    fftw_complex *signal_array = fftw_alloc_complex ( framesize );

    fftw_plan plan = fftw_plan_dft_1d (framesize, signal_array, freqs, FFTW_FORWARD, FFTW_ESTIMATE);

    for ( i = 0; i < ntables; i++ )
    {
        mwh.spectrums [ i ].num_harmonics = framesize;
        mwh.spectrums [ i ].cos_table = malloc ( sizeof (double) * framesize );
        mwh.spectrums [ i ].sin_table = malloc ( sizeof (double) * framesize );

        sf_read_double (wtfile, signal, framesize);

        for(j = 0; j < framesize; j++)
        {
            signal_array[j][0] = signal[j];
            signal_array[j][1] = 0;
        }

        fftw_execute (plan);

        for(j = 1; j < framesize; j++)
        {
            mwh.spectrums [ i ].cos_table [ j - 1 ] = freqs [ j ][ 0 ] / framesize;
            mwh.spectrums [ i ].sin_table [ j - 1 ] = freqs [ j ][ 1 ] / framesize;
        }
    }

    sf_close ( wtfile );
    fftw_destroy_plan (plan);
    fftw_free (freqs);
    fftw_free (signal_array);

    return generate_morph_wavetables (&mwh, samplerate);
}

wavetable_mipmap* load_serum_wavetable (char* path, int tables_per_octave, int oversample, long samplerate)
{
    long framesize = get_serum_wavetable_framesize ( path );
    long j;

    wavetable_harmonics wh;

    SF_INFO info = {};
    SNDFILE* wtfile = sf_open (path, SFM_READ, &info);

    fftw_complex *freqs = fftw_alloc_complex ( framesize );
    double* signal = malloc ( sizeof (double) * framesize );
    fftw_complex *signal_array = fftw_alloc_complex ( framesize );

    fftw_plan plan = fftw_plan_dft_1d (framesize, signal_array, freqs, FFTW_FORWARD, FFTW_ESTIMATE);

    wh.num_harmonics = framesize;
    wh.cos_table = malloc ( sizeof (double) * framesize );
    wh.sin_table = malloc ( sizeof (double) * framesize );

    sf_read_double (wtfile, signal, framesize);

    for(j = 0; j < framesize; j++)
    {
        signal_array[j][0] = signal[j];
        signal_array[j][1] = 0;
    }

    fftw_execute (plan);

    for(j = 1; j < framesize; j++)
    {
        wh.cos_table [ j ] = freqs [ j - 1 ][ 0 ];
        wh.sin_table [ j ] = freqs [ j - 1 ][ 1 ];
    }

    sf_close ( wtfile );
    fftw_destroy_plan (plan);
    fftw_free (freqs);
    fftw_free (signal_array);

    return generate_wavetables ( &wh, oversample, tables_per_octave );
}

wavetable_mipmap* load_wavetable(const char* path, int samplerate)
{
    uint32_t i, j;
    char magic[7];
    unsigned char* inbuf;
    uint32_t in_sample;
    uint32_t* wavetable_samples;
    uint32_t* wavetable_harmonics;
    wavetable_mipmap* wt;
    FILE* in;

    wt = malloc ( sizeof (wavetable_mipmap) );

    in = fopen(path, "r");

    fread (magic, 1, 7, in);

    if(strncmp (magic, MWT_MAGIC_NUMBER, 7) != 0)
    {
        return NULL;
    }

    /*TODO: Endianness*/
    fread (&wt->num_wavetables, 4, 1, in);

    wavetable_samples = malloc (wt->num_wavetables * 4);
    wavetable_harmonics = malloc (wt->num_wavetables * 4);
    wt->tables = calloc (wt->num_wavetables, sizeof(wavetable));

    fread (wavetable_samples, 4, wt->num_wavetables, in);
    fread (wavetable_harmonics, 4, wt->num_wavetables, in);

    inbuf = malloc (wavetable_samples[0] * 4);

    for(i = 0; i < wt->num_wavetables; i++)
    {
        wt->tables[i].num_samples = wavetable_samples[i];
        wt->tables[i].num_harmonics = wavetable_harmonics[i];
        wt->tables[i].table = malloc (wavetable_samples[i] * sizeof( double ));

        for(j = 0; j < wavetable_samples[i]; j++)
        {
            fread (inbuf, 1, 3, in);
            in_sample =  (uint32_t)(inbuf[0] << 16u) & 0xFF0000u;
            in_sample |= (uint32_t)(inbuf[1] << 8u) & 0x00FF00u;
            in_sample |= inbuf[2] & 0xFFu;

            wt->tables[i].table[j] = ((in_sample / 16777215.0f) * 2.0f) - 1.0f;
        }
    }

    fclose (in);
    free(inbuf);
    free(wavetable_samples);
    free(wavetable_harmonics);

    generate_wavetable_freqmap(wt, samplerate);

    return wt;
}

void save_morph_wavetable_harmonics ( const morph_wavetable_harmonics *mwh, const char* path )
{
    uint32_t i, num_breakpoints = 0;
    FILE *out;

    out = fopen(path, "w");

    fwrite (MWH_MAGIC_NUMBER, 1, 7, out);
    fwrite (&mwh->num_harmonics, 4, 1, out);
    fwrite (&mwh->num_wavetables, 4, 1, out);
    fwrite (&mwh->tables_per_morph, 4, 1, out);
    fwrite (&mwh->oversample, 4, 1, out);
    fwrite (&mwh->tables_per_octave, 4, 1, out);
    fwrite (&mwh->morph_style, sizeof(morph_type), 1, out);

    for(i = 0; i < mwh->num_wavetables; i++)
    {
        fwrite (&mwh->spectrums[i].num_harmonics, 4, 1, out);
        fwrite (mwh->spectrums[i].cos_table, sizeof(double), mwh->spectrums[i].num_harmonics, out);
        fwrite (mwh->spectrums[i].sin_table, sizeof(double), mwh->spectrums[i].num_harmonics, out);
    }

    if(mwh->morph_style == MORPH_ENVELOPE)
    {
        breakpoint *current_bp = mwh->morph_env->first;

        while(current_bp)
        {
            num_breakpoints++;
            current_bp = current_bp->next;
        }

        fwrite(&num_breakpoints, 4, 1, out);
        fputc ('\n', out);

        current_bp = mwh->morph_env->first;

        while ( current_bp )
        {
            fprintf ( out, "%f ", current_bp->time );
            fprintf ( out, "%f ", current_bp->value);
            fprintf ( out, "%d", current_bp->interpType);

            for ( i = 0; i < current_bp->nInterp_params; i++ )
            {
                fprintf ( out, " %f", current_bp->interp_params[i] );
            }

            fputc ( '\n', out );

            current_bp = current_bp->next;
        }
    }

    fclose(out);
}

morph_wavetable_harmonics* load_morph_wavetable_harmonics ( const char* path )
{
    uint32_t i, num_breakpoints, num_write;
    char magic_num[7];
    char buf[512];
    morph_wavetable_harmonics* out;
    FILE *in, *tmp;

    out = malloc (sizeof(morph_wavetable_harmonics));

    in = fopen(path, "r");

    fread(magic_num, 1, 7, in);

    if(strncmp (magic_num, MWH_MAGIC_NUMBER, 7) != 0)
    {
        return NULL;
    }

    fread (&out->num_harmonics, 4, 1, in);
    fread (&out->num_wavetables, 4, 1, in);
    fread (&out->tables_per_morph, 4, 1, in);
    fread (&out->oversample, 4, 1, in);
    fread (&out->tables_per_octave, 4, 1, in);
    fread (&out->morph_style, sizeof(morph_type), 1, in);

    out->spectrums = malloc (out->num_harmonics * sizeof(wavetable_harmonics));

    for(i = 0; i < out->num_wavetables; i++)
    {
        fread (&out->spectrums[i].num_harmonics, 4, 1, in);

        out->spectrums[i].cos_table = malloc (out->spectrums[i].num_harmonics * sizeof(double));
        out->spectrums[i].sin_table = malloc (out->spectrums[i].num_harmonics * sizeof(double));

        fread (out->spectrums[i].cos_table, sizeof(double), out->spectrums[i].num_harmonics, in);
        fread (out->spectrums[i].sin_table, sizeof(double), out->spectrums[i].num_harmonics, in);
    }

    if(out->morph_style == MORPH_ENVELOPE)
    {
        out->morph_env = malloc(sizeof(envelope));
        fread ( &num_breakpoints, 4, 1, in );
        fgetc ( in );

        tmp = fopen("tmp4yn8vyr32rc.bp", "w");

        while((num_write = fread (buf, 1, 512, in)) > 0)
        {
            fwrite (buf, 1, num_write, tmp);
        }

        fclose (tmp);

        load_breakpoints ("tmp4yn8vyr32rc.bp", out->morph_env);

        remove ("tmp4yn8vyr32rc.bp");
    }

    return out;
}

void free_wavetable_mipmap(wavetable_mipmap* wt)
{
    uint32_t i;

    for(i = 0; i < wt->num_wavetables; i++)
    {
        free(wt->tables[i].table);
    }

    if(wt->frequency_map)
    {
        free(wt->frequency_map);
    }

    free(wt->tables);
    free(wt);
}

void free_morph_wavetable_mipmap(morph_wavetable_mipmap* mwm)
{
    uint32_t i, j;

    for(i = 0; i < mwm->total_tables; i++)
    {
        for(j = 0; j < mwm->tables[i].num_wavetables; j++)
        {
            free(mwm->tables[i].tables[j].table);
        }

        if(mwm->tables[i].frequency_map)
        {
            free(mwm->tables[i].frequency_map);
        }
        
        free(mwm->tables[i].tables);
    }

    free(mwm->tables);
    free(mwm);
}

void generate_wavetable_freqmap(wavetable_mipmap* wt, long samplerate)
{
    long f;
    int current_mipmap;

    if (wt->frequency_map)
    {
        free (wt->frequency_map);
    }

    wt->frequency_map = calloc (samplerate / 2, sizeof(int));

    current_mipmap = 0;

    for(f = 0; f < samplerate / 2; f++)
    {
        while ((f + 1) * wt->tables[current_mipmap].num_harmonics > samplerate / 2 && current_mipmap < (wt->num_wavetables - 1))
        {
            current_mipmap++;
        }

        if ( (f + 1) * wt->tables[current_mipmap].num_harmonics > samplerate / 2
            || wt->tables[current_mipmap].num_harmonics <= 2 )
        {
            wt->frequency_map[f] = -1;
            continue;
        }

        wt->frequency_map[f] = current_mipmap;
    }
}

