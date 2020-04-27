
#include "wavetable.h"

#include <stdlib.h>

int main()
{
    int i, j;

    wavetable_harmonics square_harmonics, triangle_harmonics, sawtooth_harmonics;

    square_harmonics.num_harmonics = NUM_HARMONICS;
    square_harmonics.cos_table = calloc(NUM_HARMONICS, sizeof(double));
    square_harmonics.sin_table = calloc(NUM_HARMONICS, sizeof(double));

    for(i = 0; i < NUM_HARMONICS; i++)
    {
        if((i + 1) % 2)
        {
            square_harmonics.sin_table[ i ] = 1.0f / (double)( i + 1 );
        }
    }

    triangle_harmonics.num_harmonics = NUM_HARMONICS;
    triangle_harmonics.cos_table = calloc (NUM_HARMONICS, sizeof ( double));
    triangle_harmonics.sin_table = calloc (NUM_HARMONICS, sizeof ( double));

    j = 0;

    for(i = 0; i < NUM_HARMONICS; i++)
    {
        if((i + 1) % 2)
        {
            triangle_harmonics.sin_table[ i ] = 1.0f / (double)((i + 1) * (i + 1));

            if(j)
            {
                triangle_harmonics.sin_table[i] *= -1;
            }

            j = (j ? 0 : 1);
        }
    }

    sawtooth_harmonics.num_harmonics = NUM_HARMONICS;
    sawtooth_harmonics.cos_table = calloc (NUM_HARMONICS, sizeof(double));
    sawtooth_harmonics.sin_table = calloc (NUM_HARMONICS, sizeof(double));

    for(i = 0; i < NUM_HARMONICS; i++)
    {
        sawtooth_harmonics.sin_table[i] = 1.0f / (double)(i + 1);
    }

    wavetable_mipmap *square_mipmap, *triangle_mipmap, *sawtooth_mipmap;

    square_mipmap = generate_wavetables (&square_harmonics, OVERSAMPLE, TABLES_PER_OCTAVE);
    triangle_mipmap = generate_wavetables (&triangle_harmonics, OVERSAMPLE, TABLES_PER_OCTAVE);
    sawtooth_mipmap = generate_wavetables (&sawtooth_harmonics, OVERSAMPLE, TABLES_PER_OCTAVE);

    save_wavetable (square_mipmap, "square.mwt");
    save_wavetable (triangle_mipmap, "triangle.mwt");
    save_wavetable (sawtooth_mipmap, "sawtooth.mwt");

    free (square_harmonics.cos_table);
    free (square_harmonics.sin_table);
    free (triangle_harmonics.cos_table);
    free (triangle_harmonics.sin_table);
    free (sawtooth_harmonics.cos_table);
    free (sawtooth_harmonics.sin_table);
    free_wavetable_mipmap(square_mipmap);
    free_wavetable_mipmap(triangle_mipmap);
    free_wavetable_mipmap(sawtooth_mipmap);
}