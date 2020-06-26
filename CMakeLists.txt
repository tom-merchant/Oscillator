cmake_minimum_required(VERSION 3.15)
project(Oscillator C)

set(CMAKE_C_STANDARD 99)

link_directories(../envelope)
include_directories(../envelope)
include_directories(../AudioUtils)
link_directories(../AudioUtils/cmake-build-debug)

add_library(Oscillator SHARED oscillator.c oscillator.h wavetable.h wavetable.c)

add_executable(GenWavetables generate_wavetables.c)


target_link_libraries(Oscillator m fftw3 envelope AudioUtils sndfile)
target_link_libraries(GenWavetables m Oscillator sndfile envelope AudioUtils)