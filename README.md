# Minimal example of path tracing in F#

This is an implementation of a modified version of [smallpt](https://www.kevinbeason.com/smallpt/) in F#, which is aimed as part of the course project of UCSB 2025 winter CS 263.

This implementation is parallelized using `System.Threading.Tasks`. For more detail, please refer to the source code in `Program.fs`.

Here is a sample output of the program with 64 spp:
![Sample - 64spp](./sample.png)

