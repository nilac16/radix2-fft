
# Radix-2 decimation-in-time

The simplest FFT algorithm, implemented in C99 (with GCC builtins) using its standard complex type.
One and two-dimensional **in-place** forward and inverse transforms are provided.

The one-dimensional transform outperforms the otherwise excellent [GSL](https://www.gnu.org/software/gsl/) by a factor of two on my machine.
The two-dimensional transform suffers somewhat from poor cache utilization and associativity effects.

Here, have some images:

[inputid]:  /images/input.jpg
[outputid]: /images/output.jpg
[filterid]: /images/filtered.jpg
[inversid]: /images/inverse.jpg

| Input image  | Output FFT2   |
|            - | -             |
| ![][inputid] | ![][outputid] |

Applying a small rectangular notch filter to eliminate DC and low-frequency components yields the following results:

| Notch filter  | Inverse FFT2  |
|             - | -             |
| ![][filterid] | ![][inversid] |
