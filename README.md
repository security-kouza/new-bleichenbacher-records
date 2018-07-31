## Dependencies
### FFTW & GMP
```
sudo apt-get update && sudo apt-get install libfftw3-3 libfftw3-bin libfftw3-dev libgmp10 libgmpxx4ldbl
```

### Open MPI
  1. [Make sure that the OS has NUMA support](https://www.open-mpi.org/faq/?category=building#build-paffinity). In Ubuntu 16.04, install the following packages:
```
sudo apt install numactl hwloc libhwloc-dev libnuma1 libnuma-dev
```
  2. Follow the instruction [here](https://www.open-mpi.org/faq/?category=building#easy-build)
  3. `sudo ldconfig`

### Boost
  1. Download the source [here](https://dl.bintray.com/boostorg/release/1.65.1/source/boost_1_65_1.tar.gz) and unpack it
  2. `./bootstrap.sh --prefix=path/to/project/root --with-libraries=mpi,program_options`
  3. Open `project-config.jam` then add `using mpi ;` to the end
  4. `./b2 install` to build and install Boost.MPI

### qDSA
The qDSA source code under `qDSA/Curve25519-asm` is based on the [qDSA reference implementation](https://www.cs.ru.nl/~jrenes/software/cref-g1.tar.gz) (J. Renes) and [Ed25519](http://ed25519.cr.yp.to/software.html) (D.J.Bernstein et al.)

## Build
Run `make all` at the project root to create 3 executables: `test_fft`, `attack_mpi`, `siggen_mpi`.

## Usage of `siggen_mpi` and `attack_mpi`
### Basic invocation
```
LD_LIBRARY_PATH=./lib mpirun -np <number-of-cores> attack_mpi <options>
```

### Common options
- `--verbose`: enable verbose logging, raise the precision of progress
- `--test`: run without the real qDSA signing algorithm

### `siggen_mpi`
- `--out <filename>`: save generated signature to a file
- `--leak <num_of_bits>`: number of nonce LSBs to be leaked
- `--filter <num_of_bits>`: filter signatures by the values of h

#### Example
- Generate 2-bit biased  preprocessed qDSA signatures with the top 19-bit filtered
```
./siggen_mpi --leak 2 --filter 19 --out qdsa_a24_b2_f19
```

### `attack_mpi`
- `--in <file_prefix>`: load signature data and execute an attack
- `--fcount <number>`: specifiy number of files to be loaded
- `--known <number>`: number of known MSBs of the secret
- `--red`: perform reduction
- `--fft`: perform key recevoery

#### Example
- Quick test using pseudo-Schnorr signature over 90-bit group
``` 
OMP_NUM_THREADS=1 mpirun -np 1 -x OMP_NUM_THREADS -map-by ppr:1:core --report-bindings --display-map ./attack_mpi --in data/test90_a12_b2_f10 --test --red --fft
```

- 5-round reduction test
```
OMP_NUM_THREADS=1 mpirun -np 1 -x OMP_NUM_THREADS -map-by ppr:1:core --report-bindings --display-map ./attack_mpi --in data/test252_a15_b0_f0 --red
```

### MPI options
- `--bind-to`
- `--rank-by`
- `--map-by`
- `-x ENV_VAR_TO_PASS`

#### For debugging
- `--display-map`
- `--display-allocation`
- `--report-bindings`

### OpenMP options
- `OMP_NUM_THREADS=16`
- `GOMP_CPU_AFFINITY="0-16"`

## Further reading
### OpenMP / MPI
- [How to gain hybrid MPI-OpenMP code performance without changing a line of code a.k.a. dealing with task affinity](https://aciref.org/how-to-gain-hybrid-mpi-openmp-code-performance-without-changing-a-line-of-code-a-k-a-dealing-with-task-affinity/)
- [OpenMPI Wiki: Process Placement](https://github.com/open-mpi/ompi/wiki/ProcessPlacement)
- [mpirun man page (version 3.0.0)](https://www.open-mpi.org/doc/current/man1/mpirun.1.php)
- [FAQ: Building Open MPI](https://www.open-mpi.org/faq/?category=building)
- [FAQ: General run-time tuning](https://www.open-mpi.org/faq/?category=tuning)
- [libgomp environment variables](https://gcc.gnu.org/onlinedocs/libgomp/Environment-Variables.html#Environment-Variables)
