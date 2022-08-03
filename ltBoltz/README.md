To compile codes


## 0 - Clone Repos

Clone the amar-codes repo and initialze the submodules with:

```bash
git submodule init
```

once both submodules are initialized, update the AMReX one only with:

```bash
git submodule update submods/amrex
```

once it has cloned, proceed to step 1.

## 1 - Build LBM Code

Navigate to the ltBoltz/Build directory and run the following 3 lines:

```bash
module purge
module load openmpi/4.1.3/gcc-11.3.0-cuda-11.7
make
```

if this is the first time, it may take a while while it compiles the necessary portions of AMReX (<10 minutes). Once this is done, copy the executable to the test folder: `cp ltBoltz3d.gnu.MPI.ex ../test/`.

## 2 - Run Case

Navigate to the folder ltBoltz/test/case1 and run the following:

```bash
srun -n 10 ../ltBoltz3d.gnu.MPI.ex inputs
```

where `inputs` is the name of the input file.

