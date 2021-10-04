# Raptor utility repository

# Moved to https://github.com/seqan/raptor/tree/master/util

This contains small apps.

```bash
git clone --recurse-submodules https://github.com/eseiler/raptor_data_simulation
cd raptor_data_simulation
mkdir build
cd build
cmake ..
make -j2 install
```

There is a script at `src/simulate.sh` that simulates a dataset.<br>
Variables in upper case can be changed.<br>
`BINARY_DIR` should be the absolute path to the `build/bin` directory.<br>
`OUT_DIR` should be the absolute path to the output directory.<br>

**LENGTH % BIN_NUMBER should be 0**<br>
**(LENGTH / BIN_NUMBER) % HAPLOTYPE_COUNT should be 0**<br>
**READ_COUNT % BIN_NUMBER should be 0**<br>
The easiest way to achieve this is to set LENGTH, BIN_NUMBER. and READ_COUNT to a power of two.<br>

For
```
OUT_DIR=/some/path
BIN_NUMBER=16384
ERRORS=2
READ_LENGTHS="100 150 250"
```
, the result will look like
```
/some/path/
└── 16384
    ├── bins
    ├── info
    ├── reads_e2_100
    ├── reads_e2_150
    └── reads_e2_250
```
