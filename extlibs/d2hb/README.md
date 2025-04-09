# DHB & D2HB

This is the repository of the Dynamic Hashed Blocks (DHB) extended with the
Dynamic Doubly Hashed Blocks (D2HB) graph data structure. The work of adding the
dynamic hash table as a layer before the edge blocks also includes minor changes
in the overall library.

## Authors

### D2HB

| Name                      | E-Mail                                  | Affiliation |
|---------------------------|-----------------------------------------|-------------|
| Florian Willich           | florian.willich@informatik.hu-berlin.de | HU Berlin   |

### DHB

| Name                      | E-Mail                                  | Affiliation |
|---------------------------|-----------------------------------------|-------------|
| Alexander van der Grinten | avdgrinten@hu-berlin.de                 | HU Berlin   |
| Maria Predari             | predarim@informatik.hu-berlin.de        | HU Berlin   |
| Florian Willich           | florian.willich@informatik.hu-berlin.de | HU Berlin   |

## Publications

For a detailed discussion of DHB and a comparison with other graph data structures
please read the paper _A Fast Data Structure for Dynamic Graphs Based on 
Hash-Indexed Adjacency Blocks_ (2022) available on the [Dagstuhl Research Online 
Publication Server (DROPS)](https://drops.dagstuhl.de/opus/volltexte/2022/16545/).

# Build Library

Build the library using CMake (e.g., with the generator Ninja).

```bash
mkdir build
cd build
cmake -GNinja -DCMAKE_BUILD_TYPE=Release ..
ninja
```

# Build Tests

You must clone this repository recursively in order to obtain catch2 which is
necessary for our test environment using:

```bash
git clone --recursive https://github.com/hu-macsy/dhb.git
```

If you forgot to clone recursive, simply update with initialisation:

```bash
git submodule update --init --recursive
```

You can then build the tests by setting the CMake option `D2HB_TEST` to `On`.

```
cmake -GNinja -DCMAKE_BUILD_TYPE=Release -DD2HB_TEST=On ..
ninja
```

# Build Demo

In case you want to test some functionality using the demo application you must set `D2HB_DEMO=On` as in:

```
cmake -GNinja -DCMAKE_BUILD_TYPE=Release -DD2HB_TEST=On -DD2HB_DEMO=On ..
ninja
```

including the tests.

# Benchmarks and Experimental Evaluation

An evaluation of DHB and other dynamic graph data structures can be found at: https://github.com/hu-macsy/dhb-experiments
