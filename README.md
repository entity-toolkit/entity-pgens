# Problem Generators for the Entity code

This separate repository is intended to store non-standard user-contributed problem generators. The validity of these problem generators is not checked on a regular basis, so some of them might be outdated or use outdated API. If you find such problem generators -- please open a github issue. 

## Usage

This repository is automatically added as a submodule to the main Entity code. To use any of the problem generators with the Entity itself, simply download the submodule if not already (from the Entity root):

```sh
git submodule update --init
```

then, when configuring, specify the name of the pgen with a `pgens/` prefix. E.g.,

```sh
cmake -B build -D pgen=pgens/kelvin-helmholtz
```

The code should automatically use the pgen from this repository (which is downloaded into `extern/entity-pgens/` directory of the main code).

## Contributing a problem generator

Each problem generator is put in a specific directory that bears a descriptive and unique name. In that directory, the main file shall always named `pgen.hpp`. Technically, this is the only required file to compile the code. Together with it, however, it is a good practice to include several other files:

- `README.md` with a brief description of what the problem generator does. Within it, you can add the following:
  - input parameters with their types and default values;
  - a schematic illustration of what the setup looks like;
  - an image/plot of what the result should look like;
- `toml` file with example input parameters;
- `py` file with a plotting routine specific to that particular setup.

To add a problem generator, either create a separate branch, add it there with all the additional components listed above, and then open a pull request (alternatively, for this repository, and open PR from there).
