# Biers (Basic Inference Engine for RNA Structure)

**Biers** is a set of *MATLAB* scripts that wrap around the *RNAstructure* suite of executables to infer secondary structure models guided by chemical mapping data. Features include:

- Secondary structure modeling guided by SHAPE and DMS probing

- Secondary structure modeling from mutate-and-map (M2 and M2-seq) experiments.

- Secondary structure modeling guided by MOHCA data.

- Calculation of helix-wise confidence values using bootstrapping.

- Plotting utilities using *MATLAB* and *VARNA*.


## Installation

To install **Biers**, simply:

- Download the zip or tar file of the repository and unpack; or 
```bash
git clone https://github.com/ribokit/Biers.git
```

- In *MATLAB*, go to "**Set Path**". Then "**Add with Subfolders**" of the target `path/to/Biers/Scripts/`.

- Make sure you have set in your environment the variables `DATAPATH` (for RNAstructure) and `VARNA` (for VARNA). If you are on a Mac or Linux, put lines like the following in your `.bashrc`:

```bash
export DATAPATH=/path/to/RNAstructure/data_tables/
export VARNA=/path/to/src/VARNA.jar
```
  Note that you should have set up `DATAPATH` as part of the RNAstructure installation process.

## Usage 

Documentation (*MATLAB* tutorial) is available at https://ribokit.github.io/Biers/.

## License

Copyright &copy; of **Biers** _Source Code_ is described in [LICENSE.md](https://github.com/ribokit/Biers/blob/master/LICENSE.md).

<hr/>

Developed by **Das lab**, _Leland Stanford Junior University_.

README by [**t47**](https://t47.io/), *April 2016*. Updated, Rhiju Das, 2017.
