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
git clone https://github.com/DasLab/Biers.git
```

- In *MATLAB*, go to "**Set Path**". Then "**Add with Subfolders**" of the target `path/to/Biers/Scripts/`.

- Copy `Scripts/get_exe_dir.m.example` and `Scripts/get_varna.m.example` into `Scripts/get_exe_dir.m` and `Scripts/get_varna.m`. Edit `Scripts/get_exe_dir.m` and `Scripts/get_varna.m` following the instructions in these files.

## Usage 

Documentation (*MATLAB* tutorial) is available at https://daslab.github.io/Biers/.

## License

Copyright &copy; of **Biers** _Source Code_ is described in [LICENSE.md](https://github.com/DasLab/Biers/blob/master/LICENSE.md).

<br/>
Developed by **Das lab**, _Leland Stanford Junior University_.
<br/>
README by [**t47**](http://t47.io/), *April 2016*. Updated, Rhiju Das, 2017.
