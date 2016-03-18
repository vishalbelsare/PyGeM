# PyGeM [![Build Status](https://travis-ci.org/mathLab/PyGeM.svg)](https://travis-ci.org/mathLab/PyGeM) [![Coverage Status](https://coveralls.io/repos/github/mathLab/PyGeM/badge.svg?branch=master)](https://coveralls.io/github/mathLab/PyGeM?branch=master)
Python Geometrical Morphing.

![Python Geometrical Morphing](docs/source/_static/logo_PyGeM_small.png)

## Description

**PyGeM** is a python library using **Free Form Deformation** to parametrize and morph complex geometries.  It is ideally suited for actual industrial problems, since it allows to handle:

- Computer Aided Design files (in .iges and .stl formats)
- Mesh files (in .unv and OpenFOAM formats)

By now, it has been used with meshes with up to 14 milions of cells. Try with more and more complicated input files!

Here two applications are shown, taken from the **naval** and **automotive** engineering fields. On the other hand, the provided tutorials are related to easier geometries.

![DTMB Morphing](readme/DTMB_ffd.png)
*DTMB-5415 hull: morphing of the bulbous bow starting from an industrial .iges CAD file.*

{{< figure src="readme/DTMB_ffd.png" title="Steve Francia" >}}

![DrivAer Morphing](readme/drivAer_ffd.png)
*DrivAer model: morphing of the bumper starting from an OpenFOAM mesh file.*



If you find this collection useful, feel free to download, use it and suggest pull requests!

The official distribution is on GitHub, and you can clone the repository using

	git clone https://github.com/mathLab/PyGeM


## Documentation

**PyGeM** uses [Sphinx](http://www.sphinx-doc.org/en/stable/) for code documentation. To build the html versions of the docs simply:

```bash
> cd docs
> make html
```

The generated html can be found in `docs/build/html`. Open up the `index.html` you find there to browse

## Testing

We are using Travis CI for continuous intergration testing. You can check out the current status [here](https://travis-ci.org/mathLab/PyGeM).

To run tests locally:

```bash
> python test.py
```


## License

See the [LICENSE](LICENSE.rst) file for license rights and limitations (MIT).
