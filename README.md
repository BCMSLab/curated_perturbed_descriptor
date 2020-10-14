# curated_perturbed_descriptor
Code for the Adipocyte paper (Curated gene expression dataset of differentiating 3T3-L1 adipocytes under pharmacological and genetic perturbations)

## Setting up the docker environment

The analysis was run on a [docker](https://hub.docker.com/r/bcmslab/adiporeg_array)
image based on the the latest **rocker/verse**.
Other R packages were added to the image and were made available as an image
that can be obtained and launched on any local machine running
[docker](https://hub.docker.com/r/bcmslab/adiporeg_array).

```bash
$ docker pull bcmslab/adiporeg_array:latest
$ docker run -it bcmslab/adiporeg_array:latest bash
```

## Obtaining the source code

The source code is hosted publicly on this repository in the form of a research
compendium. This includes the scripts to reproduce the figures and tables in
this manuscript. Another public repository
[curated_perturbed_descriptor](https://github.com/BCMSLab/curated_perturbed_descriptor)
contains the scripts to download, prepare and analyze the data.

```bash
$ git clone https://github.com/BCMSLab/curated_perturbed_descriptor
```

## Runing the analysis

In the directory `curated_perturbed_descriptor`, run `make`

```bash
$ cd curated_perturbed_descriptor
$ make all
```

## Details of the R environment
The version of **R** that was used to perform this analysis is the 3.6.1
(2019-07-05) on `x86_64-pc-linux-gnu`.

## More

This manuscript was published under the title
[Curated gene expression dataset of differentiating 3T3-L1 adipocytes under pharmacological and genetic perturbations](https://www.tandfonline.com/doi/full/10.1080/21623945.2020.1829852)
