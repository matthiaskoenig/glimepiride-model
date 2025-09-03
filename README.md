[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14628700.svg)](https://doi.org/10.5281/zenodo.14628700)

# Glimepiride Model
This repository provides the glimepiride physiologically based pharmacokinetics (PBPK) model.

The model is distributed as [SBML](http://sbml.org) available from [`glimepiride_body_flat.xml`](./models/glimepiride_body_flat.xml) with 
corresponding SBML4humans model report at [https://sbml4humans.de/model_url?url=https://raw.githubusercontent.com/matthiaskoenig/glimepiride-model/main/models/glimepiride_body_flat.xml](https://sbml4humans.de/model_url?url=https://raw.githubusercontent.com/matthiaskoenig/glimepiride-model/main/models/glimepiride_body_flat.xml) and equations from [`glimepiride_body_flat.md`](./models/glimepiride_body_flat.md).

The COMBINE archive is available from [`glimepiride_model.omex`](./glimepiride_model.omex).

![model overview](./figures/glimepiride_model.png)

### Comp submodels
The liver submodel is available from [`glimepiride_liver.xml`](./models/glimepiride_liver.xml) with corresponding SBML4humans report at
[https://sbml4humans.de/model_url?url=https://raw.githubusercontent.com/matthiaskoenig/glimepiride-model/main/models/glimepiride_liver.xml](https://sbml4humans.de/model_url?url=https://raw.githubusercontent.com/matthiaskoenig/glimepiride-model/main/models/glimepiride_liver.xml) and equations from [`glimepiride_liver.md`](./models/glimepiride_liver.md).

The kidney submodel is available from [`glimepiride_kidney.xml`](./models/glimepiride_kidney.xml) with corresponding SBML4humans report at
[https://sbml4humans.de/model_url?url=https://raw.githubusercontent.com/matthiaskoenig/glimepiride-model/main/models/glimepiride_kidney.xml](https://sbml4humans.de/model_url?url=https://raw.githubusercontent.com/matthiaskoenig/glimepiride-model/main/models/glimepiride_kidney.xml) and equations from [`glimepiride_kidney.md`](./models/glimepiride_kidney.md).

The intestine submodel is available from [`glimepiride_intestine.xml`](./models/glimepiride_intestine.xml) with corresponding SBML4humans report at
[https://sbml4humans.de/model_url?url=https://raw.githubusercontent.com/matthiaskoenig/glimepiride-model/main/models/glimepiride_intestine.xml](https://sbml4humans.de/model_url?url=https://raw.githubusercontent.com/matthiaskoenig/glimepiride-model/main/models/glimepiride_intestine.xml) and equations from [`glimepiride_intestine.md`](./models/glimepiride_intestine.md).

The whole-body submodel is available from [`glimepiride_body.xml`](./models/glimepiride_body.xml) with corresponding SBML4humans report at
[https://sbml4humans.de/model_url?url=https://raw.githubusercontent.com/matthiaskoenig/glimepiride-model/main/models/glimepiride_body.xml](https://sbml4humans.de/model_url?url=https://raw.githubusercontent.com/matthiaskoenig/glimepiride-model/main/models/glimepiride_body.xml) and equations from [`glimepiride_body.md`](./models/glimepiride_body.md).

## How to cite
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14628700.svg)](https://doi.org/10.5281/zenodo.14628700)

> Elias, M., & König, M. (2025).
> *Physiologically based pharmacokinetic (PBPK) model of glimepiride.*   
> Zenodo. [https://doi.org/10.5281/zenodo.14628700](https://doi.org/10.5281/zenodo.14628700)

## License

* Source Code: [MIT](https://opensource.org/license/MIT)
* Documentation: [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)
* Models: [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.

## Run simulations
Simulations can be run within a docker container:

```bash
docker run -v "${PWD}/results:/results" -it matthiaskoenig/glimepiride:latest /bin/bash
```

Run the complete analysis:
```bash
run_glimepiride -a all -r /results
```
The results are written into the mounted `/results` folder on the host.

In case of permission issues with the mounted folder, adjust ownership and access rights with:
```bash
sudo chown $(id -u):$(id -g) -R "${PWD}/results"
sudo chmod 775 "${PWD}/results"
```

## Funding
Matthias König was supported by the Federal Ministry of Education and Research (BMBF, Germany) within LiSyM by grant number 031L0054 and ATLAS by grant number 031L0304B and by the German Research Foundation (DFG) within the Research Unit Program FOR 5151 QuaLiPerF (Quantifying Liver Perfusion-Function Relationship in Complex Resection - A Systems Medicine Approach) by grant number 436883643 and by grant number 465194077 (Priority Programme SPP 2311, Subproject SimLivA). This work was supported by the BMBF-funded de.NBI Cloud within the German Network for Bioinformatics Infrastructure (de.NBI) (031A537B, 031A533A, 031A538A, 031A533B, 031A535A, 031A537C, 031A534A, 031A532B). 

© 2024-2025 Michelle Elias & Matthias König, [Systems Medicine of the Liver](https://livermetabolism.com)
