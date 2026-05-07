<div align="center">

# OphthoImageTools

### Open ImageJ &amp; MATLAB tools for objective ophthalmic imaging analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Made with ImageJ](https://img.shields.io/badge/ImageJ-Fiji-orange)](https://imagej.net/software/fiji/)
[![Made with MATLAB](https://img.shields.io/badge/MATLAB-R2020b%2B-red)](https://www.mathworks.com/)
[![Buy Me A Coffee](https://img.shields.io/badge/donate-buymeacoffee-yellow.svg?logo=buy-me-a-coffee)](https://www.buymeacoffee.com/ophthoimagetools)
[![Editorial: AJO](https://img.shields.io/badge/editorial-AJO%202026-success)](#-citation)

*A collaborative project by clinician-scientists in uveitis and retinal imaging.*<br>
*From proprietary algorithms to open science &mdash; rethinking objectivity in uveitis.*

</div>

---

## Why this repository exists

For more than two decades, uveitis grading has relied on the SUN classification &mdash; a categorical, semi-quantitative system that was built to standardize communication, not to serve as a quantitative endpoint for modern clinical trials. Many groups have built sophisticated computational tools to fix this, but most have remained inside the labs that created them.

**OphthoImageTools** is our attempt to do the opposite: collect simple, reproducible, openly-described image-analysis workflows that any center can download, run, criticize, and improve. The macros here turn descriptive features (keratic precipitates, fundus depigmentation, macrophage-like cells, ganglion-cell topography, FAZ and vascular density) into quantitative readouts using free, widely-available platforms (ImageJ/Fiji, MATLAB).

The companion editorial in the *American Journal of Ophthalmology* explains the rationale in detail.

---

## What's inside

| Module | Platform | Purpose | Reference |
|---|---|---|---|
| [ETDRS Grid](Fiji-ImageJ/ETDRS_Grid) | ImageJ | Automatic 4-quadrant ETDRS-style overlay with particle counting per sector | Belletti 2025 |
| [Macrophage-Like Cells (MLC)](Fiji-ImageJ/Macrophage_Like_Cells_Analysis) | ImageJ | One-click full pipeline: registration, DoG/FFT denoising, MLC segmentation, 9-region ETDRS grid, NND | Pichi 2024 |
| [GCIPL 6-Sectors (Cirrus-like)](Fiji-ImageJ/GCIPL_Sectors_Cirrus) | ImageJ | Cirrus-style elliptical 6-sector overlay around the foveal center for GCIPL particle counts | Belletti 2025 |
| [Fundus Depigmentation (VKH)](Fiji-ImageJ/Fundus_Depigmentation_Analysis) | ImageJ | Quantitative ultra-widefield (Optos) measurement of fundus depigmentation in VKH | Pichi 2026 |
| [OCTA Vascular Density &amp; FAZ](Fiji-ImageJ/OCTA_Vascular_Analysis) | ImageJ | Automatic FAZ segmentation (Level Sets) + perifoveal &amp; parafoveal vascular density on en face OCTA | &mdash; |
| [Keratic Precipitates (KPs)](MATLAB/Keratic_Precipitates_Analysis) | MATLAB | AS-OCT segmentation, en face projection and volumetric quantification of KPs from Heidelberg XML exports | Pichi 2023 |
| [PubMed n8n workflow](n8n_PubMed_Workflow) | n8n | Automated PubMed search &amp; digest pipeline (importable JSON, ready to customize) | &mdash; |

Each subfolder contains a dedicated **README** with setup, dependencies, step-by-step usage and expected output.

---

## Quick start

1. **Clone or download** the repository:
   ```bash
   git clone https://github.com/OphthoImageTools/OphthoImageTools.git
   ```
   *or* click the green **Code** button above &rarr; **Download ZIP**.
2. Open the folder for the analysis you need.
3. Read its **README.md** &mdash; every macro lists the input format, parameters, and validated workflow.
4. Open the `.ijm` file in **Fiji** (Plugins &rarr; Macros &rarr; Run...) or the `.m` file in **MATLAB**.
5. Run on your dataset.

> [!TIP]
> All ImageJ macros are tested on **Fiji** (recommended distribution of ImageJ). MATLAB scripts are tested on R2020b and later, with the *Image Processing Toolbox*.

---

## Authors

|     | Affiliation |
|---|---|
| **Francesco Pichi, MD** &middot; ORCID&nbsp;[0000-0002-7357-4166](https://orcid.org/0000-0002-7357-4166) | University of Toronto, Department of Ophthalmology and Vision Science &middot; Kensington Eye Institute, Toronto, Canada |
| **Ester Carre&ntilde;o, MD, PhD** | Department of Ophthalmology, Hospital Universitario La Paz, Madrid, Spain |

Contributions, issues and pull requests are welcome &mdash; see [CONTRIBUTING](#contributing) below.

---

## Citation

If you use any of these tools in your work, please cite both the repository and the relevant primary publication.

**Repository (general use):**
> Pichi F, Carre&ntilde;o E. *OphthoImageTools: ImageJ &amp; MATLAB tools for ophthalmic imaging analysis.* GitHub repository, 2025. Available at: https://github.com/OphthoImageTools/OphthoImageTools

A machine-readable citation file ([`CITATION.cff`](CITATION.cff)) is included &mdash; GitHub will offer a one-click *"Cite this repository"* button on the right sidebar.

**Companion editorial:**
> Pichi F, Carre&ntilde;o E. *From Proprietary Algorithms to Open Science: Rethinking Objectivity in Uveitis.* American Journal of Ophthalmology, 2026.

**Primary publications per macro** are listed in each subfolder's README.

---

## Vision

We believe in **open science**. Our goal is not commercial, but to **promote reproducibility and shared standards** in ophthalmic image analysis. Methods documented well, validated across centers and devices, and available for collective review are how durable knowledge is built today &mdash; especially in a small-market field like uveitis where industry alone is unlikely to deliver the objective metrics we need.

If enough centers adopt the same tools, consensus may form among researchers first, and devices may follow.

---

## Contributing

We welcome:
- **Bug reports** and feature requests &mdash; please open an [Issue](https://github.com/OphthoImageTools/OphthoImageTools/issues).
- **New macros** for related workflows (other modalities, other diseases) &mdash; submit a Pull Request.
- **Validation data** from independent centers &mdash; reach out via Issues or email.

When contributing, please:
- Keep paths relative o parameterized (no hard-coded local paths).
- Add a short docstring at the top of each script (purpose, input, output, dependencies).
- Update the relevant README.

---

## Support the project

If these tools have been useful for your research, you can support continued development and validation here:

<a href="https://www.buymeacoffee.com/ophthoimagetools" target="_blank"><img src="https://cdn.buymeacoffee.com/buttons/v2/default-yellow.png" alt="Buy Me A Coffee" height="50"></a>

> The Buy-Me-A-Coffee handle above is currently a placeholder &mdash; it will be activated shortly. Stars &#11088; on this repository are equally appreciated.

---

## License

Released under the **MIT License** &mdash; see [LICENSE](LICENSE). You are free to use, modify, and redistribute the code, provided proper credit is given.

---

<div align="center">
<sub>Built and maintained by Francesco Pichi &amp; Ester Carre&ntilde;o &middot; Toronto &middot; Madrid</sub>
</div>
