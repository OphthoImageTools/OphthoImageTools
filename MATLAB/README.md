# MATLAB scripts

MATLAB pipelines for ophthalmic imaging analyses that go beyond what is comfortable to do in the ImageJ macro language &mdash; mostly **3D AS-OCT reconstruction and volumetry**.

## Modules

| Folder | What it does |
|---|---|
| [`Keratic_Precipitates_Analysis/`](Keratic_Precipitates_Analysis) | Segments the cornea on Heidelberg AS-OCT volume scans, projects deposits *en face*, fits an affine transform onto the localizer, and computes **per-eye precipitate volume (mm&sup3;)**, internal corneal area (mm&sup2;), and their ratio. |

## Requirements

- MATLAB **R2020b** or later.
- **Image Processing Toolbox** (uses `fitgeotrans`, `imwarp`, `alphaShape`, `surfaceArea`).
- Heidelberg Spectralis or Anterion AS-OCT volume scans **exported as XML + BMP** (the standard Heidelberg `.xml` raw export).

See the README in each subfolder for the expected input layout and step-by-step instructions.
