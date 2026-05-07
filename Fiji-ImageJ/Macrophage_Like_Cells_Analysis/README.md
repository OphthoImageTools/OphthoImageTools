# Macrophage-Like Cells (MLC) on en face OCTA &mdash; full pipeline

Single-macro pipeline that goes from **raw OCTA / OCTR frames** to a **per-region MLC count on an ETDRS grid**, with a **global nearest-neighbor distance (NND)** estimate. Replaces the older two-step macro (`Codecomplete new.ijm` + manual ETDRS step) with a unified, parameterizable workflow.

## File

| File | Role |
|---|---|
| `MLC_Full_Pipeline.ijm` | One-click pipeline. Opens a parameter dialog, registers the OCTA stack with bUnwarpJ, applies the same transform to OCTR, builds an averaged image, denoises with **DoG (Difference of Gaussians) + FFT**, thresholds with **Triangle**, runs *Analyze Particles*, lets the user click the foveal center, builds a **9-region ETDRS grid (OD/OS aware)**, computes per-region density and global NND, and saves results (CSV + figures). |

## Required Fiji plugins

- **bUnwarpJ** &mdash; for `Register Virtual Stack Slices` (*Elastic / bUnwarpJ splines*).
- **AvgNoiseRmvr** (Average Noise Remover).

Install via *Help &rarr; Update... &rarr; Manage update sites*. If a step fails with `unknown command`, a missing plugin is almost always the cause.

## Folder layout (before running)

The pipeline expects this structure for each eye, all under one **base folder**:

```
EyeFolder/
  OCTA/                    <-- raw en face OCTA frames (registration source)
  OCTR/                    <-- corresponding OCTR frames
  OutputOCTA/              <-- empty, will be filled with registered slices
  OutputOCTR/              <-- empty, will be filled with transformed slices
  Transform/               <-- empty, registration transforms saved here
  Average/                 <-- empty, average images and outputs saved here
```

All input frames should be **500&times;500 px**, covering a **6&times;6 mm** macular cube. The defaults in the dialog assume this; change *Tama&ntilde;o del campo* and *Resoluci&oacute;n* if your cube differs.

## Running

1. Open Fiji and load the macro: *Plugins &rarr; Macros &rarr; Edit...* &rarr; open `MLC_Full_Pipeline.ijm` &rarr; **Run**.
2. The **parameter dialog** opens. Set:
   - *Carpeta base* &mdash; the EyeFolder path (the macro builds the OCTA/OCTR/Transform/Average subpaths from here).
   - *Ojo* &mdash; OD or OS (drives temporal/nasal sector orientation in the ETDRS grid).
   - *Tama&ntilde;o del campo* (mm) and *Resoluci&oacute;n* (px) &mdash; defaults 6 mm / 500 px.
   - DoG sigmas, Triangle pre-blur, *Analyze Particles* size and circularity bounds.
3. The macro runs steps 2&ndash;7 automatically (registration, transform, averaging, DoG, FFT, threshold, particle detection).
4. A dialog asks you to **click the foveal center** with the Point Tool. Click once and press **OK**.
5. The macro builds the 9-region ETDRS grid (`C` central + 4 inner `IS/IT/II/IN` + 4 outer `OS/OT/OI/ON`), computes per-region counts and densities, then computes global NND.
6. Outputs saved in `Average/`:
   - `averageOCTA.png`, `averageOCTR.png`, `averageOCTRclean.png`
   - `MLCs_mask.png` &mdash; binarized MLC map.
   - ETDRS overlay PNG.
   - Per-region statistics + NND table (CSV).

## Key parameters worth tuning

- **DoG sigmas** (default 1 / 3) &mdash; suppress slow background variation. Increase sigma 2 if the macula has a strong vignette.
- **Triangle pre-blur** (default 1) &mdash; stabilizes Triangle thresholding when the histogram has a single sharp peak.
- **Particle size** (default 4&ndash;35 px&sup2;) &mdash; tune to your magnification.
- **Particle circularity** (default 0.5&ndash;1.0) &mdash; lower if vessels are being mis-segmented as elongated objects.
- **ETDRS radii** &mdash; `0.5 / 1.5 / 3.0 mm` &mdash; hard-coded to match standard ETDRS. Edit the macro if you need a different sub-foveal grid.

## What's different from the older two-step macro

| Aspect | Old (`Codecomplete new.ijm` + ETDRS step 2) | New (`MLC_Full_Pipeline.ijm`) |
|---|---|---|
| Threshold | Fixed `21&ndash;255` | **Triangle** (auto, dataset-adaptive) |
| Denoising | FFT only | **DoG + FFT** (handles slow background) |
| Sectorization | 4-quadrant grid drawn manually | **9-region ETDRS** (OD/OS aware), built automatically from a single fovea click |
| Output | ROI manager counts | CSV with per-region density (cells/mm&sup2;) + global NND |
| Parameters | Hard-coded paths | Dialog box for paths and all thresholds |

## References

- Pichi F, Neri P, Aljneibi S, Hay S, Chaudhry H, Carre&ntilde;o E. *Vitreoretinal Interface Cells Correlate In Vivo With Uveitis Activity and Decrease With Anti-Inflammatory Treatment.* **Transl Vis Sci Technol.** 2024;13(5):15.
- Pichi F, Neri P, Aljeneibi S, et&nbsp;al. *In Vivo Visualization of Macrophage-Like Cells in Patients with Uveitis by Use of En Face Swept Source Optical Coherence Tomography.* **Ocul Immunol Inflamm.** 2024;32(8):1532&ndash;1538.
