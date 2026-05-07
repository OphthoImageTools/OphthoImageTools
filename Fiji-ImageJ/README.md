# Fiji / ImageJ macros

ImageJ macro language (`.ijm`) tools for retinal &amp; uveitis imaging.

All macros are tested on **Fiji** (the recommended scientific distribution of ImageJ &mdash; download from https://imagej.net/software/fiji/).

## Modules

| Folder | What it does |
|---|---|
| [`ETDRS_Grid/`](ETDRS_Grid) | Automatic ETDRS-style 4-quadrant grid overlay with particle counting per sector. |
| [`Macrophage_Like_Cells_Analysis/`](Macrophage_Like_Cells_Analysis) | One-click full pipeline: registration of OCTA/OCTR stacks, DoG + FFT denoising, MLC segmentation, 9-region ETDRS grid, NND, CSV output. |
| [`GCIPL_Sectors_Cirrus/`](GCIPL_Sectors_Cirrus) | Cirrus-style elliptical 6-sector overlay around the foveal center for GCIPL particle counts. |
| [`Fundus_Depigmentation_Analysis/`](Fundus_Depigmentation_Analysis) | Quantitative ultra-widefield (Optos) measurement of fundus depigmentation in VKH. |
| [`OCTA_Vascular_Analysis/`](OCTA_Vascular_Analysis) | FAZ segmentation (Level Sets) + perifoveal / parafoveal vascular density on en face OCTA. |

Open the README inside each subfolder for input format, parameters, and step-by-step instructions.

## Installing a macro in Fiji

1. Open **Fiji**.
2. Go to *Plugins &rarr; Macros &rarr; Edit...* and open the `.ijm` file, **or** drag-and-drop the file onto the Fiji status bar.
3. Click **Run** (Ctrl/Cmd + R).
4. To install permanently as a menu item (and bind hotkeys like F1/F2/F3 used by `OCTA_Vascular_Analysis`): *Plugins &rarr; Macros &rarr; Install...* and select the file.

## Required Fiji plugins

Some macros depend on plugins that are **not bundled** with the default Fiji distribution. Install them once via *Help &rarr; Update... &rarr; Manage update sites*:

- **bUnwarpJ** &mdash; for elastic registration in `Register Virtual Stack Slices` (used by the MLC pipeline).
- **AvgNoiseRmvr** (Average Noise Remover) &mdash; used by the MLC pipeline for noise reduction on registered OCTA stacks.
- **Level Sets** &mdash; for Active-Contours FAZ segmentation in `OCTA_Vascular_Analysis` (a grey-threshold fallback exists if you cannot install it).

If a macro fails with `unknown command`, the missing plugin is almost always the cause.
