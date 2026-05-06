# GCIPL 6-Sector Overlay (Cirrus-like)

ImageJ macro that reproduces the **Cirrus GCIPL elliptical 6-sector grid** (inner foveal ellipse + outer macular ellipse divided into 6 sectors) on an en face image, then runs *Analyze Particles* per sector.

This is the analytic backbone of the Belletti 2025 paper that correlated MLC density with retinal ganglion-cell topography in healthy eyes &mdash; running the count on the **same anatomic grid** that Cirrus uses for GCIPL thickness.

## File

| File | Role |
|---|---|
| `Ganglion cells sectors Cirrus.ijm` | Performs FFT cleanup + threshold + particle detection on the active image, asks the user to click the **foveal center** with the Point Tool, then builds 6 elliptical sectors and counts particles inside each. |

## Workflow

1. Open the en face image (e.g. averaged OCTA slab) in Fiji.
2. Open the macro (*Plugins &rarr; Macros &rarr; Edit...*) and run it.
3. The macro:
   - Runs FFT, masks the central cross / DC, runs inverse FFT.
   - Converts to 8-bit, applies threshold `21&ndash;255`.
   - Runs *Analyze Particles* (`size=4&ndash;35, circularity=0.5&ndash;1.0`) and creates a mask renamed *MLCs*.
4. A dialog asks: *"Select foveal center with Point Tool and press OK"*. Click the **center of the fovea** with the Point Tool, then press OK.
5. The macro builds the elliptical grid:
   - **Outer ellipse**: semi-major `A = 200 px`, semi-minor `B = 167 px`.
   - **Inner foveal ellipse**: `Ai = 50 px`, `Bi = 41.5 px`.
   - 6 sectors, sector boundaries traced in polar coordinates.
   - ROIs added as `GCIPL_S1 &hellip; GCIPL_S6`.
6. The macro iterates through the 6 ROIs and runs *Analyze Particles* per sector. Read the per-sector counts from the *Summary* table.

## Calibration note

The default ellipse axes (`A=200, B=167, Ai=50, Bi=41.5`) match the pixel scaling of a **500&times;500 px** OCTA slab on Cirrus where the macular OCTA cube covers ~6&times;5 mm. **If your input dimensions differ, re-scale these constants proportionally** &mdash; otherwise the sectors will not correspond to true Cirrus boundaries.

## Reference

Belletti M, Carre&ntilde;o E, Perez Jimenez Y, Pichi F. *Relationship of macrophage-like cells and retinal ganglion cells in healthy eyes.* **Br J Ophthalmol.** 2025;109(12):1363&ndash;1369.
