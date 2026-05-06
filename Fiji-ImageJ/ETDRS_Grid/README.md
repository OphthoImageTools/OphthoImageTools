# Automatic ETDRS Grid &mdash; quadrant counting

ImageJ macro that overlays an **ETDRS-style grid** (concentric circles divided into quadrants) on en face retinal images and counts thresholded particles per sector. Designed for objective regional quantification of macrophage-like cells (MLCs), inflammatory deposits, or any other particle-like feature on en face OCTA / OCT.

## Files

| File | Role |
|---|---|
| `Codecomplete new version ETDRS.ijm` | **Step 1** &mdash; Performs FFT cleanup, thresholds the image, runs *Analyze Particles*, then asks the user to draw a line that defines the ETDRS grid extent and orientation. The macro auto-builds the concentric quadrants and adds them to the ROI Manager. |
| `Codecomplete new version ETDRS part2.ijm` | **Step 2** &mdash; Iterates through the four quadrants stored in the ROI Manager and runs *Analyze Particles* per ROI, producing a per-sector summary. |

> The two-step split is intentional: it lets the user inspect / adjust the grid before counting.

## Workflow

1. Open the en face image (e.g. averaged OCTA slab, 500&times;500 px recommended).
2. Run **Step 1** (`Codecomplete new version ETDRS.ijm`).
   - The macro performs FFT, masks the central cross + DC, applies inverse FFT and a fixed threshold (default `21&ndash;255`), then runs *Analyze Particles* with `size=4&ndash;35, circularity=0.5&ndash;1.0`.
   - Press **F2** when prompted, then draw a line through the foveal center to define the grid axes.
   - The macro builds 1 concentric ring with 4 quadrants (`R1Q1 &hellip; R1Q4`) into the ROI Manager.
3. Run **Step 2** (`Codecomplete new version ETDRS part2.ijm`) to count particles per quadrant.
4. Read the per-quadrant counts from the *Summary* table.

## Parameters you may want to change

Inside `Codecomplete new version ETDRS.ijm`:
- `nCircles` (default `1`) &mdash; number of concentric rings.
- `nQuadrants` (default `4`) &mdash; quadrants per ring.
- `setThreshold(21, 255)` &mdash; threshold range; tune to your image dynamic range.
- *Analyze Particles* `size=4-35`, `circularity=0.5-1.0` &mdash; tune to your target feature.

## Reference

Belletti M, Carre&ntilde;o E, Perez Jimenez Y, Pichi F. *Relationship of macrophage-like cells and retinal ganglion cells in healthy eyes.* **Br J Ophthalmol.** 2025;109(12):1363&ndash;1369.
