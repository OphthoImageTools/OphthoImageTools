# Keratic Precipitate (KP) volumetry on AS-OCT

MATLAB pipeline that segments the cornea on Heidelberg anterior-segment OCT volume scans, identifies keratic precipitates protruding from the endothelium, projects them *en face*, and computes their **volume in mm&sup3;** plus the **internal corneal area** they sit on (mm&sup2;).

This converts a feature that has historically been graded qualitatively (1+ to 4+, or by descriptive type) into a continuous, comparable, longitudinal readout.

## Files (run order)

| File | Role |
|---|---|
| `batchGetCorneaFromXML.m` | Parses the Heidelberg `.xml` export, loads each B-scan, segments the **external** and **internal** corneal boundaries, identifies precipitates as elevations above the endothelium, and stores everything in a single `.mat` file per eye (`cornea.mat`). |
| `corneaProcessingAndAnalysis.m` | Runs `batchGetCorneaFromXML` on a single eye-folder, builds the affine transform from B-scan space to the SLO localizer, computes the **internal corneal area** (mm&sup2;), the **precipitate volume** (mm&sup3;), and the ratio. Writes `precipitatesVolume.mat` and shows the colored elevation map overlaid on the localizer. |
| `batchCorneaProcessingAndAnalysis.m` | Runs the single-eye pipeline on **every subfolder** of a parent directory, producing one PNG per eye and a CSV results file (`results<date>.csv`) with columns `ID`, `precipitates_volume_mm3`, `area_mm2`, `ratio_PrecVol_area_mm`. |
| `batchPlotCorneaPrecipitates.m` | Re-renders the elevation map on the SLO localizer for any folder that has already been processed (i.e. has both `cornea.mat` and `precipitatesVolume.mat`). Useful for re-generating figures without re-running the segmentation. |

## Expected input layout

```
ParentDir/
  EyeFolder_001/
    export.xml             <-- Heidelberg XML
    *.bmp                  <-- B-scan images referenced by the XML
    localizer.jpg          <-- SLO localizer referenced by imagesTags(1)
  EyeFolder_002/
    ...
```

Each eye folder must contain **exactly one** `.xml` file (the script will refuse to process otherwise).

## Running

**Single eye:**
```matlab
[precVol, area, ratio] = corneaProcessingAndAnalysis('/path/to/EyeFolder_001');
```
Outputs:
- `cornea.mat` &mdash; segmented boundaries and en face projection.
- `precipitatesVolume.mat` &mdash; `precVol` (mm&sup3;), `area` (mm&sup2;), `ratio_PrecVol_area` (mm).
- A figure showing the elevation map (jet colormap, 0&ndash;0.14 mm) overlaid on the localizer con la scanned ROI outlined in green.

**Batch over a study cohort:**
```matlab
batchCorneaProcessingAndAnalysis('/path/to/ParentDir');
```
Produces:
- `results<DD-MMM-YYYY>.csv` in the parent directory.
- One PNG per eye, named after the folder ID.
- Per-eye `precipitatesVolume.mat` files (so re-running the script will skip already-processed eyes &mdash; pass `overwrite=true` in `corneaProcessingAndAnalysis` to force).

**Re-plot only:**
```matlab
batchPlotCorneaPrecipitates('/path/to/ParentDir');
```

## Method notes

- The analysis region is restricted to the **upper half** of each B-scan (`distanceFromTop = 0.5`) to focus on endothelial precipitates and exclude stromal artifacts.
- Precipitate height is converted to physical units using `scaleZ_3D` from the XML acquisition context.
- The colormap caps at **0.14 mm** with a small transparency offset (`zero_off = 0.02`) to keep low-elevation regions visible without saturating.
- The volume estimate assumes pixels are square in (x, y); spacing in the slow scan direction is computed from `y_span / (N-1)` rather than assumed.

## Reference

Pichi F, Ometto G, Invernizzi A, et&nbsp;al. *Automated quantification of uveitic keratic precipitates by use of anterior segment optical coherence tomography.* **Clin Exp Ophthalmol.** 2023;51(8):790&ndash;798.
