# OCTA Vascular Density &amp; FAZ Analyzer

ImageJ macro that, on a single en face OCTA slab, measures:

1. **Foveal avascular zone (FAZ)** &mdash; area, perimeter, and circularity (mm&sup2;, mm, dimensionless).
2. **Vascular density** &mdash; whole image, perifoveal central disc, and parafoveal ring (%).

Designed for routine perifoveal analysis without manual ROI drawing &mdash; you click the foveal center once, the macro does everything else.

## File

| File | Role |
|---|---|
| `OCTA_Vascular_Analysis.ijm` | Single-image or batch-mode pipeline. Background-subtraction (rolling ball), CLAHE, median filter; vessel segmentation with `setAutoThreshold("Otsu dark") + Open`; FAZ segmentation with **Level Sets / Active Contours** (primary) or grey-value threshold (fallback). Results to a persistent table. |

## Hotkeys (after install)

| Key | Action |
|---|---|
| **F1** | Run analysis on the active image (or pick a batch folder). |
| **F2** | Open the parameter dialog (scale, radii, FAZ method, thresholds). |
| **F3** | Export the *OCTA Results* table to CSV. |

To install permanently as menu items: *Plugins &rarr; Macros &rarr; Install...* &rarr; pick `OCTA_Vascular_Analysis.ijm`. Otherwise just *Run* once per session.

## Required Fiji plugins

- **Level Sets** &mdash; for Active-Contours FAZ segmentation. Install via *Help &rarr; Update... &rarr; Manage update sites &rarr; check Level Sets*.
- If you cannot install Level Sets, set `USAR_LEVEL_SETS = false` (F2 dialog) &mdash; the macro falls back to a grey-value threshold around the fovea.

The rest is **native ImageJ** &mdash; no MorphoLibJ dependency.

## Default scale and assumptions

- **6 mm cube, 500 px image** &rarr; `12 &micro;m / px`.
  - Change `ESCALA_UM_POR_PIXEL` (F2) if your cube is e.g. 3 mm or 8 mm.
- **Perifoveal disc**: 500 &micro;m radius (central disc).
- **Parafoveal ring**: 500&ndash;1500 &micro;m (annulus).
- **FAZ search radius**: 1500 &micro;m (large enough for ischemic FAZ in pathological eyes).
- **FAZ size limits**: 0.005&ndash;4.0 mm&sup2;.

## Workflow

1. Open the en face OCTA slab in Fiji (or in batch mode, pick a folder &mdash; the macro iterates over `.tif`, `.png`, `.jpg`, `.bmp`).
2. Press **F1** (or *Plugins &rarr; Macros &rarr; Run...* &rarr; "OCTA Vascular Density &amp; FAZ").
3. The macro:
   - Converts to 8-bit if needed.
   - Asks you to **click the foveal center** with the Point Tool (this defines the perifoveal/parafoveal grid and the FAZ search center).
   - Runs *Subtract Background* (rolling ball) &rarr; *CLAHE* &rarr; *Median*.
   - Segments vessels: `Otsu dark` &rarr; *Open* (no Fill Holes &mdash; that would erase the FAZ).
   - Computes density: total, central disc, parafoveal annulus.
   - Segments the FAZ with **Level Sets** &mdash; pick the connected region closest to the foveal click (avoids confusion with peripheral hypoperfusion).
   - Adds a row to the *OCTA Results* table with all metrics.
4. Press **F3** to export the table to CSV when you are done.

## Key parameters (F2)

| Parameter | Default | What it controls |
|---|---|---|
| `ESCALA_UM_POR_PIXEL` | `12.0` | Pixel scale (6 mm / 500 px). Critical &mdash; everything else uses it. |
| `RADIO_CENTRAL_UM` | `500` | Radius of the perifoveal central disc. |
| `RADIO_PARAFOVEAL_UM` | `1500` | Outer radius of the parafoveal ring. |
| `FAZ_RADIO_BUSQUEDA_UM` | `1500` | Region considered for FAZ detection (centered on the foveal click). |
| `FAZ_MIN_AREA_MM2` | `0.005` | Reject specks smaller than this. |
| `FAZ_MAX_AREA_MM2` | `4.0` | Allow ischemic FAZ as large as this. |
| `METODO_UMBRAL` | `Otsu` | Vessel-segmentation auto-threshold. |
| `ROLLING_BALL_PX` | `20` | Background subtraction radius. |
| `USAR_LEVEL_SETS` | `true` | If `false`, use grey-threshold fallback. |
| `LS_GREY_THRESHOLD` &hellip; | various | Tune Level Sets if FAZ contour leaks. |

## Method notes

- **No Fill Holes** on the vessel mask. That would close the FAZ and inflate density to ~100&#37;.
- The FAZ is selected as the connected region **closest to the foveal click**, not the largest &mdash; this avoids picking up peripheral hypoperfusion in ischemic eyes.
- The parafoveal density is computed as `(vasos_in_outer_disc - vasos_in_central_disc) / (annulus_area)`, so the central disc never double-counts.
- Densities are reported as **percent of pixels above the vessel threshold** within each ROI.

## Output (one row per image in *OCTA Results*)

| Column | Units |
|---|---|
| Imagen | filename |
| F&oacute;vea X / F&oacute;vea Y | px |
| Escala | &micro;m/px |
| Densidad Total | % |
| Densidad Central | % |
| Densidad Parafoveal | % |
| FAZ M&eacute;todo | "Level Sets" or "Umbral gris" |
| FAZ Encontrada | bool |
| FAZ &Aacute;rea | mm&sup2; |
| FAZ Per&iacute;metro | mm |
| FAZ Circularidad | dimensionless |

## Validation suggestions

- Compare against the device-native FAZ measurement (Cirrus, Spectralis, Heidelberg) on a small validation set.
- Run the macro twice on the same image with slightly different foveal clicks to assess intra-grader variability.
- For longitudinal use, **fix the scale** (`ESCALA_UM_POR_PIXEL`) and the foveal click strategy across visits.
