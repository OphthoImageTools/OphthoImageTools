# Fundus Depigmentation in VKH (Optos UWF)

ImageJ macro for **longitudinal quantification of fundus depigmentation** on ultra-widefield Optos imaging in Vogt&ndash;Koyanagi&ndash;Harada (VKH) disease. The macro defines a single elliptical ROI centered on the optic nerve head, isolates hypopigmented pixels with a fixed threshold, and measures their area &mdash; a comparable, longitudinally-trackable readout of disease activity.

## File

| File | Role |
|---|---|
| `VKH_Depigmentation_Optos.ijm` | Asks the user to click the **center of the optic nerve head** with the Point Tool, then to enter the elliptical ROI **width** and **height** (px), draws the ellipse, crops to it, converts to 8-bit *6 shades*, applies a fixed threshold (`80&ndash;255`), and runs *Analyze Particles* with summary output. |

## Workflow

1. Open the Optos UWF color image in Fiji (drag-and-drop or *File &rarr; Open*).
2. Open the macro (*Plugins &rarr; Macros &rarr; Edit...*) and **Run**.
3. A dialog says *"Select the center of the optic nerve (click on point tool) and press OK."* Activate the **Point Tool**, click on the **center of the optic disc**, then press OK.
4. Two number dialogs appear:
   - *Width* &mdash; horizontal axis (px) of the elliptical ROI to analyze. Default `100`.
   - *Height* &mdash; vertical axis (px). Default `50`.
5. The macro draws the ellipse centered on the optic nerve, crops to it, converts to **8-bit** with the **6 shades** LUT (compresses tones into 6 quantized bins), thresholds at `80&ndash;255`, and runs *Analyze Particles* with `summarize add`.
6. Read the **total depigmented area** from the *Summary* table.

## Why an ellipse around the optic nerve?

Depigmentation in VKH is most reproducibly tracked **outside** the macular pigment-rich zone, in the perpapillary area where the choroidal melanocyte loss is visually obvious on Optos UWF. Centering the ROI on a fixed anatomical landmark (the optic nerve head) makes the measurement **repeatable across visits** and removes ambiguity in field placement. The ellipse axes (width / height) let you adapt to gaze and image framing.

## Parameters worth tuning

- **Threshold (`80, 255`)** &mdash; calibrated for full-resolution Optos color images. If you change the input modality (Spectralis UWF, Mirante, scaled-down JPEGs), re-validate the threshold against a few control images.
- **6 shades LUT** &mdash; intentional binning that makes the threshold robust to small luminance differences across visits. Replace with another LUT if you want finer gradation.
- **Ellipse axes** &mdash; default `100 &times; 50` px. Adjust to your image scale; use le same axes at follow-up to keep the ROI comparable.

## Tips for longitudinal use

- Always re-use the **same threshold** and the **same ellipse axes** across visits for the same patient.
- Keep the optic-nerve click as central as possible; use the Optos auto-color JPEG (not the green-channel-only export) for consistent intensity.
- Track the *Total Area* (px&sup2;) over time; convert to mm&sup2; if you have a calibrated scale bar.

## Reference

Pichi F, Belletti M, Neri P, Carre&ntilde;o E. *Quantitative Ultra-Widefield Imaging Measurement of Fundus Depigmentation in Vogt&ndash;Koyanagi&ndash;Harada Disease.* **Ocul Immunol Inflamm.** 2026;34(2):334&ndash;340.
