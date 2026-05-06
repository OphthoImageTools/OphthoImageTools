# Fundus Depigmentation in VKH (Optos UWF)

ImageJ macro for **longitudinal quantification of fundus depigmentation** on ultra-widefield Optos imaging in Vogt&ndash;Koyanagi&ndash;Harada (VKH) disease.

## Status

The macro file will be added shortly. If you are reading this in the meantime, the published method (Pichi 2026) describes the workflow:

1. Export Optos UWF color images at full resolution (2600&times;2048 px or higher).
2. Convert to grayscale (or extract the green channel) to maximize depigmentation contrast.
3. Apply a fixed intensity threshold to isolate hypopigmented areas.
4. Use *Analyze Particles* with size and circularity thresholds tuned to ignore vessels and optic nerve head.
5. Compute the **percentage of fundus area** above threshold relative to the segmented retinal field.
6. Repeat at follow-up visits with the **same threshold** and field mask &rarr; longitudinal trend.

## Reference

Pichi F, Belletti M, Neri P, Carre&ntilde;o E. *Quantitative Ultra-Widefield Imaging Measurement of Fundus Depigmentation in Vogt&ndash;Koyanagi&ndash;Harada Disease.* **Ocul Immunol Inflamm.** 2026;34(2):334&ndash;340.

## To do

- [ ] Upload the `.ijm` macro file.
- [ ] Add a sample Optos image (with patient consent) and an example output overlay.
- [ ] Document threshold defaults and any preprocessing steps.
