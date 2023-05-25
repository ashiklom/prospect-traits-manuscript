# Chloroplasts to canopies: Leaf optical properties shed light on foliar trait variability at individual to global scales

Source code and ancillary files for the manuscript,
"Chloroplasts to canopies: Leaf optical properties shed light on foliar trait variability at individual to global scales"
by Alexey N. Shiklomanov and colleagues.

The manuscript is written using the [Quarto](https://quarto.org/) publishing system.
To generate PDF, HTML, and DOCX versions of the manuscript, run `quarto render` in the repository root directory;
the resulting files will be created in the `_output` directory.

## Dataset organization

- `spectra_id` -- Unique value per reflectance measurement. Associated with a single `observation_id`

- `observation_id` -- Unique value per trait. Connected to one or more `spectra_id`.
  - RTM inversions are applied to a single `observation_id`.

- `entity_id` -- Unique value per leaf. Linked to one or more `observation_id`.
