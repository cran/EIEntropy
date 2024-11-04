# EIEntropy 0.0.1.3 (2024-11-03)

### Bug Fixes

-Fixed an issue with the weights argument where, in some cases, custom weight columns were not recognized, causing equal weights to be assigned to all observations. Now, users can specify any column name as a string for weights in datahp and datahs, with the function correctly applying the specified weights for improved flexibility and usability.

# EIEntropy 0.0.1.2 (2024-11-02)

## Main changes

This update introduces significant improvements in the flexibility, accuracy, and functionality of the package. Below are the primary modifications:

-   **Correction in handling the dependent variable (`y`)**: An error in processing the dependent variable led to inaccuracies in estimations. This issue has been resolved, ensuring reliable and consistent results.

-   **Extended compatibility for the `y` variable**: The `y` variable can now be a categorical variable with multiple levels (`J > 2`). When `y` is a dummy variable, the first column represents level 1. For categorical variables, levels are ordered alphabetically.

-   **Correction in handling the error prior (`l`)**: Previously, the value `l` did not include index 3. This version corrects that omission.

-   **Enhancements in `q` assignment**: `q` can now be specified as a vector rather than only a single value per category. The function now correctly handles uniform distributions regardless of the number of categories.

-   **New default tolerance level**: The function tolerance is set to `1e-10` for improved optimization accuracy.

-   **Change in optimization method to `nlminb`**: This version uses `nlminb` instead of `optim`, improving results and eliminating the need for the `method` argument.
