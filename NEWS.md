# itsdm 0.2.0

- Convert Shapley values-based functions to usable by external models (as described in issue # 3), and add examples in function documentation and vignettes to show users how to use these functions.
- Add a function `detect_envi_change` to use Shapley values technique to analyze the potential impacts of changing environmental variables in space.
- Modify function `isotree_po` to take presence-absence dataset as well (as described in issue #7). To make this happen smoothly, another function `format_observation` is created to help the users to convert their data to fit into `itsdm` workflow.
- Reorganized reference section of the online documentation to make it user-friendly.
- Fix a few bugs in the functions.

# itsdm 0.1.3

- Fix a bug in function `print.VariableAnalysis` mentioned in issue #2: if any value is negative then it would fail. 
- As mentioned in issue #3, add a sampling step in function `plot.ShapDependence` when the number of records is larger than 1000. In order to keep the overall trend, the sampling is stratified by bins. So the points cloud can be clearer to interpret.
- Modify some text in POEvaluation plot figure.

# itsdm 0.1.2

- Fix a few bugs in README example.
- Include an independent function to calculate continuous Boyce Index (`.cont_boyce`) in `utils.R` to reduce the pool of dependencies.
- Make the start message simpler.
- Fix the duplicated printout in function `print.POEvaluation`.
- Use `inherits` function to check "try-error" in dataset functions.

# itsdm 0.1.1

- Updated lines in function `evaluate_po`, `plot.POEvaluation`, and `print.POEvaluation` related to dependency `ecospat 3.2.1`. Because now function `ecospat.boyce` supports Kendall method, `itsdm` changed to use Kendall method to calculate CBI.
- Merge the pull request made by David Cortes who is the author of package `isotree` to use more flexible way for argument passing of `isolation.forest`.
- According to David Cortes' reminder, remove argument `sample_rate` from `isotree_po`. Only use `sample_size` for sub-sampling.

# itsdm 0.1.0

This is the first release. It includes all planned features.
