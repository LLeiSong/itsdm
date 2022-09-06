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
