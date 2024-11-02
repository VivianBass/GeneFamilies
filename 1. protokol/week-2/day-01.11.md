

**Date**: 01.11.2024  
**Git Branch**: `Gene-Families-tests-Andre`

---

**Tasks**:

1. **Distance Calculations**:
- Tested `compute_exp.prof.dists.R` with new five-group input data.
- Began adapting the script for log-transformed distance calculations.

2. **Unit Testing**:
- Created unit tests for `load_data_funks.R` functions (`load_data_frame()`, `create_nested_list()`) using the `testthat` package.
- Documented test cases and scenarios for unit-testing the functions in `doc/testing_loading_functions.md`.

**Doubts and Issues**:

- **Data Alignment**: The script currently encounters mismatched gene IDs between expression profiles (FBgn IDs) and gene-groups data (FBpp IDs). Will continue with FBgn IDs for now and seek mapped data for consistent ID use.

**Next Steps**:
- Initiated and complete `roxygen2` documentation for all functions in `load_data_funks.R`.Complete `roxygen2` documentation.
