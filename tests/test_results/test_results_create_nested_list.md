
══ Testing test-create_nested_list.R ═════════════════════════════════════════════════════════════════════════════════════════════════════════════════
[ FAIL 0 | WARN 0 | SKIP 0 | PASS 0 ][ FAIL 0 | WARN 0 | SKIP 0 | PASS 1 ][ FAIL 0 | WARN 0 | SKIP 0 | PASS 2 ][ FAIL 1 | WARN 0 | SKIP 0 | PASS 2 ][ FAIL 1 | WARN 0 | SKIP 0 | PASS 3 ][ FAIL 2 | WARN 0 | SKIP 0 | PASS 3 ][ FAIL 2 | WARN 0 | SKIP 0 | PASS 4 ][ FAIL 2 | WARN 0 | SKIP 0 | PASS 5 ][ FAIL 3 | WARN 0 | SKIP 0 | PASS 5 ][ FAIL 3 | WARN 0 | SKIP 0 | PASS 6 ]

── Error ('test-create_nested_list.R:76:5'): create_nested_list handles an empty data frame ──
Error in `summarise(., gene_info = list(setNames(nested_info, paste0("(", 
    Gene_species, ", ", Gene, ")"))), .groups = "drop")`: i In argument: `gene_info = list(...)`.
Caused by error in `names(object) <- nm`:
! 'names' attribute [1] must be the same length as the vector [0]
Backtrace:
     ▆
  1. ├─global create_nested_list(df, "Ortholog") at test-create_nested_list.R:76:5
  2. │ └─... %>% deframe() at ../R/load_data_funks.R:46:5
  3. ├─tibble::deframe(.)
  4. ├─dplyr::summarise(...)
  5. ├─dplyr:::summarise.grouped_df(...)
  6. │ └─dplyr:::summarise_cols(.data, dplyr_quosures(...), by, "summarise")
  7. │   ├─base::withCallingHandlers(...)
  8. │   └─dplyr:::map(quosures, summarise_eval_one, mask = mask)
  9. │     └─base::lapply(.x, .f, ...)
 10. │       └─dplyr (local) FUN(X[[i]], ...)
 11. │         └─mask$eval_all_summarise(quo)
 12. │           └─dplyr (local) eval()
 13. ├─stats::setNames(...)
 14. └─base::.handleSimpleError(...)
 15.   └─dplyr (local) h(simpleError(msg, call))
 16.     └─dplyr (local) handler(cnd)
 17.       └─rlang::abort(message, class = error_class, parent = parent, call = error_call)

── Error ('test-create_nested_list.R:122:5'): create_nested_list throws error for missing required columns ──
Error in `group_by(., Family, Gene_species, Gene, !!sym(species_col))`: Must group by variables found in `.data`.
x Column `Gene` is not found.
Backtrace:
     ▆
  1. ├─testthat::expect_error(create_nested_list(df, "Ortholog"), "Required columns are missing.") at test-create_nested_list.R:122:5
  2. │ └─testthat:::expect_condition_matching(...)
  3. │   └─testthat:::quasi_capture(...)
  4. │     ├─testthat (local) .capture(...)
  5. │     │ └─base::withCallingHandlers(...)
  6. │     └─rlang::eval_bare(quo_get_expr(.quo), quo_get_env(.quo))
  7. ├─global create_nested_list(df, "Ortholog")
  8. │ └─... %>% deframe() at ../R/load_data_funks.R:46:5
  9. ├─tibble::deframe(.)
 10. ├─dplyr::summarise(...)
 11. ├─dplyr::group_by(., Family)
 12. ├─dplyr::summarise(...)
 13. ├─dplyr::group_by(., Family, Gene_species, Gene)
 14. ├─dplyr::summarise(...)
 15. ├─dplyr::group_by(., Family, Gene_species, Gene, !!sym(species_col))
 16. └─dplyr:::group_by.data.frame(., Family, Gene_species, Gene, !!sym(species_col))
 17.   └─dplyr::group_by_prepare(.data, ..., .add = .add, error_call = current_env())
 18.     └─rlang::abort(bullets, call = error_call)

── Failure ('test-create_nested_list.R:199:5'): create_nested_list handles non-character data types ──
`nested_list` (`actual`) not equal to `expected` (`expected`).

`actual$F1$(Sp1, G1)$SpO1` is an integer vector (1)
`expected$F1$(Sp1, G1)$SpO1` is a character vector ('1')

`actual$F1$(Sp1, G2)$SpO2` is an integer vector (2)
`expected$F1$(Sp1, G2)$SpO2` is a character vector ('2')
[ FAIL 3 | WARN 0 | SKIP 0 | PASS 6 ]
