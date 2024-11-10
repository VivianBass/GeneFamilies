
══ Testing test-calculate_exp.prof.dists.tissue.statistics.R ═════════════════════════════════════════════════════════════════════════════════════════
[ FAIL 0 | WARN 0 | SKIP 0 | PASS 0 ][ FAIL 1 | WARN 0 | SKIP 0 | PASS 0 ][ FAIL 2 | WARN 0 | SKIP 0 | PASS 0 ][ FAIL 2 | WARN 0 | SKIP 0 | PASS 1 ][ FAIL 2 | WARN 0 | SKIP 0 | PASS 2 ][ FAIL 3 | WARN 0 | SKIP 0 | PASS 2 ][ FAIL 4 | WARN 0 | SKIP 0 | PASS 2 ][ FAIL 5 | WARN 0 | SKIP 0 | PASS 2 ][ FAIL 6 | WARN 0 | SKIP 0 | PASS 2 ][ FAIL 7 | WARN 0 | SKIP 0 | PASS 2 ][ FAIL 8 | WARN 0 | SKIP 0 | PASS 2 ][ FAIL 9 | WARN 0 | SKIP 0 | PASS 2 ]

── Failure ('test-calculate_exp.prof.dists.tissue.statistics.R:26:5'): calculate_exp.prof.dists.tissue.statistics computes mean and median per tissue for each family ──
result$Mean (`actual`) not equal to `expected_means` (`expected`).

`names(actual)` is a character vector ('tissue1', 'tissue2', 'tissue1', 'tissue2')
`names(expected)` is absent

── Failure ('test-calculate_exp.prof.dists.tissue.statistics.R:27:5'): calculate_exp.prof.dists.tissue.statistics computes mean and median per tissue for each family ──
result$Median (`actual`) not equal to `expected_medians` (`expected`).

`names(actual)` is a character vector ('tissue1', 'tissue2', 'tissue1', 'tissue2')
`names(expected)` is absent

── Failure ('test-calculate_exp.prof.dists.tissue.statistics.R:49:5'): calculate_exp.prof.dists.tissue.statistics handles tissues with NA values ──
result$Mean (`actual`) not equal to `expected_means` (`expected`).

`names(actual)` is a character vector ('tissue1', 'tissue2', 'tissue1', 'tissue2')
`names(expected)` is absent

── Failure ('test-calculate_exp.prof.dists.tissue.statistics.R:50:5'): calculate_exp.prof.dists.tissue.statistics handles tissues with NA values ──
result$Median (`actual`) not equal to `expected_medians` (`expected`).

`names(actual)` is a character vector ('tissue1', 'tissue2', 'tissue1', 'tissue2')
`names(expected)` is absent

── Error ('test-calculate_exp.prof.dists.tissue.statistics.R:58:5'): calculate_exp.prof.dists.tissue.statistics returns empty tibble for empty data ──
Error in `filter_all(., all_vars(!is.na(.) & !is.infinite(.)))`: `.predicate` must match at least one column.
Backtrace:
    ▆
 1. ├─global calculate_exp.prof.dists.tissue.statistics(data) at test-calculate_exp.prof.dists.tissue.statistics.R:58:5
 2. │ └─result %>% filter_all(all_vars(!is.na(.) & !is.infinite(.))) at ../R/compute_funks.R:127:3
 3. └─dplyr::filter_all(., all_vars(!is.na(.) & !is.infinite(.)))
 4.   └─dplyr:::apply_filter_syms(.vars_predicate, syms, .tbl)
 5.     └─rlang::abort(msg, call = error_call)

── Failure ('test-calculate_exp.prof.dists.tissue.statistics.R:82:5'): calculate_exp.prof.dists.tissue.statistics filters out infinite values ──
result$Mean (`actual`) not equal to `expected_means` (`expected`).

`names(actual)` is a character vector ()
`names(expected)` is absent

  `actual`:             
`expected`: 2 5 7.5 10.5

── Failure ('test-calculate_exp.prof.dists.tissue.statistics.R:83:5'): calculate_exp.prof.dists.tissue.statistics filters out infinite values ──
result$Median (`actual`) not equal to `expected_medians` (`expected`).

`names(actual)` is a character vector ()
`names(expected)` is absent

  `actual`:             
`expected`: 2 5 7.5 10.5

── Failure ('test-calculate_exp.prof.dists.tissue.statistics.R:101:5'): calculate_exp.prof.dists.tissue.statistics handles clusters with single-value tissues ──
result$Mean (`actual`) not equal to `expected_means` (`expected`).

`names(actual)` is a character vector ('tissue1', 'tissue2', 'tissue1', 'tissue2')
`names(expected)` is absent

── Failure ('test-calculate_exp.prof.dists.tissue.statistics.R:102:5'): calculate_exp.prof.dists.tissue.statistics handles clusters with single-value tissues ──
result$Median (`actual`) not equal to `expected_medians` (`expected`).

`names(actual)` is a character vector ('tissue1', 'tissue2', 'tissue1', 'tissue2')
`names(expected)` is absent
[ FAIL 9 | WARN 0 | SKIP 0 | PASS 2 ]
