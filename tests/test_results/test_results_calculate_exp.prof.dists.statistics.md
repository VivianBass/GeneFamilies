
══ Testing test-calculate_exp.prof.dists.statistics.R ════════════════════════════════════════════════════════════════════════════════════════════════
[ FAIL 0 | WARN 0 | SKIP 0 | PASS 0 ][ FAIL 0 | WARN 0 | SKIP 0 | PASS 1 ][ FAIL 0 | WARN 0 | SKIP 0 | PASS 2 ][ FAIL 0 | WARN 0 | SKIP 0 | PASS 3 ][ FAIL 0 | WARN 0 | SKIP 0 | PASS 4 ][ FAIL 0 | WARN 0 | SKIP 0 | PASS 5 ][ FAIL 1 | WARN 0 | SKIP 0 | PASS 5 ][ FAIL 2 | WARN 0 | SKIP 0 | PASS 5 ][ FAIL 3 | WARN 0 | SKIP 0 | PASS 5 ][ FAIL 3 | WARN 0 | SKIP 0 | PASS 6 ][ FAIL 3 | WARN 0 | SKIP 0 | PASS 7 ]

── Error ('test-calculate_exp.prof.dists.statistics.R:55:5'): calculate_exp.prof.dists.statistics returns empty tibble for empty data ──
Error in `filter_all(., all_vars(!is.na(.) & !is.infinite(.)))`: `.predicate` must match at least one column.
Backtrace:
    ▆
 1. ├─global calculate_exp.prof.dists.statistics(data) at test-calculate_exp.prof.dists.statistics.R:55:5
 2. │ └─result %>% filter_all(all_vars(!is.na(.) & !is.infinite(.))) at ../R/compute_funks.R:99:3
 3. └─dplyr::filter_all(., all_vars(!is.na(.) & !is.infinite(.)))
 4.   └─dplyr:::apply_filter_syms(.vars_predicate, syms, .tbl)
 5.     └─rlang::abort(msg, call = error_call)

── Failure ('test-calculate_exp.prof.dists.statistics.R:77:5'): calculate_exp.prof.dists.statistics filters out infinite values ──
result$Mean (`actual`) not equal to `expected_means` (`expected`).

  `actual`:                                  
`expected`: 2.66666666666667 6.33333333333333

── Failure ('test-calculate_exp.prof.dists.statistics.R:78:5'): calculate_exp.prof.dists.statistics filters out infinite values ──
result$Median (`actual`) not equal to `expected_medians` (`expected`).

  `actual`:    
`expected`: 3 6
[ FAIL 3 | WARN 0 | SKIP 0 | PASS 7 ]
