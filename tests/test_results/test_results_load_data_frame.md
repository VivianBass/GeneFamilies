
══ Testing test-load_data_frame.R ════════════════════════════════════════════════════════════════════════════════════════════════════════════════════
[ FAIL 0 | WARN 0 | SKIP 0 | PASS 0 ][ FAIL 0 | WARN 0 | SKIP 0 | PASS 1 ][ FAIL 0 | WARN 0 | SKIP 0 | PASS 2 ][ FAIL 0 | WARN 0 | SKIP 0 | PASS 3 ][ FAIL 0 | WARN 0 | SKIP 0 | PASS 4 ][ FAIL 0 | WARN 0 | SKIP 0 | PASS 5 ][ FAIL 1 | WARN 0 | SKIP 0 | PASS 5 ][ FAIL 2 | WARN 0 | SKIP 0 | PASS 5 ][ FAIL 2 | WARN 0 | SKIP 0 | PASS 6 ][ FAIL 2 | WARN 0 | SKIP 0 | PASS 7 ][ FAIL 3 | WARN 0 | SKIP 0 | PASS 7 ][ FAIL 4 | WARN 0 | SKIP 0 | PASS 7 ][ FAIL 4 | WARN 1 | SKIP 0 | PASS 7 ][ FAIL 5 | WARN 1 | SKIP 0 | PASS 7 ]

── Error ('test-load_data_frame.R:50:5'): load_data_frame handles missing columns ──
Error in `scan(file = file, what = what, sep = sep, quote = quote, dec = dec, 
    nmax = nrows, skip = 0, na.strings = na.strings, quiet = TRUE, 
    fill = fill, strip.white = strip.white, blank.lines.skip = blank.lines.skip, 
    multi.line = FALSE, comment.char = comment.char, allowEscapes = allowEscapes, 
    flush = flush, encoding = encoding, skipNul = skipNul)`: line 1 did not have 5 elements
Backtrace:
    ▆
 1. ├─testthat::expect_error(load_data_frame(temp_file), "more columns than column names") at test-load_data_frame.R:50:5
 2. │ └─testthat:::expect_condition_matching(...)
 3. │   └─testthat:::quasi_capture(...)
 4. │     ├─testthat (local) .capture(...)
 5. │     │ └─base::withCallingHandlers(...)
 6. │     └─rlang::eval_bare(quo_get_expr(.quo), quo_get_env(.quo))
 7. └─global load_data_frame(temp_file)
 8.   └─utils::read.table(...) at ../R/load_data_funks.R:15:5
 9.     └─base::scan(...)

── Failure ('test-load_data_frame.R:57:5'): load_data_frame handles incorrect column names ──
`load_data_frame(temp_file)` did not throw the expected error.

── Error ('test-load_data_frame.R:73:5'): load_data_frame handles non-tab delimited file ──
Error in `scan(file = file, what = what, sep = sep, quote = quote, dec = dec, 
    nmax = nrows, skip = 0, na.strings = na.strings, quiet = TRUE, 
    fill = fill, strip.white = strip.white, blank.lines.skip = blank.lines.skip, 
    multi.line = FALSE, comment.char = comment.char, allowEscapes = allowEscapes, 
    flush = flush, encoding = encoding, skipNul = skipNul)`: line 1 did not have 5 elements
Backtrace:
    ▆
 1. └─global load_data_frame(temp_file) at test-load_data_frame.R:73:5
 2.   └─utils::read.table(...) at ../R/load_data_funks.R:15:5
 3.     └─base::scan(...)

── Failure ('test-load_data_frame.R:82:5'): load_data_frame handles extra columns gracefully ──
ncol(df) (`actual`) not equal to 5 (`expected`).

  `actual`: 6
`expected`: 5

── Warning ('test-load_data_frame.R:87:5'): load_data_frame handles file not found error ──
cannot open file 'non_existent_file.txt': No such file or directory
Backtrace:
    ▆
 1. ├─testthat::expect_error(...) at test-load_data_frame.R:87:5
 2. │ └─testthat:::expect_condition_matching(...)
 3. │   └─testthat:::quasi_capture(...)
 4. │     ├─testthat (local) .capture(...)
 5. │     │ └─base::withCallingHandlers(...)
 6. │     └─rlang::eval_bare(quo_get_expr(.quo), quo_get_env(.quo))
 7. └─global load_data_frame("non_existent_file.txt")
 8.   └─utils::read.table(...) at ../R/load_data_funks.R:15:5
 9.     └─base::file(file, "rt")

── Error ('test-load_data_frame.R:87:5'): load_data_frame handles file not found error ──
Error in `file(file, "rt")`: cannot open the connection
Backtrace:
    ▆
 1. ├─testthat::expect_error(...) at test-load_data_frame.R:87:5
 2. │ └─testthat:::expect_condition_matching(...)
 3. │   └─testthat:::quasi_capture(...)
 4. │     ├─testthat (local) .capture(...)
 5. │     │ └─base::withCallingHandlers(...)
 6. │     └─rlang::eval_bare(quo_get_expr(.quo), quo_get_env(.quo))
 7. └─global load_data_frame("non_existent_file.txt")
 8.   └─utils::read.table(...) at ../R/load_data_funks.R:15:5
 9.     └─base::file(file, "rt")
[ FAIL 5 | WARN 1 | SKIP 0 | PASS 7 ]
