✔ | F W  S  OK | Context
⠏ |          0 | create_nested_list                                                                                         ⠇ | 8        1 | create_nested_list                                                                                         ✖ | 8        1 | create_nested_list
────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
Error ('test-create_nested_list.R:21:5'): create_nested_list creates nested list for Orthologs
Error in `deframe(.)`: could not find function "deframe"
Backtrace:
    ▆
 1. └─GeneFamilies:::create_nested_list(df, "Ortholog") at test-create_nested_list.R:21:5
 2.   └─... %>% deframe() at GeneFamilies/R/load_data_funks.R:29:5

Error ('test-create_nested_list.R:47:5'): create_nested_list creates nested list for Paralogs
Error in `deframe(.)`: could not find function "deframe"
Backtrace:
    ▆
 1. └─GeneFamilies:::create_nested_list(df, "Paralog") at test-create_nested_list.R:47:5
 2.   └─... %>% deframe() at GeneFamilies/R/load_data_funks.R:29:5

Error ('test-create_nested_list.R:73:5'): create_nested_list handles an empty data frame
Error in `deframe(.)`: could not find function "deframe"
Backtrace:
    ▆
 1. └─GeneFamilies:::create_nested_list(df, "Ortholog") at test-create_nested_list.R:73:5
 2.   └─... %>% deframe() at GeneFamilies/R/load_data_funks.R:29:5

Error ('test-create_nested_list.R:91:5'): create_nested_list creates separate nested lists for multiple families
Error in `deframe(.)`: could not find function "deframe"
Backtrace:
    ▆
 1. └─GeneFamilies:::create_nested_list(df, "Ortholog") at test-create_nested_list.R:91:5
 2.   └─... %>% deframe() at GeneFamilies/R/load_data_funks.R:29:5

Error ('test-create_nested_list.R:119:5'): create_nested_list throws error for missing required columns
Error in `deframe(.)`: could not find function "deframe"
Backtrace:
    ▆
 1. ├─testthat::expect_error(create_nested_list(df, "Ortholog"), "Required columns are missing.") at test-create_nested_list.R:119:5
 2. │ └─testthat:::expect_condition_matching(...)
 3. │   └─testthat:::quasi_capture(...)
 4. │     ├─testthat (local) .capture(...)
 5. │     │ └─base::withCallingHandlers(...)
 6. │     └─rlang::eval_bare(quo_get_expr(.quo), quo_get_env(.quo))
 7. └─GeneFamilies:::create_nested_list(df, "Ortholog")
 8.   └─... %>% deframe() at GeneFamilies/R/load_data_funks.R:29:5

Error ('test-create_nested_list.R:134:5'): create_nested_list handles different header types
Error in `deframe(.)`: could not find function "deframe"
Backtrace:
    ▆
 1. └─GeneFamilies:::create_nested_list(df, "Paralog") at test-create_nested_list.R:134:5
 2.   └─... %>% deframe() at GeneFamilies/R/load_data_funks.R:29:5

Error ('test-create_nested_list.R:160:5'): create_nested_list handles duplicate genes within the same family
Error in `deframe(.)`: could not find function "deframe"
Backtrace:
    ▆
 1. └─GeneFamilies:::create_nested_list(df, "Ortholog") at test-create_nested_list.R:160:5
 2.   └─... %>% deframe() at GeneFamilies/R/load_data_funks.R:29:5

Error ('test-create_nested_list.R:185:5'): create_nested_list handles non-character data types
Error in `deframe(.)`: could not find function "deframe"
Backtrace:
    ▆
 1. └─GeneFamilies:::create_nested_list(df, "Ortholog") at test-create_nested_list.R:185:5
 2.   └─... %>% deframe() at GeneFamilies/R/load_data_funks.R:29:5
────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
⠏ |          0 | load_data_frame                                                                                            ⠙ | 4 1      7 | load_data_frame                                                                                            ✖ | 5 1      7 | load_data_frame
────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
Error ('test-load_data_frame.R:46:5'): load_data_frame handles missing columns
Error in `scan(file = file, what = what, sep = sep, quote = quote, dec = dec, 
    nmax = nrows, skip = 0, na.strings = na.strings, quiet = TRUE, 
    fill = fill, strip.white = strip.white, blank.lines.skip = blank.lines.skip, 
    multi.line = FALSE, comment.char = comment.char, allowEscapes = allowEscapes, 
    flush = flush, encoding = encoding, skipNul = skipNul)`: line 1 did not have 5 elements
Backtrace:
    ▆
 1. ├─testthat::expect_error(load_data_frame(temp_file), "more columns than column names") at test-load_data_frame.R:46:5
 2. │ └─testthat:::expect_condition_matching(...)
 3. │   └─testthat:::quasi_capture(...)
 4. │     ├─testthat (local) .capture(...)
 5. │     │ └─base::withCallingHandlers(...)
 6. │     └─rlang::eval_bare(quo_get_expr(.quo), quo_get_env(.quo))
 7. └─GeneFamilies:::load_data_frame(temp_file)
 8.   └─utils::read.table(...) at GeneFamilies/R/load_data_funks.R:9:5
 9.     └─base::scan(...)

Failure ('test-load_data_frame.R:53:5'): load_data_frame handles incorrect column names
`load_data_frame(temp_file)` did not throw the expected error.

Error ('test-load_data_frame.R:69:5'): load_data_frame handles non-tab delimited file
Error in `scan(file = file, what = what, sep = sep, quote = quote, dec = dec, 
    nmax = nrows, skip = 0, na.strings = na.strings, quiet = TRUE, 
    fill = fill, strip.white = strip.white, blank.lines.skip = blank.lines.skip, 
    multi.line = FALSE, comment.char = comment.char, allowEscapes = allowEscapes, 
    flush = flush, encoding = encoding, skipNul = skipNul)`: line 1 did not have 5 elements
Backtrace:
    ▆
 1. └─GeneFamilies:::load_data_frame(temp_file) at test-load_data_frame.R:69:5
 2.   └─utils::read.table(...) at GeneFamilies/R/load_data_funks.R:9:5
 3.     └─base::scan(...)

Failure ('test-load_data_frame.R:78:5'): load_data_frame handles extra columns gracefully
ncol(df) (`actual`) not equal to 5 (`expected`).

  `actual`: 6
`expected`: 5

Warning ('test-load_data_frame.R:83:5'): load_data_frame handles file not found error
cannot open file 'non_existent_file.txt': No such file or directory
Backtrace:
    ▆
 1. ├─testthat::expect_error(...) at test-load_data_frame.R:83:5
 2. │ └─testthat:::expect_condition_matching(...)
 3. │   └─testthat:::quasi_capture(...)
 4. │     ├─testthat (local) .capture(...)
 5. │     │ └─base::withCallingHandlers(...)
 6. │     └─rlang::eval_bare(quo_get_expr(.quo), quo_get_env(.quo))
 7. └─GeneFamilies:::load_data_frame("non_existent_file.txt")
 8.   └─utils::read.table(...) at GeneFamilies/R/load_data_funks.R:9:5
 9.     └─base::file(file, "rt")

Error ('test-load_data_frame.R:83:5'): load_data_frame handles file not found error
Error in `file(file, "rt")`: cannot open the connection
Backtrace:
    ▆
 1. ├─testthat::expect_error(...) at test-load_data_frame.R:83:5
 2. │ └─testthat:::expect_condition_matching(...)
 3. │   └─testthat:::quasi_capture(...)
 4. │     ├─testthat (local) .capture(...)
 5. │     │ └─base::withCallingHandlers(...)
 6. │     └─rlang::eval_bare(quo_get_expr(.quo), quo_get_env(.quo))
 7. └─GeneFamilies:::load_data_frame("non_existent_file.txt")
 8.   └─utils::read.table(...) at GeneFamilies/R/load_data_funks.R:9:5
 9.     └─base::file(file, "rt")
────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
Maximum number of failures exceeded; quitting at end of file.
ℹ Increase this number with (e.g.) `testthat::set_max_fails(Inf)` 

══ Results ═════════════════════════════════════════════════════════════════════════════════════════════════════════════════
Duration: 1.4 s

── Failed tests ────────────────────────────────────────────────────────────────────────────────────────────────────────────
Error ('test-create_nested_list.R:21:5'): create_nested_list creates nested list for Orthologs
Error in `deframe(.)`: could not find function "deframe"
Backtrace:
    ▆
 1. └─GeneFamilies:::create_nested_list(df, "Ortholog") at test-create_nested_list.R:21:5
 2.   └─... %>% deframe() at GeneFamilies/R/load_data_funks.R:29:5

Error ('test-create_nested_list.R:47:5'): create_nested_list creates nested list for Paralogs
Error in `deframe(.)`: could not find function "deframe"
Backtrace:
    ▆
 1. └─GeneFamilies:::create_nested_list(df, "Paralog") at test-create_nested_list.R:47:5
 2.   └─... %>% deframe() at GeneFamilies/R/load_data_funks.R:29:5

Error ('test-create_nested_list.R:73:5'): create_nested_list handles an empty data frame
Error in `deframe(.)`: could not find function "deframe"
Backtrace:
    ▆
 1. └─GeneFamilies:::create_nested_list(df, "Ortholog") at test-create_nested_list.R:73:5
 2.   └─... %>% deframe() at GeneFamilies/R/load_data_funks.R:29:5

Error ('test-create_nested_list.R:91:5'): create_nested_list creates separate nested lists for multiple families
Error in `deframe(.)`: could not find function "deframe"
Backtrace:
    ▆
 1. └─GeneFamilies:::create_nested_list(df, "Ortholog") at test-create_nested_list.R:91:5
 2.   └─... %>% deframe() at GeneFamilies/R/load_data_funks.R:29:5

Error ('test-create_nested_list.R:119:5'): create_nested_list throws error for missing required columns
Error in `deframe(.)`: could not find function "deframe"
Backtrace:
    ▆
 1. ├─testthat::expect_error(create_nested_list(df, "Ortholog"), "Required columns are missing.") at test-create_nested_list.R:119:5
 2. │ └─testthat:::expect_condition_matching(...)
 3. │   └─testthat:::quasi_capture(...)
 4. │     ├─testthat (local) .capture(...)
 5. │     │ └─base::withCallingHandlers(...)
 6. │     └─rlang::eval_bare(quo_get_expr(.quo), quo_get_env(.quo))
 7. └─GeneFamilies:::create_nested_list(df, "Ortholog")
 8.   └─... %>% deframe() at GeneFamilies/R/load_data_funks.R:29:5

Error ('test-create_nested_list.R:134:5'): create_nested_list handles different header types
Error in `deframe(.)`: could not find function "deframe"
Backtrace:
    ▆
 1. └─GeneFamilies:::create_nested_list(df, "Paralog") at test-create_nested_list.R:134:5
 2.   └─... %>% deframe() at GeneFamilies/R/load_data_funks.R:29:5

Error ('test-create_nested_list.R:160:5'): create_nested_list handles duplicate genes within the same family
Error in `deframe(.)`: could not find function "deframe"
Backtrace:
    ▆
 1. └─GeneFamilies:::create_nested_list(df, "Ortholog") at test-create_nested_list.R:160:5
 2.   └─... %>% deframe() at GeneFamilies/R/load_data_funks.R:29:5

Error ('test-create_nested_list.R:185:5'): create_nested_list handles non-character data types
Error in `deframe(.)`: could not find function "deframe"
Backtrace:
    ▆
 1. └─GeneFamilies:::create_nested_list(df, "Ortholog") at test-create_nested_list.R:185:5
 2.   └─... %>% deframe() at GeneFamilies/R/load_data_funks.R:29:5

Error ('test-load_data_frame.R:46:5'): load_data_frame handles missing columns
Error in `scan(file = file, what = what, sep = sep, quote = quote, dec = dec, 
    nmax = nrows, skip = 0, na.strings = na.strings, quiet = TRUE, 
    fill = fill, strip.white = strip.white, blank.lines.skip = blank.lines.skip, 
    multi.line = FALSE, comment.char = comment.char, allowEscapes = allowEscapes, 
    flush = flush, encoding = encoding, skipNul = skipNul)`: line 1 did not have 5 elements
Backtrace:
    ▆
 1. ├─testthat::expect_error(load_data_frame(temp_file), "more columns than column names") at test-load_data_frame.R:46:5
 2. │ └─testthat:::expect_condition_matching(...)
 3. │   └─testthat:::quasi_capture(...)
 4. │     ├─testthat (local) .capture(...)
 5. │     │ └─base::withCallingHandlers(...)
 6. │     └─rlang::eval_bare(quo_get_expr(.quo), quo_get_env(.quo))
 7. └─GeneFamilies:::load_data_frame(temp_file)
 8.   └─utils::read.table(...) at GeneFamilies/R/load_data_funks.R:9:5
 9.     └─base::scan(...)

Failure ('test-load_data_frame.R:53:5'): load_data_frame handles incorrect column names
`load_data_frame(temp_file)` did not throw the expected error.

Error ('test-load_data_frame.R:69:5'): load_data_frame handles non-tab delimited file
Error in `scan(file = file, what = what, sep = sep, quote = quote, dec = dec, 
    nmax = nrows, skip = 0, na.strings = na.strings, quiet = TRUE, 
    fill = fill, strip.white = strip.white, blank.lines.skip = blank.lines.skip, 
    multi.line = FALSE, comment.char = comment.char, allowEscapes = allowEscapes, 
    flush = flush, encoding = encoding, skipNul = skipNul)`: line 1 did not have 5 elements
Backtrace:
    ▆
 1. └─GeneFamilies:::load_data_frame(temp_file) at test-load_data_frame.R:69:5
 2.   └─utils::read.table(...) at GeneFamilies/R/load_data_funks.R:9:5
 3.     └─base::scan(...)

Failure ('test-load_data_frame.R:78:5'): load_data_frame handles extra columns gracefully
ncol(df) (`actual`) not equal to 5 (`expected`).

  `actual`: 6
`expected`: 5

Error ('test-load_data_frame.R:83:5'): load_data_frame handles file not found error
Error in `file(file, "rt")`: cannot open the connection
Backtrace:
    ▆
 1. ├─testthat::expect_error(...) at test-load_data_frame.R:83:5
 2. │ └─testthat:::expect_condition_matching(...)
 3. │   └─testthat:::quasi_capture(...)
 4. │     ├─testthat (local) .capture(...)
 5. │     │ └─base::withCallingHandlers(...)
 6. │     └─rlang::eval_bare(quo_get_expr(.quo), quo_get_env(.quo))
 7. └─GeneFamilies:::load_data_frame("non_existent_file.txt")
 8.   └─utils::read.table(...) at GeneFamilies/R/load_data_funks.R:9:5
 9.     └─base::file(file, "rt")

[ FAIL 13 | WARN 1 | SKIP 0 | PASS 8 ]
══ Terminated early ════════════════════════════════════════════════════════════════════════════════════════════════════════
