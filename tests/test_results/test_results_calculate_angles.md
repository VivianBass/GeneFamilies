
══ Testing test-calculate_angles.R ═══════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════
[ FAIL 0 | WARN 0 | SKIP 0 | PASS 0 ][ FAIL 0 | WARN 0 | SKIP 0 | PASS 1 ][ FAIL 0 | WARN 0 | SKIP 0 | PASS 2 ][ FAIL 0 | WARN 0 | SKIP 0 | PASS 3 ][ FAIL 0 | WARN 0 | SKIP 0 | PASS 4 ][ FAIL 1 | WARN 0 | SKIP 0 | PASS 4 ][ FAIL 2 | WARN 0 | SKIP 0 | PASS 4 ]

── Error ('test-calculate_angles.R:42:3'): calculate_angles handles empty gene groups ──
Error in `calculate_angles(genes, rna.seq.exp.profils, tissues)`: object 'group' not found
Backtrace:
    ▆
 1. └─global calculate_angles(genes, rna.seq.exp.profils, tissues) at test-calculate_angles.R:42:3
 2.   ├─base::warning(paste("No matching genes found for group:", group)) at ../R/angles_funks.R:48:5
 3.   └─base::paste("No matching genes found for group:", group) at ../R/angles_funks.R:48:5

── Error ('test-calculate_angles.R:56:3'): calculate_angles handles mismatched gene identifiers ──
Error in `calculate_angles(genes, rna.seq.exp.profils, tissues)`: object 'group' not found
Backtrace:
    ▆
 1. ├─testthat::expect_warning(...) at test-calculate_angles.R:56:3
 2. │ └─testthat:::expect_condition_matching(...)
 3. │   └─testthat:::quasi_capture(...)
 4. │     ├─testthat (local) .capture(...)
 5. │     │ └─base::withCallingHandlers(...)
 6. │     └─rlang::eval_bare(quo_get_expr(.quo), quo_get_env(.quo))
 7. └─global calculate_angles(genes, rna.seq.exp.profils, tissues)
 8.   ├─base::warning(paste("No matching genes found for group:", group)) at ../R/angles_funks.R:48:5
 9.   └─base::paste("No matching genes found for group:", group) at ../R/angles_funks.R:48:5
[ FAIL 2 | WARN 0 | SKIP 0 | PASS 4 ]
