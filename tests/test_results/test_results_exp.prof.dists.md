
══ Testing test-exp.prof.dists.R ═════════════════════════════════════════════════════════════════════════════════════════════════════════════════════
[ FAIL 0 | WARN 0 | SKIP 0 | PASS 0 ][ FAIL 0 | WARN 0 | SKIP 0 | PASS 1 ][ FAIL 0 | WARN 0 | SKIP 0 | PASS 2 ][ FAIL 0 | WARN 0 | SKIP 0 | PASS 3 ][ FAIL 1 | WARN 0 | SKIP 0 | PASS 3 ][ FAIL 1 | WARN 0 | SKIP 0 | PASS 4 ][ FAIL 1 | WARN 0 | SKIP 0 | PASS 5 ]

── Error ('test-exp.prof.dists.R:52:5'): exp.prof.dists handles missing columns gracefully ──
Error in `select(., all_of(tissues))`: i In argument: `all_of(tissues)`.
Caused by error in `all_of()`:
! Can't subset elements that don't exist.
x Element `tissue2` doesn't exist.
[ FAIL 1 | WARN 0 | SKIP 0 | PASS 5 ]
