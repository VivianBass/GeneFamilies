
══ Testing test-exp.prof.dists_tissue.R ══════════════════════════════════════════════════════════════════════════════════════════════════════════════
[ FAIL 0 | WARN 0 | SKIP 0 | PASS 0 ][ FAIL 0 | WARN 0 | SKIP 0 | PASS 1 ][ FAIL 0 | WARN 0 | SKIP 0 | PASS 2 ][ FAIL 0 | WARN 0 | SKIP 0 | PASS 3 ][ FAIL 0 | WARN 0 | SKIP 0 | PASS 4 ][ FAIL 0 | WARN 0 | SKIP 0 | PASS 5 ][ FAIL 0 | WARN 0 | SKIP 0 | PASS 6 ][ FAIL 0 | WARN 0 | SKIP 0 | PASS 7 ][ FAIL 1 | WARN 0 | SKIP 0 | PASS 7 ][ FAIL 1 | WARN 0 | SKIP 0 | PASS 8 ][ FAIL 1 | WARN 0 | SKIP 0 | PASS 9 ]

── Error ('test-exp.prof.dists_tissue.R:74:5'): exp.prof.dists_tissue handles missing tissue columns gracefully ──
<purrr_error_indexed/rlang_error/error/condition>
Error in `map(., ~{
    exp.profs %>% select(all_of(.x)) %>% dist(method = dist.method) %>% 
        as.vector()
})`: i In index: 2.
i With name: tissue2.
Caused by error in `select()`:
i In argument: `all_of(.x)`.
Caused by error in `all_of()`:
! Can't subset elements that don't exist.
x Element `tissue2` doesn't exist.
Backtrace:
     ▆
  1. ├─testthat::expect_error(...) at test-exp.prof.dists_tissue.R:74:5
  2. │ └─testthat:::expect_condition_matching(...)
  3. │   └─testthat:::quasi_capture(...)
  4. │     ├─testthat (local) .capture(...)
  5. │     │ └─base::withCallingHandlers(...)
  6. │     └─rlang::eval_bare(quo_get_expr(.quo), quo_get_env(.quo))
  7. ├─global exp.prof.dists_tissue(...)
  8. │ └─tissues %>% set_names() %>% ... at ../R/compute_funks.R:66:9
  9. ├─purrr::map(...)
 10. │ └─purrr:::map_("list", .x, .f, ..., .progress = .progress)
 11. │   ├─purrr:::with_indexed_errors(...)
 12. │   │ └─base::withCallingHandlers(...)
 13. │   ├─purrr:::call_with_cleanup(...)
 14. │   └─.f(.x[[i]], ...)
 15. │     └─... %>% as.vector() at ../R/compute_funks.R:69:17
 16. ├─base::as.vector(.)
 17. ├─stats::dist(., method = dist.method)
 18. │ └─base::as.matrix(x)
 19. ├─dplyr::select(., all_of(.x))
 20. ├─dplyr:::select.data.frame(., all_of(.x))
 21. │ └─tidyselect::eval_select(expr(c(...)), data = .data, error_call = error_call)
 22. │   └─tidyselect:::eval_select_impl(...)
 23. │     ├─tidyselect:::with_subscript_errors(...)
 24. │     │ └─base::withCallingHandlers(...)
 25. │     └─tidyselect:::vars_select_eval(...)
 26. │       └─tidyselect:::walk_data_tree(expr, data_mask, context_mask)
 27. │         └─tidyselect:::eval_c(expr, data_mask, context_mask)
 28. │           └─tidyselect:::reduce_sels(node, data_mask, context_mask, init = init)
 29. │             └─tidyselect:::walk_data_tree(new, data_mask, context_mask)
 30. │               └─tidyselect:::eval_context(expr, context_mask, call = error_call)
 31. │                 ├─tidyselect:::with_chained_errors(...)
 32. │                 │ └─base::withCallingHandlers(...)
 33. │                 └─rlang::eval_tidy(as_quosure(expr, env), context_mask)
 34. ├─tidyselect::all_of(.x)
 35. │ └─tidyselect:::as_indices_impl(x, vars = vars, strict = TRUE)
 36. │   └─tidyselect:::chr_as_locations(x, vars, call = call, arg = arg)
 37. │     └─vctrs::vec_as_location(...)
 38. ├─vctrs (local) `<fn>`()
 39. │ └─vctrs:::stop_subscript_oob(...)
 40. │   └─vctrs:::stop_subscript(...)
 41. │     └─rlang::abort(...)
 42. │       └─rlang:::signal_abort(cnd, .file)
 43. │         └─base::signalCondition(cnd)
 44. ├─tidyselect (local) `<fn>`(`<vctrs___>`)
 45. │ └─cli::cli_abort(c(i = msg), call = call, parent = cnd)
 46. │   └─rlang::abort(...)
 47. │     └─rlang:::signal_abort(cnd, .file)
 48. │       └─base::signalCondition(cnd)
 49. └─purrr (local) `<fn>`(`<rlng_rrr>`)
 50.   └─cli::cli_abort(...)
 51.     └─rlang::abort(...)
[ FAIL 1 | WARN 0 | SKIP 0 | PASS 9 ]
