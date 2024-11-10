
══ Testing test-validate_data.R ══════════════════════════════════════════════════════════════════════════════════════════════════════════════════════
[ FAIL 0 | WARN 0 | SKIP 0 | PASS 0 ][ FAIL 0 | WARN 1 | SKIP 0 | PASS 0 ][ FAIL 0 | WARN 2 | SKIP 0 | PASS 0 ][ FAIL 0 | WARN 3 | SKIP 0 | PASS 0 ][ FAIL 0 | WARN 4 | SKIP 0 | PASS 0 ][ FAIL 0 | WARN 5 | SKIP 0 | PASS 0 ][ FAIL 0 | WARN 6 | SKIP 0 | PASS 0 ][ FAIL 0 | WARN 7 | SKIP 0 | PASS 0 ][ FAIL 0 | WARN 8 | SKIP 0 | PASS 0 ][ FAIL 0 | WARN 9 | SKIP 0 | PASS 0 ][ FAIL 0 | WARN 10 | SKIP 0 | PASS 0 ][ FAIL 0 | WARN 11 | SKIP 0 | PASS 0 ][ FAIL 0 | WARN 12 | SKIP 0 | PASS 0 ][ FAIL 0 | WARN 13 | SKIP 0 | PASS 0 ][ FAIL 0 | WARN 14 | SKIP 0 | PASS 0 ][ FAIL 0 | WARN 15 | SKIP 0 | PASS 0 ][ FAIL 0 | WARN 16 | SKIP 0 | PASS 0 ][ FAIL 0 | WARN 17 | SKIP 0 | PASS 0 ][ FAIL 0 | WARN 18 | SKIP 0 | PASS 0 ][ FAIL 0 | WARN 19 | SKIP 0 | PASS 0 ][ FAIL 0 | WARN 20 | SKIP 0 | PASS 0 ][ FAIL 0 | WARN 21 | SKIP 0 | PASS 0 ][ FAIL 0 | WARN 22 | SKIP 0 | PASS 0 ][ FAIL 0 | WARN 23 | SKIP 0 | PASS 0 ][ FAIL 0 | WARN 24 | SKIP 0 | PASS 0 ][ FAIL 0 | WARN 25 | SKIP 0 | PASS 0 ][ FAIL 1 | WARN 25 | SKIP 0 | PASS 0 ]

── Warning ('test-validate_data.R:3:1'): (code run outside of `test_that()`) ───
data set 'codingSequences' not found
Backtrace:
     ▆
  1. └─base::library(GeneFamilies) at test-validate_data.R:3:1
  2.   ├─base::tryCatch(...)
  3.   │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
  4.   │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
  5.   │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
  6.   └─base::loadNamespace(package, lib.loc)
  7.     └─base (local) runHook(".onLoad", env, package.lib, package)
  8.       ├─base::tryCatch(fun(libname, pkgname), error = identity)
  9.       │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
 10.       │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
 11.       │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
 12.       └─GeneFamilies (local) fun(libname, pkgname)
 13.         └─utils::data("codingSequences", package = "GeneFamilies")

── Warning ('test-validate_data.R:3:1'): (code run outside of `test_that()`) ───
data set 'families' not found
Backtrace:
     ▆
  1. └─base::library(GeneFamilies) at test-validate_data.R:3:1
  2.   ├─base::tryCatch(...)
  3.   │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
  4.   │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
  5.   │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
  6.   └─base::loadNamespace(package, lib.loc)
  7.     └─base (local) runHook(".onLoad", env, package.lib, package)
  8.       ├─base::tryCatch(fun(libname, pkgname), error = identity)
  9.       │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
 10.       │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
 11.       │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
 12.       └─GeneFamilies (local) fun(libname, pkgname)
 13.         └─utils::data("families", package = "GeneFamilies")

── Warning ('test-validate_data.R:3:1'): (code run outside of `test_that()`) ───
data set 'interPro' not found
Backtrace:
     ▆
  1. └─base::library(GeneFamilies) at test-validate_data.R:3:1
  2.   ├─base::tryCatch(...)
  3.   │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
  4.   │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
  5.   │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
  6.   └─base::loadNamespace(package, lib.loc)
  7.     └─base (local) runHook(".onLoad", env, package.lib, package)
  8.       ├─base::tryCatch(fun(libname, pkgname), error = identity)
  9.       │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
 10.       │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
 11.       │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
 12.       └─GeneFamilies (local) fun(libname, pkgname)
 13.         └─utils::data("interPro", package = "GeneFamilies")

── Warning ('test-validate_data.R:3:1'): (code run outside of `test_that()`) ───
data set 'familyHumanReadableDescriptions' not found
Backtrace:
     ▆
  1. └─base::library(GeneFamilies) at test-validate_data.R:3:1
  2.   ├─base::tryCatch(...)
  3.   │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
  4.   │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
  5.   │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
  6.   └─base::loadNamespace(package, lib.loc)
  7.     └─base (local) runHook(".onLoad", env, package.lib, package)
  8.       ├─base::tryCatch(fun(libname, pkgname), error = identity)
  9.       │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
 10.       │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
 11.       │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
 12.       └─GeneFamilies (local) fun(libname, pkgname)
 13.         └─utils::data("familyHumanReadableDescriptions", package = "GeneFamilies")

── Warning ('test-validate_data.R:3:1'): (code run outside of `test_that()`) ───
data set 'orthologsTandems' not found
Backtrace:
     ▆
  1. └─base::library(GeneFamilies) at test-validate_data.R:3:1
  2.   ├─base::tryCatch(...)
  3.   │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
  4.   │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
  5.   │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
  6.   └─base::loadNamespace(package, lib.loc)
  7.     └─base (local) runHook(".onLoad", env, package.lib, package)
  8.       ├─base::tryCatch(fun(libname, pkgname), error = identity)
  9.       │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
 10.       │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
 11.       │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
 12.       └─GeneFamilies (local) fun(libname, pkgname)
 13.         └─utils::data("orthologsTandems", package = "GeneFamilies")

── Warning ('test-validate_data.R:3:1'): (code run outside of `test_that()`) ───
data set 'cafe_result' not found
Backtrace:
     ▆
  1. └─base::library(GeneFamilies) at test-validate_data.R:3:1
  2.   ├─base::tryCatch(...)
  3.   │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
  4.   │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
  5.   │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
  6.   └─base::loadNamespace(package, lib.loc)
  7.     └─base (local) runHook(".onLoad", env, package.lib, package)
  8.       ├─base::tryCatch(fun(libname, pkgname), error = identity)
  9.       │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
 10.       │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
 11.       │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
 12.       └─GeneFamilies (local) fun(libname, pkgname)
 13.         └─utils::data("cafe_result", package = "GeneFamilies")

── Warning ('test-validate_data.R:3:1'): (code run outside of `test_that()`) ───
data set 'RNA_Seq_RPKM_and_profiles' not found
Backtrace:
     ▆
  1. └─base::library(GeneFamilies) at test-validate_data.R:3:1
  2.   ├─base::tryCatch(...)
  3.   │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
  4.   │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
  5.   │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
  6.   └─base::loadNamespace(package, lib.loc)
  7.     └─base (local) runHook(".onLoad", env, package.lib, package)
  8.       ├─base::tryCatch(fun(libname, pkgname), error = identity)
  9.       │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
 10.       │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
 11.       │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
 12.       └─GeneFamilies (local) fun(libname, pkgname)
 13.         └─utils::data("RNA_Seq_RPKM_and_profiles", package = "GeneFamilies")

── Warning ('test-validate_data.R:3:1'): (code run outside of `test_that()`) ───
data set 'pairwiseKaKs' not found
Backtrace:
     ▆
  1. └─base::library(GeneFamilies) at test-validate_data.R:3:1
  2.   ├─base::tryCatch(...)
  3.   │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
  4.   │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
  5.   │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
  6.   └─base::loadNamespace(package, lib.loc)
  7.     └─base (local) runHook(".onLoad", env, package.lib, package)
  8.       ├─base::tryCatch(fun(libname, pkgname), error = identity)
  9.       │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
 10.       │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
 11.       │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
 12.       └─GeneFamilies (local) fun(libname, pkgname)
 13.         └─utils::data("pairwiseKaKs", package = "GeneFamilies")

── Warning ('test-validate_data.R:3:1'): (code run outside of `test_that()`) ───
data set 'fubar_results' not found
Backtrace:
     ▆
  1. └─base::library(GeneFamilies) at test-validate_data.R:3:1
  2.   ├─base::tryCatch(...)
  3.   │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
  4.   │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
  5.   │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
  6.   └─base::loadNamespace(package, lib.loc)
  7.     └─base (local) runHook(".onLoad", env, package.lib, package)
  8.       ├─base::tryCatch(fun(libname, pkgname), error = identity)
  9.       │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
 10.       │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
 11.       │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
 12.       └─GeneFamilies (local) fun(libname, pkgname)
 13.         └─utils::data("fubar_results", package = "GeneFamilies")

── Warning ('test-validate_data.R:3:1'): (code run outside of `test_that()`) ───
data set 'BUSTED_Results' not found
Backtrace:
     ▆
  1. └─base::library(GeneFamilies) at test-validate_data.R:3:1
  2.   ├─base::tryCatch(...)
  3.   │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
  4.   │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
  5.   │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
  6.   └─base::loadNamespace(package, lib.loc)
  7.     └─base (local) runHook(".onLoad", env, package.lib, package)
  8.       ├─base::tryCatch(fun(libname, pkgname), error = identity)
  9.       │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
 10.       │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
 11.       │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
 12.       └─GeneFamilies (local) fun(libname, pkgname)
 13.         └─utils::data("BUSTED_Results", package = "GeneFamilies")

── Warning ('test-validate_data.R:3:1'): (code run outside of `test_that()`) ───
data set 'GeneGroups' not found
Backtrace:
     ▆
  1. └─base::library(GeneFamilies) at test-validate_data.R:3:1
  2.   ├─base::tryCatch(...)
  3.   │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
  4.   │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
  5.   │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
  6.   └─base::loadNamespace(package, lib.loc)
  7.     └─base (local) runHook(".onLoad", env, package.lib, package)
  8.       ├─base::tryCatch(fun(libname, pkgname), error = identity)
  9.       │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
 10.       │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
 11.       │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
 12.       └─GeneFamilies (local) fun(libname, pkgname)
 13.         └─utils::data("GeneGroups", package = "GeneFamilies")

── Warning ('test-validate_data.R:3:1'): (code run outside of `test_that()`) ───
data set 'IprBasedEntropies' not found
Backtrace:
     ▆
  1. └─base::library(GeneFamilies) at test-validate_data.R:3:1
  2.   ├─base::tryCatch(...)
  3.   │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
  4.   │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
  5.   │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
  6.   └─base::loadNamespace(package, lib.loc)
  7.     └─base (local) runHook(".onLoad", env, package.lib, package)
  8.       ├─base::tryCatch(fun(libname, pkgname), error = identity)
  9.       │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
 10.       │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
 11.       │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
 12.       └─GeneFamilies (local) fun(libname, pkgname)
 13.         └─utils::data("IprBasedEntropies", package = "GeneFamilies")

── Warning ('test-validate_data.R:3:1'): (code run outside of `test_that()`) ───
data set 'ExpressionProfileDistances' not found
Backtrace:
     ▆
  1. └─base::library(GeneFamilies) at test-validate_data.R:3:1
  2.   ├─base::tryCatch(...)
  3.   │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
  4.   │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
  5.   │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
  6.   └─base::loadNamespace(package, lib.loc)
  7.     └─base (local) runHook(".onLoad", env, package.lib, package)
  8.       ├─base::tryCatch(fun(libname, pkgname), error = identity)
  9.       │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
 10.       │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
 11.       │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
 12.       └─GeneFamilies (local) fun(libname, pkgname)
 13.         └─utils::data("ExpressionProfileDistances", package = "GeneFamilies")

── Warning ('test-validate_data.R:3:1'): (code run outside of `test_that()`) ───
data set 'ExpressionProfileDistanceDistributions' not found
Backtrace:
     ▆
  1. └─base::library(GeneFamilies) at test-validate_data.R:3:1
  2.   ├─base::tryCatch(...)
  3.   │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
  4.   │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
  5.   │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
  6.   └─base::loadNamespace(package, lib.loc)
  7.     └─base (local) runHook(".onLoad", env, package.lib, package)
  8.       ├─base::tryCatch(fun(libname, pkgname), error = identity)
  9.       │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
 10.       │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
 11.       │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
 12.       └─GeneFamilies (local) fun(libname, pkgname)
 13.         └─utils::data("ExpressionProfileDistanceDistributions", package = "GeneFamilies")

── Warning ('test-validate_data.R:3:1'): (code run outside of `test_that()`) ───
data set 'ExpressionProfileDistancesPerTissueDistributions' not found
Backtrace:
     ▆
  1. └─base::library(GeneFamilies) at test-validate_data.R:3:1
  2.   ├─base::tryCatch(...)
  3.   │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
  4.   │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
  5.   │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
  6.   └─base::loadNamespace(package, lib.loc)
  7.     └─base (local) runHook(".onLoad", env, package.lib, package)
  8.       ├─base::tryCatch(fun(libname, pkgname), error = identity)
  9.       │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
 10.       │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
 11.       │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
 12.       └─GeneFamilies (local) fun(libname, pkgname)
 13.         └─utils::data(...)

── Warning ('test-validate_data.R:3:1'): (code run outside of `test_that()`) ───
data set 'geneCopyNumbers' not found
Backtrace:
     ▆
  1. └─base::library(GeneFamilies) at test-validate_data.R:3:1
  2.   ├─base::tryCatch(...)
  3.   │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
  4.   │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
  5.   │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
  6.   └─base::loadNamespace(package, lib.loc)
  7.     └─base (local) runHook(".onLoad", env, package.lib, package)
  8.       ├─base::tryCatch(fun(libname, pkgname), error = identity)
  9.       │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
 10.       │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
 11.       │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
 12.       └─GeneFamilies (local) fun(libname, pkgname)
 13.         └─utils::data("geneCopyNumbers", package = "GeneFamilies")

── Warning ('test-validate_data.R:3:1'): (code run outside of `test_that()`) ───
data set 'correlationExpressionCopyNumber' not found
Backtrace:
     ▆
  1. └─base::library(GeneFamilies) at test-validate_data.R:3:1
  2.   ├─base::tryCatch(...)
  3.   │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
  4.   │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
  5.   │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
  6.   └─base::loadNamespace(package, lib.loc)
  7.     └─base (local) runHook(".onLoad", env, package.lib, package)
  8.       ├─base::tryCatch(fun(libname, pkgname), error = identity)
  9.       │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
 10.       │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
 11.       │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
 12.       └─GeneFamilies (local) fun(libname, pkgname)
 13.         └─utils::data("correlationExpressionCopyNumber", package = "GeneFamilies")

── Warning ('test-validate_data.R:3:1'): (code run outside of `test_that()`) ───
data set 'mapManRootBinAnnotations' not found
Backtrace:
     ▆
  1. └─base::library(GeneFamilies) at test-validate_data.R:3:1
  2.   ├─base::tryCatch(...)
  3.   │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
  4.   │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
  5.   │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
  6.   └─base::loadNamespace(package, lib.loc)
  7.     └─base (local) runHook(".onLoad", env, package.lib, package)
  8.       ├─base::tryCatch(fun(libname, pkgname), error = identity)
  9.       │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
 10.       │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
 11.       │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
 12.       └─GeneFamilies (local) fun(libname, pkgname)
 13.         └─utils::data("mapManRootBinAnnotations", package = "GeneFamilies")

── Warning ('test-validate_data.R:3:1'): (code run outside of `test_that()`) ───
data set 'rpkmExpressionProfiles' not found
Backtrace:
     ▆
  1. └─base::library(GeneFamilies) at test-validate_data.R:3:1
  2.   ├─base::tryCatch(...)
  3.   │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
  4.   │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
  5.   │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
  6.   └─base::loadNamespace(package, lib.loc)
  7.     └─base (local) runHook(".onLoad", env, package.lib, package)
  8.       ├─base::tryCatch(fun(libname, pkgname), error = identity)
  9.       │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
 10.       │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
 11.       │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
 12.       └─GeneFamilies (local) fun(libname, pkgname)
 13.         └─utils::data("rpkmExpressionProfiles", package = "GeneFamilies")

── Warning ('test-validate_data.R:3:1'): (code run outside of `test_that()`) ───
data set 'rpkmNormalizedExpressionProfiles' not found
Backtrace:
     ▆
  1. └─base::library(GeneFamilies) at test-validate_data.R:3:1
  2.   ├─base::tryCatch(...)
  3.   │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
  4.   │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
  5.   │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
  6.   └─base::loadNamespace(package, lib.loc)
  7.     └─base (local) runHook(".onLoad", env, package.lib, package)
  8.       ├─base::tryCatch(fun(libname, pkgname), error = identity)
  9.       │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
 10.       │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
 11.       │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
 12.       └─GeneFamilies (local) fun(libname, pkgname)
 13.         └─utils::data("rpkmNormalizedExpressionProfiles", package = "GeneFamilies")

── Warning ('test-validate_data.R:3:1'): (code run outside of `test_that()`) ───
data set 'differentiallyExpressedGenes' not found
Backtrace:
     ▆
  1. └─base::library(GeneFamilies) at test-validate_data.R:3:1
  2.   ├─base::tryCatch(...)
  3.   │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
  4.   │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
  5.   │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
  6.   └─base::loadNamespace(package, lib.loc)
  7.     └─base (local) runHook(".onLoad", env, package.lib, package)
  8.       ├─base::tryCatch(fun(libname, pkgname), error = identity)
  9.       │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
 10.       │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
 11.       │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
 12.       └─GeneFamilies (local) fun(libname, pkgname)
 13.         └─utils::data("differentiallyExpressedGenes", package = "GeneFamilies")

── Warning ('test-validate_data.R:3:1'): (code run outside of `test_that()`) ───
data set 'domainBasedSubfunctionalization' not found
Backtrace:
     ▆
  1. └─base::library(GeneFamilies) at test-validate_data.R:3:1
  2.   ├─base::tryCatch(...)
  3.   │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
  4.   │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
  5.   │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
  6.   └─base::loadNamespace(package, lib.loc)
  7.     └─base (local) runHook(".onLoad", env, package.lib, package)
  8.       ├─base::tryCatch(fun(libname, pkgname), error = identity)
  9.       │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
 10.       │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
 11.       │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
 12.       └─GeneFamilies (local) fun(libname, pkgname)
 13.         └─utils::data("domainBasedSubfunctionalization", package = "GeneFamilies")

── Warning ('test-validate_data.R:3:1'): (code run outside of `test_that()`) ───
data set 'tandemsDf' not found
Backtrace:
     ▆
  1. └─base::library(GeneFamilies) at test-validate_data.R:3:1
  2.   ├─base::tryCatch(...)
  3.   │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
  4.   │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
  5.   │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
  6.   └─base::loadNamespace(package, lib.loc)
  7.     └─base (local) runHook(".onLoad", env, package.lib, package)
  8.       ├─base::tryCatch(fun(libname, pkgname), error = identity)
  9.       │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
 10.       │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
 11.       │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
 12.       └─GeneFamilies (local) fun(libname, pkgname)
 13.         └─utils::data("tandemsDf", package = "GeneFamilies")

── Warning ('test-validate_data.R:3:1'): (code run outside of `test_that()`) ───
data set 'DuplicatedKsStats' not found
Backtrace:
     ▆
  1. └─base::library(GeneFamilies) at test-validate_data.R:3:1
  2.   ├─base::tryCatch(...)
  3.   │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
  4.   │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
  5.   │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
  6.   └─base::loadNamespace(package, lib.loc)
  7.     └─base (local) runHook(".onLoad", env, package.lib, package)
  8.       ├─base::tryCatch(fun(libname, pkgname), error = identity)
  9.       │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
 10.       │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
 11.       │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
 12.       └─GeneFamilies (local) fun(libname, pkgname)
 13.         └─utils::data("DuplicatedKsStats", package = "GeneFamilies")

── Warning ('test-validate_data.R:24:1'): (code run outside of `test_that()`) ──
cannot open file 'path/to/your/validate_data_function.R': No such file or directory
Backtrace:
    ▆
 1. ├─base::source("path/to/your/validate_data_function.R") at test-validate_data.R:24:1
 2. └─base::source("path/to/your/validate_data_function.R")
 3.   └─base::file(filename, "r", encoding = encoding)

── Error ('test-validate_data.R:24:1'): (code run outside of `test_that()`) ────
Error in `file(filename, "r", encoding = encoding)`: cannot open the connection
Backtrace:
    ▆
 1. ├─base::source("path/to/your/validate_data_function.R") at test-validate_data.R:24:1
 2. └─base::source("path/to/your/validate_data_function.R")
 3.   └─base::file(filename, "r", encoding = encoding)
[ FAIL 1 | WARN 25 | SKIP 0 | PASS 0 ]
