ctl_task <-
function(job) {
  t0 <- unname(proc.time()["elapsed"])
  out <- tryCatch({
    X <- X_obs_km[job$keep]; r <- jit_test(X, job$nsim, job$seed)
    list(n = npoints(X), p = r$p, usable = r$usable, po2 = r$pobs[1], po3 = r$pobs[2], po5 = r$pobs[3],
         ps2 = r$psim[1], ps3 = r$psim[2], ps5 = r$psim[3], err = NA_character_)
  }, error = function(e) list(err = conditionMessage(e)))
  c(list(grp = job$grp, i = job$i, sec = unname(proc.time()["elapsed"]) - t0), out)
}
crow <-
function(L) do.call(rbind, lapply(L, function(x) {
  g <- function(f) if (is.null(x[[f]])) NA_real_ else as.numeric(x[[f]])
  data.frame(grp = x$grp, i = x$i, sec = x$sec, n = g("n"), p = g("p"), usable = g("usable"), po2 = g("po2"), po3 = g("po3"), po5 = g("po5"),
             ps2 = g("ps2"), ps3 = g("ps3"), ps5 = g("ps5"), err = if (is.null(x$err)) NA_character_ else as.character(x$err), stringsAsFactors = FALSE)
}))
