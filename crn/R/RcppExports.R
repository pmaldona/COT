# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

rn_init <- function(input, lsp, verbose = FALSE) {
    invisible(.Call(`_rn_rn_init`, input, lsp, verbose))
}

rn_run_cs_gen <- function(n = 20000000L) {
    invisible(.Call(`_rn_rn_run_cs_gen`, n))
}

rn_run_c_cs_gen <- function() {
    invisible(.Call(`_rn_rn_run_c_cs_gen`))
}

rn_run_c_sso_gen <- function() {
    invisible(.Call(`_rn_rn_run_c_sso_gen`))
}

rn_run_sp2p <- function(sp) {
    .Call(`_rn_rn_run_sp2p`, sp)
}

rn_run_p2sp <- function(p) {
    .Call(`_rn_rn_run_p2sp`, p)
}

rn_run_closure <- function(p) {
    .Call(`_rn_rn_run_closure`, p)
}

rn_run_is_sep <- function(p, conn = 0L) {
    .Call(`_rn_rn_run_is_sep`, p, conn)
}

rn_run_all_syn <- function(no_bass = 0L) {
    invisible(.Call(`_rn_rn_run_all_syn`, no_bass))
}

rn_run_cont_sep <- function(conn = 0L) {
    invisible(.Call(`_rn_rn_run_cont_sep`, conn))
}

rn_get <- function() {
    .Call(`_rn_rn_get`)
}

