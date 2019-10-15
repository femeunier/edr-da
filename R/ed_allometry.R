#' Calculate leaf biomass given DBH
#'
#' Allometry equation is of the form `Y = b1 * DBH ^ b2`.
#'
#' Original function in ED src is in `utils/allometry.f90`. See also
#' the `area_indices` subroutine in the same file.
#'
#' LAI is calculated as `bleaf * nplant * SLA`.
#'
#' @param dbh Cohort diameter at breast height
#' @param b1 Allometry equation base
#' @param b2 Allometry equation exponent
#' @return Leaf biomass
#' @author Alexey Shiklomanov
#' @export
size2bl <- function(dbh, b1, b2) {
  # Carbon to biomass ratio of plant tissues. Defined in
  # `memory/pft_coms`; initialized in `init/ed_params`.
  C2B <- 2.0
  b1 / C2B * dbh ^ b2
}

#' Wood allometry
#' 
#' Calculate wood area index (WAI) given DBH
#'
#' Two possibilities here, depending on `iallom`.
#' The default value is: `WAI = nplant * b1 * DBH ^ b2
#'
#' An alternative formulation is that it is always 11% of LAI: `WAI =
#' 0.11 * LAI`. But we don't use that here because it's trivial.
#' @param dbh Cohort diameter at breast height
#' @param nplant Stem density (stems m-2)
#' @param b1 Allometry equation base
#' @param b2 Allometry equation exponent
#' @return Wood area index
#' @author Alexey Shiklomanov
#' @export
wai_allometry <- function(dbh, nplant, b1, b2) {
  nplant * b1 * dbh ^ b2 
}

#' crown area allometry
#'
#' @export
#' 
dbh2ca <- function(dbh, b1, b2) {
  b1 * dbh ^ b2
}
#' cai allometry
#'
#' @export
cai_allometry <- function(dbh,nplant,b1Bl,b2Bl,sla,b1Ca,b2Ca) {
  Ca <- dbh2ca(dbh, b1Ca, b2Ca)
  Bl <- size2bl(dbh, b1Bl, b2Bl)
  loclai <- sla *Bl
  dbh2ca <- pmin(loclai, Ca)
  pmin(1.0, nplant * dbh2ca)
}
