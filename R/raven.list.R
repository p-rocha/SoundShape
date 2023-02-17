#' Sample of tables containing selection boxes from Raven Pro software
#'
#' This sample file features a list of three \code{data.frame} objects containing the selection boxes for the acoustic units of \code{\link{centralis}}, \code{\link{cuvieri}} and \code{\link{kroyeri}}. The selection was performed manually using \href{https://ravensoundsoftware.com/software/raven-pro/}{Raven Pro} software, which is commonplace in bioacoustic analysis.
#'
#' @docType data
#'
#' @usage data(raven.list)
#'
#' @keywords datasets
#'
#' @format
#' An object of the class \code{"list"} (\code{\link{base}} package).
#'
#' @details
#' This sample list was built to illustrate the usage of \code{\link{raven.to.wav}} function.
#' Each \code{data.frame} in the list represent selection boxes from \href{https://ravensoundsoftware.com/software/raven-pro/}{Raven Pro} software featuring either \code{\link{centralis}}, \code{\link{cuvieri}} or \code{\link{kroyeri}} samples. These, in turn, are acoustic recordings containing three stereotyped calls emitted by a male frog \emph{Physalaemus cuvieri}, \emph{P. centralis} or \emph{P. kroyeri}  (Amphibia, Anura, Leptodactylidae), respectively.
#'
#' @source
#' Sample data of \code{"Wave"} objects:
#' \itemize{
#'   \item{\code{\link{centralis}}: Advertisement call of \emph{Physalaemus centralis}; original recording housed at Fonoteca Neotropical Jacques Vielliard (FNJV-0031188). Recorded by Adão José Cardoso.}
#'   \item{\code{\link{cuvieri}}: Advertisement call of \emph{Physalaemus cuvieri}; original recording housed at Coleção Bioacústica da Universidade Federal de Minas Gerais (CBUFMG-00196). Recorded by Pedro Rocha.}
#'   \item{\code{\link{kroyeri}}: Advertisement call of \emph{Physalaemus kroyeri}; Original recording housed at Fonoteca Neotropical Jacques Vielliard (FNJV-0032047). Recorded by Werner Bokermann.}
#' }
#'
#' @references
#' Rocha, P. & Romano, P. (2021) The shape of sound: A new \code{R} package that crosses the bridge between Bioacoustics and Geometric Morphometrics. \emph{Methods in Ecology and Evolution, 12}(6), 1115-1121.
#'
#'
#'@examples
#'
#'
#' # Write example
#'
"raven.list"
