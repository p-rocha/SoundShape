#' Vocalization of the frog \emph{Physalaemus kroyeri}
#'
#' Recording of a series of three stereotyped calls emitted by a male frog \emph{Physalaemus kroyeri} (Amphibia, Anura, Leptodactylidae). Edited from original \code{".wav"} file for optimal sinal to noise ratio and reduced time duration.
#'
#' @docType data
#'
#' @usage data(kroyeri)
#'
#' @keywords datasets
#'
#' @details
#' Duration = 3.91 s. Sampling Frequency = 44100 Hz.
#'
#' Recorded at Ilhéus Municipality, Bahia State, Brazil, on 05 August 1972. Air temperature 24ºC.
#'
#' @format
#' An object of the class \code{"Wave"} (\code{\link{tuneR}} package).
#'
#' @source
#' Original recording housed at Fonoteca Neotropical Jacques Vielliard (FNJV-0032047). Recorded by Werner Bokermann.
#'
#' @examples
#' data(kroyeri)
#'
#' seewave::oscillo(kroyeri)
#' seewave::spectro(kroyeri)
#' threeDspectro(kroyeri,tlim=c(0, 1), flim=c(0, 4), samp.grid=FALSE, dBlevel=25)
"kroyeri"
