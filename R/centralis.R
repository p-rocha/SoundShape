#'Vocalization of the frog \emph{Physalaemus centralis}
#'
#' Recording of a series of three stereotyped calls emitted by a male frog \emph{Physalaemus centralis} (Amphibia, Anura, Leptodactylidae). Edited from original \code{".wav"} file for optimal sinal to noise ratio and reduced time duration.
#'
#' @docType data
#'
#' @usage data(centralis)
#'
#' @keywords datasets
#'
#' @details
#' Duration = 2.89 s. Sampling Frequency = 44100 Hz.
#'
#' Recorded at Formoso do Araguaia Municipality, Tocantins State, Brazil, on 9 December 1992. Air temperature 25ºC.
#'
#' @format
#' An object of the class \code{"Wave"} (\code{\link{tuneR}} package).
#'
#' @source
#' Original recording housed at Fonoteca Neotropical Jacques Vielliard (FNJV-0031188). Recorded by Adão José Cardoso.
#'
#' @examples
#' data(centralis)
#'
#' seewave::oscillo(centralis)
#' seewave::spectro(centralis)
#' threeDspectro(centralis,tlim=c(0, 0.8), flim=c(0, 4), samp.grid=FALSE, dBlevel=25)
"centralis"
