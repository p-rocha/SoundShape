#' Vocalization of the frog \emph{Physalaemus cuvieri}
#'
#' Recording of a series of three stereotyped calls emitted by a male frog \emph{Physalaemus cuvieri} (Amphibia, Anura, Leptodactylidae). Edited from original \code{".wav"} file for optimal sinal to noise ratio and reduced time duration.
#'
#' @docType data
#'
#' @usage data(cuvieri)
#'
#' @keywords datasets
#'
#' @details
#' Duration = 1.96 s. Sampling Frequency = 44100 Hz.
#'
#' Recorded at São José dos Campos Municipality, São Paulo State, Brazil, on 24 September 2013. Air temperature 22ºC.
#'
#' @format
#' An object of class \code{"Wave"}; see (\code{\link{tuneR}} package).
#'
#' @source
#' Original recording housed at Coleção Bioacústica da Universidade Federal de Minas Gerais (CBUFMG-00196). Recorded by Pedro Rocha.
#'
#' @examples
#' data(cuvieri)
#'
#' seewave::oscillo(cuvieri)
#' seewave::spectro(cuvieri)
#' threeDspectro(cuvieri, tlim=c(0, 0.5), flim=c(0, 4), samp.grid=FALSE)
"cuvieri"
