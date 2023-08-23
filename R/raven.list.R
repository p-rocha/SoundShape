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
#' This sample list was built to illustrate the usage of \code{\link{raven.to.wave}} function.
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
#' @seealso
#' \code{\link{raven.to.wave}}
#'
#' @examples
#' \donttest{
#'
#'#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#'#  Create folders on your console  #
#'#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#'
#'# Create temporary folder to store original ".wav" files containing multiple units
#'orig.wav <- file.path(base::tempdir(), "original wave")
#'if(!dir.exists(orig.wav)) dir.create(orig.wav)
#'
#'# Create temporary folder to store sample ".wav" files from original recordings
#'wav.at <- file.path(base::tempdir(), "wav samples")
#'if(!dir.exists(wav.at)) dir.create(wav.at)
#'
#'# Create temporary folder to store results
#'store.at <- file.path(base::tempdir(), "output")
#'if(!dir.exists(store.at)) dir.create(store.at)
#'
#'
#'#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#'#   Select acoustic units based on Raven Pro selections   #
#'#                                                         #
#'#              Using raven.to.wave function                #
#'#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#'
#'# Export original sample ".wav" files from SoundShape examples
#'tuneR::writeWave(centralis, extensible = TRUE,
#'                 filename = file.path(orig.wav, "centralis.wav"))
#'tuneR::writeWave(cuvieri, extensible = TRUE,
#'                 filename = file.path(orig.wav, "cuvieri.wav"))
#'tuneR::writeWave(kroyeri, extensible = TRUE,
#'                 filename = file.path(orig.wav, "kroyeri.wav"))
#'
#'# Store Raven Pro selection tables at same folder from original ".wav" files
#'for(i in 1:length(raven.list)){
#'  write.table(raven.list[i], file=file.path(orig.wav, names(raven.list)[i]),
#'                quote=FALSE, sep="\t", row.names = FALSE,
#'                col.names = colnames(raven.list[[i]]))  } # end loop
#'
#'
#'# Verify if folder has both original ".wav" files and Raven's selections
#'dir(orig.wav)
#'
#'###
#'## Start here when using your own recordings
#'
#'# Export a ".wav" sample for each selection made in Raven Pro
#'raven.to.wave(orig.wav.folder = orig.wav, wav.samples = wav.at)
#'
#'# Verify samples
#'dir(wav.at)
#'
#'
#'#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#'#  Align acoustic units using align.wave  #
#'#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#'
#'# Place sounds at the beginning of a sound window
#'align.wave(wav.at=wav.at, wav.to="Aligned",
#'           time.length = 0.8, time.perc = 0.005, dBlevel = 25)
#'
#'# Verify alignment using analysis.type = "twoDshape"
#'eigensound(analysis.type = "twoDshape", wav.at = file.path(wav.at, "Aligned"),
#'           dBlevel = 25, store.at=store.at, plot.exp=TRUE,
#'           flim=c(0, 4), tlim=c(0, 0.8), add.contour = TRUE)
#'# Go to folder specified by store.at and check jpeg files created
#'# If alignment/window dimensions are not ideal, repeat the process with new settings
#'
#' }
#'
"raven.list"
