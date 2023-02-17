#' Export sample \code{".wav"} files using selections from Raven Pro software.
#'
#' @description create one \code{".wav"} file for each selection created using \href{https://ravensoundsoftware.com/software/raven-pro/}{Raven Pro} software, which is commonplace in bioacoustical analysis. Each selection (i.e. line in table) should represent a each acoustic unit within the sample study.
#'
#' @param orig.wav filepath to the folder where original \code{".wav"} files are stored. Should be presented between quotation marks. By default: \code{orig.wav = NULL} (i.e. user must specify the filepath to \code{".wav"} files)
#'
#' @param raven.at filepath to the folder where selection tables from \href{https://ravensoundsoftware.com/software/raven-pro/}{Raven Pro} software are stored. Should be presented between quotation marks. File name should end with \code{"selections.txt"}. By default: \code{raven.at = orig.wav} (i.e. raven tables stored in the same folder as original  \code{".wav"} files)
#'
#' @param wav.samples name of the folder where new \code{".wav"} files will be stored. Should be presented between quotation marks. By default: \code{wav.samples = "wav samples"}
#'
#' @param max.dur optional argument to specify intended duration of sample \code{".wav"} files. If \code{max.dur=NULL}, function \code{raven.to.wav} will define the duration as 20% longer than the maximum duration of selections from Raven Pro software. By default: \code{max.dur=NULL}
#'
#'
#' @author
#' Pedro Rocha
#'
#' @seealso
#' \code{\link{align.wave}}
#'
#' Useful links:
#' \itemize{
#'   \item{\url{https://github.com/p-rocha/SoundShape}}
#'   \item{Report bugs at \url{https://github.com/p-rocha/SoundShape/issues}}
#'   \item{Raven Pro software \url{https://ravensoundsoftware.com/software/raven-pro/}}}
#'
#'
#' @examples
#' \donttest{
#'
#' # Write example
#'
#' }
#'
#' @references
#' MacLeod, N., Krieger, J. & Jones, K. E. (2013). Geometric morphometric approaches to acoustic signal analysis in mammalian biology. \emph{Hystrix, the Italian Journal of Mammalogy, 24}(1), 110-125.
#'
#' Rocha, P. & Romano, P. (2021) The shape of sound: A new \code{R} package that crosses the bridge between Bioacoustics and Geometric Morphometrics. \emph{Methods in Ecology and Evolution, 12}(6), 1115-1121.
#'
#'
#' @export
#'
raven.to.wave <- function(orig.wav=NULL, raven.at=orig.wav, wav.samples="wav samples", max.dur=NULL){

  if(is.null(orig.wav)) {stop("Use 'orig.wav' to specify folder path where original '.wav' files are stored")}

  # List ".wav" files
  wav.files <- dir(orig.wav, pattern=".wav")

  # List Raven tables
  raven.tables <- dir(raven.at, pattern="selections.txt")

  if(length(raven.tables)==0){
    stop("There are no tables containing selections from Raven Pro software at current folder. File path can be specified by either 'wav.at' or 'raven.at', and files must end with 'selections.txt'  files. Use 'help(raven.to.wav)' for more information")
  }

  if(length(raven.tables)!=length(wav.files)){
    warning("Number of selection tables from Raven Pro software differ from number of '.wav' files at folder specified by 'wav.at'. Some files may not have been analysed.")
  }


  # Import Raven tables
  raven.selections <- data.frame()
  for(raven in raven.tables){
    raven.temp <- read.table(raven, h=T, sep="\t", stringsAsFactors = T)

    raven.selections <- rbind(raven.selections, raven.temp)

    rm(raven.temp)
  } # import raven tables loop


  # Use delta time for maximum duration
  if(is.null(max.dur)){
    #longest acoustic unit plus 20% of longest duration
    max.dur <- max(raven.selections$End.Time..s.-raven.selections$Begin.Time..s.)*1.2  } # max.dur


  # Create folder to store sample wav files
  if(!dir.exists(file.path(orig.wav, wav.samples)))
    dir.create(file.path(orig.wav, wav.samples))


  # For each wav file, use Raven selections to create new files
  for(wav in wav.files){

    raven.temp <- read.table(grep(stringr::str_sub(wav,start=0, end = -5),
                                  raven.tables, value=T), h=T, sep="\t")

    for(i in 1:length(raven.temp$Selection[raven.temp$View=="Waveform 1"])){

      wav.temp <- tuneR::readWave(wav, units="seconds",
                                  from= raven.temp$Begin.Time..s.[
                                    raven.temp$Selection== i & raven.temp$View=="Waveform 1"],
                                  to= raven.temp$Begin.Time..s.[
                                    raven.temp$Selection== i & raven.temp$View=="Waveform 1"]+
                                    max.dur)
      tuneR::writeWave(wav.temp, extensible = T,
                       filename = file.path(orig.wav, wav.samples,
                                            paste(stringr::str_sub(wav,start=0, end = -5),
                                                  " - sample ", i, ".wav", sep="")))
      rm(wav.temp)

    } # end loop - for each selection

    rm(raven.temp)

  } # end loop - for each wav file

  rm(wav.files, raven.tables, raven.selections)

} # end function

