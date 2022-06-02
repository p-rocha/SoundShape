#' Automatic placement of calls at beggining of sound window
#'
#' @description Recreate each \code{".wav"} file on a given folder while placing calls at the beggining of sound window. New \code{".wav"} files will be stored on a new folder, which is automatically created.
#'
#' @param wav.at filepath to the folder where \code{".wav"} files are stored. Should be presented between quotation marks. By default: \code{wav.at = NULL} (i.e. user must specify the filepath to \code{".wav"} files)
#' @param wav.to name of the folder where new \code{".wav"} files will be stored. Should be presented between quotation marks. By default: \code{wav.to = "Aligned"}
#' @param time.length intended length for the time (X-axis) in seconds. Should be a value that encompasses all sounds in the study. By default: \code{time.length = 1}
#' @param time.perc slight time gap (in percentage) relative to the intended length that encompass all sounds in the study (i.e. \code{time.length}). Intervals are added before and after the minimum and maximum time coordinates (X-values) from the selected curve of relative amplitude (\code{dBlevel}). By default: \code{time.perc = 0.005} (i.e. 0.5%)
#' @param flim modifications of the frequency limits (Y-axis) to focus on acoustic units. Useful for recordings with low signal-to-noise ratio. Vector with two values in kHz. By default: \code{flim = NULL}
#' @param dBlevel absolute amplitude value to be used as relative amplitude contour, which will serve as reference for call placement. By default: \code{dBlevel = 25}
#' @param f sampling frequency of \code{".wav"} files (in Hz). By default: \code{f = 44100}
#' @param wl length of the window for the analysis. By default: \code{wl = 512}
#' @param ovlp overlap between two successive windows (in %) for increased spectrogram resolution. By default: \code{ovlp = 70}
#'
#' @author
#' Pedro Rocha
#'
#' @seealso
#' \code{\link{eigensound}}
#'
#' Useful links:
#' \itemize{
#'   \item{\url{https://github.com/p-rocha/SoundShape}}
#'   \item{Report bugs at \url{https://github.com/p-rocha/SoundShape/issues}}}
#'
#' @examples
#' \donttest{
#' library(seewave)
#' library(tuneR)
#'
#' # Create temporary folder to store ".wav" files
#' wav.at <- file.path(base::tempdir(), "align.wave")
#' if(!dir.exists(wav.at)) dir.create(wav.at)
#'
#' # Create temporary folder to store results
#' store.at <- file.path(base::tempdir(), "align.wave-output")
#' if(!dir.exists(store.at)) dir.create(store.at)
#'
#' # Select acoustic units to be analyzed
#' data(cuvieri)
#' spectro(cuvieri, flim = c(0,3)) # Visualize sound data that will be used
#'
#' # Cut acoustic units from original Wave
#' cut.cuvieri1 <- cutw(cuvieri, f=44100, from=0, to=0.5, output = "Wave")
#' cut.cuvieri2 <- cutw(cuvieri, f=44100, from=0.7, to=1.2, output = "Wave")
#' cut.cuvieri3 <- cutw(cuvieri, f=44100, from=1.4, to=1.9, output = "Wave")
#'
#' # Export ".wav" files containing selected acoustic units and store on previosly created folder
#' writeWave(cut.cuvieri1, filename = file.path(wav.at, "cut.cuvieri1.wav"), extensible = FALSE)
#' writeWave(cut.cuvieri2, filename = file.path(wav.at, "cut.cuvieri2.wav"), extensible = FALSE)
#' writeWave(cut.cuvieri3, filename = file.path(wav.at, "cut.cuvieri3.wav"), extensible = FALSE)
#'
#' # Align acoustic units selected at 1% of time lenght
#' align.wave(wav.at = wav.at, wav.to = "Aligned",
#'            time.length = 0.5, time.perc = 0.01, dBlevel = 25)
#'
#' # Verify alignment using eigensound function featuring analysis.type = "twoDshape"
#' eigensound(analysis.type = "twoDshape", wav.at = file.path(wav.at, "Aligned"), store.at = store.at,
#'            flim=c(0, 3), tlim=c(0,0.5), dBlevel = 25, plot.exp = TRUE, plot.as = "jpeg")
#' # To see jpeg files created, check folder specified by store.at
#'
#' }
#'
#' @references
#' MacLeod, N., Krieger, J. & Jones, K. E. (2013). Geometric morphometric approaches to acoustic signal analysis in mammalian biology. \emph{Hystrix, the Italian Journal of Mammalogy, 24}(1), 110-125.
#'
#' Rocha, P. & Romano, P. (\emph{in prep}) The shape of sound: A new \code{R} package that crosses the bridge between Bioacoustics and Geometric Morphometrics. \emph{Methods in Ecology and Evolution}
#'
#'
#' @export
#'
align.wave <- function(wav.at=NULL, wav.to="Aligned", time.length=1, time.perc=0.005, flim=NULL, dBlevel=25, f=44100, wl=512, ovlp=70)  {

  if(is.null(wav.at)) {stop("Use 'wav.at' to specify folder path where '.wav' files are stored")}

  # Create folder to store aligned calls
  if(!dir.exists(file.path(wav.at, wav.to))) dir.create(file.path(wav.at, wav.to))

  # Replace sounds for each ".wav" file in a folder
  for(file in list.files(wav.at, pattern = ".wav")){

    orig.wav0 <- tuneR::readWave(paste(wav.at,"/", file, sep=""))

    # Add silence to fill sound window and prevent error
    orig.wav <- seewave::addsilw(orig.wav0, f=f, at="end", d=(time.length*10), output = "Wave")

    # create spectro object
    orig.spec <- seewave::spectro(orig.wav, f=f, wl=wl, ovlp=ovlp, osc=F, grid=F, plot=F, flim=flim)

    # Acquire contours
    cont.spec <- grDevices::contourLines(x=orig.spec$time, y=orig.spec$freq, z=t(orig.spec$amp),
                              levels=seq(-dBlevel,-dBlevel,1))

    # vectors to store minimum and maximum time values
    min.spec <- numeric(length(cont.spec))
    max.spec <- numeric(length(cont.spec))

    # minimum and maximum time values among contours detected
    for(i in 1:length(min.spec)){min.spec[i] <- min(cont.spec[[i]]$x)}
    for(i in 1:length(max.spec)){max.spec[i] <- max(cont.spec[[i]]$x)}

    # minimum and maximum time values
    t.min <- min(min.spec)
    t.max <- max(max.spec)

    if(t.min==0)
    stop("Background noise is likely interfering with the alighment. Consider using 'flim' argument to focus on acoustic signals")

    if((t.min-(time.perc*time.length))<0)
      stop("Time percentage is too large. Consider a smaller value of 'time.perc'")

    # cut Wave file using minimum and maximum time values
    short.wav0 <- seewave::deletew(orig.wav, f=f, output = "Wave",
                          from = (t.max+(time.perc*time.length)), to = max(orig.spec$time))

    short.wav <- seewave::deletew(short.wav0, f=f, output = "Wave",
                         from = 0, to = (t.min-(time.perc*time.length)))


    # Add silence to fill sound window
    final.wav <- seewave::addsilw(short.wav, f=f, at="end", d=time.length, output = "Wave")

    tuneR::writeWave(final.wav, file.path(wav.at, wav.to, file), extensible = F)

    } #end loop

  } #end function






