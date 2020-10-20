#' Sound waves onto morphometric data
#'
#' @description \code{eigensound} is the main feature of \code{SoundShape} package. For each \code{".wav"} file on a given folder, the fuction will compute spectrogram data and acquire semilandmarks using a 3D representation of sound (\code{analysis.type = "threeDshape"}), allowing its users to acquire, and simultaneously store point coordinates (i.e. semilandmarkns) as an R object, and/or in TPS format – the native file format of James Rohlf’s TPS series (Rohlf, 2015).
#'
#' Moreover, \code{eigensound} also allow its user to export 2D and 3D spectrogram images (\code{plot.exp = TRUE}) that are helpful during the protocol for error verification and for illustrative purposes (see Rocha & Romano \emph{in prep}). Alternativaly, \code{eigensound} feature the option of acquiring semilandmarks as the cross-correlation between energy quantiles and a curve of relative amplitude from 2D spectrograms (\code{analysis.type = "twoDshape"}; see \code{Details} section).
#'
#' @param analysis.type type of analysis intended. If \code{analysis.type = "threeDshape"}, semilandmarks are acquired from spectrogram data using a 3D representation of sound (same as in MacLeod et al., 2013). If \code{analysis.type = "twoDshape"} and \code{add.points = TRUE}, semilandmarks are acquired using energy quantiles and a 2D curve of relative amplitude. By default: \code{analysis.type = NULL} (i.e. method must be specified before the analysis).
#' @param wav.at filepath to the folder where \code{".wav"} files are stored. Should be presented between quotation marks. By default: \code{wav.at = NULL} (i.e. user must specify the filepath to \code{".wav"} files)
#' @param store.at filepath to the folder where spectrogram plots and \code{tps} file will be stored. Should be presented between quotation marks. By default: \code{store.at = wav.at} (i.e. store outputs in the same folder as \code{".wav"} files)
#' @param dBlevel absolute amplitude value to be used as relative amplitude contour, which will serve as reference for semilandmark acquisition in both \code{analysis.type = "threeDshape"} and \code{"twoDshape"}. By default: \code{dBlevel = 25}
#' @param flim modifications of the frequency limits (Y-axis). Vector with two values in kHz. By default: \code{flim = c(0, 10)}
#' @param tlim modifications of the time limits (X-axis). Vector with two values in seconds. By default: \code{tlim = c(0, 1)}
#' @param trel only applies when \code{analysis.type = "twoDshape"}. Set the relative scale to be used on the time (X-axis); relative to the numbers displayed on the X-axis of spectrogram plots. By default: \code{trel = tlim}
#' @param x.length only applies when \code{analysis.type = "threeDshape"}. Length of sequence (i.e. number of cells per side on sound window) to be used as sampling grid coordinates on the time (X-axis).  By default: \code{x.length = 80}
#' @param y.length only applies when \code{analysis.type = "threeDshape"}. Length of sequence (i.e. number of cells per side on sound window) to be used as sampling grid coordinates on the frequency (Y-axis). By default: \code{y.length = 60}
#' @param log.scale only applies when \code{analysis.type = "threeDshape"}. A logical. If \code{TRUE}, \code{eigensound} will use a logarithmic scale on the time (X-axis), which is recommeded when the analyzed sounds present great variation on this axis (e.g. emphasize short duration sounds). If \code{FALSE}, a linear scale is used instead (same as MacLeod et al., 2013). By default: \code{log.scale = TRUE}
#' @param back.amp only applies when \code{analysis.type = "twoDshape"} and  \code{plot.exp = TRUE}. Absolute amplitude value to be used as background in the 2D spectrogram plot. By default: \code{back.amp = 35}
#' @param add.points only applies when \code{analysis.type = "twoDshape"}. A logical. If \code{TRUE}, \code{eigensound} will compute semilandmarks acquired by cross-correlation between energy quantiles (i.e. \code{EQ}) and a curve of relative amplitude (i.e. \code{dBlevel}). If \code{plot.exp = TRUE}, semilandmarks will be included in spectrogram plots. By default: \code{add.points = FALSE} (see \code{Details})
#' @param add.contour only applies when \code{analysis.type = "twoDshape"} and  \code{plot.exp = TRUE}. A logical. If \code{TRUE}, exported spectrogram plots will include the curves of relative amplitude at the level specified by \code{dBlevel}. By default: \code{add.contour = TRUE}
#' @param lwd only applies when \code{plot.exp = TRUE} and \code{add.contour = TRUE}. The line width for the curves of relative amplitude at the level specified by \code{dBlevel}. Same as in \code{\link{par}}. By default: \code{lwd = 1}
#' @param EQ only applies when \code{analysis.type = "twoDshape"} and \code{add.points = TRUE}. A vector of energy quantiles intended (with \code{0 < EQ < 1}). By default: \code{EQ = c(0.05, 0.15, 0.3, 0.5, 0.7, 0.85, 0.95)} \strong{Note:} When dealing with narrow banded calls, consider reducing the number of quantiles to prevent errors in the analysis.
#' @param mag.time only applies when \code{analysis.type = "twoDshape"}. Optional argument for magnifying the time coordinates (X-axis). This is sometimes desired for small sound windows (e.g. less than 1 s), in which the time coordinates will be on a different scale than that of frequency coordinates. In those cases, it is recommended to include \code{mag.time = 10} or \code{mag.time = 100}, depending on the lenght of sound window. By default: \code{mag.time = 1} (i.e. no magnification is performed)
#' @param f sampling frequency of \code{Wave}'s for the analysis (in Hz). By default: \code{f = 44100}
#' @param wl length of the window for the analysis. By default: \code{wl = 512}
#' @param ovlp overlap between two successive windows (in %) for increased spectrogram resolution. By default: \code{ovlp = 70}
#' @param plot.exp a logical. If \code{TRUE}, for each \code{".wav"} file on the folder indicated by \code{wav.at}, \code{eigensound} will store a spectrogram image on the folder indicated by \code{store.at}. Depending on the \code{analysis.type}, plots may consist of 2D or 3D spectrogram images. By default: \code{plot.exp = TRUE}
#' @param plot.as only applies when \code{plot.exp = TRUE}. \code{plot.as = "jpeg"} will generate compressed images for quick inspection of semilandmarks; \code{plot.as = "tiff"} or \code{"tif"} will generate uncompressed high resolution images that can be edited and used for publication. By default: \code{plot.as = "jpeg"}
#' @param plot.type only applies when \code{analysis.type = "threeDshape"} and  \code{plot.exp = TRUE}. \code{plot.type = "surface"} will produce simplified 3D sound surfaces from the calculated semilandmarks (same output employed by MacLeod et al., 2013); \code{plot.type = "points"} will produce 3D graphs with semilandmarks as points. By default: \code{plot.type = "surface"}
#' @param rotate.Xaxis only applies when \code{analysis.type = "threeDshape"} and  \code{plot.exp = TRUE}. Rotation of the X-axis. Same as \code{theta} from \code{\link{persp3D}} (\code{\link{plot3D}} package). By default: \code{rotate.Xaxis = 60}
#' @param rotate.Yaxis only applies when \code{analysis.type = "threeDshape"} and  \code{plot.exp = TRUE}. Rotation of the Y-axis. Same as \code{phi} from \code{\link{persp3D}} (\code{\link{plot3D}} package). By default: \code{rotate.Yaxis = 40}
#' @param TPS.file Desired name for the \code{tps} file containing semilandmark coordinates. Should be presented between quotation marks. \strong{Note:} Whenever \code{analysis.type = "twoDshape"}, it will only work if \code{add.points = TRUE}. By default: \code{TPS.file = NULL} (i.e. prevents \code{eigensound} from creating a \code{tps} file)
#'
#' @details
#' When \code{analysis.type = "twoDshape"} and \code{add.points = TRUE}, \code{eigensound} will compute semilandmarks acquired by cross-correlation between energy quantiles (i.e. \code{EQ}) and a curve of relative amplitude (i.e. \code{dBlevel}). However, this is often subtle and prone to incur in errors (e.g. bad alignment of acoustic units; inappropriate X and Y coordinates for the sound window; narrow banded calls). Therefore, a more robust protocol of error verification is achieved using \code{add.points = FALSE} and \code{add.contour = TRUE} (default), which allow for quick verification of acoustic units alignment and the shape of each curve of relative amplitude (specified by \code{dBlevel}).
#'
#' @note
#' In order to store the results from \code{eigensound} function and proceed with the Geometric Morphometric steps of the analysis (e.g. \code{\link{geomorph}} package; Adams et al., 2017), one can simultaneosly assign the function's output to an \code{R} object and/or store them as a \code{tps} file to be used by numerous softwares of geometric analysis of shape, such as the TPS series (Rohlf, 2015).
#'
#' Additionally, one may also export 2D or 3D plots as \code{jpeg} (compressed image) or \code{tiff} (uncompressed image) file formats, which can be edited for publication purposes.
#'
#'
#' @author
#' Pedro Rocha
#'
#' @references
#' Adams, D. C., M. L. Collyer, A. Kaliontzopoulou & Sherratt, E. (2017) \emph{Geomorph: Software for geometric morphometric analyses}. R package version 3.0.5. https://cran.r-project.org/package=geomorph.
#'
#' MacLeod, N., Krieger, J. & Jones, K. E. (2013). Geometric morphometric approaches to acoustic signal analysis in mammalian biology. \emph{Hystrix, the Italian Journal of Mammalogy, 24}(1), 110-125.
#'
#' Rocha, P. & Romano, P. (\emph{in prep}) The shape of sound: A new \code{R} package that crosses the bridge between Bioacoustics and Geometric Morphometrics.
#'
#' Rohlf, F.J. (2015) The tps series of software. \emph{Hystrix 26}, 9-12.
#'
#' @seealso
#' \code{\link{align.wave}}, \code{\link{geomorph}}, \code{\link{seewave}}
#'
#' Useful links:
#' \itemize{
#'   \item{\url{https://github.com/p-rocha/SoundShape}}
#'   \item{Report bugs at \url{https://github.com/p-rocha/SoundShape/issues}}}
#'
#' @examples
#' library(seewave)
#' library(tuneR)
#'
#' # Create temporary folder to store ".wav" files
#' wav.at <- file.path(base::tempdir(), "eigensound")
#' if(!dir.exists(wav.at)) dir.create(wav.at)
#'
#' # Create temporary folder to store results
#' store.at <- file.path(base::tempdir(), "eigensound-output")
#' if(!dir.exists(store.at)) dir.create(store.at)
#'
#' # Cut acoustic units from original Wave
#' cut.cuvieri <- cutw(cuvieri, f=44100, from=0, to=0.9, output = "Wave")
#' cut.centralis <- cutw(centralis, f=44100, from=0, to=0.9, output = "Wave")
#' cut.kroyeri <- cutw(kroyeri, f=44100, from=0.2, to=1.1, output = "Wave")
#'
#' # Export ".wav" files containing acoustic units and store on previosly created folder
#' writeWave(cut.cuvieri, filename = file.path(wav.at, "cut.cuvieri.wav"), extensible = FALSE)
#' writeWave(cut.centralis, filename = file.path(wav.at, "cut.centralis.wav"), extensible = FALSE)
#' writeWave(cut.kroyeri, filename = file.path(wav.at, "cut.kroyeri.wav"), extensible = FALSE)
#'
#'
#' \donttest{
#' # Create 2D spectrograms using analysis.type = "twoDshape"
#' eigensound(analysis.type = "twoDshape", flim=c(0, 4), tlim=c(0, 0.8),
#'            plot.exp=TRUE, wav.at = wav.at, store.at = store.at)
#'
#' # Create 3D spectrograms using analysis.type = "threeDshape" and store point coordinates
#' eig.data <- eigensound(analysis.type = "threeDshape", plot.exp=TRUE,
#'               wav.at = wav.at, store.at = store.at, flim=c(0, 4), tlim=c(0, 0.8))
#' }
#'
#' @export
#'
eigensound <- function(analysis.type = NULL, wav.at = NULL, store.at = wav.at, dBlevel=25, flim=c(0, 10), tlim = c(0, 1), trel=tlim, x.length=80, y.length=60, log.scale=TRUE, back.amp=35, add.points=FALSE, add.contour=TRUE, lwd=1, EQ= c(0.05, 0.15, 0.3, 0.5, 0.7, 0.85, 0.95), mag.time=1, f=44100, wl=512, ovlp=70, plot.exp = TRUE, plot.as = "jpeg", plot.type = "surface", rotate.Xaxis=60, rotate.Yaxis=40, TPS.file = NULL){

  # Invalid arguments
  if(is.null(wav.at)) {stop("Use 'wav.at' to specify folder path where '.wav' files are stored")}
  if(is.null(analysis.type)) {stop("Method undefined in 'analysis.type'")}
  if(!is.null(analysis.type)){
    if(analysis.type != "threeDshape" && analysis.type != "twoDshape")
      stop("Invalid analysis specified by 'analysis.type'. Choose between 'threeDshape' or 'twoDshape'")  } # invalid analysis.type

  if(plot.exp == TRUE){
    if(plot.as != "jpeg" && plot.as != "tif" && plot.as != "tiff")
      stop("Invalid image format specified by 'plot.as'. Choose between 'jpeg', 'tif' or 'tiff'") } # invalid plot.as

  if(plot.exp == TRUE && analysis.type == "threeDshape"){
    if(plot.type != "surface" && plot.type != "points")
      stop("Invalid 3D sound shape specified by 'plot.type'. Choose between 'surface' or 'points'") } # invalid plot.type

  if(analysis.type == "twoDshape"){
    if(!is.null(TPS.file) && add.points==FALSE)
      stop("Must set add.points = TRUE to acquire 2D semilandmarks and store them on a '.tps' file")}

  # Create TPS file before analysis
  if(!is.null(TPS.file)){
    TPS.path <- paste(store.at, "/", TPS.file, ".tps",  sep="")
    file.create(TPS.path, showWarnings = TRUE)
  } # create TPS file

  # threeDshape
  if(analysis.type == "threeDshape"){

    # Acquire 3D semilandmark coordinates for each ".wav" file on a given folder
    for(file in list.files(wav.at, pattern = ".wav"))
    {
      threeD <- tuneR::readWave(paste(wav.at,"/", file, sep=""))

      # Store spectrogram as R object
      e <-  seewave::spectro(threeD, f=f, wl=wl, ovlp=ovlp, flim=flim, tlim=tlim, plot=F)

      # Create sequences to use as new grid
      freq.seq <- seq(1, length(e$freq), length = y.length)

      ifelse(isTRUE(log.scale),
             time.seq <- 10^(seq(log10(1), log10(length(e$time)), length.out = x.length)),#log
             time.seq <- seq(1, length(e$time), length.out = x.length)) # linear scale

      # Subset original coordinates by new sequences
      time.sub <- e$time[time.seq]
      freq.sub <- e$freq[freq.seq]

      # Subset matrix of amplitude values using new sequences
      amp.sub <- e$amp[freq.seq, time.seq]

      # Reassign background amplitude values
      for(i in 1:length(amp.sub)){if(amp.sub[i] == -Inf |amp.sub[i] <= -dBlevel)
      {amp.sub[i] <- -dBlevel}}

      # Assign time and frequency coordinates as column and row names of amplitude matrix
      colnames(amp.sub) <- time.sub
      rownames(amp.sub) <- freq.sub

      # Transform amplitude matrix into semilandmark 3D coordinates
      ind.3D <- as.matrix(stats::setNames(reshape2::melt(t(amp.sub)), c('time', 'freq', 'amp')))

      # Store semilandmark coordinates in array
      ifelse(!exists('coord3D'),
             coord3D <- array(data=ind.3D, dim=c(dim(ind.3D), 1)),
             coord3D <- abind::abind(coord3D, ind.3D, along=3))

      # Plot semilandmarks as points or sound surface
      if(plot.exp==TRUE){
        if(plot.as == "jpeg"){grDevices::jpeg(width =5000,height = 3500, units = "px", res = 500,
                                   filename=paste(store.at,"/", sub(".wav","",file), ".jpg", sep=""))} # compressed images
        if(plot.as=="tiff"|plot.as=="tif"){grDevices::tiff(width=5000, height=3500, units="px", res=500,
                                                filename=paste(store.at, "/", sub(".wav", "", file), ".tif", sep=""))} # uncompressed images
        if(plot.type=="surface"){plot3D::persp3D(x=time.sub, y=freq.sub, z=t(amp.sub),
                                         border="black", lwd=0.1, theta=rotate.Xaxis, phi=rotate.Yaxis, resfac=1, r=3,expand=0.5, cex.axis=0.7,
                                         scale=T, axes=T, col=seewave::spectro.colors(n=100), ticktype="detailed", nticks=4,            xlab="Time (s)", ylab="Frequency (kHz)", zlab="Amplitude (dB)",
                                         main= sub(".wav", "", file), clab=expression('Amplitude dB'))}
        if(plot.type=="points"){plot3D::scatter3D(x=ind.3D[,1], y=ind.3D[,2], z=ind.3D[,3],
                                          pch=21, cex=0.5, theta=rotate.Xaxis, phi=rotate.Yaxis, resfac=1, r=3, expand=0.5, cex.axis=0.7,
                                          scale=T, axes=T, col=seewave::spectro.colors(n=100), ticktype="detailed", nticks=4,
                                          xlab="Time (s)", ylab="Frequency (kHz)", zlab="Amplitude (dB)",
                                          main= sub(".wav", "", file), clab=expression('Amplitude dB'))}

        grDevices::dev.off()  } # end threeDshape plot of each ".wav" file

      # Export semilandmark coordinates as TPS file
      if(!is.null(TPS.file)){
        lmline <- paste("LM=", dim(ind.3D)[1], sep="")
        write(lmline, TPS.path, append = TRUE)
        utils::write.table(ind.3D, TPS.path, col.names = FALSE, row.names = FALSE, append =TRUE)
        idline <- paste("ID=", sub(".wav", "", file), sep="")
        write(idline, TPS.path, append = TRUE)
        write("", TPS.path, append = TRUE)
        rm(lmline, idline)
      } #end tps file

    } #end threeDshape loop

    # Assign names to arrays' third dimension
    dimnames(coord3D)[[3]] <- sub(".wav","", list.files(wav.at, pattern = ".wav"))

  } # End threeDshape

  # twoDshape
  if(analysis.type == "twoDshape"){

    FC <- function(x) {
      sign <- ifelse(x >= 0, 1, -1)
      res <- abs(diff(sign))
      return(res)  } # FC function

    frequencies <- function (spec, EQ=NULL, f=NULL,mel=FALSE, flim=NULL)
    {
      fhz <- f
      if (!is.null(f) & mel) {f <- 2 * mel(f/2)}
      if (is.null(f)) {
        if (is.vector(spec))
          stop("'f' is missing")
        else if (is.matrix(spec))
          f <- spec[nrow(spec), 1] * 2000 * nrow(spec)/(nrow(spec) - 1) }

      if (is.matrix(spec)) {
        freq <- spec[, 1]
        freq = freq * 1000
        spec <- spec[, 2] }

      L <- length(spec)
      wl <- L * 2
      if (any(spec < 0))
        stop("The frequency spectrum to be analysed should not be in dB")
      if (!is.null(flim)) {
        if (flim[1] < 0 || flim[2] > fhz/2)
          stop("'flim' should range between 0 and f/2")
        if (mel)
          flim <- mel(flim * 1000)/1000 }
      else {flim <- c(0, f/2000)}

      spec <- spec[(flim[1] * 1000 * wl/f):(flim[2] * 1000 * wl/f)]
      amp <- spec/sum(spec)
      cumamp <- cumsum(amp)
      freq <- seq(from = flim[1] * 1000, to = flim[2] * 1000, length.out = length(spec))

      for(i in EQ) # Frequency values according to desired energy quantiles
      { if(!exists("Qlist")) {Qlist <- data.frame()}
        if(exists("Qlist")) {
          Qtemp <- data.frame("Quantiles" = paste("EQ-", i, sep=""),
                              "Frequency values" = freq[length(cumamp[cumamp <= i]) + 1])
          Qlist <- rbind(Qlist, Qtemp)
          rm(Qtemp) }  }

      results <- Qlist     } # frequencies function

    # Compute data for each file on a given folder
    for(file in list.files(wav.at, pattern = ".wav")) {
      Wav <- tuneR::readWave(paste(wav.at,"/", file, sep=""))

      if(add.points==TRUE){
        # Store spectrogram as R object
        s <- seewave::spectro(Wav,f=f,wl=wl,ovlp=ovlp,flim=flim,tlim=tlim,plot=F)

        # Selection of longest relative amplitude contour
        con <- grDevices::contourLines(x=s$time, y=s$freq, z=t(s$amp),
                                       levels=seq(-dBlevel, -dBlevel, 1))
        n <- numeric(length(con))
        for(i in 1:length(con)) n[i] <- length(con[[i]]$x)
        n.max <- which.max(n)
        con <- con[[n.max]]
        t.beg <- min(con$x)
        t.end <- max(con$x)
        f.min <- min(con$y) # Minimum frequency
        t.min <- con$x[con$y==f.min] # Time value associated
        f.max <- max(con$y) # Maximum frequency
        t.max <- con$x[con$y==f.max] # Time value associated

        # Power spectrum calculated within maximum and minimum frequency values
        pwr.spec <- seewave::spec(Wav, f=f, wl=wl, from=t.beg, to=t.end, plot=F)
        fvalues <- frequencies(pwr.spec, EQ=EQ, f=f, flim=c(f.min, f.max))
        fvalues$Frequency.values <- (fvalues$Frequency.values)/1000 # Convert to kHz

        # Crossing between frequency quantiles and relative amplitude contour
        for(i in seq_along(fvalues$Frequency.values)) {
          fc.res <- FC(con$y - fvalues$Frequency.values[i])
          t.pos <- which(fc.res==2)

          if(!exists("t1")|!exists("t2")){
            t1 <- vector(mode="integer", length=length(fvalues$Frequency.values))
            t2 <- vector(mode="integer", length=length(fvalues$Frequency.values)) }

          if(exists("t1")|exists("t2")){
            t1[i] <- min(con$x[t.pos])
            t2[i] <- max(con$x[t.pos]) }  } # end crossing EQ

        # Coordinates for plotting
        SM1 <- data.frame(cbind(fvalues, "Time"= t1))
        SM2 <- data.frame(cbind(fvalues, "Time"= t2))

        # Re-order SM1 so that points will be anticlockwise
        SM1 <- SM1[order(SM1[,2], decreasing=T), ]

        # Store SM-min and SM-max in a data frame
        SM.min <- data.frame("Quantiles" = "EQ-min",
                             "Frequency.values"= f.min, "Time"= t.min)
        SM.max <- data.frame("Quantiles" = "EQ-max",
                             "Frequency.values"= f.max, "Time"= t.max)

        # Semilandmark numbers
        SM.numbers <- as.factor(1:((length(EQ)*2) + 2)) # Twice the length of EQ plus the SM-min and SM-max

        # Complete data.frame with coordinates of x and y
        SM <- cbind(SM.numbers, data.frame(rbind(SM.max, SM1, SM.min, SM2)))

        # Magnify time coordinates
        SM[,4] <- SM[,4]*mag.time

        # Store semilandmark coordinates as array of matrices
        SM.mtx <- as.matrix(SM[,4:3])
        ifelse(!exists('coord2D'),
               coord2D <- array(data=SM.mtx, dim=c(dim(SM.mtx), 1)),
               coord2D <- abind::abind(coord2D, SM.mtx, along=3))
      } # end add.points = TRUE

      # Plot spectrogram
      if(plot.exp == TRUE){
        if(plot.as == "jpeg"){grDevices::jpeg(width =5000,height = 3500, units = "px", res = 500, filename=paste(store.at,"/", sub(".wav","",file), ".jpg", sep="")) } # compressed images
        if(plot.as=="tiff"|plot.as=="tif"){grDevices::tiff(width=5000, height=3500, units="px", res=500, filename=paste(store.at, "/", sub(".wav", "", file), ".tif", sep="")) } # uncompressed images

        if(add.points==FALSE){
          seewave::spectro(Wav, f=f, wl=wl, ovlp=ovlp, osc=F, scale=F, grid=F,
                      main= sub(".wav", "", file), trel=trel, flim=flim, tlim=tlim,
                      cont=add.contour,contlevels=seq(-dBlevel,-dBlevel,1),lwd=lwd,
                      collevels=seq(-back.amp, 0, 0.1))
        } # end contour

        if(add.points==TRUE){
          seewave::spectro(Wav, f=f, wl=wl, ovlp=ovlp, osc=F, scale=F, grid=F,
                    main= sub(".wav", "", file), trel=trel, flim=flim, tlim=tlim,
                    collevels=seq(-back.amp, 0, 0.1))
          if(add.contour==TRUE){graphics::polygon(x=con$x, y=con$y, lwd=lwd)}
          graphics::points(x=(SM$Time)/mag.time, y=SM$Frequency.values, pch=21, col="black", bg="red", cex=1.1)
          graphics::text(x=(SM$Time)/mag.time, y=SM$Frequency.values,cex=0.7, labels=SM$SM.numbers, pos=c(3, rep(2, dim(SM1)[1]), 1, rep(4, dim(SM2)[1]))) ## SM numbers
        } # end points

        grDevices::dev.off()
      } # end plot

      # Export semilandmark coordinates as TPS file
      if(!is.null(TPS.file)){
        lmline <- paste("LM=", dim(SM)[1], sep="")
        write(lmline, TPS.path, append = TRUE)
        utils::write.table(SM[,4:3], TPS.path, col.names = FALSE, row.names = FALSE, append = TRUE)
        idline <- paste("ID=", sub(".wav", "", file), sep="")
        write(idline, TPS.path, append = TRUE)
        write("", TPS.path, append = TRUE)
        rm(lmline, idline, SM)
      } # end tps file

    } # end loop

    # Assign names to arrays' third dimension
    if(add.points==TRUE){
      dimnames(coord2D)[[3]] <- sub(".wav","", list.files(wav.at, pattern = ".wav"))} # end add.poins = TRUE

  } # End twoDshape


  if(analysis.type == "twoDshape" && add.points == TRUE){ coord <- coord2D }
    else( coord <- NULL )# twoDshape results
  if(analysis.type == "threeDshape"){ coord <- coord3D } # threeDshape results

  results <- coord

} # End eigensound function
