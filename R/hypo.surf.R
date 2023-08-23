#' Hypothetical 3D sound surfaces representing a sample of sound waves
#'
#' @description
#' Using the coordinates acquired by \code{eigensound(analysis.type = "threeDshape")}, this function creates 3D plots containing hypothetical sound surfaces that represent either the mean shape configuration (consensus), or minimum and maximum deformations relative to Principal Components in a Principal Components Analysis (PCA).
#'
#' \strong{Note:} The output of \code{hypo.surf} must be interpreted along with the ordination of Principal Components (e.g. \code{\link{pca.plot}}), both featuring the same object used for \code{threeD.out} argument. By doing so, \code{hypo.surf} enhance the comprehension on how sound shape changed along the ordination plot .
#'
#' @param threeD.out the output of \code{\link{eigensound}} analysis with \code{analysis.type = "threeDshape"}. By default: \code{threeD.out = NULL} (i.e. output must be specified before ploting)
#' @param PC Principal Component intended for the plot. Alternativaly, it is also possible to create mean shape configuration (consensus) from sample \code{PC = "mean"}. By default: \code{PC = 1}
#' @param flim modifications of the frequency limits (Y-axis). Vector with two values in kHz. Should be the same employed on \code{eigensound(analysis.type="threeDshape")} By default: \code{flim = c(0, 10)}
#' @param tlim modifications of the time limits (X-axis). Vector with two values in seconds. Should be the same employed on \code{eigensound(analysis.type="threeDshape")}. By default: \code{tlim = c(0, 1)}
#' @param x.length length of sequence (i.e. number of cells per side on sound window) to be used as sampling grid coordinates on the time (X-axis). Should be the same employed on \code{eigensound(analysis.type="threeDshape")}. By default: \code{x.length = 70}
#' @param y.length length of sequence (i.e. number of cells per side on sound window) to be used as sampling grid coordinates on the frequency (Y-axis). Should be the same employed on \code{eigensound(analysis.type="threeDshape")}. By default: \code{y.length = 47}
#' @param log.scale a logical. If \code{TRUE}, \code{hypo.surf} will use a logarithmic scale on the time (X-axis), which is recommeded when the analyzed sounds present great variation on this axis (e.g. emphasize short duration sounds). If \code{FALSE}, a linear scale is used instead (same as MacLeod et al., 2013). Should be the same employed on \code{eigensound(analysis.type="threeDshape")}. By default: \code{log.scale = TRUE}
#' @param f sampling frequency of \code{".wav"} files (in Hz). Should be the same employed on \code{eigensound(analysis.type="threeDshape")}. By default: \code{f = 44100}
#' @param wl length of the window for the analysis. Should be the same employed on \code{eigensound(analysis.type="threeDshape")}. By default: \code{wl = 512}
#' @param ovlp overlap between two successive windows (in %) for increased spectrogram resolution. Should be the same employed on \code{eigensound(analysis.type="threeDshape")}. By default: \code{ovlp = 70}
#' @param plot.exp a logical. If \code{TRUE}, exports the 3D output plot containing mean shape (\code{PC = "mean"}), or minimum and maximum deformations for the desired Principal Component (e.g. \code{PC = 1}). Exported plot will be stored on the folder indicated by \code{store.at}. By default: \code{plot.exp = FALSE}
#' @param plot.as only applies when \code{plot.exp = TRUE}. \code{plot.as = "jpeg"} will generate compressed images for quick inspection; \code{plot.as = "tiff"} or \code{"tif"} will generate uncompressed high resolution images that can be edited and used for publication. By default: \code{plot.as = "jpeg"}
#' @param store.at only applies when \code{plot.exp = TRUE}. Filepath to the folder where output plots will be stored. Should be presented between quotation marks. By default: \code{store.at = NULL} (i.e. user must specify the filepath where plots of hypothetical sound surfaces will be stored)
#' @param rotate.Xaxis rotation of the X-axis. Same as \code{theta} from \code{\link{persp3D}} (\code{\link{plot3D}} package). By default: \code{rotate.Xaxis = 60}
#' @param rotate.Yaxis rotation of the Y-axis. Same as \code{phi} from \code{\link{persp3D}} (\code{\link{plot3D}} package). By default: \code{rotate.Yaxis = 40}
#' @param cex.axis similarly as in \code{\link{par}}, magnification to be used for axis values. By default: \code{cex.axis = 0.9}
#' @param cex.lab similarly as in \code{\link{par}}, magnification to be used for x and y labels. By default: \code{cex.lab = 1.2}
#' @param cex.main similarly as in \code{\link{par}}, magnification to be used for main titles. By default: \code{cex.main = 1.3}
#' @param lwd Similarly as in \code{\link{par}}, intended line width for sampling grid. By default: \code{lwd = 0.1}
#' @param xlab a character string indicating the label to be written on the *x*-axis. By default: \code{xlab="Time coordinates"}
#' @param ylab a character string indicating the label to be written on the *y*-axis. By default: \code{ylab="Frequency coordinates"}
#' @param zlab a character string indicating the label to be written on the *z*-axis. By default: \code{zlab="Amplitude"}
#' @param colkey Similarly as \code{\link{plot3D}}, a list with parameters for the color key (legend). By default: \code{colkey = list(plot = TRUE, cex.clab = 0.9, cex.axis = 0.8, side = 4, length = 0.5, width = 0.7, labels = TRUE, tick = TRUE, lty = 1, lwd = 1, lwd.ticks = 1)}. See also \code{\link{colkey}}
#'
#' @note
#' Some codes from \code{hypo.surf} were adapted from \code{plotTangentSpace} function (\code{\link{geomorph}} package version 3.1.2), which is now deprecated and replaced by current functions \code{\link{gm.prcomp}}, \code{\link{summary.gm.prcomp}} and \code{\link{plot.gm.prcomp}}. More specifically, the code chunk related to the acquisition of hypothetical point configurations from each PC (i.e. warp grids) was the same as in \code{plotTangentSpace}. However, the hypothetical configurations from \code{plotTangentSpace} were plotted along with ordination of PCs, whereas \code{hypo.surf} focuses solely on hypothetical 3D surfaces that represent minimum, maximum and mean deformations relative to each PCs.
#'
#'
#'
#'
#' @seealso
#' \code{\link{gm.prcomp}}, \code{\link{summary.gm.prcomp}}, \code{\link{plot.gm.prcomp}}, \code{\link{geomorph}}, \code{\link{eigensound}}, \code{\link{pca.plot}}
#'
#' Useful links:
#' \itemize{
#'   \item{\url{https://github.com/p-rocha/SoundShape}}
#'   \item{Report bugs at \url{https://github.com/p-rocha/SoundShape/issues}}}

#' @author
#' Pedro Rocha
#'
#'
#' @examples
#' data(eig.sample)
#'
#' # PCA using three-dimensional semilandmark coordinates
#' pca.eig.sample <- stats::prcomp(geomorph::two.d.array(eig.sample))
#'
#' # Verify names for each acoustic unit and the order in which they appear
#' dimnames(eig.sample)[[3]]
#'
#' # Create factor to use as groups in subsequent ordination plot
#' sample.gr <- factor(c(rep("centralis", 3), rep("cuvieri", 3), rep("kroyeri", 3)))
#'
#' # Clear current R plot to prevent errors
#' grDevices::dev.off()
#'
#' # Plot result of Principal Components Analysis
#' pca.plot(PCA.out = pca.eig.sample, groups = sample.gr, conv.hulls = sample.gr,
#'          main="PCA of 3D coordinates", leg=TRUE, leg.pos = "top")
#'
#' # Verify hypothetical sound surfaces using hypo.surf
#' hypo.surf(threeD.out=eig.sample, PC=1, flim=c(0, 4), tlim=c(0, 0.8), x.length=70, y.length=47)
#'
#'
#' @references
#' MacLeod, N., Krieger, J. & Jones, K. E. (2013). Geometric morphometric approaches to acoustic signal analysis in mammalian biology. \emph{Hystrix, the Italian Journal of Mammalogy, 24}(1), 110-125.
#'
#' Rocha, P. & Romano, P. (2021) The shape of sound: A new \code{R} package that crosses the bridge between Bioacoustics and Geometric Morphometrics. \emph{Methods in Ecology and Evolution, 12}(6), 1115-1121.
#'
#'
#' @export
#'
hypo.surf <- function(threeD.out=NULL, PC=1, flim=c(0, 4), tlim=c(0, 0.8), x.length=70, y.length=47, log.scale=TRUE, f=44100, wl=512, ovlp=70, plot.exp=FALSE, plot.as="jpeg", store.at = NULL, rotate.Xaxis=60, rotate.Yaxis=40, cex.axis=0.5, cex.lab=0.9, cex.main=1.1, lwd=0.1, xlab="Time coordinates", ylab="Frequency coordinates", zlab="Amplitude", colkey = list(plot = TRUE, cex.clab = 0.9, cex.axis = 0.8, side = 4, length = 0.5, width = 0.7, labels = TRUE, tick = TRUE, lty = 1, lwd = 1, lwd.ticks = 1))
{

  if(is.null(threeD.out))
  {stop('Please define threeD.out as an object containing the output of eigensound function with analysis.type="threeDshape')}

  if(is.null(PC))
  {stop('Choose a Principal Component or "mean" shape configuration to generate hypothetical semilandmark configurations. See help(hypo.surf) for details.')}

  if(plot.exp == TRUE && is.null(store.at)){
    stop("Use 'store.at' to specify folder path where plots of hypothetical sound surfaces will be stored")}

  ref = geomorph::mshape(threeD.out)

  # Use sample as reference for hypothetical surfaces
  #data("cuvieri")

  # Acquire spectrogram values
  e <- seewave::spectro(SoundShape::cuvieri, f=f, wl=wl, ovlp=ovlp, flim=flim, tlim=tlim, plot=F)

  # Create new sequences representing the sampling grid
  freq.seq <- seq(1, length(e$freq), length=y.length)
  ifelse(isTRUE(log.scale),
         time.seq <- 10^(seq(log10(1), log10(length(e$time)), length.out = x.length)),#log
         time.seq <- seq(1, length(e$time), length.out = x.length)) # linear scale

  # Subset original coordinates using new sequences
  time.sub <- e$time[time.seq]
  freq.sub <- e$freq[freq.seq]


  # Mean (consensus) shape configuration
  if(is.character(PC)){
    if(PC=="mean"){
      mean.hPC <- matrix(ref[,3], nrow=y.length, ncol=x.length, byrow = T)

      # Export plot
      if(plot.exp==TRUE){

        oldpar <- graphics::par(no.readonly = TRUE)
        base::on.exit(graphics::par(oldpar))

        if(plot.as == "jpeg"){
          grDevices::jpeg(width=4000,height=3500,units="px",res=500,
                          filename=paste(store.at,"/", "Hypothetical surfaces - Mean shape", ".jpg", sep=""))} # compressed images
        if(plot.as=="tiff"|plot.as=="tif"){
          grDevices::tiff(width=4000, height=3500, units="px", res=500,
                          filename=paste(store.at,"/", "Hypothetical surfaces - Mean shape",".tif", sep=""))} # uncompressed images
      } # end plot.exp=TRUE

      # Mean
      graphics::par(mar=c(1,2,2,2))
      plot3D::persp3D(x=time.sub, y=freq.sub, z=t(mean.hPC), theta=rotate.Xaxis, phi=rotate.Yaxis, resfac=1, r=3,
                      expand=0.5, scale=T, axes=T, ticktype="detailed", nticks=4,
                      cex.axis=cex.axis, cex.lab=cex.lab, cex.main=cex.main,
                      border="black", lwd=lwd, col=seewave::spectro.colors(n=100),
                      xlab=xlab, ylab=ylab, zlab=zlab, main="Mean shape", clab=expression('Amplitude'),
                      colkey=colkey)

      if(plot.exp==TRUE){grDevices::dev.off()}

    }# end "mean" shape

    else{stop('Invalid or misspelled PC argument. Please choose between "mean" (for mean shape configuration) or a number representing the intended Principal Component. See help(hypo.surf) for details.')}

  } # end is.character

  if(is.numeric(PC)){

    PCmax <- PC*2
    PCmin <- (PC*2)-1

    # Principal Components Analysis (PCA) using point coordinates
    PCA.out = stats::prcomp(geomorph::two.d.array(threeD.out))

    # Acquisition of hypothetical shapes (codes adapted from geomorph::plotTangentSpace)
    pcdata = PCA.out$x
    p = dim(threeD.out)[1]
    k = dim(threeD.out)[2]

    shapes <- shape.names <- NULL
    for (i in 1:ncol(pcdata)) {
      pcaxis.min <- min(pcdata[, i])
      pcaxis.max <- max(pcdata[, i])
      pc.min <- pc.max <- rep(0, dim(pcdata)[2])
      pc.min[i] <- pcaxis.min
      pc.max[i] <- pcaxis.max
      pc.min <- as.matrix(pc.min %*% (t(PCA.out$rotation))) +
        as.vector(t(ref))
      pc.max <- as.matrix(pc.max %*% (t(PCA.out$rotation))) +
        as.vector(t(ref))
      shapes <- rbind(shapes, pc.min, pc.max)
      shape.names <- c(shape.names, paste("PC", i, "min", sep = ""),
                       paste("PC", i, "max", sep = ""))
    }
    shapes <- geomorph::arrayspecs(shapes, p, k)
    shapes <- lapply(seq(dim(shapes)[3]), function(x) shapes[, , x])
    names(shapes) <- shape.names

    if(PC > (length(shapes)/2)){
      stop('Invalid number of Principal Components')}

    # Transform calculated PC into matrix of amplitude values
    min.hPC <- matrix(shapes[[PCmin]][,3], nrow=y.length, ncol=x.length, byrow = T)
    max.hPC <- matrix(shapes[[PCmax]][,3], nrow=y.length, ncol=x.length, byrow = T)


    # Hypothetical sound surfaces for desired Principal Component
    if(plot.exp==TRUE){
      if(plot.as == "jpeg"){
        grDevices::jpeg(width=8000,height=3500,units="px",res=500,
                        filename=paste(store.at,"/", "Hypothetical surfaces - PC", PC, ".jpg", sep=""))} # compressed images
      if(plot.as=="tiff"|plot.as=="tif"){
        grDevices::tiff(width=8000, height=3500, units="px", res=500,
                        filename=paste(store.at,"/", "Hypothetical surfaces - PC", PC,".tif", sep=""))} # uncompressed images
    } # end plot.exp=TRUE
    graphics::par(mfrow=c(1,2), mar=c(0,2,2,2))

    # Minimum
    plot3D::persp3D(x=time.sub, y=freq.sub, z=t(min.hPC), theta=rotate.Xaxis, phi=rotate.Yaxis, resfac=1, r=3,
                    expand=0.5, scale=T, axes=T, ticktype="detailed", nticks=4,
                    cex.axis=cex.axis, cex.lab=cex.lab, cex.main=cex.main,
                    border="black", lwd=lwd, col=seewave::spectro.colors(n=100),
                    xlab=xlab, ylab=ylab, zlab=zlab,
                    clab=expression('Amplitude'), main=paste("PC", PC, " - Minimum", sep=""),
                    colkey=colkey)

    # Maximum
    plot3D::persp3D(x=time.sub, y=freq.sub, z=t(max.hPC), theta=rotate.Xaxis, phi=rotate.Yaxis, resfac=1, r=3,
                    expand=0.5, scale=T, axes=T, ticktype="detailed", nticks=4,
                    cex.axis=cex.axis, cex.lab=cex.lab, cex.main=cex.main,
                    border="black", lwd=lwd, col=seewave::spectro.colors(n=100),
                    xlab=xlab, ylab=ylab, zlab=zlab,
                    clab=expression('Amplitude'), main=paste("PC", PC, " - Maximum", sep=""),
                    colkey=colkey)
    if(plot.exp==TRUE){grDevices::dev.off()}

  } # end is.numeric


} #end function
