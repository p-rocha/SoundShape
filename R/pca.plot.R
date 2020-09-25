#' Plot ordination of Principal Components with convex hulls
#'
#' @description
#' Ordination of Principal Components from the output of a Principal Components Analysis performed by \code{\link{prcomp}} function (\code{\link{stats}} package). Require a \code{factor} object with \code{groups}, which enable the plot to feature colored groups and convex hulls (if desired).
#'
#' @author Pedro Rocha
#'
#' @param PCA.out the output of a Principal Components Analysis performed by \code{\link{prcomp}} (\code{\link{stats}} package). By default: \code{PCA.out = NULL} (i.e. output must be specified before ploting)
#' @param groups groups to use as colors and/or convex hulls. Must be a \code{factor} object with the same length as the number of rows in the coordinates of \code{PCA.out} (i.e. \code{length(groups) == nrow(PCA.out$x)}). By default: \code{groups = NULL} (i.e. \code{factor} must be specified before ploting)
#' @param conv.hulls groups to use for convex hulls. Must be a \code{factor} object with the same length as the number of rows in the coordinates of \code{PCA.out} (i.e. \code{length(conv.hulls) == nrow(PCA.out$x)}). By default: \code{conv.hulls = NULL} (i.e. plot without convex hulls)
#' @param col.gp a \code{factor} object with the colours intended for \code{groups}. Must be the same length as the number of levels from \code{groups} (i.e. \code{length(col.gp) == length(levels(groups))}). By default: \code{col.gp = grDevices::rainbow(length(levels(groups)))} (i.e. colors defined automatically)
#' @param col.conv a \code{factor} object with the colours intended for \code{conv.hulls}. Must be the same length as the number of levels in \code{conv.hulls} (i.e. \code{length(col.conv) == length(levels(conv.hulls))}). By default: \code{col.conv = grDevices::rainbow(length(levels(conv.hulls)))} (i.e. colors defined automatically)
#' @param PCs a vector of length two with the Principal Components intended for the plot. By default: \code{PCs = c(1, 2)}
#' @param main main title of output plot. Should be presented between quotation marks. By default: \code{main = "Ordination of PCA coordinates"}
#' @param sp.as enables one to choose between ploting elements as \code{"points"} or \code{"text"}. If \code{sp.as = "text"}, also input a \code{factor} of characters to use as text (i.e. \code{sp.text}). By default: \code{sp.as = "points"}
#' @param sp.text only applies when \code{sp.as = "text"}. A \code{factor} including elements as texts intended in the plot. Has to be the same length as the number of rows in the coordinates of \code{PCA.out} (i.e. \code{length(sp.text) == nrow(PCA.out$x)}) . By default: \code{sp.text = NULL}
#' @param cross.origin A logical. If \code{TRUE}, the ordination plot will include lines crossing at the origin (i.e. x = 0, y = 0). By default: \code{cross.origin = TRUE}
#' @param lty only applies when \code{cross.origin = TRUE}. The type of lines crossing at the origin. Same as in \code{\link{par}}: "Line types can either be specified as an integer (0=blank, 1=solid (default), 2=dashed, 3=dotted, 4=dotdash, 5=longdash, 6=twodash) or as one of the character strings "blank", "solid", "dashed", "dotted", "dotdash", "longdash", or "twodash", where "blank" uses ‘invisible lines’ (i.e., does not draw them)." By default: \code{lty = "dotted"}
#' @param lwd only applies when \code{cross.origin = TRUE}. The width for lines crossing at the origin. Same as in \code{\link{par}}. By default: \code{lwd = 1}
#' @param leg a logical. If \code{TRUE}, a legend is added to plot. By default: \code{leg = TRUE}
#' @param leg.labels only applies when \code{leg = TRUE}. Must be the same length of levels in \code{groups} (i.e. \code{length(leg.labels) == length(levels(groups))}). By default: \code{leg.labels = groups} (i.e. labels will be the same as the levels from \code{groups})
#' @param leg.pos only applies when \code{leg = TRUE}. Specify legend location in the plot by using a keyword from the list: \code{"bottomright"}, \code{"bottom"}, \code{"bottomleft"}, \code{"left"}, \code{"topleft"}, \code{"top"}, \code{"topright"}, \code{"right"} and \code{"center"}. Alternatively, a single value can be provided and used for both margins, or two values, in which the first is used for X distance, and the second for Y distance. See also \code{\link{legend}} for details. By default: \code{leg.pos = "topright"}
#' @param cex.leg only applies when \code{leg = TRUE}. The magnification to be used in the legend. By default: \code{cex.leg = 1}
#' @param cex same as in \code{\link{par}}: "A numerical value giving the amount by which plotting text and symbols should be magnified relative to the default". By default: \code{cex = 1}
#' @param cex.axis same as in \code{\link{par}}. Relative to \code{cex}. The magnification to be used for axis annotation. By default: \code{cex.axix = 1}
#' @param cex.lab same as in \code{\link{par}}. Relative to \code{cex}. The magnification to be used for x and y labels. By default: \code{cex.lab = 1}
#' @param cex.main same as in \code{\link{par}}. Relative to \code{cex}. The magnification to be used for main title. By default: \code{cex.main = 1}
#'
#' @return
#' Require the output of \code{\link{prcomp}} and a vector with \code{groups} to plot. In addition, it is also possible to include convex hulls around each group (i.e. \code{conv.hulls}) and to control the colors intended for each group (i.e. \code{col.gp}) and for each convex hull (i.e. \code{col.conv}).
#'
#' @seealso
#' \code{\link{prcomp}}, \code{\link{palette}}, \code{\link{rgb}}, \code{\link{rainbow}}, \code{\link{legend}}
#'
#' Useful links:
#' \itemize{
#'   \item{\url{https://github.com/p-rocha/SoundShape}}
#'   \item{Report bugs at \url{https://github.com/p-rocha/SoundShape/issues}}}
#'
#' @references
#' Rocha, P. & Romano, P. (\emph{in prep}) The shape of sound: A new \code{R} package that crosses the bridge between Bioacoustics and Geometric Morphometrics.
#'
#' @examples
#' data(eig.sample)
#'
#' # PCA using 3D semilandmark coordinates
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
#' @export
#'
pca.plot <- function(PCA.out=NULL, groups=NULL, col.gp = grDevices::rainbow(length(levels(groups))), conv.hulls=NULL, col.conv=grDevices::rainbow(length(levels(conv.hulls))), PCs=c(1,2), main="Ordination of PCA coordinates", sp.as="points", sp.text=NULL, cross.origin=TRUE, lwd=1, lty="dotted", leg=TRUE, leg.labels=groups, leg.pos="topright", cex.leg=1, cex=1, cex.axis=1, cex.lab=1, cex.main=1){

  if(is.null(PCA.out))
  {stop('Please define PCA.out as an object containing the output of a Principal
        Components Analysis performed by prcomp function')}

  if(is.null(groups))
  {stop('Please define a vector to be used as groups')}

  # Store PC scores as new object
  PCscores <- summary(PCA.out)$importance

  # Use PC scores as labels in PCA visualization
  xlab = paste("Principal Component ", as.character(PCs[1]), " (",
               round(PCscores[2,PCs[1]]*100, 1), "%)", sep="")
  ylab = paste("Principal Component ", as.character(PCs[2])," (",
               round(PCscores[2,PCs[2]]*100, 1), "%)", sep="")

  # create a colour vector to colour the groups
  names(col.gp) = levels(groups)  # dimension labels
  col.groups = col.gp[match(groups, names(col.gp))]  # assign a colour to each specimen

  #Principal components species distribution
  graphics::plot(PCA.out$x[,PCs[1]], PCA.out$x[,PCs[2]], type='n', las=1, xlab=xlab, ylab=ylab, main=main,
       cex.axis=cex.axis, cex.lab=cex.lab, cex.main=cex.main)

  # Dotted line at origin
  if(cross.origin==TRUE){graphics::abline(h=0, v=0, lty=lty, lwd=lwd)} # end dotted line

  # Add convex hulls
  if(!is.null(conv.hulls)){ for(j in 1:nlevels(conv.hulls)) {
    # Get edge points (used to plot convex hull):
    edge_points = rownames(PCA.out$x[which(conv.hulls == levels(conv.hulls)[j]),])[
      grDevices::chull(PCA.out$x[which(conv.hulls == levels(conv.hulls)[j]), PCs])]
    # Plot convex hull as polygon:
    graphics::polygon(PCA.out$x[edge_points, PCs], col = grDevices::adjustcolor(col.conv[j], alpha.f = 0.1), border = col.conv[j])
  }}
  if(sp.as=="points"){graphics::points(PCA.out$x[,PCs[1]], PCA.out$x[,PCs[2]], pch=21, cex=cex, bg=col.groups)}
  if(sp.as=="text"){graphics::text(PCA.out$x[,PCs[1]], PCA.out$x[,PCs[2]], as.factor(sp.text), col=col.groups, cex=cex)}

  # Add legend
  if(leg==TRUE){graphics::legend(x=leg.pos, legend=unique(leg.labels), cex=cex.leg, pch=21, pt.bg=unique(col.groups))}
}
