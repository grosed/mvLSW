### Dependent libraries for functions in this package ###
#library(fields)
#library(wavethresh)

### Multivariate Locally Stationary Wavelet class ###

as.mvLSW <- function(
  x, 
  filter.number = 1, 
  family = "DaubExPhase",
  smooth.type = "all", 
  smooth.kernel = kernel("daniell", 0), 
  bias.correct = FALSE, 
  min.eig.val = -Inf,
  names = NULL
 ){

  spectrum <- x
  ##spectrum - value type and dimension
  if(!is(spectrum,"array")) stop("'spectrum' is not a numerical array.") 
  Dim <- dim(spectrum)
  if(length(Dim) != 4) stop("Dimension of 'spectrum' is not as expected.")
  if(Dim[1] == 1) stop("'spectrum' is univariate! Refer to the 'wavethresh' library.")
  if(Dim[1] != Dim[2]) stop("Dimension of 'spectrum' is not as expected.")
  if(Dim[3] != log2(Dim[4])) stop("Dimension of 'spectrum' is not as expected.")

  if(is.null(names)) names <- paste0("Series", 1:Dim[1])

  dimensions <- list(P = Dim[1], J = Dim[3], T = Dim[4])
  wavelet <- list(family = family, filter.number = filter.number)

  if(2*smooth.kernel$m + 1 >= Dim[4]) 
    stop("Kernel span is too wide for time series.")
  if(smooth.type == "all"){
    smooth.kernels <- smooth.kernel
    GCV <- NA
  }else{
    smooth.kernels <- as.list(1:Dim[3])
    for(j in 1:Dim[3]){
      smooth.kernels[[j]] <- smooth.kernel
    }
    GCV <- rep(NA, Dim[3])
  }
  
  if(is.na(min.eig.val)){
    MIN <- Inf
    for(j in 1:Dim[3]){
      for(k in 1:Dim[4]){
	    MIN <- min(MIN,eigen(spectrum[, , j, k], symmetric=TRUE)$value)
  	  }
    }
	min.eig.val <- MIN
  }

  smooth.info <- list(smooth.kernels = smooth.kernels, smooth.type = smooth.type, 
    GCV = GCV)
  correction <- list(bias.correct = bias.correct, min.eig.val = min.eig.val)
  object <- list(spectrum = spectrum,
    Information = list(names = names, dimensions = dimensions, wavelet = wavelet, 
	  smooth = smooth.info, correction = correction))
  class(object) <- append(class(object), "mvLSW")
  
  attr(object,"time") <- seq(from=0, by=1/Dim[4], len=Dim[4])

  invisible(object)
}

### Summary of a mvLSW object ###

summary.mvLSW <- function(
  object, 
  ...){ 

  if(!is.mvLSW(object)) stop("Invalid 'object' argument.")
  
  cat("== Dimensions ==\n")
  cat("P :", object$Information$dimension$P, "\n")
  cat("J :", object$Information$dimension$J, "\n")
  cat("T :", object$Information$dimension$T, "\n")
  cat("\n")
  cat("== Wavelet Transform ==\n")
  cat("Family        :", object$Information$wavelet$family, "\n")
  cat("Filter Number :", object$Information$wavelet$filter.number, "\n")
  cat("\n")
  cat("== Smoothing ==\n")
  if(object$Information$smooth$smooth.type == "all"){
    cat("Type   : all\n")
    cat("Method :", attributes(object$Information$smooth$smooth.kernels)$name)
    cat(" - GCV criterion =", object$Information$smooth$GCV, "\n")
  }else{
    cat("Type         : by.level\n")
    for(j in 1:object$Information$dimension$J){
	  kernel.name <- attributes(object$Information$smooth$smooth.kernels[[j]])$name
      cat("Method-Lev", j, ":", kernel.name)
      cat(" - GCV criterion =", object$Information$smooth$GCV[j], "\n")
    }
  }
  cat("\n")
  cat("== Sundries ==\n")
  cat("Applied Bias Correction :", object$Information$correction$bias.correct, "\n")
  cat("Minimum Eigenvalue      :", object$Information$correction$min.eig.val, "\n")
  cat("\n")
  invisible(NULL)
}

### Logical check to assess whether object is a mvLSW object and structured correctly ###

is.mvLSW <- function(
  object){

  ##General structure
  if(!("mvLSW" %in% class(object))) return(FALSE)
  if(!is.list(object)) return(FALSE)
  if(length(names(object)) != 2) return(FALSE)
  if(!("spectrum" %in% names(object))) return(FALSE)
  if(!("Information" %in% names(object))) return(FALSE)

  #Structure of Information
  if(!is.list(object$Information)) return(FALSE)
  if(length(names(object$Information)) != 5) return(FALSE)
  if(!("names" %in% names(object$Information))) return(FALSE)
  if(!("dimensions" %in% names(object$Information))) return(FALSE)
  if(!("wavelet" %in% names(object$Information))) return(FALSE)
  if(!("smooth" %in% names(object$Information))) return(FALSE)
  if(!("correction" %in% names(object$Information))) return(FALSE)

  #Structure of dimension
  if(!is.list(object$Information$dimensions)) return(FALSE)
  if(length(names(object$Information$dimensions)) != 3) return(FALSE)
  if(!("P" %in% names(object$Information$dimensions))) return(FALSE)
  if(!("J" %in% names(object$Information$dimensions))) return(FALSE)
  if(!("T" %in% names(object$Information$dimensions))) return(FALSE)
  P <- object$Information$dimensions$P
  T <- object$Information$dimensions$T
  J <- object$Information$dimensions$J
  if(!is.numeric(P) || length(P) != 1) return(FALSE)
  if(P%%1 != 0 || P <= 0) return(FALSE)
  if(!is.numeric(T) || length(T) != 1) return(FALSE)
  if(T%%1 != 0 || T <= 0) return(FALSE)
  if(!is.numeric(J) || length(J) != 1) return(FALSE)
  if(J%%1 != 0 || J <= 0) return(FALSE)
  if(J != log2(T)) return(FALSE)

  #Structure of names
  if(length(object$Information$names) != P) return(FALSE)
  if(!is.character(object$Information$names)) return(FALSE)

  #Structure of wavelet
  if(!is.list(object$Information$wavelet)) return(FALSE)
  if(length(names(object$Information$wavelet)) != 2) return(FALSE)
  if(!("family" %in% names(object$Information$wavelet))) return(FALSE)
  if(!("filter.number" %in% names(object$Information$wavelet))) return(FALSE)
  if(!is.character(object$Information$wavelet$family)) return(FALSE)
  if(length(object$Information$wavelet$family) != 1) return(FALSE)
  if(!(object$Information$wavelet$family %in% c("DaubExPhase", "DaubLeAsymm"))) return(FALSE)
  if(!is.numeric(object$Information$wavelet$filter.number)) return(FALSE)
  if(length(object$Information$wavelet$filter.number) != 1) return(FALSE)
  if(object$Information$wavelet$filter.number%%1 != 0) return(FALSE)
  if(object$Information$wavelet$filter.number < 1) return(FALSE)
  if(object$Information$wavelet$filter.number > 10) return(FALSE)

  #Structure of smooth
  if(!is.list(object$Information$smooth)) return(FALSE)
  if(length(names(object$Information$smooth)) != 3) return(FALSE)
  if(!("smooth.kernels" %in% names(object$Information$smooth))) return(FALSE)
  if(!("smooth.type" %in% names(object$Information$smooth))) return(FALSE)
  if(!("GCV" %in% names(object$Information$smooth))) return(FALSE)

  if(length(object$Information$smooth$smooth.type) != 1) return(FALSE)
  if(!is.character(object$Information$smooth$smooth.type)) return(FALSE)
  if(object$Information$smooth$smooth.type == "all"){
    if(!is.tskernel(object$Information$smooth$smooth.kernels)) return(FALSE)
    if(length(object$Information$smooth$GCV) != 1) return(FALSE)
    if(!is.na(object$Information$smooth$GCV)){
      if(!is.numeric(object$Information$smooth$GCV)) return(FALSE)
    }
  }else if(object$Information$smooth$smooth.type == "by.level"){
    if(!is.list(object$Information$smooth$smooth.kernels)) return(FALSE)
    if(length(object$Information$smooth$smooth.kernels) != J) return(FALSE)
    if(length(object$Information$smooth$GCV) != J) return(FALSE)
    if(!is.numeric(object$Information$smooth$GCV) && 
	  !all(is.na(object$Information$smooth$GCV)) ) return(FALSE)
    for(j in 1:J){
      if(!is.tskernel(object$Information$smooth$smooth.kernels[[j]])) return(FALSE)    
      if(!is.na(object$Information$smooth$GCV[j])){
        if(!is.numeric(object$Information$smooth$GCV[j])) return(FALSE)
      }
    }
  }else{
    return(FALSE)
  }

  #Structure of correction
  if(!is.list(object$Information$correction)) return(FALSE)
  if(length(names(object$Information$correction)) != 2) return(FALSE)
  if(!("bias.correct" %in% names(object$Information$correction))) return(FALSE)
  if(!("min.eig.val" %in% names(object$Information$correction))) return(FALSE)
  if(!is.logical(object$Information$correction$bias.correct)) return(FALSE)
  if(length(object$Information$correction$min.eig.val) != 1) return(FALSE)
  if(!is.numeric(object$Information$correction$min.eig.val)) return(FALSE)
  if(object$Information$correction$min.eig.val < 0 && 
    object$Information$correction$min.eig.val != -Inf) return(FALSE)

  ##spectrum
  if(!is.array(object$spectrum)) return(FALSE)
  if(anyNA(object$spectrum)) return (FALSE)
  if(!all(is.numeric(object$spectrum))) return(FALSE)
#  if(!all(is.na(object$spectrum) || is.numeric(object$spectrum))) return(FALSE)
  Dim <- dim(object$spectrum)
  if(length(Dim) != 4) return(FALSE)
  if(Dim[1] != P) return(FALSE)
  if(Dim[2] != P) return(FALSE)
  if(Dim[3] != J) return(FALSE)
  if(Dim[4] != T) return(FALSE)

  return(TRUE)
}

### Plot - for defined channel pair and level ###

single_pqj_plot <- function(
  object, 
  Int.lower = NULL, 
  Int.upper = NULL, 
  p = 1, 
  q = 1, 
  j = 1, 
  ...){

  #check object is a mvLSW object
  if(!is.mvLSW(object)) stop("Invalid 'object' argument.")
  if(is.null(attr(object,"time"))){
    Rescale <- seq(from=0, len=object$Information$dimensions$T, 
      by=1/object$Information$dimensions$T)
  }else{
    Rescale <- attr(object,"time")
	if(!any(class(Rescale)=="POSIXt")){
	  Rescale <- as.vector(Rescale)
	}
  }
  
  #check p, q & j arguments
  if(!is.numeric(p) || length(p) != 1) stop("Invalid 'p' argument.")
  if(p%%1 != 0 || p<1 || p > object$Information$dimensions$P) stop("Invalid 'p' argument.")
  if(!is.numeric(q) || length(q) != 1) stop("Invalid 'q' argument.")
  if(q%%1 != 0 || q < 1 || q > object$Information$dimensions$P) stop("Invalid 'q' argument.")
  if(!is.numeric(j) || length(j) != 1) stop("Invalid 'j' argument.")
  if(j%%1 != 0 || j < 1 || j > object$Information$dimensions$J) stop("Invalid 'j' argument.")

  #check Int.lower and Int.upper arguments
  if(!is.null(Int.lower) && !is.null(Int.upper)){
    if(!is.mvLSW(Int.lower)) stop("Invalid 'Int.lower' argument.")
    if(!is.mvLSW(Int.upper)) stop("Invalid 'Int.upper' argument.")
    Interval <- TRUE
  }else if(is.null(Int.lower) != is.null(Int.upper)){
    stop("Invalid 'Int.lower' and/or 'Int.upper' arguments.")
  }else{
    Interval <- FALSE
  }

  #check col argument supplied via ellipses
  ellipses <- list(...)
  if("col" %in% names(ellipses)){ col <- ellipses$col } else { col <- "black" }
  if(is.numeric(col)){
    if(col[1]%%1 == 0 && col[1] >= 1 && length(col) == 1){
      col <- rgb(t(col2rgb(palette())),maxColorValue=255)[(col-1) %% length(palette()) + 1]
    }
  }
  if(!is.character(col) || length(col) != 1) stop("Invalid 'col' argument.")
  if(substr(col,1,1) != "#") col <- rgb(t(col2rgb(col)), maxColorValue = 255)
  if(substr(col,1,1) != "#" || nchar(col) != 7) stop("Invalid 'col' argument.")
  if(!all(grepl("[0-9A-F]", unlist(strsplit(col, ""))[-1]))) stop("Invalid 'col' argument.")
  ellipses$col <- col

  #check ylim argument supplied via ellipses
  if(!("ylim" %in% names(ellipses))){
    ylim <- range(object$spectrum[p, q, j, ])
    if(Interval){
	  ylim <- range(ylim, Int.lower$spectrum[p, q, j, ], Int.upper$spectrum[p, q, j, ])
	}
    ellipses$ylim <- ylim
  }
  
  #check ylab, xlab and xlim arguments supplied via ellipses
  if(!("ylab" %in% names(ellipses))) ellipses$ylab <- "Spectrum"
  if(!("xlab" %in% names(ellipses))) ellipses$xlab <- "Time"#"Rescaled Time"
  if(!("xlim" %in% names(ellipses))){
    if(is.null(attr(object,"time"))){
	  ellipses$xlim <- c(0,1)
	}else{
	  ellipses$xlim <- range(Rescale)
	}
  } 

  #Create plot
  do.call(plot, c(list(x=Rescale,y=rep(ellipses$ylim[1],length(Rescale)),
    type="n"),ellipses))
  grid(col="grey60")
  if(Interval){
    #Add intervals if available
    polygon(x = c(Rescale, Rescale[object$Information$dimensions$T:1]),
	  y = c(Int.lower$spectrum[p, q, j, ], Int.upper$spectrum[p, q, j, object$Information$dimensions$T:1]),
      col = paste0(ellipses$col, "40"), border = NA)
  }
  do.call(lines, c(list(x=Rescale,y=object$spectrum[p, q, j, ]), ellipses))
}

### Plot - panel for all channel pairs at defined level ###

panel_j_plot <- function(
  object, 
  Int.lower = NULL, 
  Int.upper = NULL, 
  j = 1, 
  diag = TRUE,
  ...){

  #check object is a mvLSW object
  if(!is.mvLSW(object)) stop("Invalid 'object' argument.")

  #check j arguments
  if(!is.numeric(j) || length(j) != 1) stop("Invalid 'j' argument.")
  if(j%%1 != 0 || j < 1 || j > object$Information$dimensions$J) stop("Invalid 'j' argument.")

  #check interval arguments
  if(!is.null(Int.lower) && !is.null(Int.upper)){
    if(!is.mvLSW(Int.lower)) stop("Invalid 'Int.lower' argument.")
    if(!is.mvLSW(Int.upper)) stop("Invalid 'Int.upper' argument.")
    Interval <- TRUE
  }else if(is.null(Int.lower) != is.null(Int.upper)){
    stop("Invalid 'Int.lower' and/or 'Int.upper' arguments.")
  }else{
    Interval <- FALSE
  }
  
  #check diag arguments
  if(!is.logical(diag)) stop("Invalid 'diag' argument.")

  #Determine layout of plot window
  LAYOUT <- matrix(0, ncol = object$Information$dimensions$P, nrow = object$Information$dimensions$P)
  k <- 1
  for(p in 1:object$Information$dimensions$P){  for(q in 1:p){  LAYOUT[p, q] <- k; k <- k + 1  }  }

  par.format <- par(no.readonly = TRUE)
  LAYOUT[1, object$Information$dimensions$P] <- k
  layout(LAYOUT)
  par(mar = c(1, 1, 4, 1) / 2, oma = 2 * c(2, 1.2, 1, 1), 
    mgp = c(1.5, 0.6, 0), no.readonly = TRUE)

  #check ylim argument supplied via ellipses
  ellipses <- list(...)
  if(!("ylim" %in% names(ellipses))){
    ylim <- range(object$spectrum[, , j, ])
    if(Interval){
	  ylim <- range(ylim, Int.lower$spectrum[, , j, ], Int.upper$spectrum[, , j, ])
	}
    ellipses$ylim <- ylim
  }
  if("ylab" %in% names(ellipses)){ YLAB <- ellipses$ylab } else{ YLAB <- "Spectrum" }
  if("xlab" %in% names(ellipses)){ XLAB <- ellipses$xlab } else{ XLAB <- "Time" }
  ellipses$xlab <- ""
  ellipses$ylab <- ""
  ellipses$yaxt <- "s"
  ellipses$xaxt <- "n"
  ellipses$main <- ""
  
  #Create plots
  P <- object$Information$dimensions$P
  for(p in 1:P){
    if(p == P){ ellipses$xaxt <- "s" } else { ellipses$xaxt <- "n" }
    for(q in 1:p){
	  if(q == 1){ ellipses$yaxt <- "s" } else { ellipses$yaxt <- "n" }
      if(!diag && p==q){
        plot(NA, NA, type="n", xlim = c(-1, 1), xaxt = "n", xlab = "", 
		  ylim = c(-1, 1), ylab = "", yaxt = "n", main = "", frame=TRUE)
        text(0, 0, object$Information$names[p], font=2, cex=2)
      }else{
        do.call(single_pqj_plot, c(list(object = object, Int.lower = Int.lower, 
	      Int.upper = Int.upper, p = p, q = q, j = j), ellipses))
		if(q==1) mtext(YLAB, side=2, line=1.7)
        if(p == P) mtext(XLAB, side = 1, line=1.7)
        if(p == q) title(object$Information$names[p], font=2, cex.main=1.5)
      }
    }
  }
  
  #Add Level label
  plot(NA, NA, type = "n", xlim = c(-1, 1), xaxt = "n", xlab = "", ylim=c(-1, 1), 
    ylab = "", yaxt = "n", main = "", frame = FALSE)
  text(0, 0, paste0("Level\n", j), font = 2, cex = 2)
  par(par.format)
}

### Plot - panel for defined channel pair at all levels ###

panel_pq_plot <- function(
  object, 
  Int.lower = NULL, 
  Int.upper = NULL, 
  p = 1, 
  q = 1,
  ...){

  #check object is a mvLSW object
  if(!is.mvLSW(object)) stop("Invalid 'object' argument.")
  
  #check p & q arguments
  if(!is.numeric(p) || length(p) != 1) stop("Invalid 'p' argument.")
  if(p%%1 != 0 || p<1 || p > object$Information$dimensions$P) stop("Invalid 'p' argument.")
  if(!is.numeric(q) || length(q) != 1) stop("Invalid 'q' argument.")
  if(q%%1 != 0 || q < 1 || q > object$Information$dimensions$P) stop("Invalid 'q' argument.")

  #check Int.lower and Int.upper arguments
  if(!is.null(Int.lower) && !is.null(Int.upper)){
    if(!is.mvLSW(Int.lower)) stop("Invalid 'Int.lower' argument.")
    if(!is.mvLSW(Int.upper)) stop("Invalid 'Int.upper' argument.")
    Interval <- TRUE
  }else if(is.null(Int.lower) != is.null(Int.upper)){
    stop("Invalid 'Int.lower' and/or 'Int.upper' arguments.")
  }else{
    Interval <- FALSE
  }

  #Determine panel size
  J <- object$Information$dimensions$J
  MFROW <- c(floor(sqrt(J)), ceiling(sqrt(J)))
  if(prod(MFROW)<J) MFROW[1] <- MFROW[1]+1
  par.format <- par(no.readonly = TRUE)
  par(mfrow = MFROW)

  #check ylim argument supplied via ellipses
  ellipses <- list(...)
  if(!("ylim" %in% names(ellipses))){
    ylim <- range(object$spectrum[p, q, , ])
    if(Interval){
	  ylim <- range(ylim, Int.lower$spectrum[p, q, , ], Int.upper$spectrum[p, q, , ])
	}
    ellipses$ylim <- ylim
  }
  
  #Create panel
  for(j in 1:object$Information$dimensions$J){
    ellipses$main <- paste("Level",j)
    do.call(single_pqj_plot, c(list(object = object, Int.lower = Int.lower, 
	  Int.upper = Int.upper, p = p, q = q, j = j), ellipses))
  }

  #Fill remaining space with blanks
  j <- object$Information$dimensions$J+1
  while(j <= prod(MFROW)){
    plot(NA, NA, type = "n", xlim = c(-1, 1), xaxt = "n", xlab = "", ylim=c(-1, 1),
	  ylab = "", yaxt="n", main="", frame=FALSE)
	j <- j + 1
  }
  par(par.format)
}

### Plot - as panel_pq_plot but presented in a matrix image plot ###

matrix_j_plot <- function(
  object, 
  p = 1, 
  q = 1,
  sub = "Spectrum",
  ...){

  #check object is a mvLSW object
  if(!is.mvLSW(object)) stop("Invalid 'object' argument.")

  #check p & q arguments
  if(!is.numeric(p) || length(p) != 1) stop("Invalid 'p' argument.")
  if(p%%1 != 0 || p<1 || p > object$Information$dimensions$P) stop("Invalid 'p' argument.")
  if(!is.numeric(q) || length(q) != 1) stop("Invalid 'q' argument.")
  if(q%%1 != 0 || q < 1 || q > object$Information$dimensions$P) stop("Invalid 'q' argument.")

  #check sub argument
  if(length(sub)!=1 || any(!is.character(sub))) stop("Invalid 'sub' argument.")
  
  #check 'zlim' argument
  ellipses <- list(...)
  if(!("zlim" %in% names(ellipses))) ellipses$zlim <- range(object$spectrum[p, q, , ])

  #get subtitle
  LAB <- sub
  ellipses$sub <- ""
  
  #Create plot
  Matrix <- t(object$spectrum[p, q, , ])
  T <- object$Information$dimension$T
  if(is.null(attr(object,"time"))){
    Rescale <- seq(from = 0.5/T, len = T, by = 1/T)
  }else{
    Rescale <- as.vector(attr(object,"time"))
  }
  J <- object$Information$dimension$J

  par.format <- par(no.readonly = TRUE)
  par(mfrow=c(1, 1))
  if(!("ylim" %in% names(ellipses))) ellipses$ylim <- c(0, J) + 0.5
  if(!("xlim" %in% names(ellipses))){
    if(is.null(attr(object,"time"))){
      ellipses$xlim <- c(0, 1)
	}else{
	  ellipses$xlim <- range(attr(object,"time"))
	}
  }
  if(!("ylab" %in% names(ellipses))) ellipses$ylab <- "Level"
  if(!("xlab" %in% names(ellipses))) ellipses$xlab <- "Time"#"Rescaled Time"
  if("main" %in% names(ellipses)){ MAIN <- ellipses$main }else{ MAIN <- "" }
  ellipses$main <- ""
  do.call(image.plot, c(list(x = Rescale, y = 1:J, z = Matrix), ellipses))

  #Add title
  if(MAIN == ""){
    text <- object$Information$names[p]
    if(p != q) text <- paste0(text, " v ", object$Information$names[q])
  }else{
    text <- MAIN
  }
  text <- paste0(text, "\n", LAB)
  title(text, font = 2, cex = 1.5)
  par(par.format)
}

### Plot - Master function ###

plot.mvLSW <- function(
  x, ## Same as 'object' - set as 'x' to be compatible with 'plot' command
  style = 1, 
  info = NULL, 
  Interval = NULL,
  diag = TRUE,
  sub = "Spectrum",
  ...){
  
  if(!is.numeric(style) || length(style) != 1) stop("Invalid 'style' argument.")

  if(is.null(Interval)){
    Int.lower = NULL
    Int.upper = NULL
  }else{
    if(!is.list(Interval)) stop("Invalid 'Interval' argument.")
	if(length(Interval) != 2) stop("Invalid 'Interval' argument.")
	if(!all(c("L","U") %in% names(Interval))) stop("Invalid 'Interval' argument.")
	Int.lower <- Interval$L
	Int.upper <- Interval$U
  }

  if(style == 1){
    #single plot
    if(is.null(info)) info <- c(1, 1, 1)
    if(!is.numeric(info) || length(info) != 3) 
	  stop("Invalid 'info' argument for plotting style.")
    single_pqj_plot(object = x, Int.lower = Int.lower, Int.upper = Int.upper, 
	  p = info[1], q = info[2], j = info[3], ...)
  }else if(style == 2){
    #level panel
    if(is.null(info)) info <- c(1)
    if(!is.numeric(info) || length(info) != 1) 
	  stop("Invalid 'info' argument for plotting style.")
    panel_j_plot(object = x, Int.lower = Int.lower, Int.upper = Int.upper, 
      j = info, diag = diag, ...)
  }else if(style == 3){
    #channel pair panel
    if(is.null(info)) info <- c(1, 1)
    if(!is.numeric(info) || length(info) != 2) 
	  stop("Invalid 'info' argument for plotting style.")
    panel_pq_plot(object = x, Int.lower = Int.lower, Int.upper = Int.upper,  
	  p = info[1], q = info[2], ...)
  }else if(style == 4){
    if(is.null(info)) info <- c(1, 1)
    if(!is.numeric(info) || length(info) != 2) 
	  stop("Invalid 'info' argument for plotting style.")
    if(length(col) == 1) col <- tim.colors(64)
    matrix_j_plot(object = x, p = info[1], q = info[2], sub = sub,...)
  }else{
    stop("Invalid 'style' argument.")
  }
}

### Calculate the Evolutionary Wavelet Spectrum ###

mvEWS <- function(
  X,
  filter.number = 1, 
  family = "DaubExPhase",
  smooth = TRUE,
  type = "all", 
  kernel.name = "daniell", 
  kernel.param = floor(sqrt(nrow(X))), 
  optimize = FALSE,
  smooth.Jset = NA, 
  bias.correct = TRUE, 
  tol = 1e-10, 
  verbose = FALSE){
  
  ##Check arguments
  if(!is.matrix(X)) stop("Invalid 'X' argument")
  if(is.ts(X)){
    TIME <- time(X)
  }else if(is.zoo(X) || is.xts(X)){
    TIME <- time(X)
	X <- as.ts(X)
  }else if(is.matrix(X)){
    TIME <- 1:nrow(X)
	X <- as.ts(X)
  }else{
    stop("Invalid 'X' argument")
  }
  if(any(is.na(X)) || !is.numeric(X)) stop("Invalid 'X' argument.")
  if(ncol(X) == 1) stop("'X' is univariate! Refer to the 'wavethresh' library.")
  if(log2(nrow(X))%%1 != 0) stop("Invalid length of 'X'.")
  family <- match.arg(family, c("DaubExPhase", "DaubLeAsymm"))
  if(!is.numeric(filter.number) || length(filter.number) != 1) 
    stop("Invalid 'filter.number' argument.")
  if(filter.number%%1 != 0 || filter.number < 1 || filter.number > 10) 
    stop("Invalid 'filter.number' argument.")
  if(!is.logical(smooth) || length(smooth) != 1) stop("Invalid 'smooth' argument.")
  type <- match.arg(type, c("all", "by.level"))
  kernel.name <- match.arg(kernel.name, c("daniell", "modified.daniell", "dirichlet", "fejer"))
  if(!is.numeric(kernel.param)) stop("Invalid 'kernel.param' argument.")
  if(any(kernel.param%%1 != 0) || any(kernel.param < 1)) 
    stop("Invalid 'kernel.param' argument.")
  if(type == "all"){
    if(!is.vector(kernel.param)) 
	  stop("Invalid 'kernel.param' for type 'all'. Argument must be a vector.")
    if(kernel.name %in% c("dirichlet", "fejer") && length(kernel.param) != 2) 
	  stop("Invalid 'kernel.param' argument for kernel name.")
    if(!any(is.na(smooth.Jset))){
      if(!is.numeric(smooth.Jset) || min(smooth.Jset) < 1 || max(smooth.Jset) > log2(nrow(X)) 
        || length(smooth.Jset) > log2(nrow(X)) || any(smooth.Jset%%1 != 0)){
          stop("Invalid 'smooth.Jset' argument.")
      }else{
        smooth.Jset <- sort(unique(smooth.Jset))
      }
    }else{
      smooth.Jset <- 1:log2(nrow(X))
    }
  }else{
    if(!is.matrix(kernel.param)) 
	  stop("Invalid 'kernel.param' for type 'by.level'. Argument must be a matrix with J columns.")
    if(2^ncol(kernel.param) != nrow(X))
	  stop("Invalid 'kernel.param' for type 'by.level'. Argument must be a matrix with J columns.")
    if(kernel.name %in% c("dirichlet", "fejer") && nrow(kernel.param) != 2) 
	  stop("Invalid 'kernel.param' for 'kernel.name'. Argument must be a matrix with 2 rows.")
  }
  if(!is.logical(optimize)) stop("Invalid 'optimize' argument.")
  if(!is.logical(bias.correct) || length(bias.correct) != 1) 
    stop("Invalid 'bias.correct' argument.")
  if(length(tol) != 1)  stop("Invalid 'tol' argument.")
  if(is.na(tol)) tol <- -Inf
  if(!is.numeric(tol)) stop("Invalid 'tol' argument.")
  if(tol <= 0 && tol != -Inf) stop("Invalid 'tol' argument.")
  if(!is.logical(verbose) || length(verbose) != 1) stop("Invalid 'verbose' argument.")

  #Calculate the raw wavelet periodogram
  if(verbose) cat("Calculating the raw wavelet periodogram.\n")
  RawPeriod <- RawPeriodogram(X, filter.number = filter.number, family = family, FALSE)
  
  if(smooth){
    if(verbose) cat("Smoothing periodogram - type:", type, "\n")
    #Smooth periodogram
    if(type == "all"){ #Smooth all levels collectively
      #Make mvLSW object
      EWS <- as.mvLSW(x = RawPeriod, filter.number = filter.number, 
        family = family, smooth.type = "all", smooth.kernel = kernel("daniell", 0), 
        bias.correct = FALSE, min.eig.val = -Inf, names = colnames(X))
      EWS <- Smooth_EWS(EWS, kernel.name, kernel.param, optimize, 
        type, level = smooth.Jset)
    }else{ 
      #Smooth on a by-level basis
      #Make mvLSW object
      EWS <- as.mvLSW(x = RawPeriod, filter.number = filter.number, family = family,
        smooth.type = "by.level", smooth.kernel = kernel("daniell", 0),
        bias.correct = FALSE, min.eig.val = -Inf, names = colnames(X))
      J <- EWS$Information$dimension$J
      if(verbose) cat("- Evaluating level:")
      for(j in 1:J){
        if(verbose) cat(j)
        #Smooth level
        EWS <- Smooth_EWS(EWS, kernel.name, kernel.param = kernel.param[,j], optimize, 
          type, level = j)
      }
      if(verbose) cat("\n")
    }
  }else{ #Do not apply smoothing
    #Make mvLSW object
    EWS <- as.mvLSW(x = RawPeriod, filter.number = filter.number, family = family,
      smooth.type = "all", smooth.kernel = kernel("daniell", 0), 
      bias.correct = FALSE, names = colnames(X))
  }

  #Bias correction
  if(verbose && bias.correct) cat("Bias Correction.\n")
  if(bias.correct) EWS <- CorrectBias(EWS)

  #Adjust to positive-definite matrix
  if(verbose && tol>0) cat("Adjustment for positive definiteness.\n")
  if(tol>0) EWS <- AdjPositiveDef(EWS, tol)
  attr(EWS,"time") <- TIME
  invisible(EWS)
}
 
### Calculate the Raw Wavelet Periodogram ###

RawPeriodogram <- function(
  X,
  filter.number = 1,
  family = "DaubExPhase",
  format = TRUE){

  if(!is.matrix(X)) stop("Invalid 'X' argument")
  if(is.ts(X)){
    TIME <- time(X)
  }else if(is.zoo(X) || is.xts(X)){
    TIME <- time(X)
	X <- as.ts(X)
  }else if(is.matrix(X)){
    TIME <- 1:nrow(X)
	X <- as.ts(X)
  }else{
    stop("Invalid 'X' argument")
  }
  if(any(is.na(X)) || !is.numeric(X)) stop("Invalid 'X' argument.")
  P <- ncol(X)
  if(P == 1) stop("'X' is univariate! Refer to the 'wavethresh' library.")
  T <- nrow(X)
  J <- log2(T)
  if(J%%1 != 0) stop("Invalid length of 'X'")
  family <- match.arg(family, c("DaubExPhase", "DaubLeAsymm"))
  if(!is.numeric(filter.number) || length(filter.number) != 1) 
    stop("Invalid 'filter.number' argument.")
  if(filter.number%%1 != 0 || filter.number < 1 || filter.number > 10) 
    stop("Invalid 'filter.number' argument.")
  if(!is.logical(format)) stop("Invalid 'format' argument.")
  
  Periodogram <- array(NA, dim=c(P, P, J, T))
  #Wavelet co-efficients
  for(p in 1:P){
    #level 1 = fine // level J = coarse
    Periodogram[p, 1, , ] <- matrix(wd(as.vector(X[, p]), family = family, 
      filter.number = filter.number, type = "station", bc = "periodic")$D,
      nrow = J,ncol = T,byrow = TRUE)
  }
  #Raw periodogram: I=dd'
  min <- Inf
  for(j in 1:J){
    for(k in 1:T){
      d <- as.numeric(Periodogram[, 1, j, k])
      Periodogram[, , j, k] <- d %*% t.default(d)
	  min <- min(min,eigen(Periodogram[,,j,k], symmetric=TRUE)$val)
    }
  }
  if(format){
    RawPeriod <- as.mvLSW(x = Periodogram, filter.number = filter.number, family = family,
      smooth.type = "all", smooth.kernel = kernel("daniell", 0), 
      bias.correct = FALSE, min.eig.val = -Inf, names = colnames(X))
	  attr(RawPeriod,"time") <- TIME
    invisible(RawPeriod)
  }else{
    invisible(Periodogram)
  }
}

### Bias Correction to the Evolutionary Wavelet Spectrum Estimate ###

CorrectBias <- function(
  object){
  
  if(!is.mvLSW(object)) stop("Invalid 'object' argument.")
  if(!object$Information$correction$bias.correct){
    A <- ipndacw(-object$Information$dimensions$J, 
      family = object$Information$wavelet$family,
      filter.number = object$Information$wavelet$filter.number)
    Ainv <- solve(A)
    for(p in 1:object$Information$dimensions$P){
      for(q in 1:p){
        for(k in 1:object$Information$dimensions$T){
          object$spectrum[p, q, , k] <- object$spectrum[q, p, , k] <- Ainv %*% object$spectrum[p, q, , k]
        }
      }
    }
    object$Information$correction$bias.correct <- TRUE
    invisible(object)
  }else{
    warning("Bias correction has already been applied.")
    invisible(object)
  }
}
Smooth_EWS <- function(
  object, 
  kernel.name = "daniell", 
  kernel.param = 1, 
  optimize = FALSE, 
  type = "all", 
  level = NA){

  if(!is.mvLSW(object)) stop("Invalid 'object' argument.")
  kernel.name <- match.arg(kernel.name, c("daniell", "modified.daniell", "dirichlet", "fejer"))
  if(!is.numeric(kernel.param)) stop("Invalid 'kernel.param' argument.")
  if(any(kernel.param < 1) || any(kernel.param%%1 != 0)) stop("Invalid 'kernel.param' argument.")
  if(!is.logical(optimize)) stop("Invalid 'optimize' argument.")
  type <- match.arg(type, c("all", "by.level"))
  if(type == "by.level"){
    if(!is.numeric(level) || length(level) != 1)
      stop("Invalid 'level' argument for type 'by.level'.")
  }else{
    if(all(is.na(level))){
      level <- 1:object$Information$dimensions$J
    }else if(is.numeric(level)){
      level <- unique(level)
      if(length(level) > object$Information$dimensions$J || min(level) < 1 || 
        max(level) > object$Information$dimensions$J || any(level%%1 != 0))
          stop("Invalid 'level' argument.")
    }else{
      stop("Invalid 'level' argument.")
    }
  }
  
  if(optimize){
    options <- matrix(NA, ncol = length(kernel.param), nrow = prod(kernel.param))
    if(kernel.name %in% c("daniell", "modified.daniell") && length(kernel.param) > 1){
      kernel.param <- sort(kernel.param)
      check <- TRUE
    }else{
      check <- FALSE
    }
    val <- rep(1, length(kernel.param))
    l <- 1
    while(all(val <= kernel.param) && l <= nrow(options)){
      if(!check){
        options[l, ] <- val
        l <- l + 1
      }else if(all(diff(val) >= 0)){
        options[l, ] <- val
        l <- l + 1
      }
      val[1] <- val[1] + 1
      d <- 1
      while(d <= (length(kernel.param) - 1)){
        if(val[d] > kernel.param[d]){
          val[d + 1] <- val[d + 1] + 1
          val[d] <- 1
        }
        d <- d + 1
      }
    }
    if(l == 2){
      options <- matrix(options[1, ], nrow = 1)
    }else{
      options <- as.matrix(options[1:(l - 1), ])
    }
  }else{
    options <- matrix(kernel.param, nrow = 1)
  }

  GCV <- rep(NA, nrow(options))
  for(index in 1:nrow(options)){
    GCV[index] <- Smooth_GCV(object, kernel.name, kernel.param = as.numeric(options[index, ]), 
	  type, level, GCVonly = TRUE)
  }

  best.ind <- which(GCV == min(GCV))
  if(length(best.ind) > 1){
    best.ind <- best.ind[1]
	warning("Multiple kernels found as optimal. First occurrence is taken.")
  }
  best.kernel.param <- as.numeric(options[best.ind, ])
  object <- Smooth_GCV(object, kernel.name, kernel.param = best.kernel.param, type, 
    level, GCVonly = FALSE)
  invisible(object)
}

### Smoothining of the Evolutionary Wavelet Spectrum Estimate (calls C) ###

Smooth_GCV <- function(
  object, 
  kernel.name = "daniell", 
  kernel.param = 1, 
  type = "all", 
  level = NA, 
  GCVonly = FALSE){
  
  if(!is.mvLSW(object)) stop("Invalid 'object' argument.")
  type <- match.arg(type, c("all", "by.level"))
  if(type == "by.level"){
    if(!is.numeric(level) || length(level) != 1)
      stop("Invalid 'level' argument for type 'by.level'.")
  }else{
    if(all(is.na(level))){
      level <- 1:object$Information$dimensions$J
    }else if(is.numeric(level)){
      level <- unique(level)
      if(length(level) > object$Information$dimensions$J || min(level) < 1 || 
        max(level) > object$Information$dimensions$J || any(level%%1 != 0)){
        stop("Invalid 'level' argument.")
      }
    }else{
      stop("Invalid 'level' argument.")
    }
  }
  if(!is.logical(GCVonly)) stop("Invalid 'GCVonly' argument.")
  if(!is.numeric(kernel.param)) stop("Invalid 'kernel.param' argument.")
  if(any(kernel.param%%1 != 0) || any(kernel.param < 1)) 
    stop("Invalid 'kernel.param' argument.")
  kernel.name <- match.arg(kernel.name, c("daniell", "modified.daniell", "dirichlet", "fejer"))

  if(kernel.name %in% c("daniell", "modified.daniell")){
    kernel <- kernel(kernel.name, kernel.param)
  }else{
    if(length(kernel.param) != 2) stop("'kernel.param' is incompatible for defined kernel.")
    kernel <- kernel(kernel.name, kernel.param[1], kernel.param[2])
  }

  if(kernel$coef[1] == 1 && sum(2 * kernel$coef[-1]) < .Machine$double.eps){
    warning(paste0("Smoothing kernel ", attributes(kernel)$name, " is equivalent to the identity function."))
    #i.e. smoothed spectrum equates to the raw spectrum
    if(GCVonly){
      return(Inf)
    }else{
      if(length(level) == 1){
        object$Information$smooth$smooth.kernels[[level]] <- kernel("daniel", 0)
        object$Information$smooth$GCV[level] <- Inf
      }else{
        object$Information$smooth$smooth.kernels <- kernel("daniel", 0)
        object$Information$smooth$GCV <- Inf
      }
      invisible(object)
    }
  }
  
  if( 2 * kernel$m + 1 >= object$Information$dimensions$T ){
    stop(paste0("Smoothing kernel ", attributes(kernel)$name, " is too wide for time series."))
  }
  
  if(GCVonly){
    NCOL <- object$Information$dimensions$P
  }else{
    NCOL <- object$Information$dimensions$P * (1 + object$Information$dimensions$P) / 2
  }
  if(type == "all"){
    NCOL <- NCOL * object$Information$dimensions$J
  }
  Contribute <- rep(NA, NCOL)
  spectrum <- matrix(NA, nrow = object$Information$dimensions$T, ncol = NCOL)
  if(type == "by.level"){ Jset <- level }else{ Jset <- 1:object$Information$dimensions$J }
  l <- 1
  for(j in Jset){
    for(p in 1:object$Information$dimensions$P){
      if(GCVonly){ Q <- p }else{ Q <- 1:p }
      for(q in Q){
        spectrum[, l] <- object$spectrum[p, q, j, ]
        Contribute[l] <- as.numeric(p == q && j %in% level)
        l <- l + 1
      }
    }
  }
  
  Smoothed <- .C("SmoothEWS",
    spectrum = as.double(spectrum),
    T = as.integer(nrow(spectrum)),
    N = as.integer(ncol(spectrum)),
    M = as.integer(kernel$m),
    kernel = as.double(kernel$coef[c(kernel$m:1, 0, 1:kernel$m) + 1]),
    Contribute = as.integer(Contribute),
    eps = as.double(.Machine$double.eps),
    GCV = as.double(0.0),
    ErrorCode = as.integer(0L), PACKAGE = "mvLSW"
  )

  if(Smoothed$ErrorCode != 0L){
    cat("If you see this message then a problem has occured in executing the command.\n")
    cat("First, please check that you are supplying the arguments correctly.\n")
    cat("If the error persists, contact the package maintainer quoting error code:",Smoothed$ErrorCode,"\n")
    stop("Error in calling C when smoothing mvEWS.")
  }

  if(GCVonly){
    return(Smoothed$GCV)
  }else{
    Smoothspectrum <- matrix(Smoothed$spectrum, nrow = object$Information$dimensions$T)
    l <- 1
    for(j in Jset){
      for(p in 1:object$Information$dimensions$P){
        for(q in 1:p){
          object$spectrum[p, q, j, ] <- Smoothspectrum[, l]
          if(p != q) object$spectrum[q, p, j, ] <- Smoothspectrum[, l] 
          l <- l + 1
        }
      }
    }
    if(type == "by.level"){
      object$Information$smooth$smooth.kernels[[level]] <- kernel
      object$Information$smooth$GCV[level] <- Smoothed$GCV
    }else{
      object$Information$smooth$smooth.kernels <- kernel
      object$Information$smooth$GCV <- Smoothed$GCV
    }
    invisible(object)
  }
}

### Evaluate Multivariate Wavelet (Partial) Coherence ###

coherence <- function(
  object, 
  partial = FALSE){
  
  if(!is.mvLSW(object)) stop("Invalid 'object' argument.")
  if(object$Information$correction$min.eig.val <= 0) 
    stop("Invalid 'object' argument. Matrices must be positive definite.")
  if(!is.logical(partial) || length(partial) != 1) stop("Invalid 'partial' argument.")
  
  rho.jk <- function(Sjk, partial){
    if(partial){ M <- solve(Sjk) }else{ M<-Sjk }
	D <- diag(nrow(M))
	diag(D) <- 1 / sqrt(diag(M))
	C <- D %*% M %*% D
	if(partial){ C <- -C }
	return(C)
  }

  for(j in 1:object$Information$dimension$J){
    for(k in 1:object$Information$dimension$T){
	  object$spectrum[, , j, k] <- rho.jk(object$spectrum[, , j, k], partial)
	}
  }
  invisible(object)
}

### Sample a P-variate Multivariate Locally Stationary Process ###

rmvLSW <- function(
  Transfer = NULL, 
  Spectrum = NULL, 
  noiseFN = rnorm,
  ...){
  
  if(is.null(Transfer)){
    if(!is.mvLSW(Spectrum)) stop("Invalid 'Spectrum' argument.")
	object <- Spectrum2Transfer(Spectrum)
  }else{
    if(!is.mvLSW(Transfer)) stop("Invalid 'Transfer' argument.")
    P <- Transfer$Information$dimensions$P
    for(q in 2:P){ 
      for(p in 1:(q - 1)){
        if(any(Transfer$spectrum[p, q, , ] != 0))
          stop("Invalid 'Transfer' argument, must be lower triangular.")
      }
    }
    object <- Transfer
  }

  if(object$Information$correction$min.eig.val < 0 && 
    object$Information$correction$min.eig.val != -Inf) 
      stop("Invalid 'Transfer' or 'Spectrum' argument. 
        Matrices must be positive definite.")
  if(!is.function(noiseFN)) stop("Invalid 'noiseFN' argument.")
  if(names(formals(noiseFN))[1] != "n") stop("Invalid 'noiseFN' argument.")

  Epsilon <- noiseFN(n = prod(unlist(object$Information$dimension)), ...)
  if(!is.numeric(Epsilon)) 
    stop("Function 'noiseFN' should only return numeric values.")
  if(length(Epsilon) != prod(unlist(object$Information$dimension))) 
    stop("Function 'noiseFN' returns too many values.")
  Epsilon <- array(Epsilon,dim = object$Information$dimension)
  for(j in 1:object$Information$dimension$J){
    for(k in 1:object$Information$dimension$T){
	  Epsilon[, j, k] <- 2^j * object$spectrum[, , j, k] %*% Epsilon[, j, k]
    }
  }
  
  Transform <- function(Epsilon.p, T, wavelet){
    Univariate.Spec <- cns(T, filter.number = wavelet$filter.number, wavelet$family)
    Univariate.Spec$D <- as.numeric(t(Epsilon.p))
    return(AvBasis(convert(Univariate.Spec)))
  }
  
  Sample <- apply(Epsilon, 1, Transform, object$Information$dimension$T,
    object$Information$wavelet)
  colnames(Sample) <- object$Information$names
  
  if(!is.null(attr(object,"time"))){
    TIME <- attr(object,"time")
  }else{
    TIME <- NULL
  }
  if(any(class(TIME)=="POSIXt")){
    TSsample <- zoo(x=Sample,order.by=TIME)
  }else if(is.ts(TIME)){
    TSsample <- ts(data=Sample,start=start(TIME),frequency=frequency(TIME))
  }else{
    TSsample <- ts(data=Sample,start=1,frequency=1)
  }
  return(TSsample)
}

### Convert Spectrum Matrices to Transfer Matrices & visa versa ###

Spectrum2Transfer <- function(
  object, 
  S2V = TRUE){
  
  if(!is.mvLSW(object)) stop("Invalid 'object' argument.")
  if(!is.logical(S2V) || length(S2V) != 1) stop("Invalid 'S2V' argument.")

  my.chol <- function( A ){
    ##Evaluates a sq.root matrix for semi-positive matrix A
    if(!is.numeric(A)) stop("Matrix must be numeric!") 
    if(!is.matrix(A)) A <- as.matrix(A)
    if(ncol(A) != nrow(A)) stop("Matrix must be square.")
    if(det(A)<0) stop("Matrix must be semi-positive definite.")

    L <- diag(A)>0
    V <- matrix(0,ncol=ncol(A),nrow=nrow(A))
    if(any(L)){
      if(det(A[L,L])>0){
        V[L,L] <- chol(A[L,L])
      }else{
        B <- A[L,L]; tmp <- V[L,L]
        K <- 1:ncol(B)
        for(i in 2:ncol(B)){
          same <- which(apply(as.matrix(B[,1:(i-1)]) == B[,i],2,all))
          if(length(same)>0) K[i] <- min(same)
        }
        KK <- K == (1:ncol(B))
        tmp[KK,KK] <- chol(B[KK,KK])
        for(i in which(!KK))  tmp[K[i],i] <- sqrt(diag(B)[K[i]])
        V[L,L] <- tmp
      }
    }
    return(V)
  }

  L <- upper.tri(diag(object$Information$dimension$P))
  if(S2V){
    if(object$Information$correction$min.eig.val < 0 && 
      object$Information$correction$min.eig.val != -Inf) 
        stop("Invalid 'object' argument. Matrices must be positive definite.")

    for(j in 1:object$Information$dimension$J){
      for(k in 1:object$Information$dimension$T){
	    V <- t.default(my.chol(object$spectrum[, , j, k])) ##Fn checks for pos.def matrix
		object$spectrum[, , j, k] <- V
	  }
    }
  }else{
    for(j in 1:object$Information$dimension$J){
      for(k in 1:object$Information$dimension$T){
	    V <- object$spectrum[, , j, k]
		if(!all(V[L] == 0)) stop("There exists a transfer matrix that is not lower triangle.")
		if(!all(diag(V) >= 0)) stop("There exists a transfer matrix that has a negative diagonal element.")
	    object$spectrum[, , j, k] <- V %*% t.default(V)
      }
    }
  }
  invisible(object)
}

### Calculate the Wavelet Autocorrelations: \Psi_{j,l}(\tau) ###

PsiJLmat <- function(
  J, 
  filter.number = 1, 
  family = "DaubExPhase"){
  
  if(!is.numeric(J) || length(J) != 1) stop("Invalid 'J' argument.")
  if(J%%1 != 0 || J <= 0) stop("Invalid 'J' argument.")
  if(!is.numeric(filter.number) || length(filter.number) != 1) 
    stop("Invalid 'filter.number' argument.")
  if(filter.number%%1 != 0 || filter.number < 1 || filter.number > 10) 
    stop("Invalid 'filter.number' argument.")
  family <- match.arg(family, c("DaubExPhase", "DaubLeAsymm"))

  plus <- -1
  psi <- 1
  while( psi[1] != 0 || psi[length(psi)] != 0 ){
    plus <- plus + 1
    Jtmp <- J + plus
    blank.wd <- wd(data = rep(0, 2^Jtmp), filter.number = filter.number, family = family)
	if(plus==0){
	  replaceD <- 1
	}else if(plus == 1){
	  replaceD <- c(0,1)
	}else{
      replaceD <- c( rep(0, 2^(plus - 1)),  1, rep(0,2^(plus - 1) - 1) )
	}
    psi <- wr(putD(w = blank.wd, level = plus, v = replaceD)) * 2^(J/2)
  }
  
  L <- length(psi)
  AutoCorr <- .C("AutoCorr",
	psi = as.double(psi), 
	L = as.integer(L), 
	J = as.integer(J), 
    PsiJL = vector("double", J * J * (2 * L + 1)),
	ErrorCode = 0L, PACKAGE = "mvLSW")
  
  if(AutoCorr$ErrorCode != 0L){
    cat("If you see this message then an problem has occured in executing the command.\n")
    cat("First, please check that you are supplying the arguments correctly.\n")
    cat("If the error persists, contact the package maintainer quoting error code:",AutoCorr$ErrorCode,"\n")
    stop("Error in calling C when evaluating auto-correlation wavelet functions.")
  }

  AutoCorr_array <- array(AutoCorr$PsiJL, dim = c(2 * L + 1, J, J)) ## dim order tau, l, j
  names(dim(AutoCorr_array)) <- c("tau", "j", "l")
  invisible(AutoCorr_array)
}

### Calculate the Wavelet Autocorrelation Inner Product: A^{\lambda}_{j,l,h} ###

AutoCorrIP<- function(
  J, 
  filter.number = 1, 
  family = "DaubExPhase", 
  crop = TRUE){
  
  if(!is.numeric(J) || length(J) != 1) stop("Invalid 'J' argument.")
  if(J%%1 != 0 || J <= 0) stop("Invalid 'J' argument.")
  if(!is.numeric(filter.number) || length(filter.number) != 1) 
    stop("Invalid 'filter.number' argument.")
  if(filter.number%%1 != 0 || filter.number < 1 || filter.number > 10) 
    stop("Invalid 'filter.number' argument.")
  family <- match.arg(family, c("DaubExPhase", "DaubLeAsymm"))
  
  Psi <- PsiJLmat(J, filter.number, family) #tau, j, l

  Psi_min <- apply(Psi != 0, 2:3, function(a){ min(which(a)) })
  Psi_max <- apply(Psi != 0, 2:3, function(a){ max(which(a)) })
  L <- (dim(Psi)[1] - 1) / 2
  AutoCorr <- .C("WaveCorrInnerProd",
    Psi = as.double(Psi),
    Low = as.integer(Psi_min - 1),
    Upp = as.integer(Psi_max - 1),
    L = as.integer(L),
    J = as.integer(J),
    A = vector("double", (2 * L + 1) * J * J * J),
    ErrorCode = 0L, PACKAGE = "mvLSW"
  )
  
  if(AutoCorr$ErrorCode != 0){
    cat("If you see this message then an problem has occured in executing the command.\n")
    cat("First, please check that you are supplying the arguments correctly.\n")
    cat("If the error persists, contact the package maintainer quoting error code:",AutoCorr$ErrorCode,"\n")
    stop("Error in calling C when calculating the auto-correlation wavelet inner product.")
  }
  ACWIP <- array(AutoCorr$A, dim = c(2 * L + 1, J, J, J))
  names(dim(ACWIP)) <- c("lambda", "j", "l", "h")
  if(crop) ACWIP <- ACWIP[-2^J:2^J + L + 1, , , ]
  invisible(ACWIP)
}

### Variance of the Multivariate Evolutionary Wavelet Spectrum ###

varEWS <- function(
  object, 
  ACWIP = NULL, 
  verbose = FALSE){

  if(!is.mvLSW(object)) stop("Invalid 'object' argument.")
  if(object$Information$smooth$smooth.type != "all") 
    stop("Variance only applicable with smooth type 'all'.")
  if(object$Information$smooth$smooth.kernels$m == 0)
    warning("Identity smoothing kernel is specified.")
  if(object$Information$correction$min.eig.val <= 0) 
    stop("Invalid 'object' argument. Matrices must be positive definite.")
  if(!is.logical(verbose) || length(verbose) != 1) 
    stop("Invalid 'verbose' argument.")
  
  if(is.null(ACWIP)){
    if(verbose) cat("Calculating generalised autocorrelation wavelet inner products.\n")
    ACWIP <- AutoCorrIP(
      J = object$Information$dimensions$J, 
      filter.number = object$Information$wavelet$filter.number,
      family = object$Information$wavelet$family,
	  crop = TRUE)
  }else{
    if(!is.array(ACWIP) || !is.numeric(ACWIP)) stop("Invalid 'ACWIP' argument.")
    if(length(dim(ACWIP)) != 4) stop("Invalid 'ACWIP' argument.")
    Dim <- c(2 * object$Information$dimensions$T + 1, rep(object$Information$dimensions$J, 3))
    if(any(dim(ACWIP) != Dim)) stop("Invalid 'ACWIP' argument")
  }
  
  Amat <- matrix(NA, ncol=object$Information$dimensions$J,
    nrow=object$Information$dimensions$J)
  for(j in 1:object$Information$dimensions$J){
    Amat[, j] <- ACWIP[object$Information$dimensions$T + 1, j, j, ]
  }
  Ainv <- solve(Amat)

  PeriodVar <- object
  PeriodVar$spectrum <- array(NA, dim = dim(object$spectrum))

  if(verbose) cat("Calculating spectral variance for signal pair:\n")
  if(verbose) cat(" ", 1:object$Information$dimensions$P, "\n")
  for(p in 1:object$Information$dimensions$P){
    if(verbose) cat(p)
    for(q in 1:p){
      if(verbose) cat(" *")
      PeriodVar$spectrum[p, q, , ] <- VarSpqJ(object, ACWIP, p, q, Ainv, verbose=FALSE)
      if(p != q){
        PeriodVar$spectrum[q, p, , ] <- PeriodVar$spectrum[p, q, , ]
      }
    }
    if(verbose){
      if(p < object$Information$dimensions$P){
        for(q in (p + 1):object$Information$dimensions$P) cat(" -")
      }
	  cat("\n")
	}
  }
  invisible(PeriodVar)
}

### Variance of wavelet (cross-) spectrum S^{(p,q)}_{j,k} for all k & j ###

VarSpqJ <- function(
  object, 
  ACWIP = NULL, 
  p = 1, 
  q = 1, 
  Ainv = NULL, 
  verbose = FALSE){

  if(!is.mvLSW(object)) stop("Invalid 'object' argument.")
  if(object$Information$smooth$smooth.type != "all") 
    stop("Variance only applicable with smooth type 'all'.")
  if(object$Information$correction$min.eig.val <= 0) 
    stop("Invalid 'object' argument. Matrices must be positive definite.")
  if(!is.logical(verbose) || length(verbose) != 1) stop("Invalid 'verbose' argument.")
  P <- object$Information$dimensions$P
  J <- object$Information$dimensions$J
  T <- object$Information$dimensions$T

  if(is.null(ACWIP)){
    if(verbose) cat("Calculating generalised autocorrelation wavelet inner products.\n")
    ACWIP <- AutoCorrIP( J = J, filter.number = object$Information$wavelet$filter.number,
      family = object$Information$wavelet$family, crop=TRUE)
  }else{
    if(!is.array(ACWIP) || !is.numeric(ACWIP)) stop("Invalid 'ACWIP' argument.")
    if(any(dim(ACWIP) != c(2 * T + 1, J, J, J))) stop("Invalid 'ACWIP' argument.")
  }

  if(is.null(Ainv)){
    Amat <- matrix(NA, ncol = J, nrow = J)
    for(l in 1:J) Amat[, l] <- ACWIP[T + 1, l, l, ]
    Ainv <- solve(Amat)
  }else{
    if(!is.numeric(Ainv) || !is.matrix(Ainv)) stop("Invalid 'Ainv' argument.")
    if(nrow(Ainv) != J || nrow(Ainv) != J) stop("Invalid 'Ainv' argument.")
  }

  kernel <- object$Information$smooth$smooth.kernels
  if(kernel$m == 0){
    Wts <- 1
  }else{
    Wts <- kernel$coef[c(kernel$m:1, 0, 1:kernel$m) + 1]
  }
  if(verbose) cat("Evaluating covariance of smoothed periodogram.\n")
  SmoothCovEst <- .C("SmoothCovEst",
    SpqJ = as.numeric(object$spectrum[p, q, , ]),
    SppJ = as.numeric(object$spectrum[p, p, , ]),
    SqqJ = as.numeric(object$spectrum[q, q, , ]),
    ACWIP = as.numeric(ACWIP), #2T+1xJxJxJ (lambda x j x l x h)
    J = as.integer(J),
    T = as.integer(T),
    M = as.integer(kernel$m),
    Wts = as.double(Wts),
    SmoothCovEst = vector("double", T * J * J),
    ErrorCode = 0L, PACKAGE = "mvLSW"
  )

  if(SmoothCovEst$ErrorCode != 0L){
    cat("If you see this message then an problem has occured in executing the command.\n")
    cat("First, please check that you are supplying the arguments correctly.\n")
    cat("If the error persists, contact the package maintainer quoting error code:",SmoothCovEst$ErrorCode,"\n")
    stop("Error in calling C when evaluating the smoothed covariance estimate.")
  }
  AllCovSmoothI <- array(SmoothCovEst$SmoothCovEst,dim = c(T, J, J))

  if(verbose) cat("Applying bias correction.\n")
  VarSpqj <- matrix(NA, nrow = J, ncol = T)
  for(j in 1:J){
    for(k in 1:T){
      VarSpqj[j, k] <- t.default(Ainv[, j]) %*% AllCovSmoothI[k, , ] %*% Ainv[, j]
    }
  }
  return(VarSpqj)
}

## Approximate Quantiles
ApxQuantile <- function(
  object, 
  var = NULL, 
  prob = 0.5, 
  ...){
  
  if(!is.mvLSW(object)) stop("Invalid 'object' argument.")
  if(is.null(var)){
    var <- varEWS(object, ...)
  }else{
    if(!is.mvLSW(var)) stop("Invalid 'var' argument.")
  }
  if(length(prob) != 1 || !is.numeric(prob)) stop("invalid 'prob' argument.")
  if(prob < 0 || prob > 1) stop("Invalid 'prob' argument.")
 
  quantile <- object
  quantile$spectrum <- object$spectrum + qnorm(prob) * sqrt(var$spectrum)

  invisible(quantile)
}

##########################################

##Matrix regularization
modchol <- function(A, tol = .Machine$double.eps){
  if(!is.numeric(tol) || length(tol) != 1) stop("Invalid 'tol' argument")
  if(tol <= 0) stop("Invalid 'tol' argument")
  if(!is.numeric(A)) stop("Invalid 'A' Argument")
  if(!is.matrix(A)){
    if(length(A) == 1){
      A <- as.matrix(A)
	}else{
      stop("Invalid 'A' Argument")    
    }
  }
  if(nrow(A) != ncol(A)) stop("Invalid 'A' Argument")
  if(any(abs(t(A) - A) > tol^(1/3))) stop("Invalid 'A' Argument")
  
  n <- nrow(A)
  g <- rep(0, n)
  E <- rep(0, n)
  P <- 1:n

  tau1 <- tol^(1/3)
  tau2 <- tol^(2/3)
  
  phase1 <- TRUE
  delta <- 0
  ming <- 0
  gamma <- max(abs(diag(A)))
  taugam <- tau1 * gamma

  if(n == 1){
    delta <- tau2 * abs(A[1, 1]) - A[1, 1]
    if(delta > 0) E[1] <- delta
    if(A[1, 1] == 0) E[1] <- tau2
    A[1, 1] <- A[1, 1] + E[1]
    if(A < 0) return(NULL)
    A[1, 1] <- sqrt(A[1, 1])
    return(list(L = A, g = NA, P = P, E = E, delta = delta))
  }
  
  for(j in 1:(n - 1)){
    #PHASE1
    if(phase1){
      #Find index if maximum diagonal element A[i,i] where i>=j
      imaxd <- (j:n)[diag(A)[j:n] == max(diag(A)[j:n])]
      imaxd <- imaxd[1]

      #pivot to the top the row and column with the max diag
      if(imaxd != j){
        TMP <- A
        TMP <- t(TMP)
        TMP[lower.tri(TMP)] <- A[lower.tri(A)]
        swap <- 1:n; swap[imaxd] <- j; swap[j] <- imaxd
        TMP <- TMP[swap, swap]
        A[lower.tri(A, diag = TRUE)] <- TMP[lower.tri(TMP, diag = TRUE)]
        P <- P[swap]
      }

      #check to normal chol iteration gives pos diag element, else go to phase 2
      if(A[j,j]>0){
        jdmin <- min(A[j + 1, j + 1], diag(A)[(j + 1):n] - A[(j + 1):n, j]^2 / A[j + 1, j + 1])
        if(jdmin < taugam) phase1 <- FALSE
      }else{
        phase1 <- FALSE
      }

      if(phase1){
        #Do the normal cholesky update if still in phase 1
        if(A[j, j] < 0) return(NULL)
        A[j, j] <- sqrt(A[j, j])
        A[(j + 1):n, j] <- A[(j + 1):n, j] / A[j, j]
        for(i in (j + 1):n) A[i, (j + 1):i] <- A[i, (j + 1):i] - A[i, j]*A[(j + 1):i, j]
        if(A[n, n] < 0) return(NULL)
        if(j == (n - 1)) A[n, n] <- sqrt(A[n, n])
      }else{
        #calculate the negatives of the lower gerschgorin bounds
        for(i in j:n){
          g[i] <- -A[i, i]
          if(i != j) g[i] <- g[i] + sum(abs(A[i, j:(i - 1)]))
          if(i != n) g[i] <- g[i] + sum(abs(A[(i + 1):n, i]))
        }
      }
    }

    #PHASE2

    if(!phase1){
      if(j != (n - 1)){
        #find the minimum negative ger bound
        iming <- (j:n)[g[j:n]==min(g[j:n])]
        iming <- iming[1]
        
        #pivot
        if(iming != j){
          TMP <- A
          TMP <- t(TMP)
          TMP[lower.tri(TMP)] <- A[lower.tri(A)]
          swap <- 1:n; swap[iming] <- j; swap[j] <- iming
          TMP <- TMP[swap, swap]
          A[lower.tri(A, diag = TRUE)] <- TMP[lower.tri(TMP, diag = TRUE)]
          P <- P[swap]
          g <- g[swap]
        }
        
        # calculate delta and add to the diagonal
        normj <- sum(abs(diag(A)[(j + 1):n]))
        delta <- max(c(normj, taugam) - A[j, j], 0, delta)
        E[j] <- delta
        A[j, j] <- A[j, j] + E[j]
        
        #Update g bound estimates
        if(A[j, j] != normj){
          g[(j + 1):n] <- g[(j + 1):n] + abs(A[(j + 1):n, j]) * (normj - A[j, j])/A[j, j]
        }
        
        # Chol update
        if(A[j, j] < 0) return(NULL)
        A[j, j] <- sqrt(A[j, j])
        A[(j + 1):n, j] <- A[(j + 1):n, j] / A[j, j]
        for(i in (j + 1):n) A[i, (j + 1):i] <- A[i, (j + 1):i] - A[i, j] * A[(j + 1):i, j]
      }else{
        TMP <- A[(n - 1):n, (n - 1):n]
        TMP[2, 1] <- TMP[1, 2]
        Eig <- eigen(TMP)
        l_hi <- Eig$values[1]
        l_lo <- Eig$values[2]
        delta1 <- tau2 * max((l_hi - l_lo)/(1 - tau2), gamma) - l_lo
        delta <- max(delta, 0, delta1)
        if(delta > 0){
          A[n - 1, n - 1] <- A[n - 1, n - 1] + delta
          A[n, n] <- A[n, n] + delta
          E[n-1]  <- E[n] <- delta
        }
        if(A[n - 1, n - 1]<0) return(NULL)
        A[n - 1, n - 1] <- sqrt(A[n - 1, n - 1])
        A[n, n - 1] <- A[n, n - 1]/A[n - 1, n - 1]
        A[n, n] <- A[n, n] - A[n, n - 1]^2
        if(A[n, n] < 0) return(NULL)
        A[n, n] <- sqrt(A[n, n])
      }
    }
  }
  A[upper.tri(A)] <- 0
  return(list(L = A, g = g, P = P, E = E, delta = delta))
}

## Calculate +ve Def matrix from matrix root estimate via modchol
PosDefEst <- function(A, tol = .Machine$double.eps){
  Rt <- modchol(A, tol)
  if(is.null(Rt)) return(NULL)
  o <- order(Rt$P)
  Aest <- (Rt$L %*% t(Rt$L))[o, o]
  return(Aest)
}


AdjPositiveDef <- function (object, tol = 1e-10) 
{
    if (!is.mvLSW(object)) 
        stop("Invalid 'object' argument.")
    if (!is.numeric(tol) || length(tol) != 1) 
        stop("Invalid 'tol' argument.")
    if (tol <= 0) 
        stop("Invalid 'tol' argument.")
    POSDEF_single <- function(SqMat, tol) {
        Eig <- eigen(SqMat, symmetric=TRUE)
        if (all(Eig$values > tol)) 
            return(list(Mat = SqMat, minEval = min(Eig$values)))
        SqMatB <- try(PosDefEst(SqMat, tol), silent = TRUE)
	if (!any(class(SqMatB) == "try-error") && !is.null(SqMatB)) {
          EigB <- eigen(SqMatB, symmetric=TRUE)
          if (all(EigB$values > tol)) 
              return(list(Mat = SqMatB, minEval = min(EigB$values)))
        }
        Eig$values[Eig$values < tol] <- tol
        SqMatC <- Eig$vectors %*% (Eig$values * t.default(Eig$vectors))
        return(list(Mat = SqMatC, minEval = tol))
    }
    Min <- Inf
    for (j in 1:object$Information$dimensions$J) {
        for (k in 1:object$Information$dimensions$T) {
            L <- POSDEF_single(object$spectrum[, , j, k], tol)
            Min <- min(Min, L$minEval)
            object$spectrum[, , j, k] <- L$Mat
        }
    }
    object$Information$correction$min.eig.val <- Min
    return(object)
}

########
## print, runs summary command
print.mvLSW <- function(x, ...){
  f <- function(x, txt){
    if(is.list(x)){
      n <- length(x)
      for(i in 1:n){
        na <- names(x)[i]
        txt <- paste0(txt,"$",na)
        f(x[[i]], txt)
        p <- gregexpr("[$]",txt)[[1]]
        p <- p[length(p)]
        txt <- substr(txt,1,p-1)
      }
    }else{
      cat("                  ",txt,"\n")
    }
  }
  
  cat("Class 'mvLSW' : Multivariate Locally Stationary Wavelet Object:\n")
  cat("       ~~~~~  : List with component names\n")
  f(x,"x")
  cat("\nsummary(.):\n----------\n")
  summary.mvLSW(object = x, ...)
}

## simulate, sample a mv time series from a given mvEWS
simulate.mvLSW <- function(object, nsim = 1, seed = NULL, ...){
  if(nsim!=1) stop("Can only simulate a single mvLSW time series at a time.")
  if(!is.null(seed)) set.seed(seed)
  return(rmvLSW(Spectrum = object, ...))
}

## Approximate pointwise confidence interval
ApxCI <- function(
  object, 
  var = NULL, 
  alpha = 0.05, 
  ...){

  if(!is.mvLSW(object)) stop("Invalid 'object' argument.")
  if(is.null(var)){
    var <- varEWS(object, ...)
  }else{
    if(!is.mvLSW(var)) stop("Invalid 'var' argument.")
  }
  if(length(alpha) != 1 || !is.numeric(alpha)) stop("invalid 'alpha' argument.")
  if(alpha <= 0 || alpha > 0.5) stop("Invalid 'alpha' argument.")

  L <- ApxQuantile(object, var, prob = alpha/2)  
  U <- ApxQuantile(object, var, prob = 1 - (alpha/2))  
  invisible( list( L = L, U = U ) )
}
