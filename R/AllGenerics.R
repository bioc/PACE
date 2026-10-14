#' @rdname paceFit
#' @export
setGeneric("paceFit", function(object, ...) standardGeneric("paceFit"))

#' @rdname paceModel
#' @export
setGeneric("paceModel", function(object, ...) standardGeneric("paceModel"))

#' @rdname paceModelGaussian
#' @export
setGeneric("paceModelGaussian", function(object, ...) standardGeneric("paceModelGaussian"))

#' @rdname paceShrink
#' @export
setGeneric("paceShrink", function(object, ...) standardGeneric("paceShrink"))

#' @rdname paceDecompose
#' @export
setGeneric("paceDecompose", function(object, ...) standardGeneric("paceDecompose"))

#' @rdname paceDrivers
#' @export
setGeneric("paceDrivers", function(object, ...) standardGeneric("paceDrivers"))

#' @rdname buildNeighbourhood
#' @export
setGeneric("buildNeighbourhood", function(object, ...) standardGeneric("buildNeighbourhood"))

#' @rdname ambientField
#' @export
setGeneric("ambientField", function(object, ...) standardGeneric("ambientField"))

#' @rdname cellContamination
#' @export
setGeneric("cellContamination", function(object, ...) standardGeneric("cellContamination"))

#' @rdname neighbourSlopes
#' @export
setGeneric("neighbourSlopes", function(object, ...) standardGeneric("neighbourSlopes"))

#' @rdname varianceDecomposition
#' @export
setGeneric("varianceDecomposition", function(object, ...) standardGeneric("varianceDecomposition"))

#' @rdname topDrivers
#' @export
setGeneric("topDrivers", function(object, ...) standardGeneric("topDrivers"))

#' @rdname pairVariance
#' @export
setGeneric("pairVariance", function(object, ...) standardGeneric("pairVariance"))

#' @rdname cellTypes
#' @export
setGeneric("cellTypes", function(object, ...) standardGeneric("cellTypes"))

#' @rdname paceParams
#' @export
setGeneric("paceParams", function(object, ...) standardGeneric("paceParams"))

#' @rdname fitState
#' @export
setGeneric("fitState", function(object, ...) standardGeneric("fitState"))

#' @rdname neighbourSlopes
#' @export
setGeneric("neighbourSlopes<-", function(object, value) standardGeneric("neighbourSlopes<-"))

#' @rdname varianceDecomposition
#' @export
setGeneric("varianceDecomposition<-",
           function(object, value) standardGeneric("varianceDecomposition<-"))

#' @rdname topDrivers
#' @export
setGeneric("topDrivers<-", function(object, value) standardGeneric("topDrivers<-"))
