bhglm_prepare <- function(...) {
  utils::getFromNamespace("prepare", "BhGLM")(...)
}

bhglm_update.ptheta.group <- function(...) {
  utils::getFromNamespace("update.ptheta.group", "BhGLM")(...)
}

bhglm_update.ptheta.network <- function(...) {
  utils::getFromNamespace("update.ptheta.network", "BhGLM")(...)
}

bhglm_update.scale.p <- function(...) {
  utils::getFromNamespace("update.scale.p", "BhGLM")(...)
}

bhglm_update.theta.weights <- function(...) {
  utils::getFromNamespace("update.theta.weights", "BhGLM")(...)
}
