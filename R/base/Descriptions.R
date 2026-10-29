Descriptions <- local({

  ADJACENCY <- function(web) {
    web > 0
  }

  INTERACTIONS <- function(web) {
    adj <- web > 0
    adj > 0
  }

  WEB_SIZE <- function(web) {
    nrow(web) * ncol(web)
  }

  CONNECTANCE <- function(web) {
    adj <- web > 0
    sum(adj) / WEB_SIZE(web)
  }

  GENERALITY <- function(web) {
    rowSums(web > 0)
  }

  VULNERALITY <- function(web) {
    colSums(web > 0)
  }

  list(ADJACENCY = ADJACENCY, INTERACTIONS = INTERACTIONS, WEB_SIZE = WEB_SIZE,
       CONNECTANCE = CONNECTANCE, GENERALITY = GENERALITY, VULNERALITY = VULNERALITY)
})
