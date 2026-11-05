#' Résoudre un système d'équations linéaires
#'
#' Résout un système de la forme A %*% x = b à l'aide d'une méthode directe.
#' Trois méthodes sont proposées : élimination de Gauss avec pivot partiel,
#' décomposition LU avec pivot partiel et décomposition de Cholesky.
#'
#' @param A Matrice numérique carrée contenant les coefficients du système.
#' @param b Vecteur numérique contenant les seconds membres.
#' @param methode Méthode à utiliser : "gauss", "lu" ou "cholesky".
#' @param tol Tolérance numérique positive utilisée pour détecter les pivots
#'   nuls ou les matrices qui ne satisfont pas les conditions de Cholesky.
#'
#' @return Une liste contenant :
#' \describe{
#'   \item{solution}{Vecteur des valeurs des inconnues.}
#'   \item{methode}{Méthode utilisée.}
#'   \item{residu}{Norme euclidienne de A %*% solution - b.}
#' }
#'
#' @details
#' La méthode de Gauss transforme le système en un système triangulaire
#' supérieur, puis effectue une substitution arrière.
#'
#' La méthode LU factorise P A = L U, où P est une matrice de permutation,
#' L une matrice triangulaire inférieure et U une matrice triangulaire
#' supérieure.
#'
#' La méthode de Cholesky exige que A soit symétrique définie positive.
#' Elle factorise A = L %*% t(L). Elle n'utilise pas les équations normales.
#'
#' @examples
#' A <- matrix(c(2, 1, 1,
#'               1, 3, 2,
#'               1, 0, 0), nrow = 3, byrow = TRUE)
#' b <- c(4, 5, 1)
#' resoudre_systeme(A, b, methode = "gauss")
#'
#' # Cholesky : la matrice doit être symétrique définie positive.
#' S <- matrix(c(4, 1, 1, 3), nrow = 2)
#' resoudre_systeme(S, c(1, 2), methode = "cholesky")
#'
#' @export
resoudre_systeme <- function(A, b, methode = c("gauss", "lu", "cholesky"),
                             tol = 1e-12) {
  methode <- match.arg(methode)

  # Vérifier les données avant d'effectuer les calculs.
  .verifier_entrees(A, b, tol)

  solution <- switch(
    methode,
    gauss = .resoudre_gauss(A, b, tol),
    lu = .resoudre_lu(A, b, tol),
    cholesky = .resoudre_cholesky(A, b, tol)
  )

  # Mesurer l'écart entre le membre gauche calculé et le second membre.
  residu <- sqrt(sum((A %*% solution - b)^2))

  return(list(
    solution = as.vector(solution),
    methode = methode,
    residu = as.numeric(residu)
  ))
}


# -------------------------------------------------------------------------
# Vérification des arguments
# -------------------------------------------------------------------------

.verifier_entrees <- function(A, b, tol) {
  if (!is.matrix(A) || !is.numeric(A)) {
    stop("'A' doit être une matrice numérique.", call. = FALSE)
  }

  if (nrow(A) == 0 || ncol(A) == 0 || nrow(A) != ncol(A)) {
    stop("'A' doit être une matrice carrée non vide.", call. = FALSE)
  }

  if (any(!is.finite(A))) {
    stop("'A' ne doit contenir ni NA, ni NaN, ni Inf.", call. = FALSE)
  }

  if (!is.numeric(b) || !is.null(dim(b)) || length(b) != nrow(A)) {
    stop("'b' doit être un vecteur numérique de longueur égale à nrow(A).",
         call. = FALSE)
  }

  if (any(!is.finite(b))) {
    stop("'b' ne doit contenir ni NA, ni NaN, ni Inf.", call. = FALSE)
  }

  if (!is.numeric(tol) || length(tol) != 1 || !is.finite(tol) || tol <= 0) {
    stop("'tol' doit être un nombre strictement positif et fini.",
         call. = FALSE)
  }

  invisible(TRUE)
}


# -------------------------------------------------------------------------
# Substitution avant : résoudre L y = b
# -------------------------------------------------------------------------

.substitution_avant <- function(L, b, tol) {
  n <- length(b)
  y <- numeric(n)
  seuil <- tol * max(1, max(abs(L)))

  for (i in seq_len(n)) {
    # Les termes précédents sont déjà connus.
    somme <- if (i > 1) sum(L[i, seq_len(i - 1)] * y[seq_len(i - 1)]) else 0

    if (abs(L[i, i]) <= seuil) {
      stop("Pivot nul ou trop petit pendant la substitution avant.",
           call. = FALSE)
    }

    y[i] <- (b[i] - somme) / L[i, i]
  }

  y
}


# -------------------------------------------------------------------------
# Substitution arrière : résoudre U x = y
# -------------------------------------------------------------------------

.substitution_arriere <- function(U, y, tol) {
  n <- length(y)
  x <- numeric(n)
  seuil <- tol * max(1, max(abs(U)))

  for (i in seq.int(n, 1)) {
    # Les termes situés à droite sont déjà connus.
    somme <- if (i < n) sum(U[i, (i + 1):n] * x[(i + 1):n]) else 0

    if (abs(U[i, i]) <= seuil) {
      stop("Pivot nul ou trop petit pendant la substitution arrière.",
           call. = FALSE)
    }

    x[i] <- (y[i] - somme) / U[i, i]
  }

  x
}


# -------------------------------------------------------------------------
# Méthode de Gauss avec pivot partiel
# -------------------------------------------------------------------------

.resoudre_gauss <- function(A, b, tol) {
  n <- nrow(A)
  U <- A
  c <- b
  seuil <- tol * max(1, max(abs(A)))

  if (n > 1) {
    for (k in seq_len(n - 1)) {
      # Choisir, dans la colonne k, le plus grand pivot en valeur absolue.
      p <- which.max(abs(U[k:n, k])) + k - 1

      if (abs(U[p, k]) <= seuil) {
        stop("Le système est singulier ou numériquement mal conditionné.",
             call. = FALSE)
      }

      # Échanger les lignes de U et les valeurs correspondantes de c.
      if (p != k) {
        U[c(k, p), ] <- U[c(p, k), ]
        c[c(k, p)] <- c[c(p, k)]
      }

      # Annuler les coefficients sous le pivot.
      for (i in seq.int(k + 1, n)) {
        multiplicateur <- U[i, k] / U[k, k]
        U[i, k:n] <- U[i, k:n] - multiplicateur * U[k, k:n]
        c[i] <- c[i] - multiplicateur * c[k]
      }
    }
  }

  if (abs(U[n, n]) <= seuil) {
    stop("Le système est singulier ou numériquement mal conditionné.",
         call. = FALSE)
  }

  .substitution_arriere(U, c, tol)
}


# -------------------------------------------------------------------------
# Méthode LU avec pivot partiel : P A = L U
# -------------------------------------------------------------------------

.resoudre_lu <- function(A, b, tol) {
  n <- nrow(A)
  U <- A
  L <- diag(n)
  permutation <- seq_len(n)
  seuil <- tol * max(1, max(abs(A)))

  if (n > 1) {
    for (k in seq_len(n - 1)) {
      p <- which.max(abs(U[k:n, k])) + k - 1

      if (abs(U[p, k]) <= seuil) {
        stop("La décomposition LU a échoué : matrice singulière.",
             call. = FALSE)
      }

      if (p != k) {
        U[c(k, p), ] <- U[c(p, k), ]
        permutation[c(k, p)] <- permutation[c(p, k)]

        # Échanger uniquement les colonnes déjà calculées de L.
        if (k > 1) {
          L[c(k, p), seq_len(k - 1)] <-
            L[c(p, k), seq_len(k - 1)]
        }
      }

      for (i in seq.int(k + 1, n)) {
        multiplicateur <- U[i, k] / U[k, k]
        L[i, k] <- multiplicateur
        U[i, k:n] <- U[i, k:n] - multiplicateur * U[k, k:n]
      }
    }
  }

  if (abs(U[n, n]) <= seuil) {
    stop("La décomposition LU a échoué : matrice singulière.",
         call. = FALSE)
  }

  # P A = L U implique L y = P b, puis U x = y.
  b_permute <- b[permutation]
  y <- .substitution_avant(L, b_permute, tol)
  .substitution_arriere(U, y, tol)
}


# -------------------------------------------------------------------------
# Méthode de Cholesky : A = L %*% t(L)
# -------------------------------------------------------------------------

.resoudre_cholesky <- function(A, b, tol) {
  n <- nrow(A)
  echelle <- max(1, max(abs(A)))
  seuil <- tol * echelle

  # Une matrice de Cholesky doit être symétrique.
  if (!isTRUE(all.equal(A, t(A), tolerance = tol))) {
    stop("La méthode de Cholesky exige une matrice symétrique.",
         call. = FALSE)
  }

  L <- matrix(0, n, n)

  for (j in seq_len(n)) {
    somme_diag <- if (j > 1) sum(L[j, seq_len(j - 1)]^2) else 0
    valeur_diag <- A[j, j] - somme_diag

    if (!is.finite(valeur_diag) || valeur_diag <= seuil) {
      stop(paste(
        "La matrice n'est pas définie positive.",
        "La méthode de Cholesky ne peut pas être appliquée."
      ), call. = FALSE)
    }

    L[j, j] <- sqrt(valeur_diag)

    if (j < n) {
      for (i in seq.int(j + 1, n)) {
        somme <- if (j > 1) {
          sum(L[i, seq_len(j - 1)] * L[j, seq_len(j - 1)])
        } else {
          0
        }

        L[i, j] <- (A[i, j] - somme) / L[j, j]
      }
    }
  }

  y <- .substitution_avant(L, b, tol)
  .substitution_arriere(t(L), y, tol)
}
