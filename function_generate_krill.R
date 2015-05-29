edsu <- matrix(0, ncol = (survey_width/detected_width), nrow = (survey_length/edsu_width))

#add random krill to cells
n_cells <- round(ncol(edsu)*nrow(edsu), 0)
rows <- c(round(runif(n_cells, 1, nrow(edsu)), 0))
cols <- c(round(runif(n_cells, 1, ncol(edsu)), 0))
for (i in 1:n_cells) {
  edsu[rows[i], cols[i]] <- rexp(1, 1/20)
}

#choose random locations for 100 krill swarms
krill_col <- round(runif(100, 1, ncol(edsu)))
krill_row <- round(runif(100, 1, nrow(edsu)))

#place krill into survey area
for (i in 1:20) {
  cols <- c(krill_col[i] - 2, krill_col[i] - 1, krill_col[i], krill_col[i] + 1, krill_col[i] - 2)
  rows <- c(krill_row[i] - 2, krill_row[i] - 1, krill_row[i], krill_row[i] + 1, krill_row[i] + 2)
  if (any(cols < 1)) cols[cols < 1] <- 1
  if (any(cols > ncol(edsu))) cols[cols > ncol(edsu)] <-  ncol(edsu)
  if (any(rows < 1)) rows[rows < 1] <- 1
  if (any(rows > nrow(edsu))) rows[rows > nrow(edsu)] <- nrow(edsu)
  edsu[rows, cols] <- rnorm(25, 50)
}

for (i in 21:100) {
  cols <- c(krill_col[i] - 1, krill_col[i], krill_col[i] + 1)
  rows <- c(krill_row[i] - 1, krill_row[i], krill_row[i] + 1)
  if (any(cols < 1)) cols[cols < 1] <- 1
  if (any(cols > ncol(edsu))) cols[cols > ncol(edsu)] <-  ncol(edsu)
  if (any(rows < 1)) rows[rows < 1] <- 1
  if (any(rows > nrow(edsu))) rows[rows > nrow(edsu)] <- nrow(edsu)
  edsu[rows, cols] <- rnorm(9, 20)
}