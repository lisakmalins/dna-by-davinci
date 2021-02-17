expected_colnames <- c("abundance",
                       "number",
                       "cum_sum",
                       "cum_fraction",
                       "slope")

# Read first line only and split on any whitespace to get colnames
detected_colnames <-
  read_lines(source, n_max=1) %>%
  strsplit("\\s+") %>%
  unlist()

# Crash if columns in data are not as expected
if (!identical(detected_colnames, expected_colnames)) {
  stop(paste("Problem parsing", source,
             "\nExpected the following column names:\n",
              paste(expected_colnames, collapse=" "),
              "\nBut instead found:\n",
              paste(detected_colnames, collapse=" ")))
}
