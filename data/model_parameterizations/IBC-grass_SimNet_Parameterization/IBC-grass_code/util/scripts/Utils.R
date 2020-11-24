library(data.table)

read_data <- function(data_dir, file_type)
{
    main_dir = getwd()
    setwd(data_dir)

    files <- list.files(full.names = T)
    files <- files[which(grepl(file_type, files))]

    d <- bind_rows(map(.x = files,
                       .f = read_csv, col_names = TRUE, na = "NA", progress = TRUE))

    setwd(main_dir)
    return(d %>% as_tibble())
}

combine_data <- function(df_list, key)
{
    purrr::reduce(df_list,
                  left_join, by = key)
}

uv <- function(data) {
    r <- apply(data,
               2,
               function(x) (unique(x)))
    return(r)
}

gm_mean <- function(x, na.rm = TRUE, zero.propagate = FALSE) {
    if (any(x < 0, na.rm = TRUE)) {
        return(NaN)
    }
    if (zero.propagate) {
        if (any(x == 0, na.rm = TRUE)) {
            return(0)
        }
        exp(mean(log(x), na.rm = na.rm))
    } else {
        exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
    }
}
