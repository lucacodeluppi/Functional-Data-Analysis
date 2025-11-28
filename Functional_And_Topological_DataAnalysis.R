# 1. SET UP ENVIROMENT

# Load library
library(readr)
library(tidyverse)  # For data manipulation & visualization
library(lubridate)  # For handling datetime formats
library(ggplot2)    # For plotting
library(fda)        # For functional data analysis
library(dplyr)
library(ISOweek)
library(TDA)

# 2. LOAD DATA
data <- read_csv("Downloads/Functional/Functional_final/AUSstmfout.csv")
View(data)

# View the first few rows of the data to verify
head(data)

# Check the structure of the dataset
str(data)

# Check the shape (number of rows and columns)
dim(data)

aus <- data %>% 
  rename(
    age_0_14     = `0-14...5`,
    age_15_64    = `15-64...6`,
    age_65_74    = `65-74...7`,
    age_75_84    = `75-84...8`,
    age_85_plus  = `85+...9`,
    total        = `Total...10`,
    prop_0_14    = `0-14...11`,
    prop_15_64   = `15-64...12`,
    prop_65_74   = `65-74...13`,
    prop_75_84   = `75-84...14`,
    prop_85_plus = `85+...15`,
    total_prop   = `Total...16`
  )
str(aus)

# Check for missing values
missing_values <- colSums(is.na(aus))
print(missing_values)  # Displays missing values per column

# Remove duplicates (if any)
aus <- aus %>% distinct()

# Summary statistics
summary(aus)

aus <- aus %>%
  mutate(
    week_str = sprintf("%s-W%02d-1", Year, Week),
    date = ISOweek2date(week_str)
  )

# Visualizing Australia Weekly Mortality for both sex
ggplot(aus %>% filter(Sex == "b"), aes(x = date)) +
  geom_line(aes(y = total), color = "blue") +
  labs(
    title = "Australia Weekly Mortality (Both Sexes)",
    x = "Time",
    y = "Total Deaths"
  ) +
  theme_minimal()

# Distributions of weekly death
ggplot(aus %>% filter(Sex == "b"), aes(x = total)) +
  geom_histogram(fill = "blue", bins = 50, alpha = 0.6) +
  labs(
    title = "Distribution of Weekly Deaths (Australia)",
    x = "Deaths per Week",
    y = "Count of Weeks"
  ) +
  theme_minimal()

# Comparison between Male-Female-both
ggplot(aus, aes(x = date, y = total, color = Sex)) +
  geom_line() +
  theme_minimal()


# Check for outliers in total deaths
boxplot(aus$total, main="Boxplot of Total deaths", col="lightblue", horizontal=TRUE)

# 3. B SPLINE
# Save functional objects (one each year)
yearly_fd_list <- list()

# B-spline basis parameters
norder <- 6      # ordinal B-spline
nbasis <- 25     # number of basis (we have 52–53 week, so 25 should be fine)

# List of years in the dataset
years <- sort(unique(data$Year))

# from 2015 to 2025
for (yr in years) {
  
  # Filter data for years and both sexes ("b")
  df_year <- data %>%
    filter(Year == yr, Sex == "b") %>%
    arrange(Week) %>%
    dplyr::select(Week, `Total...10`)
  
  # Time index = week of the year
  time_index <- df_year$Week
  
  # Response variable = total weekly deaths
  deaths <- df_year$`Total...10`
  
  # Function domain: from first to last observed week
  rangeval <- c(min(time_index), max(time_index)) 
  
  # B-spline basis
  basis_spline <- create.bspline.basis(rangeval = rangeval, # c (1, 53)
                                       nbasis   = nbasis,   # 25
                                       norder   = norder)   # 6
  
  # Functional smoothing object (penalized)
  fdParobj <- fdPar(basis_spline, Lfdobj = 2, lambda = 0.1)
  
  # Smoothed function estimation
  smooth_fd <- smooth.basis(time_index, deaths, fdParobj)$fd
  yearly_fd_list[[as.character(yr)]] <- smooth_fd
  
  # Smoothed function plot for that year
  plot(
    smooth_fd,
    xlab = "Week",
    ylab = "Smoothed weekly deaths (Total)",
    main = paste("B-spline Smoothing - Australia Mortality - Year", yr)
  )
}

# 4. PERFORMING FUNCTIONAL PCA (FPCA)
# Extract Year and Week from datetime
# Create a list to store weekly functional data objects
weekly_fd_list <- list()

# Define B-spline parameters (over age, not time)
norder <- 4   # cubic B-splines
nbasis <- 5   # we have 5 age points so 5 basis functions is reasonable

# Define the "age" grid (e.g. midpoints of the age groups)
age_grid <- c(7, 40, 70, 80, 90)   # midpoints for 0–14, 15–64, 65–74, 75–84, 85+

# Restrict to both sexes, or choose "m"/"f" if you prefer
aus_both <- aus %>% filter(Sex == "b")

# Loop over years
for (yr in unique(aus_both$Year)) {
  
  df_year <- aus_both %>%
    filter(Year == yr)
  
  # Loop over weeks in that year
  for (wk in unique(df_year$Week)) {
    
    df_week <- df_year %>%
      filter(Week == wk)
    
    # If data for that year-week exists
    if (nrow(df_week) == 1) {
      
      # Vector of deaths by age group for that year-week
      y_values <- c(
        df_week$age_0_14,
        df_week$age_15_64,
        df_week$age_65_74,
        df_week$age_75_84,
        df_week$age_85_plus
      )
      
      # Range of the age domain
      rangeval <- range(age_grid)
      
      # B-spline basis over age
      basis_spline <- create.bspline.basis(rangeval, nbasis, norder)
      
      # Functional parameter object (smoothing over age)
      fdParobj <- fdPar(basis_spline, Lfdobj = 2, lambda = 0.1)
      
      # Smooth functional data: age → deaths(age)
      smooth_fd <- smooth.basis(age_grid, y_values, fdParobj)$fd
      
      # Store functional object, indexed by "Year_WeekX"
      key <- paste0(yr, "_Week", wk)
      weekly_fd_list[[key]] <- smooth_fd
      
      plot(
        smooth_fd,
        xlab = "Age",
        ylab = "Smoothed deaths by age",
        main = paste("Age-mortality curve -", yr, "Week", wk)
      )
    }
  }
}


# Convert list of fd objects into a coefficient matrix
fd_matrix <- do.call(
  cbind,
  lapply(weekly_fd_list, function(fd) fd$coefs)
)

# Basis object (same for all weekly curves)
fd_basis <- weekly_fd_list[[1]]$basis

# Build final fd object
weekly_fd_data <- fd(coef = fd_matrix, basisobj = fd_basis)

# PC1: overall level of age-specific mortality

# PC2: contrast between old vs younger ages

# PC3: shifts or curvature changes across intermediate ages

# FPCA with first 3 principal components
pca_results <- pca.fd(weekly_fd_data, nharm = 3)

par(mfrow = c(1,3))

plot(pca_results$harmonics[1],
     main = "1st Principal Component\n(Age-Mortality Variation)")

plot(pca_results$harmonics[2],
     main = "2nd Principal Component\n(Age-Mortality Variation)")

plot(pca_results$harmonics[3],
     main = "3rd Principal Component\n(Age-Mortality Variation)")

par(mfrow = c(1,1))

cat("Variance Explained:\n")
print(pca_results$varprop)


# 5. - FUNCTIONAL REGRESSION (FRegress)

# For each Year–Week, build two functions of age:
# – 𝑌(age) deaths by age (levels)
# – X (age) = proportion by age
## ================================
## Functional regression over WEEKS
## Response: total weekly deaths
## Predictor: proportion 85+ (prop_85_plus)
## ================================

week_range <- range(aus_both$Week)
week_range

# B-spline on weekly domain
norder_week <- 6      # B-spline order
nbasis_week <- 25     # number of basis

basis_week <- create.bspline.basis(
  rangeval = week_range,
  nbasis   = nbasis_week,
  norder   = norder_week
)

# Smoothing object 
fdPar_week <- fdPar(basis_week, Lfdobj = 2, lambda = 0.1)
deaths_fd_list   <- list()  # Y(t): total deaths
prop85_fd_list   <- list()  # X(t): prop_85_plus
year_ids         <- c()

years <- sort(unique(aus_both$Year))

for (yr in years) {
  
  df_year <- aus_both %>%
    filter(Year == yr) %>%
    arrange(Week)
  
  # time grid
  t_grid <- df_year$Week
  
  # response: total deaths
  y_vals <- df_year$total
  
  # predictor: proportion 85+
  x_vals <- df_year$prop_85_plus
  
  # smoothing su settimana
  deaths_fd <- smooth.basis(
    argvals = t_grid,
    y       = y_vals,
    fdParobj = fdPar_week
  )$fd
  
  prop85_fd <- smooth.basis(
    argvals = t_grid,
    y       = x_vals,
    fdParobj = fdPar_week
  )$fd
  
  deaths_fd_list[[as.character(yr)]] <- deaths_fd
  prop85_fd_list[[as.character(yr)]] <- prop85_fd
  year_ids <- c(year_ids, yr)
}

# Combine curves
deaths_fd_combined <- fd(
  coef = do.call(cbind, lapply(deaths_fd_list, function(fdobj) fdobj$coefs)),
  basisobj = basis_week
)

prop85_fd_combined <- fd(
  coef = do.call(cbind, lapply(prop85_fd_list, function(fdobj) fdobj$coefs)),
  basisobj = basis_week
)

nyears <- length(year_ids)

# -------------------------------
# Functional regression setup:
# Y(t) = beta0(t) + beta1(t) * X(t) + error
# -------------------------------

xfdlist <- list(
  const   = rep(1, nyears),     
  prop85  = prop85_fd_combined    
)

# Beta
betabasis_const <- create.constant.basis(week_range)
betafd_const    <- fd(0, betabasis_const)
betafdPar_const <- fdPar(betafd_const)

# Beta for prop85: same B-spline base on week
betabasis_prop  <- basis_week
nbasis_prop     <- betabasis_prop$nbasis
betafd_prop     <- fd(rep(0, nbasis_prop), betabasis_prop)
betafdPar_prop  <- fdPar(betafd_prop)

betalist <- list(
  const  = betafdPar_const,
  prop85 = betafdPar_prop
)

fregress_time <- fRegress(
  y = deaths_fd_combined,    
  xfdlist = xfdlist,           
  betalist = betalist          
)

# -------------------------------
# Plot beta(t) on weeks
# -------------------------------
plot(
  fregress_time$betaestlist$prop85,
  xlab = "Week",
  ylab = expression(beta[1](t) ~ "for prop85(t)"),
  main = "Functional regression coefficient over weeks"
)

# 5.1 ACTUAL VS FITTED
# fitted curves from functional regression
fitted_fd <- fregress_time$yhatfdobj
nyears <- length(year_ids)

for (i in seq_len(nyears)) {
  
  # observed curve for i-esimo year
  obs_i <- fd(
    coef     = deaths_fd_combined$coefs[, i, drop = FALSE],
    basisobj = deaths_fd_combined$basis
  )
  
  # fitted curve for the same year
  fit_i <- fd(
    coef     = fitted_fd$coefs[, i, drop = FALSE],
    basisobj = fitted_fd$basis
  )
  
  plot(
    obs_i,
    col  = "black",
    lwd  = 2,
    xlab = "Week",
    ylab = "Total deaths",
    main = paste("Observed vs fitted mortality curve - Year", year_ids[i])
  )
  
  lines(fit_i, col = "red", lwd = 2)
  
  legend(
    "topleft",
    legend = c("Observed", "Fitted"),
    col    = c("black", "red"),
    lwd    = 2,
    bty    = "n"
  )
}


# 6. – FUNCTIONAL ANOVA (FANOVA) 
# Defining group yearly based [PRE-POST COVID]
# Define groups based on year, from names like "2015_Week1" 
group_labels <- sapply(names(weekly_fd_list), function(id) {
  year <- as.numeric(substr(id, 1, 4))
  if (year <= 2019) {
    return("Group_A")  # Pre-pandemic (2015–2019)
  } else {
    return("Group_B")  # Pandemic period (2020+)
  }
})

# Split curves into groups 
group_A_fds <- weekly_fd_list[group_labels == "Group_A"]
group_B_fds <- weekly_fd_list[group_labels == "Group_B"]

# Proceed only if both groups have data 
if (length(group_A_fds) > 0 & length(group_B_fds) > 0) {
  
  fd_matrix_A <- do.call(cbind, lapply(group_A_fds, function(fdobj) fdobj$coefs))
  fd_matrix_B <- do.call(cbind, lapply(group_B_fds, function(fdobj) fdobj$coefs))
  
  basis_obj <- weekly_fd_list[[1]]$basis
  
  groupA_fd <- fd(coef = fd_matrix_A, basisobj = basis_obj)
  groupB_fd <- fd(coef = fd_matrix_B, basisobj = basis_obj)
  
  # --- 1) Mean curves (come prima, parte descrittiva) ---
  mean_A <- mean.fd(groupA_fd)
  mean_B <- mean.fd(groupB_fd)
  
  ylim_range <- range(c(mean_A$coefs, mean_B$coefs), na.rm = TRUE)
  
  plot(
    mean_A,
    col  = "blue",
    lwd  = 2,
    ylab = "Mortality (age profile)",
    xlab = "Age",
    main = "Functional ANOVA: Mean Age-Mortality by Group",
    ylim = ylim_range
  )
  lines(mean_B, col = "red", lwd = 2)
  legend(
    "topleft",
    legend = c("Group A (<=2019)", "Group B (>=2020)"),
    col    = c("blue", "red"),
    lwd    = 2
  )
  
  ## FANOVA with Fperm.fd (test + p-value)
  
  # Number of curves per group
  nA <- dim(groupA_fd$coefs)[2]
  nB <- dim(groupB_fd$coefs)[2]
  n  <- nA + nB
  
  # fd object
  Y_fd <- fd(
    coef     = cbind(groupA_fd$coefs, groupB_fd$coefs),
    basisobj = basis_obj
  )
  
  # Dummy: 0 = A (pre-2020), 1 = B (2020+)
  group_indicator <- c(rep(0, nA), rep(1, nB))
  
  xfdlist <- list(
    const = rep(1, n),           
    group = group_indicator   
  )
  
  rangeval <- basis_obj$rangeval
  
  # Beta for the intercept: constant on t
  betabasis_const <- create.constant.basis(rangeval)
  betafd_const    <- fd(0, betabasis_const)
  betafdPar_const <- fdPar(betafd_const)
  
  # Beta for the group effect:
  betabasis_group <- create.constant.basis(rangeval)
  betafd_group    <- fd(0, betabasis_group)
  betafdPar_group <- fdPar(betafd_group)
  
  betalist <- list(
    const = betafdPar_const,
    group = betafdPar_group
  )
  
  
  # Permutation F-test (FANOVA) 
  set.seed(123)
  fanova_res <- Fperm.fd(
    yfdPar  = Y_fd,
    xfdlist = xfdlist,
    betalist = betalist,
    nperm   = 999,
    plotres = TRUE
  )
  
  # p-value: H0 = 0 no difference between Group A and Group B
  cat("Global FANOVA p-value (group effect):", fanova_res$pval, "\n")
  
} else {
  print("Error: One of the groups has no functional data!")
}


# 7. TDA ON WEEKLY MORTALITY
# Total weekly death in 2020, both sexes ("b")

# Time series extraction: total weekly deaths in 2020, both sexes
df_sample <- aus %>%
  dplyr::filter(Year == 2020, Sex == "b") %>%
  dplyr::arrange(Week) %>%
  dplyr::select(Week, total)

# Time series: total deaths per week
ts_data <- df_sample$total
ts_scaled <- as.numeric(scale(ts_data))
# detrending
ts_detr <- ts_scaled - stats::filter(ts_scaled, rep(1/5, 5), sides = 2)
ts_detr[is.na(ts_detr)] <- ts_scaled[is.na(ts_detr)]

# Sliding window embedding to build point cloud
sliding_window <- function(ts, dim = 3, delay = 1) {
  n <- length(ts)
  m <- n - (dim - 1) * delay
  if (m <= 0) stop("Time series too short for embedding.")
  embed_matrix <- matrix(NA, nrow = m, ncol = dim)
  for (i in 1:dim) {
    embed_matrix[, i] <- ts[(1:m) + (i - 1) * delay]
  }
  return(embed_matrix)
}
# sliding window embedding
point_cloud <- sliding_window(ts_detr, dim = 5, delay = 2)

# Point cloud embedded plot
pairs(point_cloud,
      main = "Embedded Point Cloud (Weekly Mortality – Sliding Window)")

diag <- ripsDiag(
  X            = point_cloud,
  maxdimension = 1,
  maxscale     = 1,
  dist         = "euclidean"
)

# Persistence diagram
plot(diag[["diagram"]],
     main = "Persistence Diagram (H0 and H1) – Weekly Mortality (2020, both sexes)")

# Plot Persistence Barcode (H0 e H1)
plot(diag[["diagram"]],
     barcode = TRUE,
     main = "Persistence Barcode (H0 and H1) – Weekly Mortality")

# Function to extract a point cloud from a window of weeks
extract_window_cloud <- function(aus_data, year_val, week_min, week_max) {
  df_win <- aus_data %>%
    dplyr::filter(Year == year_val,
                  Sex == "b",
                  Week >= week_min,
                  Week <= week_max) %>%
    dplyr::arrange(Week)
  
  ts <- df_win$total
  
  if (length(ts) < 5) {
    stop("Time window too short for embedding. Increase the week range.")
  }
  
  sliding_window(ts, dim = 3, delay = 1)
}

# Window 1: weeks 5–15
pc1_df <- aus %>%
  dplyr::filter(Year == 2020, Sex == "b", Week >= 5, Week <= 15) %>%
  dplyr::arrange(Week)

# Window 2: weeks 30–40
pc2_df <- aus %>%
  dplyr::filter(Year == 2020, Sex == "b", Week >= 30, Week <= 40) %>%
  dplyr::arrange(Week)

cat("Serie 1 (Week 5–15):\n")
print(pc1_df$total)

cat("Serie 2 (Week 30–40):\n")
print(pc2_df$total)

print(all.equal(pc1_df$total, pc2_df$total))

# Extract point clouds for two week windows
pc1 <- extract_window_cloud(aus, year_val = 2020, week_min = 5,  week_max = 15)
pc2 <- extract_window_cloud(aus, year_val = 2020, week_min = 30, week_max = 40)


pairs(pc1, main = "Point cloud – Weeks 5–15")
pairs(pc2, main = "Point cloud – Weeks 30–40")
maxscale1 <- max(dist(pc1))
maxscale2 <- max(dist(pc2))
maxscale  <- max(maxscale1, maxscale2)

d1 <- ripsDiag(
  X            = pc1,
  maxdimension = 1,
  maxscale     = maxscale,
  dist         = "euclidean",
  library      = "GUDHI"
)

d2 <- ripsDiag(
  X            = pc2,
  maxdimension = 1,
  maxscale     = maxscale,
  dist         = "euclidean",
  library      = "GUDHI"
)

par(mfrow = c(1,2))
plot(d1[["diagram"]], main = "Persistence Diagram – Weeks 5–15")
plot(d2[["diagram"]], main = "Persistence Diagram – Weeks 30–40")
par(mfrow = c(1,1))

# Bottleneck distance
bdist <- bottleneck(d1[["diagram"]], d2[["diagram"]])
cat("Bottleneck distance:", round(bdist, 4), "\n")
