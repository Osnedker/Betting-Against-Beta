# Indlæs nødvendige pakker
library(dplyr)
library(zoo)
library(matrixStats)
library(lubridate)

# Indlæs data
Marketindlæs <- read.csv2(file.choose())
FFM3indlæs <- read.csv2(file.choose())

Market <- data.frame(lapply(Marketindlæs, function(x) as.numeric(as.character(x))))
FFM3 <- data.frame(lapply(FFM3indlæs, function(x) as.numeric(as.character(x))))

Market[Market == -99.99] <- NA

Date <- Market[1]
View(FFM3)

# Tager log 
Market <- log(Market[-1]/100+1)

# Tager log 
FFM3 <- log(FFM3/100+1)

# Tilføjer market return 
Market$MarketReturn <- FFM3[,2] + FFM3[,5]

# Udregn volatilitet (rolling)
Volatilities <- rollapply(
  Market,                   
  width = 250,               
  FUN = function(x) {
    valid_data <- na.omit(x)  # Fjern NA'er
    if (length(valid_data) < 120) return(NA)  # Kræver mindst 120 ikke-NA observationer
    sd(valid_data)  # Beregn standardafvigelsen
  }, 
  fill = NA,                
  align = "right"           
)

# Hent navne på alle kolonner undtagen kolonne 50
cols <- colnames(Market)

# Rolling 3 days
Market_3d <- rollapply(Market, width = 3, FUN = sum, align = "right", na.rm = TRUE)/3

# Beregn rolling correlation mellem hver kolonne og kolonne 50
Correlation <- do.call(cbind, lapply(cols, function(col) {
  rollapply(
    data = Market_3d[, c(col,"MarketReturn"), drop = FALSE],  
    width = 1250, 
    FUN = function(x) {
      valid_data <- na.omit(x)  # Fjern NA'er
      if (nrow(valid_data) < 750) return(NA)  # Kræver mindst 36 ikke-NA observationer
      cor(valid_data[, 1], valid_data[, 2], use = "pairwise.complete.obs")  # Beregn korrelation
    }, 
    by.column = FALSE, 
    fill = NA, 
    align = "right"
  )
}))

# Sæt kolonnenavne korrekt
colnames(Correlation) <- cols

str(Volatilities)
str(Correlation)

Volatilities <- Volatilities[-c(1,2), ]

# 1. Gem de originale betaer uden shrinkage
Beta <- Correlation * Volatilities / matrix(
  Volatilities[,"MarketReturn"], 
  nrow = nrow(Correlation), 
  ncol = ncol(Correlation), 
  byrow = FALSE
)

# Vi antager at du allerede har Beta_TS, Market og FFM3 som i din nuværende kode

# Funktion til at beregne shrinkede betaer
compute_beta_shrink <- function(Beta_TS, w, Beta_XS = 1) {
  return(w * Beta_TS + (1 - w) * Beta_XS)
}

# Shrinkage weights du vil teste
w_values <- c(0.4, 0.6, 0.8, 0.95, 1)

# Gem resultater for hver w
results_list <- list()
cumulative_returns_list <- list()

# Sørg for at Market har YearMonth kolonne
Market$Date <- as.Date(as.character(Date[,]), format = "%Y%m%d")
Market$YearMonth <- format(Market$Date, "%Y-%m")

for (w in w_values) {
  cat("\n--- Shrinkage weight w =", w, "---\n")
  
  # Beregn shrinkede betaer
  Beta_shrink <- compute_beta_shrink(Beta, w = w)
  
  # Tilføj dato og rens kolonner
  Beta_shrink <- as.data.frame(Beta_shrink)
  Beta_shrink$Date <- Date[3:nrow(Date), ]
  Beta_shrink <- Beta_shrink[, c(ncol(Beta_shrink), 1:(ncol(Beta_shrink)-1))]
  Beta_shrink <- Beta_shrink[, !colnames(Beta_shrink) %in% c("MarketReturn", "YearMonth")]
  
  Beta_shrink$Date <- as.Date(as.character(Beta_shrink$Date), format="%Y%m%d")
  Beta_shrink <- Beta_shrink %>%
    mutate(YearMonth = format(Date, "%Y-%m"))
  
  # Find sidste dato i hver måned
  last_days <- Beta_shrink %>%
    group_by(Year = year(Date), Month = month(Date)) %>%
    filter(Date == max(Date)) %>%
    ungroup()
  
  clean_cols <- setdiff(colnames(last_days), c("Date", "Year", "Month", "YearMonth"))
  
  z <- apply(last_days[, clean_cols], 1, function(x) {
    rank(x, ties.method = "average", na.last = NA)
  })
  
  z_bar <- sapply(z, function(row) mean(row, na.rm = TRUE))
  
  z_diff <- mapply(function(z_values, z_mean) {
    if (is.null(z_values)) return(NULL)
    return(z_values - z_mean)
  }, z, z_bar, SIMPLIFY = FALSE)
  
  k <- lapply(z_diff, function(x) {
    if (!is.null(x) && length(x) > 0) {
      return(2 / sum(abs(x), na.rm = TRUE))
    } else {
      return(NULL)
    }
  })
  
  w_H <- mapply(function(x, k_val) {
    if (!is.null(x) && !is.null(k_val)) pmax(x, 0) * k_val else NULL
  }, z_diff, k, SIMPLIFY = FALSE)
  
  w_L <- mapply(function(x, k_val) {
    if (!is.null(x) && !is.null(k_val)) pmax(-x, 0) * k_val else NULL
  }, z_diff, k, SIMPLIFY = FALSE)
  
  last_days <- last_days %>% mutate(YearMonth = format(Date, "%Y-%m"))
  
  # Tilføjer date til marked
  Market <- as.data.frame(Market)
  Market$Date <- Date[1:nrow(Date),]
  Market <- Market[, c(ncol(Market), 1:(ncol(Market)-1))]
  
  Market <- Market %>%
    mutate(Date = as.Date(as.character(Date), format="%Y%m%d"),
           YearMonth = format(Date, "%Y-%m"))  # YYYY-MM format
  
  # Opret også YYYYMM-kolonne i last_days for at matche w_H
  last_days <- last_days %>%
    mutate(YearMonth = format(Date, "%Y-%m"))  
  
  w_H_months <- unique(last_days$YearMonth)[1:length(w_H)]
  w_L_months <- unique(last_days$YearMonth)[1:length(w_L)]
  
  weights_lookup_H <- setNames(w_H, w_H_months)
  weights_lookup_L <- setNames(w_L, w_L_months)
  
  r_H_next <- sapply(1:(nrow(Market)-1), function(i) {
    current_month <- Market$YearMonth[i]  # Find YYYY-MM fra Market
    
    # Hvis vi ikke har vægte for denne måned, returnér NA
    if (!(current_month %in% names(weights_lookup_H))) {
      print(paste("Fejl! Måneden", current_month, "findes ikke i weights_lookup"))
      return(NA)
    }
    
    # Beregn vægtet afkast for næste dag
    return(sum(weights_lookup_H[[current_month]] * 
                 Market[i+1, names(weights_lookup_H[[current_month]]), drop = FALSE], 
               na.rm = TRUE))
  })
  
  # Beregn r_L_next baseret på YYYYMM matchning
  r_L_next <- sapply(1:(nrow(Market)-1), function(i) {
    current_month <- Market$YearMonth[i]  # Find YYYY-MM fra Market
    
    # Hvis vi ikke har vægte for denne måned, returnér NA
    if (!(current_month %in% names(weights_lookup_L))) {
      print(paste("Fejl! Måneden", current_month, "findes ikke i weights_lookup"))
      return(NA)
    }
    
    # Beregn vægtet afkast for næste dag
    return(sum(weights_lookup_L[[current_month]] * 
                 Market[i+1, names(weights_lookup_L[[current_month]]), drop = FALSE], 
               na.rm = TRUE))
  })
  
  Beta_H <- sapply(1:(nrow(Beta_shrink)-1), function(i) {
    current_month <- Beta_shrink$YearMonth[i]
    weights <- weights_lookup_H[[current_month]]
    if (is.null(weights)) return(NA)
    
    valid_assets <- intersect(names(weights), colnames(Beta_shrink))
    if (length(valid_assets) == 0) return(NA)
    
    weights_vec <- weights[valid_assets]
    beta_vec <- as.numeric(Beta_shrink[i, valid_assets, drop = FALSE])
    
    if (length(beta_vec) != length(weights_vec)) return(NA)
    if (any(is.na(beta_vec)) || any(is.na(weights_vec))) return(NA)
    
    return(sum(weights_vec * beta_vec, na.rm = TRUE))
  })
  
  Beta_L <- sapply(1:(nrow(Beta_shrink)-1), function(i) {
    current_month <- Beta_shrink$YearMonth[i]
    weights <- weights_lookup_L[[current_month]]
    if (is.null(weights)) return(NA)
    
    valid_assets <- intersect(names(weights), colnames(Beta_shrink))
    if (length(valid_assets) == 0) return(NA)
    
    weights_vec <- weights[valid_assets]
    beta_vec <- as.numeric(Beta_shrink[i, valid_assets, drop = FALSE])
    
    if (length(beta_vec) != length(weights_vec)) return(NA)
    if (any(is.na(beta_vec)) || any(is.na(weights_vec))) return(NA)
    
    return(sum(weights_vec * beta_vec, na.rm = TRUE))
  })
  
  rf <- FFM3[,5]
  min_length <- min(length(r_H_next), length(r_L_next), length(Beta_H), length(Beta_L), length(rf))
  r_H_next <- r_H_next[3:(min_length+2)]
  r_L_next <- r_L_next[3:(min_length+2)]
  rf <- rf[4:(min_length+3)]
  
  Beta_H[1:1250] <- NaN
  Beta_L[1:1250] <- NaN
  
  r_BAB <- (1 / Beta_L) * (r_L_next - rf) - (1 / Beta_H) * (r_H_next - rf)
  
  cat("Mean Beta_H:", mean(Beta_H, na.rm = TRUE), "\n")
  cat("Mean Beta_L:", mean(Beta_L, na.rm = TRUE), "\n")
  cat("Mean r_H_next:", mean(r_H_next, na.rm = TRUE), "\n")
  cat("Mean r_L_next:", mean(r_L_next, na.rm = TRUE), "\n")
  cat("Mean rf:", mean(rf, na.rm = TRUE), "\n")
  cat("Sum af vægt_h_60:", sum(w_H[[60]], na.rm = TRUE), "\n")
  cat("Sum af vægt_l_60:", sum(w_L[[60]], na.rm = TRUE), "\n")
  cat("Mean BAB:", mean(r_BAB, na.rm = TRUE),"\n")
  cat("Original BAB månedligseret afkast:", round((((1+mean(r_BAB, na.rm = TRUE))^(250/12))-1)*100, 10))
  
  mean_BAB <- mean(r_BAB, na.rm = TRUE)
  vol_BAB <- sd(r_BAB, na.rm = TRUE) * sqrt(250)
  sharpe <- (((1+mean_BAB)^250)-1) / vol_BAB
  
  results_list[[paste0("w_", w)]] <- list(
    return = ((1+mean_BAB)^250)-1,
    volatility = vol_BAB,
    sharpe = sharpe
  )
  
  cum_ret <- cumprod(1 + na.omit(r_BAB))
  cumulative_returns_list[[paste0("w_", w)]] <- cum_ret
}

print(do.call(rbind, results_list))

library(ggplot2)
library(tidyr)

cum_df <- do.call(cbind, cumulative_returns_list)
cum_df <- as.data.frame(cum_df)
cum_df$Index <- 1:nrow(cum_df)
cum_df_long <- pivot_longer(cum_df, -Index, names_to = "Shrinkage", values_to = "CumulativeReturn")

ggplot(cum_df_long, aes(x = Index, y = CumulativeReturn, color = Shrinkage)) +
  geom_line() +
  labs(title = "Kumulative afkast for forskellige shrinkage weights",
       x = "Tid (dage)", y = "Kumulativt afkast") +
  theme_minimal()
