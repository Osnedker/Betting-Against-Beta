# Indlæs nødvendige pakker
library(dplyr)
library(zoo)
library(matrixStats)
library(lubridate)

davidson_test <- function(A, B, Trim, order) {
  v <- sort(unique(c(A, B)))

  T_val <- ceiling(length(v) * Trim / 100)
  if (length(v) <= 2*T_val) {
    return(list(pvalue = NA, tvalue = NA))
  }
  Z <- v[(T_val + 1):(length(v) - T_val)]
  
  ZZ <- matrix(rep(Z, each = length(A)), nrow = length(A))
  Da <- pmax(ZZ - matrix(rep(A, length(Z)), ncol = length(Z)), 0)^(order - 1) / (order - 1)
  Db <- pmax(ZZ - matrix(rep(B, length(Z)), ncol = length(Z)), 0)^(order - 1) / (order - 1)
  DA <- colMeans(Da)
  DB <- colMeans(Db)
  
  t_vals <- sqrt(length(A)) * (DB - DA) / apply(Db - Da, 2, sd)
  
  tvalue <- min(t_vals, na.rm = TRUE)
  
  pvalue <- 2 * pnorm(tvalue)
  
  return(list(pvalue = pvalue, tvalue = tvalue))
}

# Indlæs data
Marketindlæs <- read.csv2(file.choose())
FFM3indlæs <- read.csv2(file.choose())

Market <- data.frame(lapply(Marketindlæs, function(x) as.numeric(as.character(x))))
FFM3 <- data.frame(lapply(FFM3indlæs, function(x) as.numeric(as.character(x))))

FFM3 <- FFM3[-1,]
Market <- Market[-1,]

Market[Market == -99.99] <- NA

Date <- Market[1]

# Tager log 
Market <- log(Market[-1]+1)

# Tager log 
FFM3 <- log(FFM3+1)

# Tilføjer market return 
Market$MarketReturn <- FFM3[,2]

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

Volatilities <- Volatilities[-c(1,2), ]

# --------- BEREGN VASICEK-SHRINKAGE VÆGTE ---------

# Gem de originale betaer uden shrinkage
Beta <- Correlation * Volatilities / matrix(
  Volatilities[,"MarketReturn"], 
  nrow = nrow(Correlation), 
  ncol = ncol(Correlation), 
  byrow = FALSE
)

# Beregn tværsnitsvariansen af betaer for hver dato
sigma2_XS <- apply(Beta, 1, function(row) var(row, na.rm = TRUE))

# Rullende varians (250-dages vindue) for hver aktie
sigma2_TS <- rollapply(
  Beta,
  width = 250,
  FUN = function(x) apply(x, 2, var, na.rm = TRUE),
  by.column = FALSE,
  fill = NA,
  align = "right"
)

# Justér sigma2_XS så den passer i længde med sigma2_TS
sigma2_XS_trimmed <- sigma2_XS[(length(sigma2_XS) - nrow(sigma2_TS) + 1):length(sigma2_XS)]

w_matrix <- matrix(NA, nrow = nrow(sigma2_TS), ncol = ncol(sigma2_TS))
colnames(w_matrix) <- colnames(sigma2_TS)

for (i in 1:nrow(sigma2_TS)) {
  for (j in 1:ncol(sigma2_TS)) {
    if (!is.na(sigma2_TS[i, j]) && !is.na(sigma2_XS_trimmed[i])) {
      w_matrix[i, j] <- 1 - (sigma2_TS[i, j] / (sigma2_TS[i, j] + sigma2_XS_trimmed[i]))
    }
  }
}

w <- mean(w_matrix, na.rm = TRUE)
Beta_XS <- 1

Beta_shrink <- w * Beta + (1-w) * Beta_XS

# Tilføjer Date til Beta shrink
Beta_shrink <- as.data.frame(Beta_shrink)
Date <- as.data.frame(Date)
Beta_shrink$Date <- Date[3:nrow(Date),]
Beta_shrink <- Beta_shrink[, c(ncol(Beta_shrink), 1:(ncol(Beta_shrink)-1))]
Beta_shrink <- Beta_shrink[, !colnames(Beta_shrink) %in% "MarketReturn"]

# Find den sidste dag i hver måned
Beta_shrink$Date <- as.Date(as.character(Beta_shrink$Date), format="%Y%m%d")

# Find den sidste dag i hver måned
last_days <- Beta_shrink %>%
  group_by(Year = year(Date), Month = month(Date)) %>%
  filter(Date == max(Date)) %>%
  ungroup()

# Find vægte 
z <- apply(last_days[, -c(1, (ncol(last_days)-1), ncol(last_days))], 1, function(x) {
  rank(x, ties.method = "average", na.last = NA)
})

z_bar <- sapply(z, function(row) mean(row, na.rm = TRUE))

z_diff <- mapply(function(z_values, z_mean) {
  if (is.null(z_values)) return(NULL)  # Håndter tomme lister
  return(z_values - z_mean)  # Beregn forskellen
}, z, z_bar, SIMPLIFY = FALSE)  # Returner som liste

# Definer k
k <- lapply(z_diff, function(x) {
  if (!is.null(x) && length(x) > 0) {
    return(2 / sum(abs(x), na.rm = TRUE))  # Normaliseringskonstant k
  } else {
    return(NULL)  # Bevar NULL hvis der ikke er data
  }
})

# Anvend k til at skalere vægtene
# Beregn w_H fra z_diff
w_H <- lapply(z_diff, function(x) {
  if (!is.null(x)) pmax(x, 0) else NULL  # Beholder kun positive værdier
})

w_H <- mapply(function(w, k_val) {
  if (!is.null(w) && !is.null(k_val)) return(w * k_val) else return(NULL)
}, w_H, k, SIMPLIFY = FALSE)

# Beregn w_L fra z_diff
w_L <- lapply(z_diff, function(x) {
  if (!is.null(x)) pmax(-x, 0) else NULL  # Bevarer kun negative værdier som positive
})

w_L <- mapply(function(w, k_val) {
  if (!is.null(w) && !is.null(k_val)) return(w * k_val) else return(NULL)
}, w_L, k, SIMPLIFY = FALSE)

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

# Opret en mapping fra YearMonth til w_H (brug kun YYYYMM som nøgle)
w_H_months <- unique(last_days$YearMonth)[1:length(w_H)]  # Sikrer korrekt mapping
w_L_months <- unique(last_days$YearMonth)[1:length(w_L)]

weights_lookup_H <- setNames(w_H, w_H_months)  
weights_lookup_L <- setNames(w_L, w_L_months)  

# ------------- STOKASTISK DOMINANS -------------
perform_monthly_sd_tests <- function(Market, sd_order = 3, trim_pct = 5) {

  dates <- unique(format(Market$Date, "%Y-%m"))
  
  market_dominates <- matrix(NA, nrow = length(dates), 
                             ncol = ncol(Market)-2)  
  asset_dominates <- matrix(NA, nrow = length(dates), 
                            ncol = ncol(Market)-2)
  
  colnames(market_dominates) <- colnames(Market)[2:(ncol(Market)-1)]
  colnames(asset_dominates) <- colnames(Market)[2:(ncol(Market)-1)]
  
  for (i in 1:length(dates)) {
    month_data <- Market[format(Market$Date, "%Y-%m") == dates[i], ]
    if (nrow(month_data) < 15) next 
    
    market_returns <- month_data$MarketReturn
    
    for (j in 2:(ncol(Market)-1)) {
      asset_name <- colnames(Market)[j]
      asset_returns <- month_data[[asset_name]]

      valid_indices <- !is.na(asset_returns) & !is.na(market_returns)
      if (sum(valid_indices) < 15) next
      
      asset_returns <- asset_returns[valid_indices]
      market_returns_filtered <- market_returns[valid_indices]

      md_test <- davidson_test(market_returns_filtered, asset_returns, trim_pct, sd_order)
      market_dominates[i, j-1] <- md_test$tvalue

      ad_test <- davidson_test(asset_returns, market_returns_filtered, trim_pct, sd_order)
      asset_dominates[i, j-1] <- ad_test$tvalue
    }
    
    if (i %% 10 == 0) {
      cat("Processed", i, "of", length(dates), "months for SD tests\n")
    }
  }
  
  return(list(
    dates = dates,
    market_dominates = market_dominates,
    asset_dominates = asset_dominates
  ))
}

identify_dominated_assets <- function(sd_results, threshold = -2.576) { 
  dates <- sd_results$dates
  market_dominates <- sd_results$market_dominates
  asset_dominates <- sd_results$asset_dominates

  dominated_by_market <- list()
  dominating_market <- list()

  for (i in 1:(length(dates)-1)) {  
    dom_by_market <- rep(FALSE, ncol(market_dominates))
    dom_market <- rep(FALSE, ncol(asset_dominates))
    
    for (j in 1:ncol(market_dominates)) {
      if (!is.na(market_dominates[i, j])) {
        if (market_dominates[i, j] < threshold) {
          dom_by_market[j] <- TRUE
        }
      }
      
      if (!is.na(asset_dominates[i, j])) {
        if (asset_dominates[i, j] < threshold) {
          dom_market[j] <- TRUE
        }
      }
    }
    
    dominated_by_market[[dates[i + 1]]] <- colnames(market_dominates)[dom_by_market]
    dominating_market[[dates[i + 1]]] <- colnames(asset_dominates)[dom_market]
  }
  
  return(list(
    dominated_by_market = dominated_by_market,
    dominating_market = dominating_market
  ))
}

modify_weights_with_sd <- function(weights_H, weights_L, dominated_assets, month) {
  
  if (!(month %in% names(dominated_assets$dominated_by_market)) || 
      !(month %in% names(dominated_assets$dominating_market))) {
    return(list(w_H = weights_H, w_L = weights_L))
  }
  
  dominated <- dominated_assets$dominated_by_market[[month]]
  dominating <- dominated_assets$dominating_market[[month]]
  
  new_weights_H <- weights_H
  new_weights_L <- weights_L
  
  for (asset in dominated) {
    if (asset %in% names(new_weights_L)) {
      new_weights_L[asset] <- 0
    }
  }
  
  for (asset in dominating) {
    if (asset %in% names(new_weights_H)) {
      new_weights_H[asset] <- 0
    }
  }
  
  new_weights_H[is.na(new_weights_H)] <- 0
  new_weights_L[is.na(new_weights_L)] <- 0
  
  if (sum(new_weights_L, na.rm = TRUE) > 0) {
    new_weights_L <- new_weights_L / sum(new_weights_L, na.rm = TRUE)
  }
  if (sum(new_weights_H, na.rm = TRUE) > 0) {
    new_weights_H <- new_weights_H / sum(new_weights_H, na.rm = TRUE)
  }
  
  return(list(w_H = new_weights_H, w_L = new_weights_L))
}

cat("Computing stochastic dominance tests...\n")
sd_results <- perform_monthly_sd_tests(Market, sd_order = 3, trim_pct = 5)

cat("Identifying persistently dominated assets...\n")
dominated_assets <- identify_dominated_assets(sd_results, threshold = -2.576)  # Using -1.96 for two-sided 5% test

weights_lookup_H_sd <- weights_lookup_H
weights_lookup_L_sd <- weights_lookup_L

cat("Modifying BAB weights based on stochastic dominance...\n")
for (month in names(weights_lookup_H)) {
  modified_weights <- modify_weights_with_sd(
    weights_lookup_H[[month]], 
    weights_lookup_L[[month]], 
    dominated_assets, 
    month
  )
  weights_lookup_H_sd[[month]] <- modified_weights$w_H
  weights_lookup_L_sd[[month]] <- modified_weights$w_L
}


# ------------- ORIGINAL CODE -------------
# Beregn r_H_next baseret på YYYYMM matchning
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

# Calculate r_H_next and r_L_next with SD-enhanced weights
r_H_next_sd <- sapply(1:(nrow(Market)-1), function(i) {
  current_month <- Market$YearMonth[i]
  
  if (!(current_month %in% names(weights_lookup_H_sd))) {
    return(NA)
  }
  
  return(sum(weights_lookup_H_sd[[current_month]] * 
               Market[i+1, names(weights_lookup_H_sd[[current_month]]), drop = FALSE], 
             na.rm = TRUE))
})

r_L_next_sd <- sapply(1:(nrow(Market)-1), function(i) {
  current_month <- Market$YearMonth[i]
  
  if (!(current_month %in% names(weights_lookup_L_sd))) {
    return(NA)
  }
  
  return(sum(weights_lookup_L_sd[[current_month]] * 
               Market[i+1, names(weights_lookup_L_sd[[current_month]]), drop = FALSE], 
             na.rm = TRUE))
})


# Tilføjer date til Beta_shrink
Beta_shrink <- as.data.frame(Beta_shrink)
Date <- as.data.frame(Date)
Beta_shrink$Date <- Date[3:nrow(Date),]
Beta_shrink <- Beta_shrink %>%
  mutate(Date = as.Date(as.character(Date), format="%Y%m%d"),
         YearMonth = format(Date, "%Y-%m"))  # YYYY-MM format

#Shrink betaer
Beta_H <- sapply(1:(nrow(Beta_shrink)-1), function(i) {
  current_month <- Beta_shrink$YearMonth[i]  # Find YYYY-MM fra Market
  
  # Hvis vi ikke har vægte for denne måned, returnér NA
  if (!(current_month %in% names(weights_lookup_H))) {
    print(paste("Fejl! Måneden", current_month, "findes ikke i weights_lookup"))
    return(NA)
  }
  
  # Beregn vægtet afkast for næste dag
  return(sum(weights_lookup_H[[current_month]] * 
               Beta_shrink[i, names(weights_lookup_H[[current_month]]), drop = FALSE], 
             na.rm = TRUE))
})

Beta_L <- sapply(1:(nrow(Beta_shrink)-1), function(i) {
  current_month <- Beta_shrink$YearMonth[i]  # Find YYYY-MM fra Market
  
  # Hvis vi ikke har vægte for denne måned, returnér NA
  if (!(current_month %in% names(weights_lookup_L))) {
    print(paste("Fejl! Måneden", current_month, "findes ikke i weights_lookup"))
    return(NA)
  }
  
  # Beregn vægtet afkast for næste dag
  return(sum(weights_lookup_L[[current_month]] * 
               Beta_shrink[i, names(weights_lookup_L[[current_month]]), drop = FALSE], 
             na.rm = TRUE))
})

Beta_H_SD <- sapply(1:(nrow(Beta_shrink) - 1), function(i) {
  current_month <- Beta_shrink$YearMonth[i]  # Find YYYY-MM from Beta
  
  if (!(current_month %in% names(weights_lookup_H_sd))) {
    return(NA)
    
  } else {
    return(sum(weights_lookup_H_sd[[current_month]] * 
                 Beta_shrink[i, names(weights_lookup_H_sd[[current_month]]), drop = FALSE], 
               na.rm = TRUE))
  }
})

Beta_L_SD <- sapply(1:(nrow(Beta_shrink) - 1), function(i) {
  current_month <- Beta_shrink$YearMonth[i]  # Find YYYY-MM from Beta
  
  if (!(current_month %in% names(weights_lookup_L_sd))) {
    return(NA)
    
  } else {
    return(sum(weights_lookup_L_sd[[current_month]] * 
                 Beta_shrink[i, names(weights_lookup_L_sd[[current_month]]), drop = FALSE], 
               na.rm = TRUE))
  }
})

## Udregner BAB-faktor
# Antag at rf er en vektor med 1181 risikofrie renter
# Træk den risikofrie rente ud
rf <- FFM3[,5]

# Find den mindste fælles længde
min_length <- min(length(r_H_next), length(r_L_next), 
                  length(r_H_next_sd), length(r_L_next_sd), 
                  length(Beta_H), length(Beta_L), 
                  length(Beta_H_SD), length(Beta_L_SD), 
                  length(rf))

# Trunkér alle variabler til samme længde
r_H_next <- r_H_next[3:(min_length+2)]
r_L_next <- r_L_next[3:(min_length+2)]
r_H_next_sd <- r_H_next_sd[3:(min_length+2)]
r_L_next_sd <- r_L_next_sd[3:(min_length+2)]
rf <- rf[4:(min_length+3)]

Beta_H[1:1250] <- NaN
Beta_L[1:1250] <- NaN
Beta_H_SD[1:1250] <- NaN
Beta_L_SD[1:1250] <- NaN


# Beregn BAB-faktoren (original)
r_BAB <- (1 / Beta_L) * (r_L_next - rf) - (1 / Beta_H) * (r_H_next - rf)

# Beregn SD-enhanced BAB-faktoren
r_BAB_sd <- (1 / Beta_L_SD) * (r_L_next_sd - rf) - (1 / Beta_H_SD) * (r_H_next_sd - rf)

# -------- LÅNEOMKOSTNINGSPARAMETRE ----------
transaction_cost_rate <- 0.001        # 10 bp
borrowing_cost_rate <- 0.02           # 2% 

calculate_turnover <- function(prev_weights, curr_weights, returns_next) {
  if (is.null(prev_weights) || is.null(curr_weights)) return(NA)
  
  shared_assets <- intersect(names(prev_weights), names(curr_weights))
  if (length(shared_assets) == 0) return(NA)
  
  returns_vec <- returns_next[shared_assets]
  if (any(is.na(returns_vec))) return(NA)
  
  adjusted_prev <- prev_weights[shared_assets] * (1 + returns_vec)
  adjusted_prev <- adjusted_prev / sum(adjusted_prev, na.rm = TRUE)
  
  turnover <- sum(abs(curr_weights[shared_assets] - adjusted_prev), na.rm = TRUE) / 2
  return(turnover)
}

transaction_costs_full <- rep(0, length(r_BAB))
borrowing_costs_full <- rep(0, length(r_BAB))

valid_months <- c()
month_lookup <- c()

months <- names(weights_lookup_H)

for (i in 2:length(months)) {
  month_prev <- months[i - 1]
  month_curr <- months[i]
  
  if (!(month_curr %in% Market$YearMonth)) next
  idx_curr <- which(Market$YearMonth == month_curr)[1] + 1
  if (idx_curr > nrow(Market)) next
  
  returns_next <- Market[idx_curr, ]
  returns_next <- returns_next[!(names(returns_next) %in% c("Date", "YearMonth", "MarketReturn"))]
  returns_next <- as.numeric(returns_next)
  names(returns_next) <- colnames(Market)[!(colnames(Market) %in% c("Date", "YearMonth", "MarketReturn"))]
  
  t_H <- calculate_turnover(weights_lookup_H[[month_prev]], weights_lookup_H[[month_curr]], returns_next)
  t_L <- calculate_turnover(weights_lookup_L[[month_prev]], weights_lookup_L[[month_curr]], returns_next)
  if (is.na(t_H) || is.na(t_L)) next
  
  transaction_cost <- (t_H + t_L) * transaction_cost_rate
  borrowing_cost <- sum(weights_lookup_H[[month_curr]], na.rm = TRUE) * (borrowing_cost_rate / 12)
  
  rebalance_day_index <- which(Market$YearMonth == month_curr)[1]
  transaction_costs_full[rebalance_day_index] <- transaction_cost
  borrowing_costs_full[rebalance_day_index] <- borrowing_cost
}  

transaction_costs_full_sd <- rep(0, length(r_BAB))
borrowing_costs_full_sd <- rep(0, length(r_BAB))

valid_months_sd <- c()
month_lookup_sd <- c()

months_sd <- names(weights_lookup_H_sd)

for (i in 2:length(months_sd)) {
  month_prev <- months_sd[i - 1]
  month_curr <- months_sd[i]
  
  if (!(month_curr %in% Market$YearMonth)) next
  idx_curr <- which(Market$YearMonth == month_curr)[1] + 1
  if (idx_curr > nrow(Market)) next
  
  returns_next <- Market[idx_curr, ]
  returns_next <- returns_next[!(names(returns_next) %in% c("Date", "YearMonth", "MarketReturn"))]
  returns_next <- as.numeric(returns_next)
  names(returns_next) <- colnames(Market)[!(colnames(Market) %in% c("Date", "YearMonth", "MarketReturn"))]
  
  t_H_sd <- calculate_turnover(weights_lookup_H_sd[[month_prev]], weights_lookup_H_sd[[month_curr]], returns_next)
  t_L_sd <- calculate_turnover(weights_lookup_L_sd[[month_prev]], weights_lookup_L_sd[[month_curr]], returns_next)
  if (is.na(t_H_sd) || is.na(t_L_sd)) next
  
  transaction_cost_sd <- (t_H_sd + t_L_sd) * transaction_cost_rate
  borrowing_cost_sd <- sum(weights_lookup_H_sd[[month_curr]], na.rm = TRUE) * (borrowing_cost_rate / 12)
  
  rebalance_day_index_sd <- which(Market$YearMonth == month_curr)[1]
  transaction_costs_full_sd[rebalance_day_index_sd] <- transaction_cost_sd
  borrowing_costs_full_sd[rebalance_day_index_sd] <- borrowing_cost_sd
  
  valid_months_sd <- c(valid_months_sd, month_curr)
  month_lookup_sd <- c(month_lookup_sd, rebalance_day_index_sd)
}

# Beregn BAB-faktoren (original) minus transaction cost
r_BAB_net <- (1 / Beta_L) * (r_L_next - rf) - (1 / Beta_H) * (r_H_next - rf) - transaction_costs_full - borrowing_costs_full

# Beregn SD-enhanced BAB-faktoren minus transaction cost
r_BAB_sd_net <- (1 / Beta_L_SD) * (r_L_next_sd - rf) - (1 / Beta_H_SD) * (r_H_next_sd - rf)- transaction_costs_full_sd - borrowing_costs_full_sd

# Tjek begge BAB faktorer
par(mfrow = c(2, 1))
plot(r_BAB, type = "l", main = "Original BAB Factor", col = "blue")
plot(r_BAB_sd, type = "l", main = "SD-Enhanced BAB Factor", col = "green")
par(mfrow = c(1, 1))

# Sammenlign performance
cat("\nPerformance sammenligning:\n")
cat("Original BAB månedligseret afkast:", round((((1+mean(r_BAB, na.rm = TRUE))^(250/12))-1)*100, 10), "%\n")
cat("SD-Enhanced BAB månedligseret afkast:", round((((1+mean(r_BAB_sd, na.rm = TRUE))^(250/12))-1)*100, 10), "%\n")

# Transaktionsomkostninger
cat("Original BAB månedligseret afkast minus transaktionsomkostninger:", round((((1+mean(r_BAB_net, na.rm = TRUE))^(250/12))-1)*100, 10), "%\n")
cat("SD-Enhanced BAB månedligseret afkast minus transaktionsomkostninger:", round((((1+mean(r_BAB_sd_net, na.rm = TRUE))^(250/12))-1)*100, 10), "%\n")


# ------------- FORTSÆT MED ORIGINAL ALPHA OG BETA BEREGNINGER -------------
Market_excess_return <- FFM3[,2] - FFM3[,5]
Market_excess_return <- Market_excess_return[4:(min_length+3)]

capm_model <- lm(r_BAB ~ Market_excess_return)

capm_alpha <- coef(capm_model)[1]*100*250/12
print(paste("CAPM Alpha (Original BAB):", round(capm_alpha, 2), "%"))

# Udtræk realiseret beta (hældningskoefficienten)
Beta_BAB_real <- coef(capm_model)[2]  # Beta er den anden koefficient
print(paste("Realiseret Beta for BAB:", round(Beta_BAB_real, 4)))

capm_model_sd <- lm(r_BAB_sd ~ Market_excess_return)

capm_alpha_sd <- coef(capm_model_sd)[1]*100*250/12
print(paste("CAPM Alpha (SD-Enhanced BAB):", round(capm_alpha_sd, 2), "%"))

# Udtræk realiseret beta (hældningskoefficienten)
Beta_SD_real <- coef(capm_model_sd)[2]  # Beta er den anden koefficient
print(paste("Realiseret Beta for SD:", round(Beta_SD_real, 4)))

###############################TRANS###############################

capm_model_trans <- lm(r_BAB_net ~ Market_excess_return)

capm_alpha_trans <- coef(capm_model_trans)[1]*100*250/12
print(paste("CAPM Alpha (Original BAB):", round(capm_alpha_trans, 2), "%"))

# Udtræk realiseret beta (hældningskoefficienten) med trans
Beta_BAB_real_trans <- coef(capm_model_trans)[2]  # Beta er den anden koefficient
print(paste("Realiseret Beta for BAB med trans:", round(Beta_BAB_real_trans, 4)))

capm_model_sd_trans <- lm(r_BAB_sd_net ~ Market_excess_return)

capm_alpha_sd_trans <- coef(capm_model_sd_trans)[1]*100*250/12
print(paste("CAPM Alpha (SD-Enhanced BAB):", round(capm_alpha_sd_trans, 2), "%"))

# Udtræk realiseret beta (hældningskoefficienten)
Beta_SD_real_trans <- coef(capm_model_sd_trans)[2]  # Beta er den anden koefficient
print(paste("Realiseret Beta for BAB SD med trans:", round(Beta_SD_real_trans, 4)))

SMB <- FFM3[, 3]
HML <- FFM3[, 4]

SMB <- SMB[3:(min_length+2)]
HML <- HML[3:(min_length+2)]

# Run the Fama-French 3-factor regression for original BAB
ff3_model <- lm(r_BAB ~ Market_excess_return + SMB + HML)

ff3_alpha <- coef(ff3_model)[1]*100*250/12
print(paste("Fama-French 3-Factor Alpha (Original BAB):", round(ff3_alpha, 2), "%"))

ff3_model_sd <- lm(r_BAB_sd ~ Market_excess_return + SMB + HML)

ff3_alpha_sd <- coef(ff3_model_sd)[1]*100*250/12
print(paste("Fama-French 3-Factor Alpha (SD-Enhanced BAB):", round(ff3_alpha_sd, 2), "%"))

#############################TRANS#############################

ff3_model_trans <- lm(r_BAB_net ~ Market_excess_return + SMB + HML)

ff3_alpha_trans <- coef(ff3_model_trans)[1]*100*250/12
print(paste("Fama-French 3-Factor Alpha (Original BAB) med trans:", round(ff3_alpha_trans, 2), "%"))

ff3_model_sd_trans <- lm(r_BAB_sd_net ~ Market_excess_return + SMB + HML)

ff3_alpha_sd_trans <- coef(ff3_model_sd_trans)[1]*100*250/12
print(paste("Fama-French 3-Factor Alpha (SD-Enhanced BAB) med trans:", round(ff3_alpha_sd_trans, 2), "%"))


#######################################Vola og Sharpe#######################################
# Beregn daglig volatilitet for original BAB
volatility_BAB <- sd(r_BAB, na.rm = TRUE)
volatility_BAB_sd <- sd(r_BAB_sd, na.rm = TRUE)
volatility_BAB_trans <- sd(r_BAB_net, na.rm = TRUE)
volatility_BAB_sd_trans <- sd(r_BAB_sd_net, na.rm = TRUE)

# Konverter til årlig volatilitet og angiv i procent
volatility_BAB_annual <- volatility_BAB * sqrt(250)  
print(paste("Årlig volatilitet for BAB:", round(volatility_BAB_annual*100, 2), "%"))

volatility_BAB_sd_annual <- volatility_BAB_sd * sqrt(250)  
print(paste("Årlig volatilitet for SD-Enhanced BAB:", round(volatility_BAB_sd_annual*100, 2), "%"))

volatility_BAB_annual_trans <- volatility_BAB_trans * sqrt(250)  
print(paste("Årlig volatilitet for BAB med trans:", round(volatility_BAB_annual_trans*100, 2), "%"))

volatility_BAB_sd_annual_trans <- volatility_BAB_sd_trans * sqrt(250)  
print(paste("Årlig volatilitet for SD-Enhanced BAB med trans:", round(volatility_BAB_sd_annual_trans*100, 2), "%"))


# Beregn gennemsnitligt dagligt afkast and daglig risikofri rente
mean_BAB <- mean(r_BAB, na.rm = TRUE)
mean_BAB_sd <- mean(r_BAB_sd, na.rm = TRUE)
mean_BAB_trans <- mean(r_BAB_net, na.rm = TRUE)
mean_BAB_sd_trans <- mean(r_BAB_sd_net, na.rm = TRUE)

# Konverter til årlige afkast (gang med 250) og omregn til procent
mean_BAB_annual <- ((1+mean_BAB)^250)-1
mean_BAB_sd_annual <- ((1+mean_BAB_sd)^250)-1
mean_BAB_annual_trans <- ((1+mean_BAB_trans)^250)-1
mean_BAB_sd_annual_trans <- ((1+mean_BAB_sd_trans)^250)-1

# Beregn årlig Sharpe-ratio
sharpe_BAB <- mean_BAB_annual / volatility_BAB_annual
sharpe_BAB_sd <- mean_BAB_sd_annual / volatility_BAB_sd_annual
sharpe_BAB_trans <- mean_BAB_annual_trans / volatility_BAB_annual_trans
sharpe_BAB_sd_trans <- mean_BAB_sd_annual_trans / volatility_BAB_sd_annual_trans
print(paste("Årlig Sharpe Ratio for BAB:", round(sharpe_BAB, 4)))
print(paste("Årlig Sharpe Ratio for SD-Enhanced BAB:", round(sharpe_BAB_sd, 4)))
print(paste("Årlig Sharpe Ratio for BAB med trans:", round(sharpe_BAB_trans, 4)))
print(paste("Årlig Sharpe Ratio for SD-Enhanced BAB med trans:", round(sharpe_BAB_sd_trans, 4)))

# Opsummerende statistik tabel
summary_table <- data.frame(
  Metric = c("Annualized Return (%)", "Annualized Volatility (%)", "Sharpe Ratio", 
             "CAPM Alpha (%)", "FF3 Alpha (%)"),
  Original_BAB = c(round(mean_BAB_annual*100, 2), 
                   round(volatility_BAB_annual*100, 2), 
                   round(sharpe_BAB, 2),
                   round(capm_alpha, 2),
                   round(ff3_alpha, 2)),
  SD_Enhanced_BAB = c(round(mean_BAB_sd_annual*100, 2), 
                      round(volatility_BAB_sd_annual*100, 2), 
                      round(sharpe_BAB_sd, 2),
                      round(capm_alpha_sd, 2),
                      round(ff3_alpha_sd, 2)),
  Original_BAB_trans = c(round(mean_BAB_annual_trans*100, 2), 
                   round(volatility_BAB_annual_trans*100, 2), 
                   round(sharpe_BAB_trans, 2),
                   round(capm_alpha_trans, 2),
                   round(ff3_alpha_trans, 2)),
  SD_Enhanced_BAB_trans = c(round(mean_BAB_sd_annual_trans*100, 2), 
                      round(volatility_BAB_sd_annual_trans*100, 2), 
                      round(sharpe_BAB_sd_trans, 2),
                      round(capm_alpha_sd_trans, 2),
                      round(ff3_alpha_sd_trans, 2))
)

print(summary_table)

# Lav en sammenlignende plot af de kumulerede afkast
r_BAB_clean <- na.omit(r_BAB)
r_BAB_SD_clean <- na.omit(r_BAB_sd)
r_BAB_clean_trans <- na.omit(r_BAB_net)
r_BAB_SD_clean_trans <- na.omit(r_BAB_sd_net)

cumulative_BAB <- cumprod(1 + r_BAB_clean)
cumulative_BAB_sd <- cumprod(1 + r_BAB_SD_clean)
cumulative_BAB_trans <- cumprod(1 + r_BAB_clean_trans)
cumulative_BAB_sd_trans <- cumprod(1 + r_BAB_SD_clean_trans)

plot(cumulative_BAB, type="l", col="blue", xlab="Tid", ylab="Kumulative afkast", ylim = c(0,3.5),
     main="Kumulative afkast for BAB og SDBAB")
lines(cumulative_BAB_sd, col="green")
lines(cumulative_BAB_trans, col="yellow")
lines(cumulative_BAB_sd_trans, col="red")
legend("topleft", legend=c("BAB", "SDBAB", "BAB-T", "SDBAB-T"), 
       col=c("blue", "green", "yellow", "red"), lty=1)

n_obs_BAB <- sum(!is.na(r_BAB))
se_return_BAB <- sd(r_BAB, na.rm = TRUE) / sqrt(n_obs_BAB)
t_stat_return_BAB <- mean(r_BAB, na.rm = TRUE) / se_return_BAB

capm_summary <- summary(capm_model)
t_stat_capm_alpha_BAB <- capm_summary$coefficients[1, 3]  # t-stat for intercept

ff3_summary <- summary(ff3_model)
t_stat_ff3_alpha_BAB <- ff3_summary$coefficients[1, 3]  # t-stat for intercept

# For SD-enhanced BAB
n_obs_BAB_sd <- sum(!is.na(r_BAB_sd))
se_return_BAB_sd <- sd(r_BAB_sd, na.rm = TRUE) / sqrt(n_obs_BAB_sd)
t_stat_return_BAB_sd <- mean(r_BAB_sd, na.rm = TRUE) / se_return_BAB_sd

capm_summary_sd <- summary(capm_model_sd)
t_stat_capm_alpha_BAB_sd <- capm_summary_sd$coefficients[1, 3]  # t-stat for intercept

ff3_summary_sd <- summary(ff3_model_sd)
t_stat_ff3_alpha_BAB_sd <- ff3_summary_sd$coefficients[1, 3]  # t-stat for intercept

##################################TRANS##################################

n_obs_BAB_trans <- sum(!is.na(r_BAB_net))
se_return_BAB_trans <- sd(r_BAB_net, na.rm = TRUE) / sqrt(n_obs_BAB_trans)
t_stat_return_BAB_trans <- mean(r_BAB_net, na.rm = TRUE) / se_return_BAB_trans

capm_summary_trans <- summary(capm_model_trans)
t_stat_capm_alpha_BAB_trans <- capm_summary_trans$coefficients[1, 3]  

ff3_summary_trans <- summary(ff3_model_trans)
t_stat_ff3_alpha_BAB_trans <- ff3_summary_trans$coefficients[1, 3]  

n_obs_BAB_sd_trans <- sum(!is.na(r_BAB_sd_net))
se_return_BAB_sd_trans <- sd(r_BAB_sd_net, na.rm = TRUE) / sqrt(n_obs_BAB_sd_trans)
t_stat_return_BAB_sd_trans <- mean(r_BAB_sd_net, na.rm = TRUE) / se_return_BAB_sd_trans

capm_summary_sd_trans <- summary(capm_model_sd_trans)
t_stat_capm_alpha_BAB_sd_trans <- capm_summary_sd_trans$coefficients[1, 3]  

ff3_summary_sd_trans <- summary(ff3_model_sd_trans)
t_stat_ff3_alpha_BAB_sd_trans <- ff3_summary_sd_trans$coefficients[1, 3] 

add_significance <- function(t_stat) {
  if (abs(t_stat) > 2.58) return(" ***")  # 1% level
  else if (abs(t_stat) > 1.96) return(" **")  # 5% level
  else if (abs(t_stat) > 1.65) return(" *")  # 10% level
  else return("")
}

monthly_BAB <- (((1+mean(r_BAB, na.rm = TRUE))^(250/12))-1)*100
monthly_BAB_sd <- (((1+mean(r_BAB_sd, na.rm = TRUE))^(250/12))-1)*100
monthly_BAB_trans <- (((1+mean(r_BAB_net, na.rm = TRUE))^(250/12))-1)*100
monthly_BAB_sd_trans <- (((1+mean(r_BAB_sd_net, na.rm = TRUE))^(250/12))-1)*100

summary_table_enhanced <- data.frame(
  Metric = c("Monthly Return (%)", "t-statistic", 
             "CAPM Alpha (%)", "t-statistic", 
             "FF3 Alpha (%)", "t-statistic"),
  Original_BAB = c(
    paste0(round(monthly_BAB, 2), add_significance(t_stat_return_BAB)), 
    round(t_stat_return_BAB, 2),
    paste0(round(capm_alpha, 2), add_significance(t_stat_capm_alpha_BAB)),
    round(t_stat_capm_alpha_BAB, 2),
    paste0(round(ff3_alpha, 2), add_significance(t_stat_ff3_alpha_BAB)),
    round(t_stat_ff3_alpha_BAB, 2)
  ),
  SD_Enhanced_BAB = c(
    paste0(round(monthly_BAB_sd, 2), add_significance(t_stat_return_BAB_sd)), 
    round(t_stat_return_BAB_sd, 2),
    paste0(round(capm_alpha_sd, 2), add_significance(t_stat_capm_alpha_BAB_sd)),
    round(t_stat_capm_alpha_BAB_sd, 2),
    paste0(round(ff3_alpha_sd, 2), add_significance(t_stat_ff3_alpha_BAB_sd)),
    round(t_stat_ff3_alpha_BAB_sd, 2)
  ),
  Original_BAB_trans = c(
    paste0(round(monthly_BAB_trans, 2), add_significance(t_stat_return_BAB_trans)), 
    round(t_stat_return_BAB_trans, 2),
    paste0(round(capm_alpha_trans, 2), add_significance(t_stat_capm_alpha_BAB_trans)),
    round(t_stat_capm_alpha_BAB_trans, 2),
    paste0(round(ff3_alpha_trans, 2), add_significance(t_stat_ff3_alpha_BAB_trans)),
    round(t_stat_ff3_alpha_BAB_trans, 2)
  ),
  SD_Enhanced_BAB_trans = c(
    paste0(round(monthly_BAB_sd_trans, 2), add_significance(t_stat_return_BAB_sd_trans)), 
    round(t_stat_return_BAB_sd_trans, 2),
    paste0(round(capm_alpha_sd_trans, 2), add_significance(t_stat_capm_alpha_BAB_sd_trans)),
    round(t_stat_capm_alpha_BAB_sd_trans, 2),
    paste0(round(ff3_alpha_sd_trans, 2), add_significance(t_stat_ff3_alpha_BAB_sd_trans)),
    round(t_stat_ff3_alpha_BAB_sd_trans, 2)
  )
)

print(summary_table_enhanced)


#P1 til 10 porteføljer
# Antal porteføljer
num_portfolios <- 10

# Find kvantilgrænser for opdeling
portfolio_breaks <- lapply(z, function(rank_values) {
  if (!is.null(rank_values)) {
    quantile(rank_values, probs = seq(0, 1, length.out = num_portfolios + 1), na.rm = TRUE)
  } else {
    return(NULL)
  }
})

# Initialiser en liste til hver portefølje
portfolio_weights <- vector("list", num_portfolios)
names(portfolio_weights) <- paste0("w_P", 1:num_portfolios)

# Loop gennem hver portefølje og tildel vægte
for (p in 1:num_portfolios) {
  portfolio_weights[[p]] <- mapply(function(z_values, breaks) {
    if (!is.null(z_values) && !is.null(breaks)) {
      # Initialiser en vektor med 0 for alle aktiver
      portfolio_vector <- setNames(rep(0, length(z_values)), names(z_values))
      
      # Justeret betingelse for at sikre, at det højeste rangnummer inkluderes
      if (p == num_portfolios) {
        in_portfolio <- (z_values >= breaks[p]) & (z_values <= breaks[p + 1])
      } else {
        in_portfolio <- (z_values >= breaks[p]) & (z_values < breaks[p + 1])
      }
      
      # Opdater porteføljeværdier
      portfolio_vector[in_portfolio] <- z_values[in_portfolio]
      
      return(portfolio_vector)
    } else {
      return(NULL)
    }
  }, z, portfolio_breaks, SIMPLIFY = FALSE)
}

# Omdan porteføljerne til ligevægtede vægte
portfolio_weights <- lapply(portfolio_weights, function(portfolio) {
  lapply(portfolio, function(weights) {
    if (!is.null(weights)) {
      # Find aktiver i porteføljen (dem der ikke er 0)
      active_assets <- weights[weights > 0]
      
      # Antal aktiver i porteføljen
      num_assets <- length(active_assets)
      
      # Hvis der er aktiver i porteføljen, sæt vægten til 1 / num_assets
      if (num_assets > 0) {
        weights[names(active_assets)] <- 1 / num_assets
      }
      
      return(weights)
    } else {
      return(NULL)
    }
  })
})

Market <- data.frame(lapply(Marketindlæs, function(x) as.numeric(as.character(x))))
FFM3 <- data.frame(lapply(FFM3indlæs, function(x) as.numeric(as.character(x))))

Market[Market == -99.99] <- NA

Date <- Market[1]

# Tilføjer market return 
Market$MarketReturn <- FFM3[,2]

# Tager log 
Market <- log(Market[-1]+1)
FFM3 <- log(FFM3+1)

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

# Opret en mapping fra YearMonth til w_H (brug kun YYYYMM som nøgle)
w_P1_months <- unique(last_days$YearMonth)[1:length(portfolio_weights$w_P1)]
w_P2_months <- unique(last_days$YearMonth)[1:length(portfolio_weights$w_P2)]
w_P3_months <- unique(last_days$YearMonth)[1:length(portfolio_weights$w_P3)]
w_P4_months <- unique(last_days$YearMonth)[1:length(portfolio_weights$w_P4)]
w_P5_months <- unique(last_days$YearMonth)[1:length(portfolio_weights$w_P5)]
w_P6_months <- unique(last_days$YearMonth)[1:length(portfolio_weights$w_P6)]
w_P7_months <- unique(last_days$YearMonth)[1:length(portfolio_weights$w_P7)]
w_P8_months <- unique(last_days$YearMonth)[1:length(portfolio_weights$w_P8)]
w_P9_months <- unique(last_days$YearMonth)[1:length(portfolio_weights$w_P9)]
w_P10_months <- unique(last_days$YearMonth)[1:length(portfolio_weights$w_P10)]

weights_lookup_P1 <- setNames(portfolio_weights$w_P1, w_P1_months)  
weights_lookup_P2 <- setNames(portfolio_weights$w_P2, w_P2_months)  
weights_lookup_P3 <- setNames(portfolio_weights$w_P3, w_P3_months)  
weights_lookup_P4 <- setNames(portfolio_weights$w_P4, w_P4_months)  
weights_lookup_P5 <- setNames(portfolio_weights$w_P5, w_P5_months)  
weights_lookup_P6 <- setNames(portfolio_weights$w_P6, w_P6_months)  
weights_lookup_P7 <- setNames(portfolio_weights$w_P7, w_P7_months)  
weights_lookup_P8 <- setNames(portfolio_weights$w_P8, w_P8_months)  
weights_lookup_P9 <- setNames(portfolio_weights$w_P9, w_P9_months)  
weights_lookup_P10 <- setNames(portfolio_weights$w_P10, w_P10_months)

# Liste over porteføljenavne
portfolio_names <- paste0("P", 1:10)

# Liste over weights_lookup for hver portefølje
weights_lookups <- list(
  weights_lookup_P1, weights_lookup_P2, weights_lookup_P3, weights_lookup_P4, weights_lookup_P5, 
  weights_lookup_P6, weights_lookup_P7, weights_lookup_P8, weights_lookup_P9, weights_lookup_P10
)
names(weights_lookups) <- portfolio_names

# Funktion til at beregne r_PX_next for en given portefølje
compute_r_P_next <- function(weights_lookup) {
  sapply(1:(nrow(Market)-1), function(i) {
    current_month <- Market$YearMonth[i]  # Find YYYY-MM fra Market
    
    # Hvis vi ikke har vægte for denne måned, returnér NA
    if (!(current_month %in% names(weights_lookup))) {
      print(paste("Fejl! Måneden", current_month, "findes ikke i weights_lookup"))
      return(NA)
    }
    
    # Beregn vægtet afkast for næste dag
    return(sum(weights_lookup[[current_month]] * 
                 Market[i+1, names(weights_lookup[[current_month]]), drop = FALSE], 
               na.rm = TRUE))
  })
}

# Beregn r_PX_next for alle porteføljer og gem som individuelle variabler
for (p in portfolio_names) {
  assign(paste0("r_", p, "_next"), compute_r_P_next(weights_lookups[[p]]))
}

# Tilføjer date til Beta Shrink
Beta_shrink <- as.data.frame(Beta_shrink)
Date <- as.data.frame(Date)
Beta_shrink$Date <- Date[4:nrow(Date),]
Beta_shrink <- Beta_shrink %>%
  mutate(Date = as.Date(as.character(Date), format="%Y%m%d"),
         YearMonth = format(Date, "%Y-%m"))  # YYYY-MM format

# Funktion til at beregne Beta_PX for en given portefølje
compute_Beta_P <- function(weights_lookup) {
  sapply(1:(nrow(Beta_shrink)-1), function(i) {
    current_month <- Beta_shrink$YearMonth[i]  # Find YYYY-MM fra Beta
    
    # Hvis vi ikke har vægte for denne måned, returnér NA
    if (!(current_month %in% names(weights_lookup))) {
      print(paste("Fejl! Måneden", current_month, "findes ikke i weights_lookup"))
      return(NA)
    }
    
    # Beregn vægtet beta
    return(sum(weights_lookup[[current_month]] * 
                 Beta_shrink[i, names(weights_lookup[[current_month]]), drop = FALSE], 
               na.rm = TRUE))
  })
}

# Beregn Beta_PX for alle porteføljer og gem som individuelle variabler
for (p in portfolio_names) {
  assign(paste0("Beta_", p), compute_Beta_P(weights_lookups[[p]]))
}

# Tilføjer date til Beta_shrink
Beta_shrink <- as.data.frame(Beta_shrink)
Date <- as.data.frame(Date)
Beta_shrink$Date <- Date[4:nrow(Date),]
Beta_shrink <- Beta_shrink %>%
  mutate(Date = as.Date(as.character(Date), format="%Y%m%d"),
         YearMonth = format(Date, "%Y-%m"))  # YYYY-MM format

rf <- FFM3[,5]
rf <- rf[2:(length(r_P1_next)+1)]
månedligrf <- (((1+mean(rf, na.rm = TRUE))^(250/12))-1)*100


###################################### Alpha Beregning ######################################
# Definer porteføljenavne
portfolio_names <- paste0("P", 1:10)

# Opret en liste til at gemme alpha-værdierne
capm_alphas <- numeric(length(portfolio_names))
ff3_alphas <- numeric(length(portfolio_names))
real_betas <- numeric(length(portfolio_names))

# Definer markedsafkast og faktorer
Market_excess_return <- FFM3[, 2][2:(length(r_P1_next)+1)]
SMB <- FFM3[, 3][2:(length(r_P1_next)+1)]
HML <- FFM3[, 4][2:(length(r_P1_next)+1)]

# Loop gennem hver portefølje og beregn alpha
for (i in seq_along(portfolio_names)) {
  portfolio_return <- get(paste0("r_", portfolio_names[i], "_next"))  # Hent r_PX_next
  
  # CAPM Regression
  capm_model <- lm(portfolio_return ~ Market_excess_return)
  capm_alphas[i] <- coef(capm_model)[1] * 100 * 250 / 12  # Justering til månedlig alpha
  real_betas[i] <- coef(capm_model)[2] # Justering til månedlig alpha
  
  # Fama-French 3-Factor Regression
  ff3_model <- lm(portfolio_return ~ Market_excess_return + SMB + HML)
  ff3_alphas[i] <- coef(ff3_model)[1] * 100 * 250 / 12  # Justering til årlig alpha
}

# Saml resultater i en data frame
alpha_results <- data.frame(
  Portfolio = portfolio_names,
  CAPM_Alpha = capm_alphas,
  FF3_Alpha = ff3_alphas,
  Real_Beta = real_betas
)

# Print resultaterne
print(alpha_results)

####################################SHARPE####################################

portfolio_names <- paste0("P", 1:10)

portfolio_volatility <- numeric(length(portfolio_names))
sharpe_ratios <- numeric(length(portfolio_names))

for (i in seq_along(portfolio_names)) {
  portfolio_return <- get(paste0("r_", portfolio_names[i], "_next"))  # Get portfolio returns
  
  portfolio_volatility[i] <- sd(portfolio_return, na.rm = TRUE) * sqrt(250)  # Annualized volatility
  
  mean_excess_return <- mean(portfolio_return, na.rm = TRUE)
  
  sharpe_ratios[i] <- mean_excess_return / portfolio_volatility[i]
}

sharpe_results <- data.frame(
  Portfolio = portfolio_names,
  Volatility = portfolio_volatility*100,
  Sharpe_Ratio = sharpe_ratios*100
)

print(sharpe_results)

portfolio_return_tstat <- numeric(10)
capm_alpha_tstat <- numeric(10)
ff3_alpha_tstat <- numeric(10)


add_significance <- function(t_stat) {
  if (abs(t_stat) > 2.58) return(" ***")  # 1%
  else if (abs(t_stat) > 1.96) return(" **")  # 5% 
  else if (abs(t_stat) > 1.65) return(" *")  # 10% 
  else return("")
}

for (i in 1:10) {
  portfolio_return <- get(paste0("r_P", i, "_next"))
  
  n_obs <- sum(!is.na(portfolio_return))
  se <- sd(portfolio_return, na.rm = TRUE) / sqrt(n_obs)
  portfolio_return_tstat[i] <- mean(portfolio_return, na.rm = TRUE) / se
  
  capm_model <- lm(portfolio_return ~ Market_excess_return)
  capm_summary <- summary(capm_model)
  capm_alpha_tstat[i] <- capm_summary$coefficients[1, 3]
  
  ff3_model <- lm(portfolio_return ~ Market_excess_return + SMB + HML)
  ff3_summary <- summary(ff3_model)
  ff3_alpha_tstat[i] <- ff3_summary$coefficients[1, 3]
}

monthly_excess_returns <- sapply(1:10, function(i) {
  r <- get(paste0("r_P", i, "_next"))
  monthly_return <- (((1+mean(r, na.rm = TRUE))^(250/12))-1)*100 - månedligrf
  return(paste0(round(monthly_return, 2), add_significance(portfolio_return_tstat[i])))
})

capm_alphas_formatted <- sapply(1:10, function(i) {
  paste0(round(capm_alphas[i], 2), add_significance(capm_alpha_tstat[i]))
})

ff3_alphas_formatted <- sapply(1:10, function(i) {
  paste0(round(ff3_alphas[i], 2), add_significance(ff3_alpha_tstat[i]))
})

portfolio_table <- data.frame(
  Portfolio = paste0("P", 1:10),
  Beta = sapply(1:10, function(i) round(mean(get(paste0("Beta_P", i)), na.rm = TRUE), 2)),
  Monthly_Excess_Return = monthly_excess_returns,
  Return_tstat = round(portfolio_return_tstat, 2),
  CAPM_Alpha = capm_alphas_formatted,
  CAPM_tstat = round(capm_alpha_tstat, 2),
  FF3_Alpha = ff3_alphas_formatted,
  FF3_tstat = round(ff3_alpha_tstat, 2)
)

print(portfolio_table)

#------------------------------ Ted Spread (Proposition 3) ------------------------------
# BetaSpread (ligger allerede i CSV)
BetaSpread <- (Beta_H - Beta_L) / (Beta_H*Beta_L)
BetaSpreadSD <- (Beta_H_SD - Beta_L_SD) /(Beta_H_SD*Beta_L_SD)

# Definér funktionen
export_r_BAB_to_excel <- function(r_BAB_vector, date_vector, filename = "r_BAB_export.xlsx") {
  # Lav et data.frame med dato og afkast
  export_df <- data.frame(Date = date_vector, r_BAB = r_BAB_vector)
  
  # Gem til Excel
  write_xlsx(export_df, path = filename)
  
  cat("r_BAB er nu eksporteret til", filename, "\n")
}

export_r_BAB_to_excel(r_BAB_vector = r_BAB, Date[4:4151,], filename = "/Users/OliverSnedker_1/Desktop/Kandidat/Speciale/BAB 49.xlsx")
export_r_BAB_to_excel(r_BAB_vector = r_BAB_sd, Date[4:4151,], filename = "/Users/OliverSnedker_1/Desktop/Kandidat/Speciale/BAB_SD 49.xlsx")
export_r_BAB_to_excel(r_BAB_vector = BetaSpread, Date[4:4151,], filename = "/Users/OliverSnedker_1/Desktop/Kandidat/Speciale/BS 49.xlsx")
export_r_BAB_to_excel(r_BAB_vector = BetaSpreadSD, Date[4:4151,], filename = "/Users/OliverSnedker_1/Desktop/Kandidat/Speciale/BS_SD 49.xlsx")

# Indlæs TED-data
TED <- read.csv(file.choose(), sep = ";")
TED$MarketReturn <- as.numeric(TED$MarketReturn)
View(TED)

model1 <- lm(R_BAB ~ Ted_lag + Ted_diff, data = TED)
summary(model1)
model2 <- lm(R_BAB_SD ~ Ted_lag + Ted_diff, data = TED)
summary(model2)

## Tilføjer kontrol variable
model3 <- lm(R_BAB ~ Ted_lag + Ted_diff + BS + BAB_lagged + MarketReturn, data = TED)
summary(model3)

model4 <- lm(R_BAB_SD ~ Ted_lag + Ted_diff + BS_SD + BAB_SD_lagged + MarketReturn, data = TED)
summary(model4)

## Tests
coeftest(model1, vcov. = NeweyWest(model1, prewhite = FALSE))
coeftest(model2, vcov. = NeweyWest(model2, prewhite = FALSE))
coeftest(model3, vcov. = NeweyWest(model3, prewhite = FALSE))
coeftest(model4, vcov. = NeweyWest(model4, prewhite = FALSE))

#--------------------------- Break even transaktionsomkostningsanalyse ----------------------------------
calculate_net_returns <- function(r_BAB, Beta_L, Beta_H, r_L_next, r_H_next, rf, 
                                  original_transaction_costs, borrowing_costs,
                                  transaction_cost_multiplier) {
  scaled_transaction_costs <- original_transaction_costs * transaction_cost_multiplier
  
  r_BAB_net <- (1 / Beta_L) * (r_L_next - rf) - (1 / Beta_H) * (r_H_next - rf) - 
    scaled_transaction_costs - borrowing_costs
  
  return(((1 + mean(r_BAB_net, na.rm = TRUE))^(250/12)) - 1)
}

find_breakeven_tc <- function(r_BAB, Beta_L, Beta_H, r_L_next, r_H_next, rf, 
                              original_transaction_costs, borrowing_costs,
                              base_tc_rate = 0.001,
                              tolerance = 1e-6, max_iter = 100) {
  
  lower <- 0
  upper <- 100
  
  zero_tc_return <- calculate_net_returns(r_BAB, Beta_L, Beta_H, r_L_next, r_H_next, rf,
                                          original_transaction_costs, borrowing_costs, 0)
  if (zero_tc_return <= 0) {
    return(0)  
  }
  
  high_tc_return <- calculate_net_returns(r_BAB, Beta_L, Beta_H, r_L_next, r_H_next, rf,
                                          original_transaction_costs, borrowing_costs, upper)
  if (high_tc_return > 0) {
    return(Inf) 
  }
  
  iter <- 0
  while ((upper - lower) > tolerance && iter < max_iter) {
    mid <- (upper + lower) / 2
    mid_return <- calculate_net_returns(r_BAB, Beta_L, Beta_H, r_L_next, r_H_next, rf,
                                        original_transaction_costs, borrowing_costs, mid)
    
    if (abs(mid_return) < tolerance) {
      return(mid * base_tc_rate) 
    } else if (mid_return > 0) {
      lower <- mid
    } else {
      upper <- mid
    }
    
    iter <- iter + 1
  }
  
  return(((upper + lower) / 2) * base_tc_rate)
}

plot_tc_sensitivity <- function(r_BAB, r_BAB_sd, Beta_L, Beta_H, Beta_L_SD, Beta_H_SD,
                                r_L_next, r_H_next, r_L_next_sd, r_H_next_sd, rf,
                                transaction_costs_full, borrowing_costs_full,
                                transaction_costs_full_sd, borrowing_costs_full_sd,
                                base_tc_rate = 0.001) {
  
  tc_multipliers <- seq(0, 100, by = 0.2)
  
  returns_original <- sapply(tc_multipliers, function(mult) {
    return <- calculate_net_returns(r_BAB, Beta_L, Beta_H, r_L_next, r_H_next, rf,
                                    transaction_costs_full, borrowing_costs_full, mult)
    return(return * 100)  
  })
  
  returns_sd <- sapply(tc_multipliers, function(mult) {
    return <- calculate_net_returns(r_BAB_sd, Beta_L_SD, Beta_H_SD, r_L_next_sd, r_H_next_sd, rf,
                                    transaction_costs_full_sd, borrowing_costs_full_sd, mult)
    return(return * 100)
  })
  
  crosses_zero_orig <- which(diff(sign(returns_original)) != 0)
  if (length(crosses_zero_orig) > 0) {
    idx <- crosses_zero_orig[1]
    x1 <- tc_multipliers[idx]
    x2 <- tc_multipliers[idx + 1]
    y1 <- returns_original[idx]
    y2 <- returns_original[idx + 1]
    breakeven_orig <- x1 + (0 - y1) * (x2 - x1) / (y2 - y1)
  } else {
    breakeven_orig <- max(tc_multipliers)
  }
  
  crosses_zero_sd <- which(diff(sign(returns_sd)) != 0)
  if (length(crosses_zero_sd) > 0) {
    idx <- crosses_zero_sd[1]
    x1 <- tc_multipliers[idx]
    x2 <- tc_multipliers[idx + 1]
    y1 <- returns_sd[idx]
    y2 <- returns_sd[idx + 1]
    breakeven_sd <- x1 + (0 - y1) * (x2 - x1) / (y2 - y1)
  } else {
    breakeven_sd <- max(tc_multipliers)
  }
  
  min_breakeven <- min(breakeven_sd, breakeven_orig) * 0.8
  max_breakeven <- max(breakeven_sd, breakeven_orig) * 1.2
  
  min_x_index <- which.min(abs(tc_multipliers - 0))  
  max_x_index <- which.min(abs(tc_multipliers - max_breakeven))
  
  min_x_index <- max(1, min_x_index)
  max_x_index <- min(length(tc_multipliers), max_x_index)
  
  y_values_in_range <- c(returns_original[min_x_index:max_x_index], 
                         returns_sd[min_x_index:max_x_index])
  y_min <- min(y_values_in_range, na.rm = TRUE) - 0.2
  y_max <- max(y_values_in_range, na.rm = TRUE) + 0.2
  
  plot(tc_multipliers * base_tc_rate * 10000, returns_original, type = "l", col = "blue",
       xlab = "Transaktionsomkostninger (basispoint)", ylab = "Årligt afkast (%)",
       main = "Følsomhedstest af BAB afkast over for transaktionsomkostninger",
       xlim = c(0, max_breakeven * base_tc_rate * 10000 * 1.05),  # Add 5% buffer
       ylim = c(y_min, y_max))
  
  lines(tc_multipliers * base_tc_rate * 10000, returns_sd, col = "green")
  
  abline(h = 0, lty = 2)
  abline(v = breakeven_orig * base_tc_rate * 10000, col = "blue", lty = 2)
  abline(v = breakeven_sd * base_tc_rate * 10000, col = "green", lty = 2)
  
  legend("topright", legend = c("Original BAB", "SD-forbedret BAB"),
         col = c("blue", "green"), lty = 1, cex = 0.8)
  
  return(list(
    breakeven_orig = breakeven_orig * base_tc_rate,
    breakeven_sd = breakeven_sd * base_tc_rate
  ))
}

cat("\n----- BREAK-EVEN TRANSAKTIONSOMKOSTNINGSANALYSE -----\n")

breakeven_tc_original <- find_breakeven_tc(
  r_BAB, Beta_L, Beta_H, r_L_next, r_H_next, rf,
  transaction_costs_full, borrowing_costs_full,
  base_tc_rate = 0.001
)

breakeven_tc_sd <- find_breakeven_tc(
  r_BAB_sd, Beta_L_SD, Beta_H_SD, r_L_next_sd, r_H_next_sd, rf,
  transaction_costs_full_sd, borrowing_costs_full_sd,
  base_tc_rate = 0.001
)

cat("Break-even transaction cost for Original BAB:", round(breakeven_tc_original * 10000, 2), "basis points\n")
cat("Break-even transaction cost for SD-Enhanced BAB:", round(breakeven_tc_sd * 10000, 2), "basis points\n")

par(mfrow = c(1, 1))
breakeven_results <- plot_tc_sensitivity(
  r_BAB, r_BAB_sd, Beta_L, Beta_H, Beta_L_SD, Beta_H_SD,
  r_L_next, r_H_next, r_L_next_sd, r_H_next_sd, rf,
  transaction_costs_full, borrowing_costs_full,
  transaction_costs_full_sd, borrowing_costs_full_sd
)

improvement_ratio <- breakeven_tc_sd / breakeven_tc_original
cat("\nImprovement in break-even transaction costs: ", 
    round((improvement_ratio - 1) * 100, 2), "% higher for SD-Enhanced BAB\n")

annual_turnover_original <- 1 / breakeven_tc_original
annual_turnover_sd <- 1 / breakeven_tc_sd

cat("\nBreak-even annual turnover for Original BAB:", round(annual_turnover_original, 2), "times\n")
cat("Break-even annual turnover for SD-Enhanced BAB:", round(annual_turnover_sd, 2), "times\n")

#--------------------------- Break even låneomkostningsanalyse --------------------------------------
calculate_net_returns_borrowing <- function(r_BAB, Beta_L, Beta_H, r_L_next, r_H_next, rf, 
                                            transaction_costs, original_borrowing_costs,
                                            borrowing_cost_multiplier) {
  scaled_borrowing_costs <- original_borrowing_costs * borrowing_cost_multiplier
  
  r_BAB_net <- (1 / Beta_L) * (r_L_next - rf) - (1 / Beta_H) * (r_H_next - rf) - 
    transaction_costs - scaled_borrowing_costs
  
  return(((1 + mean(r_BAB_net, na.rm = TRUE))^(250/12)) - 1)
}

find_breakeven_borrowing <- function(r_BAB, Beta_L, Beta_H, r_L_next, r_H_next, rf, 
                                     transaction_costs, original_borrowing_costs,
                                     base_borrowing_rate = 0.02,
                                     tolerance = 1e-6, max_iter = 100) {
  
  lower <- 0
  upper <- 30 
  
  zero_bc_return <- calculate_net_returns_borrowing(r_BAB, Beta_L, Beta_H, r_L_next, r_H_next, rf,
                                                    transaction_costs, original_borrowing_costs, 0)
  if (zero_bc_return <= 0) {
    return(0)  
  }
  
  high_bc_return <- calculate_net_returns_borrowing(r_BAB, Beta_L, Beta_H, r_L_next, r_H_next, rf,
                                                    transaction_costs, original_borrowing_costs, upper)
  if (high_bc_return > 0) {
    return(Inf) 
  }
  
  iter <- 0
  while ((upper - lower) > tolerance && iter < max_iter) {
    mid <- (upper + lower) / 2
    mid_return <- calculate_net_returns_borrowing(r_BAB, Beta_L, Beta_H, r_L_next, r_H_next, rf,
                                                  transaction_costs, original_borrowing_costs, mid)
    
    if (abs(mid_return) < tolerance) {
      return(mid * base_borrowing_rate) 
    } else if (mid_return > 0) {
      lower <- mid
    } else {
      upper <- mid
    }
    
    iter <- iter + 1
  }
  
  return(((upper + lower) / 2) * base_borrowing_rate)  
}

plot_borrowing_sensitivity <- function(r_BAB, r_BAB_sd, Beta_L, Beta_H, Beta_L_SD, Beta_H_SD,
                                       r_L_next, r_H_next, r_L_next_sd, r_H_next_sd, rf,
                                       transaction_costs_full, borrowing_costs_full,
                                       transaction_costs_full_sd, borrowing_costs_full_sd,
                                       base_borrowing_rate = 0.02) {
  bc_multipliers <- seq(0, 15, by = 0.1)
  
  returns_original <- sapply(bc_multipliers, function(mult) {
    return <- calculate_net_returns_borrowing(r_BAB, Beta_L, Beta_H, r_L_next, r_H_next, rf,
                                              transaction_costs_full, borrowing_costs_full, mult)
    return(return * 100)  
  })
  
  returns_sd <- sapply(bc_multipliers, function(mult) {
    return <- calculate_net_returns_borrowing(r_BAB_sd, Beta_L_SD, Beta_H_SD, r_L_next_sd, r_H_next_sd, rf,
                                              transaction_costs_full_sd, borrowing_costs_full_sd, mult)
    return(return * 100)  
  })
  
  crosses_zero_orig <- which(diff(sign(returns_original)) != 0)
  if (length(crosses_zero_orig) > 0) {
    idx <- crosses_zero_orig[1]
    x1 <- bc_multipliers[idx]
    x2 <- bc_multipliers[idx + 1]
    y1 <- returns_original[idx]
    y2 <- returns_original[idx + 1]
    breakeven_orig <- x1 + (0 - y1) * (x2 - x1) / (y2 - y1)
  } else {
    breakeven_orig <- max(bc_multipliers)
  }
  
  crosses_zero_sd <- which(diff(sign(returns_sd)) != 0)
  if (length(crosses_zero_sd) > 0) {
    idx <- crosses_zero_sd[1]
    x1 <- bc_multipliers[idx]
    x2 <- bc_multipliers[idx + 1]
    y1 <- returns_sd[idx]
    y2 <- returns_sd[idx + 1]
    breakeven_sd <- x1 + (0 - y1) * (x2 - x1) / (y2 - y1)
  } else {
    breakeven_sd <- max(bc_multipliers)
  }
  
  min_breakeven <- min(breakeven_sd, breakeven_orig) * 0.8
  max_breakeven <- max(breakeven_sd, breakeven_orig) * 1.2
  
  plot(bc_multipliers * base_borrowing_rate * 10000, returns_original, type = "l", col = "blue",
       xlab = "Årlige låneomkostninger (basispoint)", ylab = "Årligt afkast (%)",
       main = "Følsomhedstest af BAB-afkast over for låneomkostninger",
       xlim = c(0, max_breakeven * base_borrowing_rate * 10000),
       ylim = range(c(returns_original[1:which.min(abs(bc_multipliers - max_breakeven))], 
                      returns_sd[1:which.min(abs(bc_multipliers - max_breakeven))]), na.rm = TRUE))
  
  lines(bc_multipliers * base_borrowing_rate * 10000, returns_sd, col = "green")
  
  abline(h = 0, lty = 2)
  abline(v = breakeven_orig * base_borrowing_rate * 10000, col = "blue", lty = 2)
  abline(v = breakeven_sd * base_borrowing_rate * 10000, col = "green", lty = 2)
  
  legend("topright", legend = c("Original BAB", "SD-forbedret BAB"),
         col = c("blue", "green"), lty = 1, cex = 0.8)
  
  return(list(
    breakeven_orig = breakeven_orig * base_borrowing_rate,
    breakeven_sd = breakeven_sd * base_borrowing_rate
  ))
}

cat("\n----- ANALYSE AF BREAK-EVEN LÅNEOMKOSTNINGER -----\n")

breakeven_bc_original <- find_breakeven_borrowing(
  r_BAB, Beta_L, Beta_H, r_L_next, r_H_next, rf,
  transaction_costs_full, borrowing_costs_full,
  base_borrowing_rate = 0.02
)

breakeven_bc_sd <- find_breakeven_borrowing(
  r_BAB_sd, Beta_L_SD, Beta_H_SD, r_L_next_sd, r_H_next_sd, rf,
  transaction_costs_full_sd, borrowing_costs_full_sd,
  base_borrowing_rate = 0.02
)

cat("Break-even låneomkostning for Original BAB:", round(breakeven_bc_original * 10000, 2), "basispoint årligt\n")
cat("Break-even låneomkostning for SD-forbedret BAB:", round(breakeven_bc_sd * 10000, 2), "basispoint årligt\n")

par(mfrow = c(1, 1))
breakeven_results <- plot_borrowing_sensitivity(
  r_BAB, r_BAB_sd, Beta_L, Beta_H, Beta_L_SD, Beta_H_SD,
  r_L_next, r_H_next, r_L_next_sd, r_H_next_sd, rf,
  transaction_costs_full, borrowing_costs_full,
  transaction_costs_full_sd, borrowing_costs_full_sd
)

if (is.finite(breakeven_bc_original) && is.finite(breakeven_bc_sd)) {
  improvement_ratio <- breakeven_bc_sd / breakeven_bc_original
  if (improvement_ratio > 1) {
    cat("\nForbedring i break-even låneomkostninger: ", 
        round((improvement_ratio - 1) * 100, 2), "% højere for SD-forbedret BAB\n")
  } else {
    cat("\nOriginal BAB kan modstå ", 
        round((1/improvement_ratio - 1) * 100, 2), "% højere låneomkostninger end SD-forbedret BAB\n")
  }
} else {
  cat("\nKan ikke beregne forbedringsforhold: en eller begge strategier har ubegrænsede break-even låneomkostninger\n")
}



