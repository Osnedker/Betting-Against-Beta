# Indlæs nødvendige pakker
library(dplyr)
library(zoo)
library(matrixStats)
library(lubridate)
library(writexl)
library(lmtest)
library(sandwich)

# Davidson stokastisk dominans test funktion
davidson_test <- function(A, B, Trim, order) {
  # H null: A dominerer ikke B
  # Trim: % af hale som trimmes, f.eks., 5 for 5%
  # order: orden af stokastisk dominans (heltal >= 2)
  
  v <- sort(unique(c(A, B)))
  
  # Trim øvre og nedre grænse af data
  T_val <- ceiling(length(v) * Trim / 100)
  if (length(v) <= 2*T_val) {
    return(list(pvalue = NA, tvalue = NA))  # NA hvis der ikke er nok data
  }
  Z <- v[(T_val + 1):(length(v) - T_val)]
  
  ZZ <- matrix(rep(Z, each = length(A)), nrow = length(A))
  Da <- pmax(ZZ - matrix(rep(A, length(Z)), ncol = length(Z)), 0)^(order - 1) / (order - 1)
  Db <- pmax(ZZ - matrix(rep(B, length(Z)), ncol = length(Z)), 0)^(order - 1) / (order - 1)
  DA <- colMeans(Da)
  DB <- colMeans(Db)
  
  # Udregn t-værdier
  t_vals <- sqrt(length(A)) * (DB - DA) / apply(Db - Da, 2, sd)
  
  tvalue <- min(t_vals, na.rm = TRUE)
  
  # For to-sidet test har vi brug for double p-værdi eller bruge en anden udregning
  # Bruger to-sidet p-værdi
  pvalue <- 2 * pnorm(tvalue)  # To-sidet test p-værdi
  
  return(list(pvalue = pvalue, tvalue = tvalue))
}

# Indlæs data
Marketindlæs <- read.csv2(file.choose())
FFM3indlæs <- read.csv2(file.choose())

Market <- data.frame(lapply(Marketindlæs, function(x) as.numeric(as.character(x))))
FFM3 <- data.frame(lapply(FFM3indlæs, function(x) as.numeric(as.character(x))))

Market[Market == -99.99] <- NA

Date <- Market[1]

# Tager log 
Market <- log(Market[-1]/100+1)

# Tager log 
FFM3 <- log(FFM3/100+1)

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

# Rullende 3 dags afkast
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

# Udregn Vasicek-shrinkage vægte w_i
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

Beta_shrink$Date <- as.Date(as.character(Beta_shrink$Date), format="%Y%m%d")

last_days <- Beta_shrink %>%
  group_by(Year = year(Date), Month = month(Date)) %>%
  filter(Date == max(Date)) %>%
  ungroup()

z <- apply(last_days[, -c(1, (ncol(last_days)-1), ncol(last_days))], 1, function(x) {
  rank(x, ties.method = "average", na.last = NA)
})

z_bar <- sapply(z, function(row) mean(row, na.rm = TRUE))

z_diff <- mapply(function(z_values, z_mean) {
  if (is.null(z_values)) return(NULL)  # Håndter tomme lister
  return(z_values - z_mean)  # Beregn forskellen
}, z, z_bar, SIMPLIFY = FALSE)  # Returner som liste

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

# ------------- Stokastisk dominans implementation -------------
# Funktion der laver stokastisk dominans test for hver måned 
perform_monthly_sd_tests <- function(Market, sd_order = 3, trim_pct = 5) {
  # Lav en data frame til resultaterne 
  dates <- unique(format(Market$Date, "%Y-%m"))
  
  # Igangsæt matricer til at gemme t-værdier 
  # En for markeds dominerede aktiver, en for aktiv dominerende aktiver
  market_dominates <- matrix(NA, nrow = length(dates), 
                             ncol = ncol(Market)-2)  # Fjerner Date og YearMonth
  asset_dominates <- matrix(NA, nrow = length(dates), 
                            ncol = ncol(Market)-2)
  
  colnames(market_dominates) <- colnames(Market)[2:(ncol(Market)-1)]
  colnames(asset_dominates) <- colnames(Market)[2:(ncol(Market)-1)]
  
  # For hver måned 
  for (i in 1:length(dates)) {
    # Få månedsdata
    month_data <- Market[format(Market$Date, "%Y-%m") == dates[i], ]
    if (nrow(month_data) < 15) next  # Spring måneder med ikke nok data over 
    
    market_returns <- month_data$MarketReturn
    
    # For hvert aktiv 
    for (j in 2:(ncol(Market)-1)) {
      asset_name <- colnames(Market)[j]
      asset_returns <- month_data[[asset_name]]
      
      # Spring over hvis der er for mange NA'er 
      valid_indices <- !is.na(asset_returns) & !is.na(market_returns)
      if (sum(valid_indices) < 15) next
      
      asset_returns <- asset_returns[valid_indices]
      market_returns_filtered <- market_returns[valid_indices]
      
      # Test om markedet dominerer aktivet (3. orden SD)
      md_test <- davidson_test(market_returns_filtered, asset_returns, trim_pct, sd_order)
      market_dominates[i, j-1] <- md_test$tvalue
      
      # Test om aktivet dominerer markedet (3. orden SD)
      ad_test <- davidson_test(asset_returns, market_returns_filtered, trim_pct, sd_order)
      asset_dominates[i, j-1] <- ad_test$tvalue
    }
    
    # Udskriv hvor langt koden er kommet 
    if (i %% 10 == 0) {
      cat("Gennemgået", i, "af", length(dates), "måneder for SD test\n")
    }
  }
  
  return(list(
    dates = dates,
    market_dominates = market_dominates,
    asset_dominates = asset_dominates
  ))
}

identify_dominated_assets <- function(sd_results, threshold = -2.58) {  
  dates <- sd_results$dates
  market_dominates <- sd_results$market_dominates
  asset_dominates <- sd_results$asset_dominates
  
  # Igangsæt lister til at gemme resultater 
  dominated_by_market <- list()
  dominating_market <- list()
  
  # For hver måned 
  for (i in 1:(length(dates)-1)) {  # Tjek hver måned pånær den sidste
    # Igangsæt vektorer for nuværende måned 
    dom_by_market <- rep(FALSE, ncol(market_dominates))
    dom_market <- rep(FALSE, ncol(asset_dominates))
    
    # Tjek om markedet dominerer eller aktivet dominerer for den næste måned 
    for (j in 1:ncol(market_dominates)) {
      # Tjek om markedet dominerer aktivet i nuværende måned 
      if (!is.na(market_dominates[i, j])) {
        if (market_dominates[i, j] < threshold) {
          dom_by_market[j] <- TRUE
        }
      }
      
      # Tjek om aktivet dominerer markedet i nuværende måned 
      if (!is.na(asset_dominates[i, j])) {
        if (asset_dominates[i, j] < threshold) {
          dom_market[j] <- TRUE
        }
      }
    }
    
    # Gem resultater for næste måned 
    dominated_by_market[[dates[i + 1]]] <- colnames(market_dominates)[dom_by_market]
    dominating_market[[dates[i + 1]]] <- colnames(asset_dominates)[dom_market]
  }
  
  return(list(
    dominated_by_market = dominated_by_market,
    dominating_market = dominating_market
  ))
}

# Tilpas BAB portefølje vægte baseret på stokastisk dominans
modify_weights_with_sd <- function(weights_H, weights_L, dominated_assets, month) {
  if (!(month %in% names(dominated_assets$dominated_by_market)) || 
      !(month %in% names(dominated_assets$dominating_market))) {
    return(list(w_H = weights_H, w_L = weights_L))
  }
  
  # Lad aktivet blive udelukket for nuværende måned
  dominated <- dominated_assets$dominated_by_market[[month]]
  dominating <- dominated_assets$dominating_market[[month]]
  
  # Laver kopier af vægtene 
  new_weights_H <- weights_H
  new_weights_L <- weights_L
  
  # Fjern de dominerede aktiver fra det lange ben (lav beta)
  for (asset in dominated) {
    if (asset %in% names(new_weights_L)) {
      new_weights_L[asset] <- 0
    }
  }
  
  # Fjern de dominerede aktiver fra det korte bet (høj beta)
  for (asset in dominating) {
    if (asset %in% names(new_weights_H)) {
      new_weights_H[asset] <- 0
    }
  }
  
  # Sæt NA værdier til 0 
  new_weights_H[is.na(new_weights_H)] <- 0
  new_weights_L[is.na(new_weights_L)] <- 0
  
  # Renormaliser vægtene hvis der fjernes nogle aktiver 
  if (sum(new_weights_L, na.rm = TRUE) > 0) {
    new_weights_L <- new_weights_L / sum(new_weights_L, na.rm = TRUE)
  }
  if (sum(new_weights_H, na.rm = TRUE) > 0) {
    new_weights_H <- new_weights_H / sum(new_weights_H, na.rm = TRUE)
  }
  
  return(list(w_H = new_weights_H, w_L = new_weights_L))
}

# Udfør SD test
cat("Udfører stokastisk dominans test...\n")
sd_results <- perform_monthly_sd_tests(Market, sd_order = 3, trim_pct = 5)

# Identificer dominerede aktiver
cat("Indentificerer dominerede aktiver...\n")
dominated_assets <- identify_dominated_assets(sd_results, threshold = -2.58)

# Igangsæt lister for modificerede vægte
weights_lookup_H_sd <- weights_lookup_H
weights_lookup_L_sd <- weights_lookup_L

# Modificer vægtene baseret på stokastisk dominans resultaterne 
cat("Modificerer BAB vægte baseret på stokastisk dominans...\n")
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


# ------------- Tilbage til original kode -------------
# Beregn r_H_next baseret på YYYYMM matchning
r_H_next <- sapply(1:(nrow(Market)-1), function(i) {
  current_month <- Market$YearMonth[i]  # Find YYYY-MM fra Market
  
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
  
  if (!(current_month %in% names(weights_lookup_L))) {
    print(paste("Fejl! Måneden", current_month, "findes ikke i weights_lookup"))
    return(NA)
  }
  
  # Beregn vægtet afkast for næste dag
  return(sum(weights_lookup_L[[current_month]] * 
               Market[i+1, names(weights_lookup_L[[current_month]]), drop = FALSE], 
             na.rm = TRUE))
})

# Udregn r_H_next og r_L_next med SD vægte
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

Beta_H <- sapply(1:(nrow(Beta_shrink)-1), function(i) {
  current_month <- Beta_shrink$YearMonth[i]  # Find YYYY-MM fra Market
  
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
  
  if (!(current_month %in% names(weights_lookup_L))) {
    print(paste("Fejl! Måneden", current_month, "findes ikke i weights_lookup"))
    return(NA)
  }
  
  # Beregn vægtet afkast for næste dag
  return(sum(weights_lookup_L[[current_month]] * 
               Beta_shrink[i, names(weights_lookup_L[[current_month]]), drop = FALSE], 
             na.rm = TRUE))
})

# Udregn Beta_H_SD med ændrede vægte, ellers bruge originale 
Beta_H_SD <- sapply(1:(nrow(Beta_shrink) - 1), function(i) {
  current_month <- Beta_shrink$YearMonth[i]  # Find YYYY-MM fra Beta
  
  # Tjek om vi har SD-justerede vægte for nuværende måned 
  if (!(current_month %in% names(weights_lookup_H_sd))) {
    return(NA)
    
  } else {
    return(sum(weights_lookup_H_sd[[current_month]] * 
                 Beta_shrink[i, names(weights_lookup_H_sd[[current_month]]), drop = FALSE], 
               na.rm = TRUE))
  }
})

# Udregn Beta_L_SD med ændrede vægte, ellers bruge originale 
Beta_L_SD <- sapply(1:(nrow(Beta_shrink) - 1), function(i) {
  current_month <- Beta_shrink$YearMonth[i]  # Find YYYY-MM fra Beta
  
  # Tjek om vi har SD-justerede vægte for nuværende måned 
  if (!(current_month %in% names(weights_lookup_L_sd))) {
    return(NA)
    
  } else {
    return(sum(weights_lookup_L_sd[[current_month]] * 
                 Beta_shrink[i, names(weights_lookup_L_sd[[current_month]]), drop = FALSE], 
               na.rm = TRUE))
  }
})

#------------------- Udregner BAB-faktor ---------------------------
# Træk den risikofrie rente ud
rf <- FFM3[,5]

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

r_BAB_sd[is.infinite(r_BAB_sd)] <- NA # Givet at der er INF

# ------------------------- Transaktionsomkostninger ---------------------------------
transaction_cost_rate <- 0.001        # 10 bps
borrowing_cost_rate <- 0.02           # 2% årligt

# Funktion til udregning af overskud
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

# Igangsæt lister til resultater 
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
  
  # Udregn overskud 
  t_H <- calculate_turnover(weights_lookup_H[[month_prev]], weights_lookup_H[[month_curr]], returns_next)
  t_L <- calculate_turnover(weights_lookup_L[[month_prev]], weights_lookup_L[[month_curr]], returns_next)
  if (is.na(t_H) || is.na(t_L)) next
  
  # Udregn omkostninger 
  transaction_cost <- (t_H + t_L) * transaction_cost_rate
  borrowing_cost <- sum(weights_lookup_H[[month_curr]], na.rm = TRUE) * (borrowing_cost_rate / 12)
  
  # Gem under forskellige indeks
  rebalance_day_index <- which(Market$YearMonth == month_curr)[1]
  transaction_costs_full[rebalance_day_index] <- transaction_cost
  borrowing_costs_full[rebalance_day_index] <- borrowing_cost
}  

# Igangsæt SD omkostnings vektorer 
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
  
  # Udregn overskud for SD porteføljer 
  t_H_sd <- calculate_turnover(weights_lookup_H_sd[[month_prev]], weights_lookup_H_sd[[month_curr]], returns_next)
  t_L_sd <- calculate_turnover(weights_lookup_L_sd[[month_prev]], weights_lookup_L_sd[[month_curr]], returns_next)
  if (is.na(t_H_sd) || is.na(t_L_sd)) next
  
  # Udregn omkostninger for SD porteføljer
  transaction_cost_sd <- (t_H_sd + t_L_sd) * transaction_cost_rate
  borrowing_cost_sd <- sum(weights_lookup_H_sd[[month_curr]], na.rm = TRUE) * (borrowing_cost_rate / 12)
  
  # Gem resultater under de rigtige indeks
  rebalance_day_index_sd <- which(Market$YearMonth == month_curr)[1]
  transaction_costs_full_sd[rebalance_day_index_sd] <- transaction_cost_sd
  borrowing_costs_full_sd[rebalance_day_index_sd] <- borrowing_cost_sd
  
  valid_months_sd <- c(valid_months_sd, month_curr)
  month_lookup_sd <- c(month_lookup_sd, rebalance_day_index_sd)
}


# Beregn BAB-faktoren (original) fratrukket transaction cost
r_BAB_net <- (1 / Beta_L) * (r_L_next - rf) - (1 / Beta_H) * (r_H_next - rf) - transaction_costs_full - borrowing_costs_full

# Beregn SD-enhanced BAB-faktoren fratrukket transaction cost
r_BAB_sd_net <- (1 / Beta_L_SD) * (r_L_next_sd - rf) - (1 / Beta_H_SD) * (r_H_next_sd - rf)- transaction_costs_full_sd - borrowing_costs_full_sd


# Tjek begge BAB faktorer
par(mfrow = c(2, 1))
plot(r_BAB, type = "l", main = "Original BAB Factor", col = "blue")
plot(r_BAB_sd, type = "l", main = "SD-Enhanced BAB Factor", col = "green")
par(mfrow = c(1, 1))

cat("\nPerformance sammenligning:\n")
cat("Original BAB månedligseret afkast:", round((((1+mean(r_BAB, na.rm = TRUE))^(250/12))-1)*100, 10), "%\n")
cat("SD-Enhanced BAB månedligseret afkast:", round((((1+mean(r_BAB_sd, na.rm = TRUE))^(250/12))-1)*100, 10), "%\n")

# Transaktionsomkostninger
cat("Original BAB månedligseret afkast minus transaktionsomkostninger:", round((((1+mean(r_BAB_net, na.rm = TRUE))^(250/12))-1)*100, 10), "%\n")
cat("SD-Enhanced BAB månedligseret afkast minus transaktionsomkostninger:", round((((1+mean(r_BAB_sd_net, na.rm = TRUE))^(250/12))-1)*100, 10), "%\n")


# ------------- Fortsæt med original Alpha og Beta beregninger -------------
Market_excess_return <- FFM3[,2]
Market_excess_return <- Market_excess_return[4:(min_length+3)]

# Lav CAPM regression for original BAB
capm_model <- lm(r_BAB ~ Market_excess_return)

capm_alpha <- coef(capm_model)[1]*100*250/12
print(paste("CAPM Alpha (Original BAB):", round(capm_alpha, 2), "%"))

# Udtræk realiseret beta (hældningskoefficienten)
Beta_BAB_real <- coef(capm_model)[2]  # Beta er den anden koefficient
print(paste("Realiseret Beta for BAB:", round(Beta_BAB_real, 4)))

# Lav CAPM regression for SDBAB
capm_model_sd <- lm(r_BAB_sd ~ Market_excess_return)

capm_alpha_sd <- coef(capm_model_sd)[1]*100*250/12
print(paste("CAPM Alpha (SD-Enhanced BAB):", round(capm_alpha_sd, 2), "%"))

# Udtræk realiseret beta (hældningskoefficienten)
Beta_SD_real <- coef(capm_model_sd)[2]  # Beta er den anden koefficient
print(paste("Realiseret Beta for SD:", round(Beta_SD_real, 4)))

#------------------------ Transaktionsomkostninger ------------------------------------

# Lav CAPM regression for BAB med transaktionsomkostninger
capm_model_trans <- lm(r_BAB_net ~ Market_excess_return)

capm_alpha_trans <- coef(capm_model_trans)[1]*100*250/12
print(paste("CAPM Alpha (Original BAB):", round(capm_alpha_trans, 2), "%"))

# Udtræk realiseret beta (hældningskoefficienten) med trans
Beta_BAB_real_trans <- coef(capm_model_trans)[2]  # Beta er den anden koefficient
print(paste("Realiseret Beta for BAB med trans:", round(Beta_BAB_real_trans, 4)))

# Lav CAPM regression for SDBAB med transaktionsomkostninger 
capm_model_sd_trans <- lm(r_BAB_sd_net ~ Market_excess_return)

capm_alpha_sd_trans <- coef(capm_model_sd_trans)[1]*100*250/12
print(paste("CAPM Alpha (SD-Enhanced BAB):", round(capm_alpha_sd_trans, 2), "%"))

# Udtræk realiseret beta (hældningskoefficienten)
Beta_SD_real_trans <- coef(capm_model_sd_trans)[2]  # Beta er den anden koefficient
print(paste("Realiseret Beta for BAB SD med trans:", round(Beta_SD_real_trans, 4)))

# Træk SMB og HML ud fra Fama-French datasæt 
SMB <- FFM3[, 3]
HML <- FFM3[, 4]

# Tilpas længder med BAB faktor
SMB <- SMB[3:(min_length+2)]
HML <- HML[3:(min_length+2)]

# Lav 3-faktor regression for original BAB 
ff3_model <- lm(r_BAB ~ Market_excess_return + SMB + HML)

ff3_alpha <- coef(ff3_model)[1]*100*250/12
print(paste("Fama-French 3-Factor Alpha (Original BAB):", round(ff3_alpha, 2), "%"))

# Lav 3-faktor regression for SDBAB 
ff3_model_sd <- lm(r_BAB_sd ~ Market_excess_return + SMB + HML)

ff3_alpha_sd <- coef(ff3_model_sd)[1]*100*250/12
print(paste("Fama-French 3-Factor Alpha (SD-Enhanced BAB):", round(ff3_alpha_sd, 2), "%"))

#---------------------------- Transaktionsomkostninger ----------------------------------

# Lav 3-faktor regression for BAB med transaktionsomkostninger 
ff3_model_trans <- lm(r_BAB_net ~ Market_excess_return + SMB + HML)

ff3_alpha_trans <- coef(ff3_model_trans)[1]*100*250/12
print(paste("Fama-French 3-Factor Alpha (Original BAB) med trans:", round(ff3_alpha_trans, 2), "%"))

# Lav 3-faktor regression for SDBAB med transaktionsomkostninger 
ff3_model_sd_trans <- lm(r_BAB_sd_net ~ Market_excess_return + SMB + HML)

ff3_alpha_sd_trans <- coef(ff3_model_sd_trans)[1]*100*250/12
print(paste("Fama-French 3-Factor Alpha (SD-Enhanced BAB) med trans:", round(ff3_alpha_sd_trans, 2), "%"))


#------------------------------- Volatilitet og Sharpe ratio ------------------------------------------
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

# Beregn gennemsnitligt dagligt afkast og daglig risikofri rente
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

plot(cumulative_BAB, type="l", col="blue", xlab="Time", ylab="Cumulative Return", ylim = c(0,6.5),
     main="Cumulative Returns: Original vs SD-Enhanced BAB vs Original med trans vs SD med trans")
lines(cumulative_BAB_sd, col="green")
lines(cumulative_BAB_trans, col="yellow")
lines(cumulative_BAB_sd_trans, col="black")
legend("topleft", legend=c("Original BAB", "SD-Enhanced BAB", "Original med trans", "SD med trans"), 
       col=c("blue", "green", "yellow", "black"), lty=1)

# Tilføj t-tests for afkast og alpha'er
# For original BAB
# Udregn SE for afkast (dagligt) 
n_obs_BAB <- sum(!is.na(r_BAB))
se_return_BAB <- sd(r_BAB, na.rm = TRUE) / sqrt(n_obs_BAB)
t_stat_return_BAB <- mean(r_BAB, na.rm = TRUE) / se_return_BAB

# Udregn SE og t-stat for CAPM alpha 
capm_summary <- summary(capm_model)
t_stat_capm_alpha_BAB <- capm_summary$coefficients[1, 3]  # t-stat for intercept

# Udregn SE og t-stat for FF3 alpha 
ff3_summary <- summary(ff3_model)
t_stat_ff3_alpha_BAB <- ff3_summary$coefficients[1, 3]  # t-stat for intercept

# For SD BAB
n_obs_BAB_sd <- sum(!is.na(r_BAB_sd))
se_return_BAB_sd <- sd(r_BAB_sd, na.rm = TRUE) / sqrt(n_obs_BAB_sd)
t_stat_return_BAB_sd <- mean(r_BAB_sd, na.rm = TRUE) / se_return_BAB_sd

# Udregn SE og t-stat for CAPM alpha (SDBAB)
capm_summary_sd <- summary(capm_model_sd)
t_stat_capm_alpha_BAB_sd <- capm_summary_sd$coefficients[1, 3]  # t-stat for intercept

# Udregn SE og t-stat for FF3 alpha (SDBAB)
ff3_summary_sd <- summary(ff3_model_sd)
t_stat_ff3_alpha_BAB_sd <- ff3_summary_sd$coefficients[1, 3]  # t-stat for intercept

#------------------------------ Transaktionsomkostninger -------------------------------------
n_obs_BAB_trans <- sum(!is.na(r_BAB_net))
se_return_BAB_trans <- sd(r_BAB_net, na.rm = TRUE) / sqrt(n_obs_BAB_trans)
t_stat_return_BAB_trans <- mean(r_BAB_net, na.rm = TRUE) / se_return_BAB_trans

# Udregn SE og t-stat for CAPM alpha (BAB m. transaktionsomkostninger)
capm_summary_trans <- summary(capm_model_trans)
t_stat_capm_alpha_BAB_trans <- capm_summary_trans$coefficients[1, 3]  # t-stat for intercept

# Udregn SE og t-stat for FF3 alpha (BAB m. transaktionsomkostninger)
ff3_summary_trans <- summary(ff3_model_trans)
t_stat_ff3_alpha_BAB_trans <- ff3_summary_trans$coefficients[1, 3]  # t-stat for intercept

# For SD BAB
n_obs_BAB_sd_trans <- sum(!is.na(r_BAB_sd_net))
se_return_BAB_sd_trans <- sd(r_BAB_sd_net, na.rm = TRUE) / sqrt(n_obs_BAB_sd_trans)
t_stat_return_BAB_sd_trans <- mean(r_BAB_sd_net, na.rm = TRUE) / se_return_BAB_sd_trans

# Udregn SE og t-stat for CAPM alpha (SDBAB m. transaktionsomkostninger)
capm_summary_sd_trans <- summary(capm_model_sd_trans)
t_stat_capm_alpha_BAB_sd_trans <- capm_summary_sd_trans$coefficients[1, 3]  # t-stat for intercept

# Udregn SE og t-stat for FF3 alpha (SDBAB m. transaktionsomkostninger)
ff3_summary_sd_trans <- summary(ff3_model_sd_trans)
t_stat_ff3_alpha_BAB_sd_trans <- ff3_summary_sd_trans$coefficients[1, 3]  # t-stat for intercept

# Funktion der tilføjer signifikansniveauer 
add_significance <- function(t_stat) {
  if (abs(t_stat) > 2.58) return(" ***")  # 1% niveau
  else if (abs(t_stat) > 1.96) return(" **")  # 5% niveau
  else if (abs(t_stat) > 1.65) return(" *")  # 10% niveau
  else return("")
}

# Udregn månedlig afkast 
monthly_BAB <- (((1+mean(r_BAB, na.rm = TRUE))^(250/12))-1)*100
monthly_BAB_sd <- (((1+mean(r_BAB_sd, na.rm = TRUE))^(250/12))-1)*100
monthly_BAB_trans <- (((1+mean(r_BAB_net, na.rm = TRUE))^(250/12))-1)*100
monthly_BAB_sd_trans <- (((1+mean(r_BAB_sd_net, na.rm = TRUE))^(250/12))-1)*100

# Saml resultater i tabel med t-stat og signifikansniveauer 
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

# Udskriv tabel
print(summary_table_enhanced)

#--------------------------- Bear / Bull marked ---------------------
Market_Return <- FFM3[,2] + FFM3[,5]

# Lav en markedspris indeks af afkast 
market_price_index <- cumprod(1 + Market_Return)

# Igangsæt variable til at holde øje med stigninger og fald 
peak <- market_price_index[1]
peak_idx <- 1
trough <- market_price_index[1]
trough_idx <- 1
in_bear <- FALSE
threshold <- 0.2

# Lav en markedsbetingelsesvektor
market_condition <- rep(NA, length(market_price_index))

# Identificer bull og bear markeder 
for (i in 1:length(market_price_index)) {
  if (!in_bear) {  # I bull marked
    if (market_price_index[i] > peak) {
      # Ny stigning 
      peak <- market_price_index[i]
      peak_idx <- i
    } else if (market_price_index[i] < peak * (1 - threshold)) {
      # Gået ind i bear marked (faldet 20% fra toppen)
      in_bear <- TRUE
      market_condition[peak_idx:i] <- "Bear"
      trough <- market_price_index[i]
      trough_idx <- i
    }
  } else {  # I bear marked
    if (market_price_index[i] < trough) {
      # Ny tilbage gang
      trough <- market_price_index[i]
      trough_idx <- i
    } else if (market_price_index[i] > trough * (1 + threshold)) {
      # Gået ind i bull marked (20% stigning fra fald)
      in_bear <- FALSE
      market_condition[trough_idx:i] <- "Bull"
      peak <- market_price_index[i]
      peak_idx <- i
    }
  }
  
  if (is.na(market_condition[i])) {
    market_condition[i] <- ifelse(in_bear, "Bear", "Bull")
  }
}

# Sikre at markedsbetingelserne har samme længde som BAB 
# Antager at BAB afkast starter på samme tid som markedsafkast
market_condition_for_bab <- market_condition

# Deler BAB afkast ud ift. markedssituation
r_BAB_bull <- r_BAB[market_condition_for_bab == "Bull"]
r_BAB_bear <- r_BAB[market_condition_for_bab == "Bear"]
r_BAB_sd_bull <- r_BAB_sd[market_condition_for_bab == "Bull"]
r_BAB_sd_bear <- r_BAB_sd[market_condition_for_bab == "Bear"]

# Bull marked statstikker 
mean_BAB_bull <- mean(r_BAB_bull, na.rm = TRUE)
sd_BAB_bull <- sd(r_BAB_bull, na.rm = TRUE)
sharpe_BAB_bull <- mean_BAB_bull / sd_BAB_bull * sqrt(250)

mean_BAB_sd_bull <- mean(r_BAB_sd_bull, na.rm = TRUE)
sd_BAB_sd_bull <- sd(r_BAB_sd_bull, na.rm = TRUE)
sharpe_BAB_sd_bull <- mean_BAB_sd_bull / sd_BAB_sd_bull * sqrt(250)

# Bear marked statistikker 
mean_BAB_bear <- mean(r_BAB_bear, na.rm = TRUE)
sd_BAB_bear <- sd(r_BAB_bear, na.rm = TRUE)
sharpe_BAB_bear <- mean_BAB_bear / sd_BAB_bear * sqrt(250)

mean_BAB_sd_bear <- mean(r_BAB_sd_bear, na.rm = TRUE)
sd_BAB_sd_bear <- sd(r_BAB_sd_bear, na.rm = TRUE)
sharpe_BAB_sd_bear <- mean_BAB_sd_bear / sd_BAB_sd_bear * sqrt(250)

# Lav CAPM regressioner for hver markedsbetingelse 
# Bull marked regressioner 
market_bull <- Market_Return[market_condition_for_bab == "Bull"]
capm_BAB_bull <- lm(r_BAB_bull ~ market_bull)
capm_BAB_sd_bull <- lm(r_BAB_sd_bull ~ market_bull)

# Bear marked regressioner
market_bear <- Market_Return[market_condition_for_bab == "Bear"]
capm_BAB_bear <- lm(r_BAB_bear ~ market_bear)
capm_BAB_sd_bear <- lm(r_BAB_sd_bear ~ market_bear)

# Træk alpha'er ud (årlig)
alpha_BAB_bull <- coef(capm_BAB_bull)[1] * 250 * 100
alpha_BAB_sd_bull <- coef(capm_BAB_sd_bull)[1] * 250 * 100
alpha_BAB_bear <- coef(capm_BAB_bear)[1] * 250 * 100
alpha_BAB_sd_bear <- coef(capm_BAB_sd_bear)[1] * 250 * 100

# Opsumer resultater i tabel
bull_bear_summary <- data.frame(
  Market_Condition = c("Bull", "Bull", "Bear", "Bear"),
  Strategy = c("Original BAB", "SD-Enhanced BAB", "Original BAB", "SD-Enhanced BAB"),
  Annualized_Return = c(mean_BAB_bull * 250 * 100, mean_BAB_sd_bull * 250 * 100, 
                        mean_BAB_bear * 250 * 100, mean_BAB_sd_bear * 250 * 100),
  Annualized_Volatility = c(sd_BAB_bull * sqrt(250) * 100, sd_BAB_sd_bull * sqrt(250) * 100,
                            sd_BAB_bear * sqrt(250) * 100, sd_BAB_sd_bear * sqrt(250) * 100),
  Sharpe_Ratio = c(sharpe_BAB_bull, sharpe_BAB_sd_bull, sharpe_BAB_bear, sharpe_BAB_sd_bear),
  CAPM_Alpha = c(alpha_BAB_bull, alpha_BAB_sd_bull, alpha_BAB_bear, alpha_BAB_sd_bear)
)

# Udskriv resultater 
print("Performance of BAB Strategies in Bull and Bear Markets")
print(bull_bear_summary)

# Lav plots som sammenligner resultater 
# Udregn kumulative resultater
cum_BAB_bull <- cumprod(1 + na.omit(r_BAB_bull))
cum_BAB_sd_bull <- cumprod(1 + na.omit(r_BAB_sd_bull))
cum_BAB_bear <- cumprod(1 + na.omit(r_BAB_bear))
cum_BAB_sd_bear <- cumprod(1 + na.omit(r_BAB_sd_bear))

par(mfrow=c(2,1))
plot(cum_BAB_bull, type="l", col="blue", main="Bull Market Performance", 
     xlab="Time", ylab="Cumulative Return")
lines(cum_BAB_sd_bull, col="green")
legend("topleft", legend=c("Original BAB", "SD-Enhanced BAB"), 
       col=c("blue", "green"), lty=1)

plot(cum_BAB_bear, type="l", col="blue", main="Bear Market Performance", 
     xlab="Time", ylab="Cumulative Return")
lines(cum_BAB_sd_bear, col="green")
legend("topleft", legend=c("Original BAB", "SD-Enhanced BAB"), 
       col=c("blue", "green"), lty=1)
par(mfrow=c(1,1))

# Udskriv opsummerende statistikker
cat("Total observations:", length(market_condition), "\n")
cat("Bull market periods:", sum(market_condition == "Bull"), "days (", 
    round(sum(market_condition == "Bull")/length(market_condition)*100, 1), "%)\n")
cat("Bear market periods:", sum(market_condition == "Bear"), "days (", 
    round(sum(market_condition == "Bear")/length(market_condition)*100, 1), "%)\n")

bull_vs_bear_BAB <- t.test(r_BAB_bull, r_BAB_bear)
bull_vs_bear_BAB_sd <- t.test(r_BAB_sd_bull, r_BAB_sd_bear)

cat("\nStatistical Tests - Original BAB Bull vs. Bear:\n")
print(bull_vs_bear_BAB)

cat("\nStatistical Tests - SD-Enhanced BAB Bull vs. Bear:\n")
print(bull_vs_bear_BAB_sd)

# Udregning markedsafkast i bull og bear perioder 
market_bull <- Market_Return[market_condition == "Bull"]
market_bear <- Market_Return[market_condition == "Bear"]

# Udregn markeds statistikker 
mean_market_bull <- mean(market_bull, na.rm = TRUE)
sd_market_bull <- sd(market_bull, na.rm = TRUE)
sharpe_market_bull <- mean_market_bull / sd_market_bull * sqrt(250)

mean_market_bear <- mean(market_bear, na.rm = TRUE)
sd_market_bear <- sd(market_bear, na.rm = TRUE)
sharpe_market_bear <- mean_market_bear / sd_market_bear * sqrt(250)

# Tilføj marked til tabellen
market_summary <- data.frame(
  Market_Condition = c("Bull", "Bear"),
  Strategy = c("Market", "Market"),
  Annualized_Return = c(mean_market_bull * 250 * 100, mean_market_bear * 250 * 100),
  Annualized_Volatility = c(sd_market_bull * sqrt(250) * 100, sd_market_bear * sqrt(250) * 100),
  Sharpe_Ratio = c(sharpe_market_bull, sharpe_market_bear),
  CAPM_Alpha = c(NA, NA)  # Market alpha relative to itself is always 0
)

full_summary <- rbind(bull_bear_summary, market_summary)

bull_outperformance_original <- (mean_BAB_bull * 250) - (mean_market_bull * 250)
bull_outperformance_sd <- (mean_BAB_sd_bull * 250) - (mean_market_bull * 250)
bear_outperformance_original <- (mean_BAB_bear * 250) - (mean_market_bear * 250)
bear_outperformance_sd <- (mean_BAB_sd_bear * 250) - (mean_market_bear * 250)

# Udregn op/ned marked ratio
up_capture_original <- (mean_BAB_bull / mean_market_bull) * 100
up_capture_sd <- (mean_BAB_sd_bull / mean_market_bull) * 100
down_capture_original <- (mean_BAB_bear / mean_market_bear) * 100
down_capture_sd <- (mean_BAB_sd_bear / mean_market_bear) * 100

# Udskriv statistikker 
cat("\nRelative Performance to Market:\n")
cat("Bull Markets - BAB vs Market Outperformance (pp):", round(bull_outperformance_original * 100, 2), "%\n")
cat("Bull Markets - SD-BAB vs Market Outperformance (pp):", round(bull_outperformance_sd * 100, 2), "%\n")
cat("Bear Markets - BAB vs Market Outperformance (pp):", round(bear_outperformance_original * 100, 2), "%\n")
cat("Bear Markets - SD-BAB vs Market Outperformance (pp):", round(bear_outperformance_sd * 100, 2), "%\n")

cat("\nUp/Down Market Capture Ratios:\n")
cat("Original BAB - Up Market Capture:", round(up_capture_original, 2), "%\n")
cat("SD-Enhanced BAB - Up Market Capture:", round(up_capture_sd, 2), "%\n")
cat("Original BAB - Down Market Capture:", round(down_capture_original, 2), "%\n")
cat("SD-Enhanced BAB - Down Market Capture:", round(down_capture_sd, 2), "%\n")

# Korrelation af markedsregime 
cor_bull_original <- cor(r_BAB_bull, market_bull, use="pairwise.complete.obs")
cor_bull_sd <- cor(r_BAB_sd_bull, market_bull, use="pairwise.complete.obs")
cor_bear_original <- cor(r_BAB_bear, market_bear, use="pairwise.complete.obs")
cor_bear_sd <- cor(r_BAB_sd_bear, market_bear, use="pairwise.complete.obs")

cat("\nCorrelations with Market:\n")
cat("Bull Markets - Original BAB:", round(cor_bull_original, 3), "\n")
cat("Bull Markets - SD-Enhanced BAB:", round(cor_bull_sd, 3), "\n")
cat("Bear Markets - Original BAB:", round(cor_bear_original, 3), "\n")
cat("Bear Markets - SD-Enhanced BAB:", round(cor_bear_sd, 3), "\n")

par(mfrow=c(2,1))

# Bull markeds sammenligning
cum_market_bull <- cumprod(1 + na.omit(market_bull))
plot(cum_market_bull, type="l", col="red", 
     main="Bull Market: BAB vs Market", 
     xlab="Time", ylab="Cumulative Return")
lines(cum_BAB_bull, col="blue")
lines(cum_BAB_sd_bull, col="green")
legend("topleft", legend=c("Market", "Original BAB", "SD-Enhanced BAB"), 
       col=c("red", "blue", "green"), lty=1)

# Bear markeds sammenligning
cum_market_bear <- cumprod(1 + na.omit(market_bear))
plot(cum_market_bear, type="l", col="red", 
     main="Bear Market: BAB vs Market", 
     xlab="Time", ylab="Cumulative Return")
lines(cum_BAB_bear, col="blue")
lines(cum_BAB_sd_bear, col="green")
legend("topleft", legend=c("Market", "Original BAB", "SD-Enhanced BAB"), 
       col=c("red", "blue", "green"), lty=1)

par(mfrow=c(1,1))

Market_Return <- FFM3[,2] + FFM3[,5]

# Lav en markedspris ideks fra afkast 
market_price_index <- cumprod(1 + Market_Return)

# Igangsæt variabler der søger efter stigninger og fald 
peak <- market_price_index[1]
peak_idx <- 1
trough <- market_price_index[1]
trough_idx <- 1
in_bear <- FALSE
threshold <- 0.2

# Lav en markedsbetingelsesvektor 
market_condition <- rep(NA, length(market_price_index))

for (i in 1:length(market_price_index)) {
  if (!in_bear) {  # I bull marked
    if (market_price_index[i] > peak) {
      # Ny stigning 
      peak <- market_price_index[i]
      peak_idx <- i
    } else if (market_price_index[i] < peak * (1 - threshold)) {
      # Gået ind i bear marked (20% fald fra top) 
      in_bear <- TRUE
      market_condition[peak_idx:i] <- "Bear"
      trough <- market_price_index[i]
      trough_idx <- i
    }
  } else {  # I bear marked
    if (market_price_index[i] < trough) {
      # Nyt fald
      trough <- market_price_index[i]
      trough_idx <- i
    } else if (market_price_index[i] > trough * (1 + threshold)) {
      # Gået ind i bull marked (20% stigning fra bund)
      in_bear <- FALSE
      market_condition[trough_idx:i] <- "Bull"
      peak <- market_price_index[i]
      peak_idx <- i
    }
  }
  
  if (is.na(market_condition[i])) {
    market_condition[i] <- ifelse(in_bear, "Bear", "Bull")
  }
}

# Sikrer at markedsbetingelsen har samme længde som BAB
market_condition_for_bab <- market_condition

# Uddeler BAB afkast ift. markedet (Bull/bear)
r_BAB_bull <- r_BAB[market_condition_for_bab == "Bull"]
r_BAB_bear <- r_BAB[market_condition_for_bab == "Bear"]
r_BAB_sd_bull <- r_BAB_sd[market_condition_for_bab == "Bull"]
r_BAB_sd_bear <- r_BAB_sd[market_condition_for_bab == "Bear"]

# Udregn daglig gns. afkast og standard afvigelse
mean_daily_BAB_bull <- mean(r_BAB_bull, na.rm = TRUE)
sd_daily_BAB_bull <- sd(r_BAB_bull, na.rm = TRUE)
mean_daily_BAB_sd_bull <- mean(r_BAB_sd_bull, na.rm = TRUE)
sd_daily_BAB_sd_bull <- sd(r_BAB_sd_bull, na.rm = TRUE)
mean_daily_BAB_bear <- mean(r_BAB_bear, na.rm = TRUE)
sd_daily_BAB_bear <- sd(r_BAB_bear, na.rm = TRUE)
mean_daily_BAB_sd_bear <- mean(r_BAB_sd_bear, na.rm = TRUE)
sd_daily_BAB_sd_bear <- sd(r_BAB_sd_bear, na.rm = TRUE)

# Konverter daglig afkast til månedlig
mean_monthly_BAB_bull <- ((1 + mean_daily_BAB_bull)^21 - 1) * 100
sd_annual_BAB_bull <- sd_daily_BAB_bull * sqrt(250) * 100
sharpe_BAB_bull <- (mean_daily_BAB_bull * 250) / (sd_daily_BAB_bull * sqrt(250))

mean_monthly_BAB_sd_bull <- ((1 + mean_daily_BAB_sd_bull)^21 - 1) * 100
sd_annual_BAB_sd_bull <- sd_daily_BAB_sd_bull * sqrt(250) * 100
sharpe_BAB_sd_bull <- (mean_daily_BAB_sd_bull * 250) / (sd_daily_BAB_sd_bull * sqrt(250))

mean_monthly_BAB_bear <- ((1 + mean_daily_BAB_bear)^21 - 1) * 100
sd_annual_BAB_bear <- sd_daily_BAB_bear * sqrt(250) * 100
sharpe_BAB_bear <- (mean_daily_BAB_bear * 250) / (sd_daily_BAB_bear * sqrt(250))

mean_monthly_BAB_sd_bear <- ((1 + mean_daily_BAB_sd_bear)^21 - 1) * 100
sd_annual_BAB_sd_bear <- sd_daily_BAB_sd_bear * sqrt(250) * 100
sharpe_BAB_sd_bear <- (mean_daily_BAB_sd_bear * 250) / (sd_daily_BAB_sd_bear * sqrt(250))

# Lav CAPM regressioner for hver markedsstatus
market_bull <- Market_Return[market_condition_for_bab == "Bull"]
capm_BAB_bull <- lm(r_BAB_bull ~ market_bull)
capm_BAB_sd_bull <- lm(r_BAB_sd_bull ~ market_bull)

market_bear <- Market_Return[market_condition_for_bab == "Bear"]
capm_BAB_bear <- lm(r_BAB_bear ~ market_bear)
capm_BAB_sd_bear <- lm(r_BAB_sd_bear ~ market_bear)

# Træk månedlige alpha'er ud
alpha_monthly_BAB_bull <- coef(capm_BAB_bull)[1] * 21 * 100
alpha_monthly_BAB_sd_bull <- coef(capm_BAB_sd_bull)[1] * 21 * 100
alpha_monthly_BAB_bear <- coef(capm_BAB_bear)[1] * 21 * 100
alpha_monthly_BAB_sd_bear <- coef(capm_BAB_sd_bear)[1] * 21 * 100

# Lav en tabel med månedlig afkast og årlig volatilitet 
bull_bear_summary <- data.frame(
  Market_Condition = c("Bull", "Bull", "Bear", "Bear"),
  Strategy = c("Original BAB", "SD-Enhanced BAB", "Original BAB", "SD-Enhanced BAB"),
  Monthly_Return = c(mean_monthly_BAB_bull, mean_monthly_BAB_sd_bull, 
                     mean_monthly_BAB_bear, mean_monthly_BAB_sd_bear),
  Annual_Volatility = c(sd_annual_BAB_bull, sd_annual_BAB_sd_bull,
                        sd_annual_BAB_bear, sd_annual_BAB_sd_bear),
  Sharpe_Ratio = c(sharpe_BAB_bull, sharpe_BAB_sd_bull, sharpe_BAB_bear, sharpe_BAB_sd_bear),
  CAPM_Alpha_Monthly = c(alpha_monthly_BAB_bull, alpha_monthly_BAB_sd_bull, 
                         alpha_monthly_BAB_bear, alpha_monthly_BAB_sd_bear)
)

# Udskriv tabel
print("Performance of BAB Strategies in Bull and Bear Markets")
print(bull_bear_summary)

# Lav plots
# Udregn kumulative afkast
cum_BAB_bull <- cumprod(1 + na.omit(r_BAB_bull))
cum_BAB_sd_bull <- cumprod(1 + na.omit(r_BAB_sd_bull))
cum_BAB_bear <- cumprod(1 + na.omit(r_BAB_bear))
cum_BAB_sd_bear <- cumprod(1 + na.omit(r_BAB_sd_bear))

par(mfrow=c(2,1))
plot(cum_BAB_bull, type="l", col="blue", main="Bull Market Performance", 
     xlab="Time", ylab="Cumulative Return")
lines(cum_BAB_sd_bull, col="green")
legend("topleft", legend=c("Original BAB", "SD-Enhanced BAB"), 
       col=c("blue", "green"), lty=1)

plot(cum_BAB_bear, type="l", col="blue", main="Bear Market Performance", 
     xlab="Time", ylab="Cumulative Return")
lines(cum_BAB_sd_bear, col="green")
legend("topleft", legend=c("Original BAB", "SD-Enhanced BAB"), 
       col=c("blue", "green"), lty=1)
par(mfrow=c(1,1))

# Udskriv statistikker
cat("Total observations:", length(market_condition), "\n")
cat("Bull market periods:", sum(market_condition == "Bull"), "days (", 
    round(sum(market_condition == "Bull")/length(market_condition)*100, 1), "%)\n")
cat("Bear market periods:", sum(market_condition == "Bear"), "days (", 
    round(sum(market_condition == "Bear")/length(market_condition)*100, 1), "%)\n")

bull_vs_bear_BAB <- t.test(r_BAB_bull, r_BAB_bear)
bull_vs_bear_BAB_sd <- t.test(r_BAB_sd_bull, r_BAB_sd_bear)

cat("\nStatistical Tests - Original BAB Bull vs. Bear:\n")
print(bull_vs_bear_BAB)

cat("\nStatistical Tests - SD-Enhanced BAB Bull vs. Bear:\n")
print(bull_vs_bear_BAB_sd)

# Udregn markedsafkast i bull / bear perioder
market_bull <- Market_Return[market_condition == "Bull"]
market_bear <- Market_Return[market_condition == "Bear"]

# Udregn daglig markeds statistiskker 
mean_daily_market_bull <- mean(market_bull, na.rm = TRUE)
sd_daily_market_bull <- sd(market_bull, na.rm = TRUE)
mean_daily_market_bear <- mean(market_bear, na.rm = TRUE)
sd_daily_market_bear <- sd(market_bear, na.rm = TRUE)

# Konverter til månedlig afkast men beholder volatilitet som årlig
mean_monthly_market_bull <- ((1 + mean_daily_market_bull)^21 - 1) * 100
sd_annual_market_bull <- sd_daily_market_bull * sqrt(250) * 100
sharpe_market_bull <- (mean_daily_market_bull * 250) / (sd_daily_market_bull * sqrt(250))

mean_monthly_market_bear <- ((1 + mean_daily_market_bear)^21 - 1) * 100
sd_annual_market_bear <- sd_daily_market_bear * sqrt(250) * 100
sharpe_market_bear <- (mean_daily_market_bear * 250) / (sd_daily_market_bear * sqrt(250))

# Tilføj marked til tabel
market_summary <- data.frame(
  Market_Condition = c("Bull", "Bear"),
  Strategy = c("Market", "Market"),
  Monthly_Return = c(mean_monthly_market_bull, mean_monthly_market_bear),
  Annual_Volatility = c(sd_annual_market_bull, sd_annual_market_bear),
  Sharpe_Ratio = c(sharpe_market_bull, sharpe_market_bear),
  CAPM_Alpha_Monthly = c(NA, NA)  # Market alpha relative to itself is always 0
)

full_summary <- rbind(bull_bear_summary, market_summary)

bull_outperformance_original <- mean_monthly_BAB_bull - mean_monthly_market_bull
bull_outperformance_sd <- mean_monthly_BAB_sd_bull - mean_monthly_market_bull
bear_outperformance_original <- mean_monthly_BAB_bear - mean_monthly_market_bear
bear_outperformance_sd <- mean_monthly_BAB_sd_bear - mean_monthly_market_bear

# Udregn op/ned marked ratio
up_capture_original <- (mean_daily_BAB_bull / mean_daily_market_bull) * 100
up_capture_sd <- (mean_daily_BAB_sd_bull / mean_daily_market_bull) * 100
down_capture_original <- (mean_daily_BAB_bear / mean_daily_market_bear) * 100
down_capture_sd <- (mean_daily_BAB_sd_bear / mean_daily_market_bear) * 100

# Udskriv resultater
cat("\nRelative Performance to Market (Monthly):\n")
cat("Bull Markets - BAB vs Market Outperformance (pp):", round(bull_outperformance_original, 2), "%\n")
cat("Bull Markets - SD-BAB vs Market Outperformance (pp):", round(bull_outperformance_sd, 2), "%\n")
cat("Bear Markets - BAB vs Market Outperformance (pp):", round(bear_outperformance_original, 2), "%\n")
cat("Bear Markets - SD-BAB vs Market Outperformance (pp):", round(bear_outperformance_sd, 2), "%\n")

cat("\nUp/Down Market Capture Ratios:\n")
cat("Original BAB - Up Market Capture:", round(up_capture_original, 2), "%\n")
cat("SD-Enhanced BAB - Up Market Capture:", round(up_capture_sd, 2), "%\n")
cat("Original BAB - Down Market Capture:", round(down_capture_original, 2), "%\n")
cat("SD-Enhanced BAB - Down Market Capture:", round(down_capture_sd, 2), "%\n")

# Korrelationer af markedsregime
cor_bull_original <- cor(r_BAB_bull, market_bull, use="pairwise.complete.obs")
cor_bull_sd <- cor(r_BAB_sd_bull, market_bull, use="pairwise.complete.obs")
cor_bear_original <- cor(r_BAB_bear, market_bear, use="pairwise.complete.obs")
cor_bear_sd <- cor(r_BAB_sd_bear, market_bear, use="pairwise.complete.obs")

cat("\nCorrelations with Market:\n")
cat("Bull Markets - Original BAB:", round(cor_bull_original, 3), "\n")
cat("Bull Markets - SD-Enhanced BAB:", round(cor_bull_sd, 3), "\n")
cat("Bear Markets - Original BAB:", round(cor_bear_original, 3), "\n")
cat("Bear Markets - SD-Enhanced BAB:", round(cor_bear_sd, 3), "\n")

par(mfrow=c(2,1))

# Bull marked
cum_market_bull <- cumprod(1 + na.omit(market_bull))
plot(cum_market_bull, type="l", col="red", 
     main="Bull Market: BAB vs Market", 
     xlab="Time", ylab="Cumulative Return")
lines(cum_BAB_bull, col="blue")
lines(cum_BAB_sd_bull, col="green")
legend("topleft", legend=c("Market", "Original BAB", "SD-Enhanced BAB"), 
       col=c("red", "blue", "green"), lty=1)

# Bear marked
cum_market_bear <- cumprod(1 + na.omit(market_bear))
plot(cum_market_bear, type="l", col="red", 
     main="Bear Market: BAB vs Market", 
     xlab="Time", ylab="Cumulative Return")
lines(cum_BAB_bear, col="blue")
lines(cum_BAB_sd_bear, col="green")
legend("topleft", legend=c("Market", "Original BAB", "SD-Enhanced BAB"), 
       col=c("red", "blue", "green"), lty=1)

par(mfrow=c(1,1))

# Udskriv hele tabel
print("Complete Performance Summary (Monthly Returns, Annual Volatility):")
print(full_summary)

#------------------------ P1 til 10 porteføljer --------------------------------
# Antal porteføljer
num_portfolios <- 10

portfolio_breaks <- lapply(z, function(rank_values) {
  if (!is.null(rank_values)) {
    quantile(rank_values, probs = seq(0, 1, length.out = num_portfolios + 1), na.rm = TRUE)
  } else {
    return(NULL)
  }
})

portfolio_weights <- vector("list", num_portfolios)
names(portfolio_weights) <- paste0("w_P", 1:num_portfolios)

for (p in 1:num_portfolios) {
  portfolio_weights[[p]] <- mapply(function(z_values, breaks) {
    if (!is.null(z_values) && !is.null(breaks)) {
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

Market$MarketReturn <- FFM3[,2] + FFM3[,5]

# Tager log 
Market <- log(Market[-1]/100+1)
FFM3 <- log(FFM3/100+1)

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

# Tilføjer date til Beta
Beta <- as.data.frame(Beta)
Date <- as.data.frame(Date)
Beta$Date <- Date[3:nrow(Date),]
Beta <- Beta %>%
  mutate(Date = as.Date(as.character(Date), format="%Y%m%d"),
         YearMonth = format(Date, "%Y-%m"))  # YYYY-MM format

# Tilføjer date til Beta Shrink
Beta_shrink <- as.data.frame(Beta_shrink)
Date <- as.data.frame(Date)
Beta_shrink$Date <- Date[3:nrow(Date),]
Beta_shrink <- Beta_shrink %>%
  mutate(Date = as.Date(as.character(Date), format="%Y%m%d"),
         YearMonth = format(Date, "%Y-%m"))  # YYYY-MM format

# Funktion til at beregne Beta_PX for en given portefølje
compute_Beta_P <- function(weights_lookup) {
  sapply(1:(nrow(Beta_shrink)-1), function(i) {
    current_month <- Beta_shrink$YearMonth[i]  # Find YYYY-MM fra Beta
    
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

# Funktion til at beregne Beta_PX for en given portefølje uden shrink
compute_Beta_P_uden_shrink <- function(weights_lookup) {
  sapply(1:(nrow(Beta)-1), function(i) {
    current_month <- Beta$YearMonth[i]  # Find YYYY-MM fra Beta
    
    if (!(current_month %in% names(weights_lookup))) {
      print(paste("Fejl! Måneden", current_month, "findes ikke i weights_lookup"))
      return(NA)
    }
    
    # Beregn vægtet beta
    return(sum(weights_lookup[[current_month]] * 
                 Beta[i, names(weights_lookup[[current_month]]), drop = FALSE], 
               na.rm = TRUE))
  })
}

# Beregn Beta_PX for alle porteføljer og gem som individuelle variabler uden shrink
for (p in portfolio_names) {
  assign(paste0("Beta_uden_shrink_", p), compute_Beta_P_uden_shrink(weights_lookups[[p]]))
}

# Tilføjer date til Beta
Beta <- as.data.frame(Beta)
Date <- as.data.frame(Date)
Beta$Date <- Date[3:nrow(Date),]
Beta <- Beta %>%
  mutate(Date = as.Date(as.character(Date), format="%Y%m%d"),
         YearMonth = format(Date, "%Y-%m"))  # YYYY-MM format

# Tilføjer date til Beta_shrink
Beta_shrink <- as.data.frame(Beta_shrink)
Date <- as.data.frame(Date)
Beta_shrink$Date <- Date[3:nrow(Date),]
Beta_shrink <- Beta_shrink %>%
  mutate(Date = as.Date(as.character(Date), format="%Y%m%d"),
         YearMonth = format(Date, "%Y-%m"))  # YYYY-MM format

rf <- FFM3[,5]
rf <- rf[2:(min_length+1)]
månedligrf <- (((1+mean(rf, na.rm = TRUE))^(250/12))-1)*100

#--------------------------- Alpha beregning -----------------------------------
portfolio_names <- paste0("P", 1:10)

# Opret en liste til at gemme alpha-værdierne
capm_alphas <- numeric(length(portfolio_names))
ff3_alphas <- numeric(length(portfolio_names))
real_betas <- numeric(length(portfolio_names))

Market_excess_return <- FFM3[, 2][2:(min_length+3)]
SMB <- FFM3[, 3][2:(min_length+3)]
HML <- FFM3[, 4][2:(min_length+3)]

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

#----------------------------- Sharpe ratio ----------------------------
portfolio_names <- paste0("P", 1:10)

# Igangsæt vektorer for resultater 
portfolio_volatility <- numeric(length(portfolio_names))
sharpe_ratios <- numeric(length(portfolio_names))

# Kør igennem hver portefølje
for (i in seq_along(portfolio_names)) {
  # Afkast
  portfolio_return <- get(paste0("r_", portfolio_names[i], "_next"))
  
  # Volatilitet 
  portfolio_volatility[i] <- sd(portfolio_return, na.rm = TRUE) * sqrt(250)
  
  # Merafkast
  mean_excess_return <- mean(portfolio_return, na.rm = TRUE)
  
  # Sharpe ratio
  sharpe_ratios[i] <- mean_excess_return / portfolio_volatility[i]
}

# Gem resultater i data frame
sharpe_results <- data.frame(
  Portfolio = portfolio_names,
  Volatility = portfolio_volatility*100,
  Sharpe_Ratio = sharpe_ratios*100
)

# Print resultater
print(sharpe_results)

# Tilføj t-test for afkast og alpha'er
# Igangsæt vektorer til at gemme t-stat
portfolio_return_tstat <- numeric(10)
capm_alpha_tstat <- numeric(10)
ff3_alpha_tstat <- numeric(10)

# Funktion der tilføjer signifikansniveau
add_significance <- function(t_stat) {
  if (abs(t_stat) > 2.58) return(" ***")  # 1% nivea
  else if (abs(t_stat) > 1.96) return(" **")  # 5% niveau
  else if (abs(t_stat) > 1.65) return(" *")  # 10% niveau
  else return("")
}

# Udregn t-stat for hver portefølje
for (i in 1:10) {
  # Afkast
  portfolio_return <- get(paste0("r_P", i, "_next"))
  
  # t-stat for gns. afkast
  n_obs <- sum(!is.na(portfolio_return))
  se <- sd(portfolio_return, na.rm = TRUE) / sqrt(n_obs)
  portfolio_return_tstat[i] <- mean(portfolio_return, na.rm = TRUE) / se
  
  # t-stat for alpha'er
  capm_model <- lm(portfolio_return ~ Market_excess_return)
  capm_summary <- summary(capm_model)
  capm_alpha_tstat[i] <- capm_summary$coefficients[1, 3]  # t-stat for intercept
  
  ff3_model <- lm(portfolio_return ~ Market_excess_return + SMB + HML)
  ff3_summary <- summary(ff3_model)
  ff3_alpha_tstat[i] <- ff3_summary$coefficients[1, 3]  # t-stat for intercept
}

# Udregn månedlig merafkast med signifikans
monthly_excess_returns <- sapply(1:10, function(i) {
  r <- get(paste0("r_P", i, "_next"))
  monthly_return <- (((1+mean(r, na.rm = TRUE))^(250/12))-1)*100 - månedligrf
  return(paste0(round(monthly_return, 2), add_significance(portfolio_return_tstat[i])))
})

# CAPM alpha'er med signifikans
capm_alphas_formatted <- sapply(1:10, function(i) {
  paste0(round(capm_alphas[i], 2), add_significance(capm_alpha_tstat[i]))
})

# FF3 alpha'er med signifikans
ff3_alphas_formatted <- sapply(1:10, function(i) {
  paste0(round(ff3_alphas[i], 2), add_significance(ff3_alpha_tstat[i]))
})

# Lav tabel
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

# Print tabel
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

export_r_BAB_to_excel(r_BAB_vector = r_BAB, Date[4:22693,], filename = "/Users/OliverSnedker_1/Desktop/Kandidat/Speciale/BAB 49.xlsx")
export_r_BAB_to_excel(r_BAB_vector = r_BAB_sd, Date[4:22693,], filename = "/Users/OliverSnedker_1/Desktop/Kandidat/Speciale/BAB_SD 49.xlsx")
export_r_BAB_to_excel(r_BAB_vector = BetaSpread, Date[4:22693,], filename = "/Users/OliverSnedker_1/Desktop/Kandidat/Speciale/BS 49.xlsx")
export_r_BAB_to_excel(r_BAB_vector = BetaSpreadSD, Date[4:22693,], filename = "/Users/OliverSnedker_1/Desktop/Kandidat/Speciale/BS_SD 49.xlsx")

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
  
  # For original BAB
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
  
  # For SD-forbedret BAB
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
       xlim = c(0, max_breakeven * base_tc_rate * 10000 * 1.05),
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

cat("Break-even transaktionsomkostning for BAB:", round(breakeven_tc_original * 10000, 2), "basis point\n")
cat("Break-even transaktionsomkostning for SDBAB:", round(breakeven_tc_sd * 10000, 2), "basis point\n")

par(mfrow = c(1, 1))
breakeven_results <- plot_tc_sensitivity(
  r_BAB, r_BAB_sd, Beta_L, Beta_H, Beta_L_SD, Beta_H_SD,
  r_L_next, r_H_next, r_L_next_sd, r_H_next_sd, rf,
  transaction_costs_full, borrowing_costs_full,
  transaction_costs_full_sd, borrowing_costs_full_sd
)

improvement_ratio <- breakeven_tc_sd / breakeven_tc_original
cat("\nForbedring i break-even transaktionsomkostninger: ", 
    round((improvement_ratio - 1) * 100, 2), "% højere for SDBAB\n")

annual_turnover_original <- 1 / breakeven_tc_original
annual_turnover_sd <- 1 / breakeven_tc_sd

cat("\nBreak-even årlig turnover for BAB :", round(annual_turnover_original, 2), "gange\n")
cat("Break-even årlig turnover for SDBAB:", round(annual_turnover_sd, 2), "gange\n")

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
  
  # For original BAB
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
  
  # For SD-forbedret BAB
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
