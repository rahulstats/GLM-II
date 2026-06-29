
K <- 20

############################################################
# Convert nu vector to dominance matrix
############################################################

nu_to_matrix <- function(nu, K) {
  
  if(length(nu) != choose(K, 2))
    stop("Length of nu must equal choose(K,2).")
  
  M <- matrix(0, K, K)
  
  idx <- 1
  
  for(i in 1:(K - 1)) {
    for(j in (i + 1):K) {
      
      if(nu[idx] > 0) {
        M[i, j] <- 1
      } else if(nu[idx] < 0) {
        M[j, i] <- 1
      }
      
      idx <- idx + 1
    }
  }
  
  diag(M) <- 0
  
  M
}

############################################################
# Recursive tie breaking
############################################################

recursive_scores <- function(M,
                             teams,
                             level = 1,
                             result) {
  
  scores <- rowSums(M[teams, teams, drop = FALSE])
  
  level_name <- paste0("Level_", level)
  
  if(!(level_name %in% names(result)))
    result[[level_name]] <- NA
  
  result[teams, level_name] <- scores
  
  tied_groups <- split(teams, scores)
  
  for(g in tied_groups) {
    
    if(length(g) > 1) {
      
      sub_scores <- rowSums(M[g, g, drop = FALSE])
      
      # recurse only if scores are not all identical
      if(length(unique(sub_scores)) > 1) {
        
        result <- recursive_scores(
          M = M,
          teams = g,
          level = level + 1,
          result = result
        )
        
      }
    }
  }
  
  result
}

############################################################
# Main routine
############################################################

dominance_ranking <- function(nu, K = 20) {
  
  M <- nu_to_matrix(nu, K)
  
  result <- data.frame(
    Team = 1:K
  )
  
  result <- recursive_scores(
    M = M,
    teams = 1:K,
    level = 1,
    result = result
  )
  
  ##########################################################
  # Compute final dominance score
  #
  # D = L1 + L2/K + L3/K^2 + ...
  ##########################################################
  
  level_cols <- grep("^Level_", names(result),
                     value = TRUE)
  
  result[level_cols] <- lapply(
    result[level_cols],
    function(x) ifelse(is.na(x), 0, x)
  )
  
  result$DominanceScore <- 0
  
  for(l in seq_along(level_cols)) {
    
    result$DominanceScore <-
      result$DominanceScore +
      result[[level_cols[l]]] / K^(l - 1)
  }
  
  ##########################################################
  # Lexicographic ranking
  ##########################################################
  
  rank_order <- do.call(
    order,
    c(
      lapply(level_cols,
             function(x) -result[[x]])
    )
  )
  
  rank_table <- result[rank_order, ]
  
  rank_table$Rank <- seq_len(nrow(rank_table))
  
  ##########################################################
  # Return rows ordered by Team
  ##########################################################
  
  final_table <- rank_table[
    order(rank_table$Team),
  ]
  
  final_table <- final_table[
    ,
    c("Team",
      "Rank",
      "DominanceScore",
      level_cols)
  ]
  
  rownames(final_table) <- NULL
  
  final_table
}

############################################################
# EPL dominance ranks
############################################################


final_table <- dominance_ranking(nu_hat, K)

print(final_table)

############################################################
# Teams sorted by final rank
############################################################

cat("\nTeams ordered by dominance:\n")

print(
  final_table[
    order(final_table$Rank),
    c("Rank", "Team", "DominanceScore")
  ]
)