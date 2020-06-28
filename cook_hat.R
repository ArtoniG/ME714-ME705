library(tidyverse)
library(gridExtra)

cook_hat <- function(fit.model){
  cook_hat_df <- tibble(cooks = cooks.distance(fit.model),
                        hat_measure = hatvalues(fit.model),
                        index = names(hatvalues(fit.model)) %>% as.numeric())
  cut <- sum(cook_hat_df$hat_measure)*2/nrow(cook_hat_df)
  hat_fit <- cook_hat_df %>% ggplot(aes(x = index, y = hat_measure)) + geom_point() +
    geom_hline(yintercept = cut, linetype = "dashed") +
    labs(x="Índice", y="Medida h") + theme_bw()
  cook_fit <- cook_hat_df %>% ggplot(aes(x = index, y = cooks)) + geom_point() +
    labs(x="Índice", y="Distância de Cook") + theme_bw()
  grid.arrange(hat_fit,cook_fit, ncol = 2)
}
