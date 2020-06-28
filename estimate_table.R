library(generics)
library(magrittr)
estimate_tibble <- function(fit.model){
  summary(fit.model) %>% tidy() %>%
    mutate(lower = estimate - std.error*1.96, upper = estimate + std.error*1.96) %>%
    mutate_if(is.numeric, funs(round(.,3))) %>%
    mutate(IC = paste("[",lower,";",upper,"]"),
           p.value = if_else(p.value < 0.001, "< 0.001", as.character(p.value))) %>%
    dplyr::select(1:3,IC,statistic, p.value) %>%
    set_colnames(c("Termo", "Estimativa", "Erro Padrão", "IC (95%)", "Estatística t", "p-valor"))
}