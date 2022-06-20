#set the path for your comparison folder
setwd("path_new_directory_data_comparison/")

fragmentation <- function(data){
#I. Dynamic Programming for HL
    data <- fragment_HL(
    probe = data,
    cores = 60,
    pen = pen_HL[[1]][[9]],
    pen_out = pen_HL[[1]][[10]]
    )

#II. Dynamic Programming for intensity
    data <- fragment_inty(
    probe = data,
    cores = 60,
    pen = pen_int[[1]][[9]],
    pen_out = pen_int[[1]][[10]]
    )
    return(data)
}

df_comb <- fragmentation(data=df_comb)
save(df_comb, file="df_comb_se.rda")
