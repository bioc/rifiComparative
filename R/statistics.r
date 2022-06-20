##############################Statistics###############################
#check if each segment is significant using t-test checking distance from 0
#t-test for Half-life
#select HL fragments
setwd("path_new_directory_data_comparison/")
statistics <- function(data){

    frag_HL <-
    unique(data[,"HL_comb_fragment"][
        grep(paste0("Dc_\\d+$"), df_comb[,"HL_comb_fragment"])])
    
    frag_int <-
    unique(df_comb[",intensity_comb_fragment"][
        grep(paste0("I_\\d+$"), df_comb[,"intensity_comb_fragment"])])

    data <-
        t_test_function(
            data = data,
            par = "HL",
            par1 = "half_life",
            cdt1 = "sc",
            cdt2 = "fe",
            frag = frag_HL
        )
    p_adjusted <-
        as.numeric(p.adjust(as.numeric(data[ ,"p_value_distance_HL"]),
                            method = "fdr"))
    data <-
        tibble::add_column(data, formatC(p_adjusted, format = "e", digits = 2),
                           .after = ncol(data)-1)
    
    colnames(data)[ncol(data)] <- "p_adjusted_HL"
    
    data <-
        t_test_function(
            data = data,
            par = "intensity",
            par1 = "intensity",
            cdt1 = "sc",
            cdt2 = "fe",
            frag = frag_int
        )

        p_adjusted <-
            as.numeric(p.adjust(as.numeric(data[ ,"p_value_distance_intensity"]),
                                method = "fdr"))
        data <-
            tibble::add_column(data, formatC(p_adjusted, format = "e", digits = 2),
                               .after = ncol(data)-1)
        
        colnames(data)[ncol(data)] <- "p_adjusted_intensity"
        return(data)
        }



df_comb <- statistics(data = df_comb)

#remove duplicated intensity fragments to save it on table
df_comb_uniq <-
    df_comb[!duplicated(df_comb$intensity_comb_fragment), ]

#save the data
write_xlsx(df_comb_uniq, "df_comb_uniq_se.xlsx")
save(df_comb, file="df_comb_se.rda")
