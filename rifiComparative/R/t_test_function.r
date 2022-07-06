t_test_function <- function(data, par, par1, frag_HL, frag_int){
    
    data[,paste0("p_value_distance_", par)] <- NA
    
    if(par == "HL"){
        frag <- frag_HL
    }else{
        frag <- frag_int
    }
    for (i in seq_along(frag)) {
        par.cd1 <- data[which(data[,paste0(par,"_comb_fragment")] == frag[i]),
                        paste0(par1, ".cdt1")]
        par.cd2 <- data[which(data[,paste0(par,"_comb_fragment")] == frag[i]),
                        paste0(par1, ".cdt2")]
        if (length(par.cd1) < 3 | length(par.cd2) < 3 ) {
            next ()
        } 
        t_h <-
            t.test(par.cd1,
                   par.cd2,
                   alternative = "two.sided",
                   var.equal = FALSE)
        #fill in the dataframe with the p_value 
        data[which(data[,paste0(par,"_comb_fragment")] == frag[i]), 
             paste0("p_value_distance_", par)] <- t_h[[3]]
    }
    return(data)
}
