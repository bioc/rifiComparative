#' An example SummarizedExperiment from Synechosystis PCC 6803 first condition
#' obtained from rifi_statistics and used as input for rifiComparative
#'
#' @format A rowRanges of SummarizedExperiment with 500 rows and 50 variables:
#' \describe{
#'   \item{seqnames:}{The sequence name chromosome}
#'   \item{start:}{The bin/probe start position}
#'   \item{end:}{The bin/probe end position}
#'   \item{width:}{The bin/probe length}
#'   \item{strand:}{The strand specific}
#'   \item{position:}{The bin/probe specific position}
#'   \item{ID:}{The bin/probe specific ID}
#'   \item{FLT:}{The bin/probe flag for background level}
#'   \item{intensity:}{The relative intensity at time point 0}
#'   \item{probe_TI:}{An internal value to determine which fitting model is
#'   applied}
#'   \item{flag:}{Information on which fitting model is applied}
#'   \item{position_segment:}{The position based segment}
#'   \item{delay:}{The delay value of the bin/probe}
#'   \item{half_life:}{The half-life of the bin/probe}
#'   \item{TI_termination_factor:}{The termination factor of the bin/probe}
#'   \item{delay_fragment:}{The delay fragment the bin belongs to}
#'   \item{velocity_fragment:}{The velocity value of the respective delay
#'   fragment}
#'   \item{intercept:}{The vintercept of fit through the respective delay
#'   fragment}
#'   \item{slope:}{The slope of the fit through the respective delay fragment}
#'   \item{HL_fragment:}{The half-life fragment the bin belongs to}
#'   \item{HL_mean_fragment:}{The mean half-life value of the respective
#'   half-life fragment}
#'   \item{intensity_fragment:}{The intensity fragment the bin belongs to}
#'   \item{intensity_mean_fragment:}{The mean intensity value of the respective
#'   intensity fragment}
#'   \item{TU:}{The overarching transcription unit}
#'   \item{TI_termination_fragment:}{The TI fragment the bin belongs to}
#'   \item{TI_mean_termination_factor:}{The mean termination factor of the
#'   respective TI fragment}
#'   \item{seg_ID:}{The combined ID of the fragment}
#'   \item{pausing_site:}{presence of pausing site indicated by +/-}
#'   \item{iTSS_I:}{presence of iTSS_I indicated by +/-}
#'   \item{ps_ts_fragment:}{The fragments involved in pausing site or iTSS_I}
#'   \item{event_ps_itss_p_value_Ttest:}{p_value of pausing site or iTSS_I}#'   
#'   \item{delay_frg_slope:}{the slope value of the respective delay fragment}
#'   \item{p_value_slope:}{p_value of the slope}
#'   \item{velocity_ratio:}{Integer, ratio of velocity between 2 delay fragments}
#'   \item{event_duration:}{Integer, the duration between two delay fragments}
#'   \item{event_position:}{Integer, the position middle between 2 fragments with an event}
#'   \item{FC_HL:}{Integer, the fold change value of 2 HL fragments}
#'   \item{FC_fragment_HL:}{Integer, the fold change value of 2 intensity fragments}
#'   \item{p_value_HL:}{p_value of the fold change of HL fragments}
#'   \item{FC_intensity:}{Integer, the fold change value of 2 intensity fragments}
#'   \item{FC_fragment_intensity:}{String, fragments involved in fold change 
#'   between 2 intensity fragments}
#'   \item{p_value_intensity:}{p_value of the fold change of intensity fragments}
#'   \item{FC_HL_intensity:}{ratio of fold change between 2 half-life fragments 
#'   and fold change between 2 intensity fragments}
#'   \item{FC_HL_intensity_fragment:}{fragments involved on ratio of fold change
#'   between 2 half-life fragments and fold change between 2 intensity fragments}
#'   \item{FC_HL_adapted:}{Integer, the fold change of half-life/ fold change
#'   of intensity, position of the half-life fragment is adapted to intensity 
#'   fragment}
#'   \item{synthesis_ratio:}{Integer, the value correspomding to synthesis rate}
#'   \item{synthesis_ratio_event:}{String, the event assigned by synthesis rate
#'   either Termination or iTSS}
#'   \item{p_value_Manova:}{p_value of the variance between two fold-changes, 
#'   HL and intensity}
#'   \item{p_value_TI:}{p_value of TI fragment} 
#'   \item{TI_fragments_p_value:}{p_value of 2 TI fragments}
#' }
#'
#' @source \url{https://github.com/CyanolabFreiburg/rifiComparative}
#'
#' @usage data(stats_se_cdt1)  
#'
#'
"stats_se_cdt1"


#' An example SummarizedExperiment from Synechosystis PCC 6803 second condition
#' obtained from rifi_statistics and used as input for rifiComparative
#'
#' @format A rowRanges of SummarizedExperiment with 500 rows and 50 variables:
#' \describe{
#'   \item{seqnames:}{The sequence name chromosome}
#'   \item{start:}{The bin/probe start position}
#'   \item{end:}{The bin/probe end position}
#'   \item{width:}{The bin/probe length}
#'   \item{strand:}{The strand specific}
#'   \item{position:}{The bin/probe specific position}
#'   \item{ID:}{The bin/probe specific ID}
#'   \item{FLT:}{The bin/probe flag for background level}
#'   \item{intensity:}{The relative intensity at time point 0}
#'   \item{probe_TI:}{An internal value to determine which fitting model is
#'   applied}
#'   \item{flag:}{Information on which fitting model is applied}
#'   \item{position_segment:}{The position based segment}
#'   \item{delay:}{The delay value of the bin/probe}
#'   \item{half_life:}{The half-life of the bin/probe}
#'   \item{TI_termination_factor:}{The termination factor of the bin/probe}
#'   \item{delay_fragment:}{The delay fragment the bin belongs to}
#'   \item{velocity_fragment:}{The velocity value of the respective delay
#'   fragment}
#'   \item{intercept:}{The vintercept of fit through the respective delay
#'   fragment}
#'   \item{slope:}{The slope of the fit through the respective delay fragment}
#'   \item{HL_fragment:}{The half-life fragment the bin belongs to}
#'   \item{HL_mean_fragment:}{The mean half-life value of the respective
#'   half-life fragment}
#'   \item{intensity_fragment:}{The intensity fragment the bin belongs to}
#'   \item{intensity_mean_fragment:}{The mean intensity value of the respective
#'   intensity fragment}
#'   \item{TU:}{The overarching transcription unit}
#'   \item{TI_termination_fragment:}{The TI fragment the bin belongs to}
#'   \item{TI_mean_termination_factor:}{The mean termination factor of the
#'   respective TI fragment}
#'   \item{seg_ID:}{The combined ID of the fragment}
#'   \item{pausing_site:}{presence of pausing site indicated by +/-}
#'   \item{iTSS_I:}{presence of iTSS_I indicated by +/-}
#'   \item{ps_ts_fragment:}{The fragments involved in pausing site or iTSS_I}
#'   \item{event_ps_itss_p_value_Ttest:}{p_value of pausing site or iTSS_I}
#'   \item{p_value_slope:}{p_value of the slope}
#'   \item{delay_frg_slope:}{the slope value of the respective delay fragment}
#'   \item{velocity_ratio:}{Integer, ratio of velocity between 2 delay fragments}
#'   \item{event_duration:}{Integer, the duration between two delay fragments}
#'   \item{event_position:}{Integer, the position middle between 2 fragments 
#'   with an event}
#'   \item{FC_HL:}{Integer, the fold change value of 2 HL fragments}
#'   \item{FC_fragment_HL:}{Integer, the fold change value of 2 intensity fragments}
#'   \item{p_value_HL:}{p_value of the fold change of HL fragments}
#'   \item{FC_intensity:}{Integer, the fold change value of 2 intensity fragments}
#'   \item{FC_fragment_intensity:}{String, fragments involved in fold change 
#'   between 2 intensity fragments}
#'   \item{p_value_intensity:}{p_value of the fold change of intensity fragments}
#'   \item{FC_HL_intensity:}{ratio of fold change between 2 half-life fragments 
#'   and fold change between 2 intensity fragments}
#'    \item{FC_HL_intensity_fragment:}{fragments involved on ratio of fold 
#'    change between 2 half-life fragments and fold change between 2 intensity 
#'    fragments}
#'   \item{FC_HL_adapted:}{Integer, the fold change of half-life/ fold change
#'   of intensity, position of the half-life fragment is adapted to intensity 
#'   fragment}
#'   \item{synthesis_ratio:}{Integer, the value correspomding to synthesis rate}
#'   \item{synthesis_ratio_event:}{String, the event assigned by synthesis rate 
#'   either  Termination or iTSS}
#'   \item{p_value_Manova:}{p_value of the variance between two fold-changes, 
#'   HL and intensity}
#'   \item{p_value_TI:}{p_value of TI fragment} 
#'   \item{TI_fragments_p_value:}{p_value of 2 TI fragments} 
#' }
#'
#' @source \url{https://github.com/CyanolabFreiburg/rifiComparative}
#'
#' @usage data(stats_se_cdt2)
#'
#'
"stats_se_cdt2"


#' An example data frame from Synechosystis PCC 6803 differential probes
#' expression obtained from limma package and only interesting variables were
#' selected. The data frame was used entirely.
#'
#' @format A data frame of differential_expression with 55508 rows and 4 variables:
#' \describe{
#'   \item{position:}{The bin/probe specific position}
#'   \item{strand:}{The strand specific}
#'   \item{logFC_int:}{The bin/probe differential expression}
#'   \item{P.Value:}{The bin/probe p_value adjusted}
#' }
#'
#' @source \url{https://github.com/CyanolabFreiburg/rifiComparative}
#'
#' @usage data(differential_expression)
#'
#'
"differential_expression"

#' The result of loading_fun for stats_se_cdt1 example data
#' Two data frame containing the output of loading_fun as first element of a
#' list.
#'
#' @format A data frame with 500 rows and 49 variables:
#' \describe{
#'   \item{strand:}{The strand specific}
#'   \item{position:}{The bin/probe specific position}
#'   \item{ID:}{The bin/probe specific ID}
#'   \item{FLT:}{The bin/probe flag for background level}
#'   \item{intensity:}{The relative intensity at time point 0}
#'   \item{probe_TI:}{An internal value to determine which fitting model is
#'   applied}
#'   \item{flag:}{Information on which fitting model is applied}
#'   \item{position_segment:}{The position based segment}
#'   \item{delay:}{The delay value of the bin/probe}
#'   \item{half_life:}{The half-life of the bin/probe}
#'   \item{TI_termination_factor:}{The termination factor of the bin/probe}
#'   \item{delay_fragment:}{The delay fragment the bin belongs to}
#'   \item{velocity_fragment:}{The velocity value of the respective delay
#'   fragment}
#'   \item{intercept:}{The vintercept of fit through the respective delay
#'   fragment}
#'   \item{slope:}{The slope of the fit through the respective delay fragment}
#'   \item{HL_fragment:}{The half-life fragment the bin belongs to}
#'   \item{HL_mean_fragment:}{The mean half-life value of the respective
#'   half-life fragment}
#'   \item{intensity_fragment:}{The intensity fragment the bin belongs to}
#'   \item{intensity_mean_fragment:}{The mean intensity value of the respective
#'   intensity fragment}
#'   \item{TU:}{The overarching transcription unit}
#'   \item{TI_termination_fragment:}{The TI fragment the bin belongs to}
#'   \item{TI_mean_termination_factor:}{The mean termination factor of the
#'   respective TI fragment}
#'   \item{seg_ID:}{The combined ID of the fragment}
#'   \item{pausing_site:}{presence of pausing site indicated by +/-}
#'   \item{iTSS_I:}{presence of iTSS_I indicated by +/-}
#'   \item{ps_ts_fragment:}{The fragments involved in pausing site or iTSS_I}
#'   \item{event_ps_itss_p_value_Ttest:}{p_value of pausing site or iTSS_I}
#'   \item{p_value_slope:}{p_value of the slope}
#'   \item{delay_frg_slope:}{the slope value of the respective delay fragment}
#'   \item{velocity_ratio:}{Integer, ratio of velocity between 2 delay fragments}
#'   \item{event_duration:}{Integer, the duration between two delay fragments}
#'   \item{event_position:}{Integer, the position middle between 2 fragments with an event}
#'   \item{FC_HL:}{Integer, the fold change value of 2 HL fragments}
#'   \item{FC_fragment_HL:}{Integer, the fold change value of 2 intensity fragments}
#'   \item{p_value_HL:}{p_value of the fold change of HL fragments}
#'   \item{FC_intensity:}{Integer, the fold change value of 2 intensity fragments}
#'   \item{FC_fragment_intensity:}{String, fragments involved in fold change 
#'   between 2 intensity fragments}
#'   \item{p_value_intensity:}{p_value of the fold change of intensity fragments}
#'   \item{FC_HL_intensity:}{ratio of fold change between 2 half-life fragments 
#'   and fold change between 2 intensity fragments}
#'   \item{FC_HL_intensity_fragment:}{fragments involved on ratio of fold 
#'   change between 2 half-life fragments and fold change between 2 intensity 
#'   fragments}
#'   \item{FC_HL_adapted:}{Integer, the fold change of half-life/ fold change of
#'   intensity, position of the half-life fragment is adapted to intensity fragment}
#'   \item{synthesis_ratio:}{Integer, the value correspomding to synthesis rate}
#'   \item{synthesis_ratio_event:}{String, the event assigned by synthesis rate either  Termination or iTSS}
#'   \item{p_value_Manova:}{p_value of the variance between two fold-changes, 
#'   HL and intensity}
#'   \item{p_value_TI:}{p_value of TI fragment}
#'   \item{TI_fragments_p_value:}{p_value of 2 TI fragments} 
#'   \item{cdt:}{The condition assigned to the experiment here cdt1}
#'   \item{logFC_int:}{The bin/probe log2 fold change of intensity at time 0}
#'   \item{P.Value:}{The bin/probe p_value adjusted }
#' }
#'
#' @source \url{https://github.com/CyanolabFreiburg/rifiComparative}
#'
#' @usage data(inp_s)  
#'
#'
"inp_s"


#' The result of loading_fun for stats_se_cdt2 example data
#' Two data frame containing the output of loading_fun as second element of a
#' list.
#'
#' @format A data frame with 500 rows and 49 variables:
#' \describe{
#'   \item{strand:}{The strand specific}
#'   \item{position:}{The bin/probe specific position}
#'   \item{ID:}{The bin/probe specific ID}
#'   \item{FLT:}{The bin/probe flag for background level}
#'   \item{intensity:}{The relative intensity at time point 0}
#'   \item{probe_TI:}{An internal value to determine which fitting model is
#'   applied}
#'   \item{flag:}{Information on which fitting model is applied}
#'   \item{position_segment:}{The position based segment}
#'   \item{delay:}{The delay value of the bin/probe}
#'   \item{half_life:}{The half-life of the bin/probe}
#'   \item{TI_termination_factor:}{The termination factor of the bin/probe}
#'   \item{delay_fragment:}{The delay fragment the bin belongs to}
#'   \item{velocity_fragment:}{The velocity value of the respective delay
#'   fragment}
#'   \item{intercept:}{The vintercept of fit through the respective delay
#'   fragment}
#'   \item{slope:}{The slope of the fit through the respective delay fragment}
#'   \item{HL_fragment:}{The half-life fragment the bin belongs to}
#'   \item{HL_mean_fragment:}{The mean half-life value of the respective
#'   half-life fragment}
#'   \item{intensity_fragment:}{The intensity fragment the bin belongs to}
#'   \item{intensity_mean_fragment:}{The mean intensity value of the respective
#'   intensity fragment}
#'   \item{TU:}{The overarching transcription unit}
#'   \item{TI_termination_fragment:}{The TI fragment the bin belongs to}
#'   \item{TI_mean_termination_factor:}{The mean termination factor of the
#'   respective TI fragment}
#'   \item{seg_ID:}{The combined ID of the fragment}
#'   \item{pausing_site:}{presence of pausing site indicated by +/-}
#'   \item{iTSS_I:}{presence of iTSS_I indicated by +/-}
#'   \item{ps_ts_fragment:}{The fragments involved in pausing site or iTSS_I}
#'   \item{event_ps_itss_p_value_Ttest:}{p_value of pausing site or iTSS_I}
#'   \item{p_value_slope:}{p_value of the slope}
#'   \item{delay_frg_slope:}{the slope value of the respective delay fragment}
#'   \item{velocity_ratio:}{Integer, ratio of velocity between 2 delay fragments}
#'   \item{event_duration:}{Integer, the duration between two delay fragments}
#'   \item{event_position:}{Integer, the position middle between 2 fragments with an event}
#'   \item{FC_HL:}{Integer, the fold change value of 2 HL fragments}
#'   \item{FC_fragment_HL:}{Integer, the fold change value of 2 intensity fragments}
#'   \item{p_value_HL:}{p_value of the fold change of HL fragments}
#'   \item{FC_intensity:}{Integer, the fold change value of 2 intensity fragments}
#'   \item{FC_fragment_intensity:}{String, fragments involved in fold change 
#'   between 2 intensity fragments}
#'   \item{p_value_intensity:}{p_value of the fold change of intensity fragments}
#'   \item{FC_HL_intensity:}{ratio of fold change between 2 half-life fragments 
#'   and fold change between 2 intensity fragments}
#'    \item{FC_HL_intensity_fragment:}{fragments involved on ratio of fold 
#'    change between 2 half-life fragments and fold change between 2 intensity 
#'    fragments}
#'   \item{FC_HL_adapted:}{Integer, the fold change of half-life/ fold change of
#'   intensity, position of the half-life fragment is adapted to intensity 
#'   fragment}
#'   \item{synthesis_ratio:}{Integer, the value correspomding to synthesis rate}
#'   \item{synthesis_ratio_event:}{String, the event assigned by synthesis 
#'   rate either  Termination or iTSS}
#'   \item{p_value_Manova:}{p_value of the variance between two fold-changes,
#'    HL and intensity}
#'   \item{p_value_TI:}{p_value of TI fragment} 
#'   \item{TI_fragments_p_value:}{p_value of 2 TI fragments} 
#'   \item{cdt:}{The condition assigned to the experiment here cdt2}
#'   \item{logFC_int:}{The bin/probe log2 fold change of intensity at time 0}
#'   \item{P.Value:}{The bin/probe p_value adjusted }
#' }
#'
#' @source \url{https://github.com/CyanolabFreiburg/rifiComparative}
#'
#' @usage data(inp_f)  
#'
#'
"inp_f"


#' The result of joining_by_row for inp_s and inp_f example data
#' A data frame containing the output of joining_by_row as a data frame
#' 
#' @format A data frame with 600 rows and 49 variables:
#' \describe{
#'   \item{strand:}{The strand specific}
#'   \item{position:}{The bin/probe specific position}
#'   \item{ID:}{The bin/probe specific ID}
#'   \item{FLT:}{The bin/probe flag for background level}
#'   \item{intensity:}{The relative intensity at time point 0}
#'   \item{probe_TI:}{An internal value to determine which fitting model is
#'   applied}
#'   \item{flag:}{Information on which fitting model is applied}
#'   \item{position_segment:}{The position based segment}
#'   \item{delay:}{The delay value of the bin/probe}
#'   \item{half_life:}{The half-life of the bin/probe}
#'   \item{TI_termination_factor:}{The termination factor of the bin/probe}
#'   \item{delay_fragment:}{The delay fragment the bin belongs to}
#'   \item{velocity_fragment:}{The velocity value of the respective delay
#'   fragment}
#'   \item{intercept:}{The vintercept of fit through the respective delay
#'   fragment}
#'   \item{slope:}{The slope of the fit through the respective delay fragment}
#'   \item{HL_fragment:}{The half-life fragment the bin belongs to}
#'   \item{HL_mean_fragment:}{The mean half-life value of the respective
#'   half-life fragment}
#'   \item{intensity_fragment:}{The intensity fragment the bin belongs to}
#'   \item{intensity_mean_fragment:}{The mean intensity value of the respective
#'   intensity fragment}
#'   \item{TU:}{The overarching transcription unit}
#'   \item{TI_termination_fragment:}{The TI fragment the bin belongs to}
#'   \item{TI_mean_termination_factor:}{The mean termination factor of the
#'   respective TI fragment}
#'   \item{seg_ID:}{The combined ID of the fragment}
#'   \item{pausing_site:}{presence of pausing site indicated by +/-}
#'   \item{iTSS_I:}{presence of iTSS_I indicated by +/-}
#'   \item{ps_ts_fragment:}{The fragments involved in pausing site or iTSS_I}
#'   \item{event_ps_itss_p_value_Ttest:}{p_value of pausing site or iTSS_I}
#'   \item{p_value_slope:}{p_value of the slope}
#'   \item{delay_frg_slope:}{the slope value of the respective delay fragment}
#'   \item{velocity_ratio:}{Integer, ratio of velocity between 2 delay fragments}
#'   \item{event_duration:}{Integer, the duration between two delay fragments}
#'   \item{event_position:}{Integer, the position middle between 2 fragments 
#'   with an event}
#'   \item{FC_fragment_HL:}{Integer, the fold change value of 2 intensity fragments}
#'   \item{FC_HL:}{Integer, the fold change value of 2 HL fragments}#'   
#'   \item{p_value_HL:}{p_value of the fold change of HL fragments}
#'   \item{FC_intensity:}{Integer, the fold change value of 2 intensity 
#'   fragments}
#'   \item{FC_fragment_intensity:}{String, fragments involved in fold change 
#'   between 2 intensity fragments}
#'   \item{p_value_intensity:}{p_value of the fold change of intensity fragments}
#'   \item{FC_HL_intensity:}{ratio of fold change between 2 half-life fragments 
#'   and fold change between 2 intensity fragments}
#'    \item{FC_HL_intensity_fragment:}{fragments involved on ratio of fold 
#'    change between 2 half-life fragments and fold change between 2 intensity 
#'    fragments}
#'   \item{FC_HL_adapted:}{Integer, the fold change of half-life/ fold change of
#'   intensity, position of the half-life fragment is adapted to intensity 
#'   fragment}
#'   \item{synthesis_ratio:}{Integer, the value correspomding to synthesis rate}
#'   \item{synthesis_ratio_event:}{String, the event assigned by synthesis 
#'   rate either  Termination or iTSS}
#'   \item{p_value_Manova:}{p_value of the variance between two fold-changes, 
#'   HL and intensity}
#'   \item{p_value_TI:}{p_value of TI fragment} 
#'   \item{TI_fragments_p_value:}{p_value of 2 TI fragments} 
#'   \item{cdt:}{The condition assigned to the experiment here cdt2}
#'   \item{logFC_int:}{The bin/probe log2 fold change of intensity at time 0}
#'   \item{P.Value:}{The bin/probe p_value adjusted }
#' }
#'
#' @source \url{https://github.com/CyanolabFreiburg/rifi}
#'
#' @usage data(data_combined_minimal)
#'
"data_combined_minimal"


#' The result of joining_by_column for data_combined_minimal example data
#' A data frame containing the output of joining_by_row as a data frame
#'
#' @format A data frame with 300 rows and 18 variables:
#' \describe{
#'   \item{strand:}{The strand specific}
#'   \item{position:}{The bin/probe specific position}
#'   \item{ID:}{The bin/probe specific ID}
#'   \item{intensity.cdt1:}{The relative intensity at time point 0 for 
#'   condition 1}
#'   \item{position_segment:}{The position based segment}
#'   \item{half_life.cdt1:}{The half-life of the bin/probe condition 1}
#'   \item{TI_termination_factor.cdt1:}{The termination factor of the bin/probe
#'   condition 1}
#'   \item{HL_fragment.cdt1:}{The half-life fragment the bin belongs to  
#'   condition 1}
#'   \item{intensity_fragment.cdt1:}{The intensity fragment the bin belongs to
#'   condition 1}
#'   \item{TI_termination_fragment.cdt1:}{The TI fragment the bin belongs to 
#'   condition 1}
#'   \item{logFC_int:}{The bin/probe log2 fold change of intensity at time 0}
#'   \item{P.Value:}{The bin/probe p_value adjusted }
#'   \item{intensity.cdt2:}{The relative intensity at time point 0  
#'   condition 2}
#'   \item{half_life.cdt2:}{The half-life of the bin/probe condition 2}
#'   \item{TI_termination_factor.cdt2:}{The termination factor of the bin/probe
#'   condition 2}
#'   \item{HL_fragment.cdt2:}{The half-life fragment the bin belongs to  
#'   condition 2}
#'   \item{intensity_fragment.cdt2:}{The intensity fragment the bin belongs to
#'   condition 2}
#'   \item{TI_termination_fragment.cdt2:}{The TI fragment the bin belongs to 
#'   condition 2}
#' }
#'
#' @source \url{https://github.com/CyanolabFreiburg/rifiComparative}
#'
#' @usage data(df_comb_minimal)  
#'
#'
"df_comb_minimal"

#' The result of penalties for df_comb_minimal example data
#' A data frame containing the output of penalties as a data frame
#'
#' @format A data frame with 300 rows and 20 variables:
#' \describe{
#'   \item{strand:}{The strand specific}
#'   \item{position:}{The bin/probe specific position}
#'   \item{ID:}{The bin/probe specific ID}
#'   \item{intensity.cdt1:}{The relative intensity at time point 0 for 
#'   condition 1}
#'   \item{position_segment:}{The position based segment}
#'   \item{half_life.cdt1:}{The half-life of the bin/probe condition 1}
#'   \item{TI_termination_factor.cdt1:}{The termination factor of the bin/probe
#'   condition 1}
#'   \item{HL_fragment.cdt1:}{The half-life fragment the bin belongs to  
#'   condition 1}
#'   \item{intensity_fragment.cdt1:}{The intensity fragment the bin belongs to
#'   condition 1}
#'   \item{TI_termination_fragment.cdt1:}{The TI fragment the bin belongs to 
#'   condition 1}
#'   \item{logFC_int:}{The bin/probe log2 fold change of intensity at time 0}
#'   \item{P.Value:}{The bin/probe p_value adjusted }
#'   \item{intensity.cdt2:}{The relative intensity at time point 0  
#'   condition 2}
#'   \item{half_life.cdt2:}{The half-life of the bin/probe condition 2}
#'   \item{TI_termination_factor.cdt2:}{The termination factor of the bin/probe
#'   condition 2}
#'   \item{HL_fragment.cdt2:}{The half-life fragment the bin belongs to  
#'   condition 2}
#'   \item{intensity_fragment.cdt2:}{The intensity fragment the bin belongs to
#'   condition 2}
#'   \item{TI_termination_fragment.cdt2:}{The TI fragment the bin belongs to 
#'   condition 2}
#'   \item{distance_HL:}{The bin/probe difference of half-life from both
#'   conditions}
#'   \item{distance_int:}{The bin/probe log2 fold change of intensity at time 0}
#' }
#'
#' @source \url{https://github.com/CyanolabFreiburg/rifiComparative}
#'
#' @usage data(penalties_df)  
#'
#'
"penalties_df"


#' The result of penalties for df_comb_minimal example data.
#' A list containing the output from penalties including the logbook and two
#' penalty objects.
#' 
#' @format A list with 5 items:
#' \describe{
#'   \item{pen_obj_HL:}{A list with 4 items:
#'     \describe{
#'       \item{logbook:}{The logbook vector containing half-life penalty
#'       information}
#'       \item{HL_penalties:}{a vetor with the half-life penalty and half-life
#'       outlier penalty}
#'       \item{correct:}{a matrix of the correct splits}
#'       \item{wrong:}{a matrix of the incorrect splits}
#'     }
#'   }
#' }
#'
#' @source \url{https://github.com/CyanolabFreiburg/rifi}
#'
#' @usage data(pen_HL)
#'
"pen_HL"

#' The result of penalties for df_comb_minimal example data.
#' A list containing the output from penalties including the logbook and two
#' penalty objects.
#' 
#' @format A list with 5 items:
#' \describe{
#'   \item{pen_int:}{A list with 4 items:
#'     \describe{
#'       \item{logbook:}{The logbook vector containing intensity penalty
#'       information}
#'       \item{int_penalties:}{a vector with the intensity penalty and
#'       intensity outlier penalty}
#'       \item{correct:}{a matrix of the correct splits}
#'       \item{wrong:}{a matrix of the incorrect splits}
#'     }
#'   }
#' }
#'
#' @source \url{https://github.com/CyanolabFreiburg/rifi}
#'
#' @usage data(pen_int)
#'
"pen_int"


#' The result of fragmentation for df_comb_minimal example data
#' A data frame containing the output of fragmentation as a data frame
#'
#' @format A data frame with 500 rows and 24 variables:
#' \describe{
#'   \item{strand:}{The strand specific}
#'   \item{position:}{The bin/probe specific position}
#'   \item{ID:}{The bin/probe specific ID}
#'   \item{intensity.cdt1:}{The relative intensity at time point 0 for 
#'   condition 1}
#'   \item{position_segment:}{The position based segment}
#'   \item{half_life.cdt1:}{The half-life of the bin/probe condition 1}
#'   \item{TI_termination_factor.cdt1:}{The termination factor of the bin/probe
#'   condition 1}
#'   \item{HL_fragment.cdt1:}{The half-life fragment the bin belongs to  
#'   condition 1}
#'   \item{intensity_fragment.cdt1:}{The intensity fragment the bin belongs to
#'   condition 1}
#'   \item{TI_termination_fragment.cdt1:}{The TI fragment the bin belongs to 
#'   condition 1}
#'   \item{logFC_int:}{The bin/probe log2 fold change of intensity at time 0}
#'   \item{P.Value:}{The bin/probe p_value adjusted }
#'   \item{intensity.cdt2:}{The relative intensity at time point 0  
#'   condition 2}
#'   \item{half_life.cdt2:}{The half-life of the bin/probe condition 2}
#'   \item{TI_termination_factor.cdt2:}{The termination factor of the bin/probe
#'   condition 2}
#'   \item{HL_fragment.cdt2:}{The half-life fragment the bin belongs to  
#'   condition 2}
#'   \item{intensity_fragment.cdt2:}{The intensity fragment the bin belongs to
#'   condition 2}
#'   \item{TI_termination_fragment.cdt2:}{The TI fragment the bin belongs to 
#'   condition 2}
#'   \item{distance_HL:}{The bin/probe difference of half-life from both
#'   conditions}
#'   \item{distance_int:}{The bin/probe log2 fold change of intensity at time 0}
#'   \item{HL_comb_fragment:}{The half-life fragment the bin belongs to both
#'   conditions}
#'   \item{HL_mean_comb_fragment:}{The half-life mean of the fragment the bin
#'   belongs to both conditions}
#'   \item{intensity_comb_fragment:}{The intensity fragment the bin belongs to 
#'   both conditions}
#'   \item{intensity_mean_comb_fragment:}{The intensity mean of the fragment the
#'   bin belongs to both conditions}
#' }
#'
#' @source \url{https://github.com/CyanolabFreiburg/rifiComparative}
#'
#' @usage data(fragment_int)
#'
#'
"fragment_int"

#' The result of statistics for fragment_int example data
#' A data frame containing the output of statistics as a data frame
#'
#' @format A data frame with 500 rows and 26 variables:
#' \describe{
#'   \item{strand:}{The strand specific}
#'   \item{position:}{The bin/probe specific position}
#'   \item{ID:}{The bin/probe specific ID}
#'   \item{intensity.cdt1:}{The relative intensity at time point 0 for 
#'   condition 1}
#'   \item{position_segment:}{The position based segment}
#'   \item{half_life.cdt1:}{The half-life of the bin/probe condition 1}
#'   \item{TI_termination_factor.cdt1:}{The termination factor of the bin/probe
#'   condition 1}
#'   \item{HL_fragment.cdt1:}{The half-life fragment the bin belongs to  
#'   condition 1}
#'   \item{intensity_fragment.cdt1:}{The intensity fragment the bin belongs to
#'   condition 1}
#'   \item{TI_termination_fragment.cdt1:}{The TI fragment the bin belongs to 
#'   condition 1}
#'   \item{logFC_int:}{The bin/probe log2 fold change of intensity at time 0}
#'   \item{P.Value:}{The bin/probe p_value adjusted }
#'   \item{intensity.cdt2:}{The relative intensity at time point 0  
#'   condition 2}
#'   \item{half_life.cdt2:}{The half-life of the bin/probe condition 2}
#'   \item{TI_termination_factor.cdt2:}{The termination factor of the bin/probe
#'   condition 2}
#'   \item{HL_fragment.cdt2:}{The half-life fragment the bin belongs to  
#'   condition 2}
#'   \item{intensity_fragment.cdt2:}{The intensity fragment the bin belongs to
#'   condition 2}
#'   \item{TI_termination_fragment.cdt2:}{The TI fragment the bin belongs to 
#'   condition 2}
#'   \item{distance_HL:}{The bin/probe difference of half-life from both
#'   conditions}
#'   \item{distance_int:}{The bin/probe log2 fold change of intensity at time 0}
#'   \item{HL_comb_fragment:}{The half-life fragment the bin belongs to both
#'   conditions}
#'   \item{HL_mean_comb_fragment:}{The half-life mean of the fragment the bin
#'   belongs to both conditions}
#'   \item{intensity_comb_fragment:}{The intensity fragment the bin belongs to 
#'   both conditions}
#'   \item{intensity_mean_comb_fragment:}{The intensity mean of the fragment the
#'   bin belongs to both conditions}
#'   \item{p_value_distance_HL:}{The p_value adjusted of the half-life fragment 
#'   the bin belongs to both conditions}
#'   \item{p_value_distance_intensity:}{The p_value adjusted of the intensity
#'   fragment the bin belongs to both conditions}
#' }
#'
#' @source \url{https://github.com/CyanolabFreiburg/rifiComparative}
#'
#' @usage data(stats_df_comb_minimal)
#'
#'
"stats_df_comb_minimal"


#' The result of adjusting_HLToInt for stats_df_comb_minimal and annotation
#' example data
#' A data frame containing the output of adjusting_HLToInt as a data frame
#'
#' @format A data frame with 52 rows and 15 variables:
#' \describe{
#'   \item{position:}{The bin/probe specific position}
#'   \item{region:}{the region from the gff file}
#'   \item{gene:}{the annotated gene name}
#'   \item{locus_tag:}{the annotated locus tag}
#'   \item{strand:}{The strand specific}
#'   \item{fragment_HL:}{The half-life fragment the bin belongs}
#'   \item{fragment_int:}{The intensity fragment the bin belongs}
#'   \item{position_frg_int:}{The position of the first fragment and the last 
#'   position of the last fragment}
#'   \item{mean_HL_fragment:}{The mean half-life value of the respective
#'   half-life fragments}
#'   \item{mean_int_fragment:}{The mean intensity value of the respective
#'   intensity fragments}
#'   \item{log2FC(decay_rate):}{log2FC(decay(condition1)/decay(condition2))}
#'   \item{Log2FC(HL)-Log2FC(int):}{log2FC(decay_rate/intensity)}
#'   \item{log2FC(synthesis_rate):}{log2FC(decay_rate) + 
#'      log2FC(intensity)}
#'   \item{intensity_FC:}{log2FC(mean(intensity(condition1))/mean(
#'      intensity(condition2)))}
#'   \item{p_value:}{indicated by "*" means at least one fragment 
#'      either HL fragment or intensity fragment has a significant p_value}
#' }
#'
#' @source \url{https://github.com/CyanolabFreiburg/rifiComparative}
#'
#' @usage data(df_mean_minimal)  
#'
#'
"df_mean_minimal"


#' The result of gff3_preprocessing of gff3 file
#' A list containing all necessary information from a gff file for 
#' adjusting_HLToInt and visualization. 
#'
#' @format A list with 2 items:
#' \describe{
#'   \item{data annotation:}{a data frame with 5853 rows and 6 variables
#'     \describe{
#'       \item{region:}{the region from the gff file}
#'       \item{start:}{the start of the annotation}
#'       \item{end:}{the end of the annotation}
#'       \item{strand:}{the strand of the annotation}
#'       \item{gene:}{the annotated gene name}
#'       \item{locus_tag:}{the annotated locus tag}
#'     }
#'   }
#'   \item{genome length:}{a numeric vector containing the length of the genome}
#' }
#' @source \url{https://github.com/CyanolabFreiburg/rifiComparative}
#'
#' @usage data(annot_g)  
#'
#'
"annot_g"