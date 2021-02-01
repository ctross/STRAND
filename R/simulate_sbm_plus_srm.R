#' A function to simulate single layer SRM data, with sender-receiver effects and dyadic reciprocity
#'
#' This function appends header data to automatically entered data batch-wise. The end result is a "csv" with the same structure as that
#' which is exported by the standard enter_data function, but with one extra column for game ID.
#' @param 
#' path Full path to main folder.
#' @param 
#' results A cell from the results list exported by auto_enter_all. 
#' @param 
#' HHID Household ID of focal individual/respondent.
#' @param 
#' RID ID of researcher.
#' @param 
#' day Day of interview.
#' @param 
#' month Month of interview.
#' @param 
#' year Year of interview.
#' @param 
#' name Name of focal individual/respondent.
#' @param 
#' ID Unique ID of focal individual/respondent.
#' @param 
#' game ID for case/game/question.
#' @param 
#' order Order of frames/panels of photos as presented to the respodent: e.g., with 4 frames, "ABCD", "CDBA", etc. are legal entries.
#' @param 
#' seed A seed for the random number generator to sort the order of photos in the array. This should match the seed used to make the survey.
#' @export
#' @examples
#' \dontrun{
#'annotate_batch_data(path = path, results=Game_all6, HHID="JKF", RID="CR", day=11, month=3, year=2020, 
#'              name = "Walter W.", PID="CVD", game="LikertData", order="AB", seed = 1)
#'                    }

simulate_sbm_plus_srm_network = function( N_id = 50,
                                          N_groups = 3L, 
                                          group_probs = c(0.2, 0.5, 0.3),
                                          in_block = 0.5, 
                                          out_block = 0.1,
                                          group_structure=NA,

                                          sr_mu = rep(0,2),
                                          dr_mu = rep(0,2),

                                          s_sigma = 1, 
                                          r_sigma = 1.5,
                                          sr_rho = 0.3,
                                       
                                          dr_sigma = 2.1,
                                          dr_rho = 0.6
                                         )
{

# sample ppl into groups
groups = sample( 1:N_groups , size=N_id , replace=TRUE , prob=group_probs )

# define interaction matrix across groups
if(is.na(group_structure)){
B = diag(N_groups)
for ( i in 1:length(B) ) if ( B[i]==0 ) B[i] = out_block
for ( i in 1:length(B) ) if ( B[i]==1 ) B[i] = in_block
} else{
B = group_structure  
}

SR_rho = diag(2)
SR_rho[2,1] = SR_rho[1,2] = sr_rho

DR_rho = diag(2)
DR_rho[2,1] = DR_rho[1,2] = dr_rho

sr = rmvnorm2( N_id , Mu=c(sr_mu), sigma=c(s_sigma, r_sigma), Rho=SR_rho )

y1_true = matrix(0, N_id, N_id)

DR = array(NA, c(N_id,N_id))

for ( i in 1:(N_id-1) ) {
    for ( j in (i+1):N_id ) {
          dr = rmvnorm2(1, Mu=c(dr_mu), sigma=c(dr_sigma, dr_sigma), Rho=DR_rho)
          DR[i,j] = dr[1]
          DR[j,i] = dr[2]

          p_ij = inv_logit( logit(B[ groups[i] , groups[j] ]) + sr[i,1] + sr[j,2] + dr[1])
          y1_true[i,j] = rbern( 1 , p_ij )

          p_ji = inv_logit( logit(B[ groups[j] , groups[i] ]) + sr[j,1] + sr[i,2] + dr[2])
          y1_true[j,i] = rbern( 1 , p_ji )

        }#j
    }#i

return(list(Net = y1_true, SR=sr, DR=DR))
}

set.seed(8675309)
G = simulate_sbm_plus_srm_network()

