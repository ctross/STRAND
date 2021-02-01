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

simulate_multiplex_srm_network = function(N_id = 50,
                                          N_layers=3,
                                          theta=c(-3.8, -4, -3.5),
                                          sr_mu = rep(0,6),
                                          dr_mu = rep(0,6),

                                          s_sigma = c(1, 0.3, 0.7), 
                                          r_sigma = c(1.5, 1, 0.5),
                                          sr_rho = SR_Rho,
                                       
                                          dr_sigma = c(2, 1.75, 1.2),
                                          dr_rho = DR_Rho
                                         )
{

sr = rmvnorm2( N_id , Mu=c(sr_mu), sigma=c(s_sigma, r_sigma), Rho=sr_rho )

y1_true = array(0, c(N_id, N_id, N_layers))

DR = array(NA, c(N_id, N_id, N_layers))

for ( i in 1:(N_id-1) ) {
    for ( j in (i+1):N_id ) {

         dr = rmvnorm2(1, Mu=c(dr_mu), sigma=c(dr_sigma, dr_sigma), Rho=dr_rho)
         for(k in 1:N_layers){
          DR[i,j,k] = dr[k]
          DR[j,i,k] = dr[k+N_layers]
          }
         
         for(k in 1:N_layers){ 
         ########## Layer 1
         p_ij = inv_logit( theta[k] + sr[i,k] + sr[j,N_layers+k] + dr[k])
         y1_true[i,j,k] = rbern( 1 , p_ij )

         p_ji = inv_logit( theta[k] + sr[j,k] + sr[i,N_layers+k] + dr[N_layers+k])
         y1_true[j,i,k] = rbern( 1 , p_ji )
         }

        }#j
    }#i

return(list(Net = y1_true, SR=sr, DR=DR))
}

set.seed(8675309)
G = simulate_multiplex_srm_network()

