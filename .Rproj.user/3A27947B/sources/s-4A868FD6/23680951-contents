cirt <- function(data, 
                 family,
                 n_nodes = 121,
                 stand_errors = TRUE,
                 constraints = list(
                   fix_disps = NULL, fix_alphas = NULL,
                   same_disps = FALSE, same_alphas = FALSE
                 ),
                 control = list(
                   thres = Inf, prob = 0, init_disp_one = TRUE,
                   maxiter = 1000, 
                   convtol = 1e-5, ctol_maxstep = 1e-8,
                   m_method = "nleqslv", convcrit = "marglik"
                 )) {
  # TODO checks and data prep
  
  if (family == "cmp") {
    start_values <- get_start_values(
      data = data, init_disp_one = control$init_disp_one
    )
    
    fit <- run_newem(
      data = data, 
      init_params = start_values, 
      n_nodes = n_nodes, 
      fix_disps = constraints$fix_disps,
      fix_alphas = constraints$fix_alphas,
      same_disps = constraints$same_disps, 
      same_alphas = constraints$same_alphas,
      thres = control$thres,
      prob = control$prob,
      maxiter = control$maxiter, 
      convtol = control$convtol, 
      ctol_maxstep = control$ctol_maxstep,
      m_method = control$m_method, 
      convcrit = control$convcrit
    )
  
    
    +ü+üpwp
  
  if (stand_errors) {
    
  }
  
}