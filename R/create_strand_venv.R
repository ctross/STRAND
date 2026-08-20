#' A function to create a Python venv for STRAND to use for NumPyro models
#'
#' This is a simple helper function to get Python running smoothly
#'
#' @param rebuild The user needs to run this once as TRUE to rebuild. Inside of STRAND we will load the pre-built venv.
#' @param verbose If TRUE, it will print details of GPU detections.
#' @param mcmc_parameters A parameter list to control NUTS sampling. NumPyro uses the same control set as Stan.
#' @export
#'

create_strand_venv = function(rebuild = FALSE, 
                              verbose = TRUE,
                              mcmc_parameters
                              ){
  mcmc_parameters = merge_mcmc_parameters(mcmc_parameters)
  if(!mcmc_parameters$jax_device_type %in% c("cuda", "gpu", "cpu")){stop("JAX device type must be 'cpu', 'gpu', or 'cuda'.")} 

  if(rebuild == TRUE | reticulate::virtualenv_exists("strand_venv")==FALSE){
   reticulate::virtualenv_create("strand_venv", force = TRUE)
   reticulate::use_virtualenv("strand_venv")
   reticulate::py_config()

   if(mcmc_parameters$jax_device_type == "cpu"){
    # Install CPU toolchain 
    reticulate::py_install(mcmc_parameters$cpu_versions)
   }

   if(mcmc_parameters$jax_device_type %in% c("gpu","cuda")){
    # Check OS
    user_os = Sys.info()[["sysname"]]
     if(user_os != 'Linux'){
      stop("No GPU support for non-Linux devices.")
      }

    # Check for GPU
     if(suppressWarnings(system("nvidia-smi", ignore.stdout = TRUE, ignore.stderr = TRUE)) != 0){
      stop("No NVIDIA GPU detected.")
     }

    # Install GPU toolchain 
    reticulate::py_install(mcmc_parameters$gpu_versions)

    # Deal with missing cuSPARSE
    if(!is.null(mcmc_parameters$cusparse_path)){
     if(!file.exists(mcmc_parameters$cusparse_path)){
      stop("cuSPARSE library not found at: ", mcmc_parameters$cusparse_path)
     }

      ctypes = reticulate::import("ctypes")
      ctypes$CDLL(mcmc_parameters$cusparse_path)
      message("cuSPARSE library loaded successfully")
    }
   }
  
  if(verbose == TRUE){
    if(mcmc_parameters$jax_device_type %in% c("gpu","cuda")){
    jax = reticulate::import("jax")
    devices = jax$devices()
    device_strs = sapply(seq_along(devices), function(i) reticulate::py_str(devices[[i]]))
    gpu_count = sum(grepl("Cuda|Gpu|cuda|gpu", device_strs))
    
    message("Found ", gpu_count, " GPU(s): ", paste(device_strs, collapse=", "))
    
    if(gpu_count == 0){
      stop("No GPUs detected! Ensure that your local toolchain is up-to-date.")
     }
    }

    if(mcmc_parameters$jax_device_type == "cpu"){
    jax = reticulate::import("jax")
    devices = jax$devices()
    device_strs = sapply(seq_along(devices), function(i) reticulate::py_str(devices[[i]]))
    cpu_count = sum(grepl("Cpu|cpu|CPU", device_strs))
    
    message("Found ", cpu_count, " CPU(s): ", paste(device_strs, collapse=", "))
    
    if(cpu_count == 0){
      stop("No CPUs detected! How did you find yourself in this position?")
     }
    }
   }

   stop("This looks like you first time running a NumPyro model. Your toolchain has now been built, but you need to restart R before changes take effect. Try to run your model in a fresh R session.")
  }
 
 if(rebuild == FALSE){
   reticulate::use_virtualenv("strand_venv")
   reticulate::py_config()
   
    # Deal with missing cuSPARSE
    if(!is.null(mcmc_parameters$cusparse_path)){
     if(!file.exists(mcmc_parameters$cusparse_path)){
      stop("cuSPARSE library not found at: ", mcmc_parameters$cusparse_path)
     }

      ctypes = reticulate::import("ctypes")
      ctypes$CDLL(mcmc_parameters$cusparse_path)
      message("cuSPARSE library loaded successfully")
    }
 }

}



