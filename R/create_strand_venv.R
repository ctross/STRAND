#' A function to create a Python venv for STRAND to use for NumPyro models
#'
#' This is a simple helper function to get Python running smoothly
#'
#' @param back_end A string: either "cpu" for basic models, or "gpu" for Linux systems with NVIDIA GPUs.
#' @param rebuild The user needs to run this once as TRUE to rebuild. Inside of STRAND we will load the pre-built venv.
#' @param verbose If TRUE, it will print details of GPU detections.
#' @param cpu_versions Libraries to install for CPU build.
#' @param gpu_versions Libraries to install for GPU build. Only cuda13 is tested.
#' @export
#'

create_strand_venv = function(back_end = "cpu", 
                              rebuild = FALSE, 
                              verbose = TRUE,
                              cpu_versions = c("numpy>=1.27","numpyro>=0.18", "jax>=0.7", "jaxlib>=0.7"), 
                              gpu_versions = c("numpy>=1.27","numpyro>=0.18", "jax[cuda13]>=0.7", "jaxlib>=0.7")
                              ){
  if(rebuild == TRUE | reticulate::virtualenv_exists("strand_venv")==FALSE){
   reticulate::virtualenv_create("strand_venv", force = TRUE)
   reticulate::use_virtualenv("strand_venv")
   reticulate::py_config()

   if(back_end == "cpu"){
    # Install CPU toolchain 
    reticulate::py_install(cpu_versions)
   }

   if(back_end == "gpu"){
    # Check OS
    user_os = Sys.info()[["sysname"]]
     if(user_os != 'Linux'){
      stop("No GPU support for non-Linux devices.")
      }

    # Check for GPU
     if(suppressWarnings(system("nvidia-smi", ignore.stdout = TRUE, ignore.stderr = TRUE)) != 0){
      stop("No NVIDIA GPU detected")
     }

    # Install GPU toolchain 
    reticulate::py_install(gpu_versions)
   }
  
  if(verbose == TRUE){
    if(back_end == "gpu"){
    jax = reticulate::import("jax")
    devices = jax$devices()
    device_strs = sapply(seq_along(devices), function(i) reticulate::py_str(devices[[i]]))
    gpu_count = sum(grepl("Cuda|Gpu|cuda|gpu", device_strs))
    
    message("Found ", gpu_count, " GPU(s): ", paste(device_strs, collapse=", "))
    
    if(gpu_count == 0){
      stop("No GPUs detected! Ensure that your local toolchain is up-to-date.")
     }
    }

    if(back_end == "cpu"){
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
 }

}
