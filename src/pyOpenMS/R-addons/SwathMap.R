
SwathMap$$set("public","getSpectrumPtr",function() {
  tryCatch({
    obj <- private$$py_obj$$getSpectrumPtr()
    ans <- eval(parse(text = paste0(class_to_wrap(obj),"$$","new(obj)")))
    return(ans)
           }, error = function(e) { "Did not find suitable conversion to Python object" })
})

SwathMap$$set("public","setSpectrumPtr",function(arg) {
  if(!(is.R6(arg) && class(arg)[1] %in% c("SpectrumAccessOpenMS","SpectrumAccessOpenMSCached","SpectrumAccessOpenMSInMemory","SpectrumAccessQuadMZTransforming") )) { stop("Need to provide suitable ISpectrumAccess-derived child class") }
  private$$py_obj$$setSpectrumPtr(r_to_py(arg))
  invisible()
})