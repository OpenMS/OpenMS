
# void getAAFrequencies(Map[String, size_t]) nogil except + # wrap-ignore
AASequence$$set("public","getAAFrequencies", function(map) {
  if( !(is.environment(map) && identical(parent.env(map), asNamespace("collections")) && strsplit(capture.output(map$$print())," ")[[1]][1] == "dict"
      && all(sapply(map$$keys(),function(k)  is_scalar_character(k) ))
      && all(sapply(map$$values(),function(v) is_scalar_integer(v) )))
  ) { stop("arg map wrong type")}
    k <- map$$keys()
    v <- lapply(map$$values(),as.integer)
    mmap <- py_dict(k,v)
    private$$py_obj$$getAAFrequencies(mmap)
    invisible()
} )

