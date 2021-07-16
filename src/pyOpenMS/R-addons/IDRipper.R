
# Input types:
# ripped   :  dict object ( see Collections::dict() ) of key -> R6 object of <String> & value -> list(list(<ProteinIdentification>),list(<PeptideIdentification>))
# proteins :  list(<ProteinIdentification>)
# peptides :  list(<PeptideIdentification>)
IDRipper$$set("public","rip",function(ripped,proteins,peptides) {
    if( !(is.environment(ripped) && identical(parent.env(ripped), asNamespace("collections")) && strsplit(capture.output(ripped$$print())," ")[[1]][1] == "dict"
      && all(sapply(ripped$$keys(),function(k) is.R6(k) && class(k)[1] == "String" ))
      && all(sapply(ripped$$values(),function(v) len(v) == 2 && all(sapply(v[[1]],function(p) is.R6(p) && class(p)[1]=="ProteinIdentification")) && all(sapply(v[[2]],function(p) is.R6(p) && class(p)[1]=="PeptideIdentification")) )))
  ) { stop("arg ripped wrong type")}
  if(!(is_list(proteins) && all(sapply(proteins,function(p) is.R6(p) && class(p)[1] == "ProteinIdentification"))) ) { stop("wrong arg proteins") }
  if(!(is_list(peptides) && all(sapply(peptides,function(p) is.R6(p) && class(p)[1] == "PeptideIdentification"))) ) { stop("wrong arg peptides") }
  d <- py_dict(lapply(ripped$$keys(), function(k) r_to_py(k)),modify_depth(ripped$$values(), 3, function(v) r_to_py(v)))
  private$$py_obj$$rip(d,lapply(proteins,function(p) r_to_py(p)),lapply(peptides,function(p) r_to_py(p)))
  k <- lapply(py_to_r(py_builtin$$list(d$$keys())),function(k) String$$new(k))
  v <- modify_depth( py_to_r(py_builtin$$list(d$$values())), 3, function(v) eval(parse(text = paste(class_to_wrap(v),"$$","new(v)"))) )
  out_ripped <- collections::dict(v,k)
  tryCatch({
      eval.parent(substitute(ripped <- out_ripped))
      invisible()
           }, error = function(e){ invisible() })
} )

