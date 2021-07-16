
String$$set("public","initialize",function(in_0){
  if(missing(in_0)){
    private$$py_obj <- Pymod$$String()
  } else {
    if(is_scalar_character(in_0)){
      private$$py_obj <- Pymod$$String(in_0)
    } else if (is.R6(in_0) && class(in_0)[1] == "String") {
       private$$py_obj <- Pymod$$String(r_to_py(in_0))
    } else if ("python.builtin.object" %in% class(in_0) && class_to_wrap(in_0) == "String") {
       private$$py_obj <- in_0
    } else {
       stop(paste("invalid argument","can only be character or <String> r6 object",sep = " "))
    }
  }
},overwrite=TRUE)

String$$set("public","toString",function(){
  private$$py_obj$$toString()
})