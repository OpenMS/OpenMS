
MzXMLFile$$set("public","transform", function(path,transformer) {
    if(!((is.R6(path) && class(path)[1]=="String") || is_scalar_character(path))) { stop("arg path wrong type") }
    if(!( is.R6(transformer) && class(transformer)[1] %in% c("CachedSwathFileConsumer","MSDataAggregatingConsumer","MSDataCachedConsumer","MSDataSqlConsumer","MSDataStoringConsumer","MzMLSwathFileConsumer","NoopMSDataWritingConsumer","PlainMSDataWritingConsumer","RegularSwathFileConsumer")  )) { stop("arg transformer wrong type") }
    private$$py_obj$$transform(r_to_py(path),r_to_py(transformer))
    invisible()
} )
