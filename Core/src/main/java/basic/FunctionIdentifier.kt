package basic

data class FunctionIdentifier<T>(val type: Class<T>, val number: Int)

data class FunctionSignature(val valueT: Class<*>, val gradientT: Class<*>, val hessianT: Class<*>)


//HOW TO COMPILE KOTLIN FIRST IN MAVEN?
