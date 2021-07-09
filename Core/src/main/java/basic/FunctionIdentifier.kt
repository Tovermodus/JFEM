package basic

data class FunctionIdentifier<valueT, gradientT, hessianT>(val vT: Class<valueT>,
                                                              val gT: Class<gradientT>,
                                                              val hT: Class<hessianT>,
                                                              val number: Int)