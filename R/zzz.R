
.onLoad<- function(lib, pkg) 
{    
    packageStartupMessage("\n") 
    packageStartupMessage("Hola!\n") 
    packageStartupMessage("welcome to inevRsion package. \n \n \n")
    packageStartupMessage("type: manual() for full manual \n      vignette(\"inveRsion\") for a quick start \n")
}



.onUnload <- function(libpath)
{
    library.dynam.unload("inveRsion",libpath)

     cat("\n Hasta luego! \n")

}


manual<-function()
{
  pdf<-system.file("doc/Manual.pdf",package="inveRsion")

   if (.Platform$OS.type == "windows") {
     shell.exec(pdf)
    }else{
    system(paste(shQuote(getOption("pdfviewer")), shQuote(pdf)),
            wait = FALSE)
    }          
}
