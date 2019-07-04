suscess <- require(ggplot2)

if(suscess) {
    
    library(ggplot2)
    
    args <- commandArgs(trailingOnly=TRUE)
    data <- read.csv(args[1])

    products <- c()
    l <- length(data[,'PRODUTO'])

    for(i in c(1 : l) ) {

        products <- c(products, as.vector(rep(data[i,'PRODUTO'], 6)))

    }

    conditions <- rep(c("EQUALS 0","0% ~ 0.1%","0.1% ~ 1%","1% ~ 10%","bigger than 10%","NA"), l)

    values <- c()

    for(i in c(1 : l) ) {

        values <- c(values, as.numeric(data[i, 2:7]))

    }

    data=data.frame(products,conditions,values)

    print(data)

    p <- ggplot(data, aes(fill=conditions, y=values, x=products)) +
        geom_bar( stat="identity", position="fill")
    
    plot(p)

} else {

    install.packages('ggplot2')

}
    