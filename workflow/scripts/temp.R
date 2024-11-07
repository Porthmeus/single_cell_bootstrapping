require(data.table)
pops <- round(seq(30,1000, length.out =10))
draws <- seq(0.01,0.5, length.out = 10)
n <- 2:100

calcJacc <- function(x){
    jaccs <-c()
    for(i in 1:length(x)){
        x1 <- x[[i]]
        for(j in 2:length(x)){
            if(j > i){
                x2 <- x[[j]]
                jaccs <- c(jaccs,length(intersect(x1,x2))/length(union(x1,x2)))
            }
        }
    }
    return(jaccs)
}

sims <- length(pops)*length(draws)*length(n)

sim_pops <- rep(0, sims)
sim_draws <- rep(0, sims)
sim_n <- rep(0, sims)
jacc <- rep(0,sims)
y <- 1
for(pop in pops){
    for(draw in draws){
        for(nn in n){
            draw <- round(pop*draw)
            samples <- lapply(1:nn, sample, x=1:pop, size=draw)
            jacc_mean <- mean(calcJacc(samples))
            sim_pops[y] <- pop
            sim_draws[y] <- draw
            sim_n[y] <- nn
            jacc[y] <- jacc_mean
            y <- y+1
        }
    }
}
                

N <- 10000
n <- 50

a <- sapply(1:1000,function(x){
    samples <- lapply(1:2, sample, x=1:N, size=n)
    calcJacc(samples)})
t.test(a, mu = n/(2*N))

n/(2*N)

