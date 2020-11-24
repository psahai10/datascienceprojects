a <- 3
class(a)
b <- 3.14
class(b)
a <- as.integer(a)
class(a)

vec_a <- cbind(1,2,3,4,5)
vec_b <- cbind(6,7,8,9,10)

matrix <- rbind(vec_a, vec_b)

class(matrix)

list1 <- list(1,2)

list2 <- list(c(1,2), c(3,4))

list3 <- list(3, c(1, 2), "lists are amazing!")

list <- list(list1, list2, list3)


seq(10)
seq_len(10)
seq_along(10)

paste("hello", 'amigo', sep='--')
paste0('hello', 'amigo')

paste0(c('Mary','Joseph', 'Jesus'), collapse=', and ')

       