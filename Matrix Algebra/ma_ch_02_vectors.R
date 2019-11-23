# Vectors


a <- c(1, 2, 3)
a

b <- c(4, 5, 6)
b

c <- matrix(c(3, 2, 1), 3, 1, byrow = T)
d <- matrix(c(6, 5, 4), nrow = 3, ncol = 1, byrow = T)

c
d

class(b)

b <- as.matrix(b)
b

class(b)

u <- matrix(c(3, 2, 1), 1, 3,
            byrow = T)

v <- matrix(c(6, 5, 4),
            ncol = 1, nrow = 3, byrow = T)


# Matrices

U <- matrix(c(1, 2, 3, 4), 2, 2, byrow = T)
V <- matrix(c(5, 6, 7, 8), 2, 2, byrow = T)
W <- matrix(c(2, 2, 3, 5), 2, 2, byrow = T)
A <- matrix(c(1, 2, 3, 4, 5, 6), 2, 3, byrow = T)
B <- matrix(c(1, 2, 3, 4, 5, 6), 3, 2, byrow = T)

# Example Operations

U; V; W;

A; B

U + V; 2 * U

A %*% B;B %*% A

U %*% V; V %*% U

U %*% A; U %*% t(B)

U %*% W; W %*% U

V %*% W; W %*% V

