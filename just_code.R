eg_dat     <- runif(n = 48, min = -4, max = 8);
eg_mat     <- matrix(data = eg_dat, ncol = 4);
eg_mat[,2] <- eg_mat[,1] + runif(n = 12, min = -2, max = 2);
rownames(eg_mat) <- c("Sample_1", "Sample_2",  "Sample_3",  "Sample_4",
                      "Sample_5", "Sample_6",  "Sample_7",  "Sample_8",
                      "Sample_9", "Sample_10", "Sample_11", "Sample_12");
colnames(eg_mat) <- c("Variable_1", "Variable_2", "Variable_3", "Variable_4");


par(mfrow = c(1, 2), mar = c(5, 5, 1, 0.5), lwd = 2)
hist(x = eg_mat[,1], cex.lab = 1.25, cex.axis = 1.25, col = "blue", 
     xlab = "Variable 1", main = "");
par(mar = c(5, 5, 1, 1));
plot(x = eg_mat[,1], y = eg_mat[,2], ylim = c(-8, 10), xlim = c(-8, 10), 
     asp = 1, cex.lab = 1.25, cex.axis = 1.25, pch = 20, cex = 4, 
     xlab = "Variable 1", ylab = "Variable 2");
text(x = eg_mat[,1], y = eg_mat[,2], labels = 1:12, col = "white", 
     cex = 0.8);



# Plotting the actual data
par(mfrow = c(1, 2), mar = c(5, 5, 1, 0.5), lwd = 2)
plot(x = eg_mat[,1], y = eg_mat[,2], ylim = c(-9, 11), xlim = c(-9, 11), 
     asp = 1, cex.lab = 1.25, cex.axis = 1.25, pch = 20, cex = 4, 
     xlab = "Variable 1", ylab = "Variable 2");
text(x = eg_mat[,1], y = eg_mat[,2], labels = 1:12, col = "white", 
     cex = 0.8);
mod <- lm(eg_mat[,2] ~ eg_mat[,1]);
b0  <- as.numeric(mod$coefficients[1]);
b1  <- as.numeric(mod$coefficients[2]);
abline(a = b0, b = b1, lwd = 2, col = "red");
points(x = eg_mat[,1], y = eg_mat[,2], pch = 20, cex = 4);
text(x = eg_mat[,1], y = eg_mat[,2], labels = 1:12, col = "white", 
     cex = 0.8);
# Principal component analysis
pc_eg <- prcomp(x = eg_mat[,1:2]);
par(mar = c(5, 5, 1, 1));
plot(x = pc_eg$x[,1], y = pc_eg$x[,2], ylim = c(-4, 4), xlim = c(-10, 10), 
     asp = 1, cex.lab = 1.25, cex.axis = 1.25, pch = 20, cex = 4, 
     xlab = "Principal Component 1", ylab = "Principal Component 2");
text(x = pc_eg$x[,1], y = pc_eg$x[,2], labels = 1:12, col = "white", 
     cex = 0.8);
abline(h = 0, lwd = 2, col = "red");
points(x = pc_eg$x[,1], y = pc_eg$x[,2], cex = 4, pch = 20);
text(x = pc_eg$x[,1], y = pc_eg$x[,2], labels = 1:12, col = "white", 
     cex = 0.8);



# No correlation PCA
A_comp <- NULL;
A_dat  <- rnorm(n = 400*400, mean = 0, sd = 0.4);
A_mat  <- matrix(data = A_dat, nrow = 400);
A0_e   <- eigen(A_mat)$values;
A0_r   <- Re(A0_e) + 8;
A0_i   <- Im(A0_e) + 8;
A0_meg <- cbind(A0_r, A0_i);
par(mfrow = c(1, 2), mar = c(5, 5, 1, 0.5), lwd = 2);
plot(x = A0_meg[,1], y = A0_meg[,2], ylim = c(-2, 18), xlim = c(-2, 18), 
     asp = 1, cex.lab = 1.25, cex.axis = 1.25, pch = 20, cex = 1, 
     xlab = "Variable 1", ylab = "Variable 2");
A0_pc <- prcomp(x = A0_meg[,1:2]);
par(mar = c(5, 5, 1, 1));
plot(x = A0_pc$x[,1], y = A0_pc$x[,2], ylim = c(-10, 10), xlim = c(-10, 10), 
     asp = 1, cex.lab = 1.25, cex.axis = 1.25, pch = 20, cex = 1, 
     xlab = "Principal Component 1", ylab = "Principal Component 2");


# Total correlation
ecrr1  <-  runif(n = 48, min = 2, max = 16);
ecrr2  <-  ecrr1 + 0.4;
ecorr  <- cbind(ecrr1, ecrr2);
par(mfrow = c(1, 2), mar = c(5, 5, 1, 0.5), lwd = 2);
plot(x = ecorr[,1], y = ecorr[,2], ylim = c(-2, 22), xlim = c(-2, 22), 
     asp = 1, cex.lab = 1.25, cex.axis = 1.25, pch = 20, cex = 1, 
     xlab = "Variable 1", ylab = "Variable 2");
ecorr <- prcomp(x = ecorr[,1:2]);
par(mar = c(5, 5, 1, 1));
plot(x = ecorr$x[,1], y = ecorr$x[,2], ylim = c(-12, 12), xlim = c(-12, 12), 
     asp = 1, cex.lab = 1.25, cex.axis = 1.25, pch = 20, cex = 1, 
     xlab = "Principal Component 1", ylab = "Principal Component 2");


# Fig wasps
fig_wasps    <- read.csv("data/wing_loadings.csv", header = TRUE);
fig_cols     <- colnames(fig_wasps);
print(fig_cols);

# Heterandrium
het1 <- fig_wasps[fig_wasps[,1] == "Het1",] # Species in panel 'f' above
het2 <- fig_wasps[fig_wasps[,1] == "Het2",] # Species in panel 'e' above
het  <- rbind(het1, het2);



dat <- het[,-(1:4)]; # Remove columns 1 through 4
dat <- as.matrix(dat);
rownames(dat) <- NULL;
head(dat);
pairs(x = dat, gap = 0, cex.labels = 0.5);


# PCA
dat <- scale(dat);
pca_dat <- prcomp(dat);
print(pca_dat);
par(mar = c(5, 5, 1, 1));
plot(x = pca_dat$x[,1], y = pca_dat$x[,2], asp = 1, cex.lab = 1.25, 
     cex.axis = 1.25, xlab = "PC1", ylab = "PC2");
par(mar = c(5, 5, 1, 1));
plot(x = pca_dat$x[,1], y = pca_dat$x[,3], asp = 1, cex.lab = 1.25, 
     cex.axis = 1.25, xlab = "PC1", ylab = "PC3");


par(mar = c(5, 5, 1, 1));
plot(x = pca_dat$x[,1], y = pca_dat$x[,2], asp = 1, cex.lab = 1.25, 
     cex.axis = 1.25, xlab = "PC1", ylab = "PC2");
h1_rw <- which(het[,1] == "Het1"); # Recall the original data set
h2_rw <- which(het[,1] == "Het2"); # Recall the original data set
points(x = pca_dat$x[h1_rw,1], y = pca_dat$x[h1_rw,2], pch = 20, col = "red");
points(x = pca_dat$x[h2_rw,1], y = pca_dat$x[h2_rw,2], pch = 20, col = "blue");
legend("topleft", fill = c("red", "blue"), cex = 1.5, 
       legend = c( expression(paste(italic("Heterandrium "), 1)), 
                   expression(paste(italic("Heterandrium "), 2))));


# Screeplot
pca_variance      <- (pca_dat$sdev) * (pca_dat$sdev); # Get variance
print(pca_variance);
screeplot(pca_dat, npcs = 11, main = "", type = "lines", cex.lab = 1.5);

# Barplot
pca_pr_explained  <- pca_variance / sum(pca_variance);
print(pca_pr_explained);
pc_names <- paste("PC", 1:length(pca_pr_explained), sep = "");
barplot(height = pca_pr_explained * 100, names = pc_names, cex.names = 0.8,
        ylab = "Per cent of total variation explained", cex.lab = 1.25);

# Biplot
biplot(pca_dat, cex = 0.8, asp = 1);
print(pca_dat$rotation[,1:2]); # Note: 11 total columns; one for each PC
# Note that the arrows are scaled by a factor of ca 4.5, hence the scaling below
points(x = 4.5 * pca_dat$sdev[1] * pca_dat$rotation[,1],
       y = 4.5 * pca_dat$sdev[2] * pca_dat$rotation[,2],
       col = "blue", pch = 20, cex = 1.5);

# Recreate PCA
plot(x = 0, y = 0, type = "n", xlim = c(-6.2, 4.2), ylim = c(-4, 3), asp = 1,
     xlab = "PC1", ylab = "PC2", cex.axis = 1.25, cex.lab = 1.25);
for(i in 1:dim(dat)[1]){ # Take the sum of the measurements times PC loadings
    pt_PC1 <- sum(dat[i,] * pca_dat$rotation[,1]); # measurements * loadings 1
    pt_PC2 <- sum(dat[i,] * pca_dat$rotation[,2]); # measurements * loadings 2
    points(x = pt_PC1, y = pt_PC2, cex = 4, col = "black", pch = 20);
    text(x = pt_PC1, y = pt_PC2, col = "white", labels = i, cex = 0.8);
}


# Matrix algebra
print(eg_mat);
eg_m <- eg_mat[,1:2]; # Now we have an R object with just two columns
prcomp(eg_m);
V <- cov(eg_m); # Variance covariance matrix of two simulated measurements
print(V);
eigen(V); # Eigenvalues and eigenvectors of V

# A vector
u <- c(2, 1);
plot(x = 0, y = 0, type = "n", xlim = c(-1, 3), ylim = c(-1, 3), xlab = "",
     ylab = "");
arrows(x0 = 0, y0 = 0, x1 = u[1], y1 = u[2], length = 0.1, lwd = 3);
abline(h = 0, lty = "dotted", lwd = 0.8);
abline(v = 0, lty = "dotted", lwd = 0.8);


x  <- matrix(data = c(2, 1), ncol = 1); 
Vx <- V %*% x;

U <- eigen(V)$vectors;
print(U);

U_inv <- solve(U);
print(U_inv);

# Get back to covariance matrix with eigenvalues and vectors
U %*% U_inv; # Identity matrix, with a bit of a rounding error.
L       <- matrix(data = 0, nrow = 2, ncol = 2);
L[1, 1] <- eigen(V)$values[1];
L[2, 2] <- eigen(V)$values[2];

U %*% L %*% U_inv;

print(V); print(U); print(U_inv); print(L);






