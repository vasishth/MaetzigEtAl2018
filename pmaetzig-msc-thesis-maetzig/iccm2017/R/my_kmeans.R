my_kmeans <- function(data, k_max, iter.max=10, nstart=1, maintitle="") {
  kmeans_results <- vector(mode="list", length=k_max)
  
  for (i in 1:k_max) {
    kmeans_results[[i]] <- kmeans(data, centers=i, iter.max=iter.max)
  }
  
  kmeans_totss <- rep(0, k_max)
  kmeans_totss[1] <- kmeans_results[[1]]$tot.withinss
  
  for (k in 1:(k_max-1)) {
    kmeans_totss[k+1] <- kmeans_results[[k]]$tot.withinss - kmeans_results[[k+1]]$tot.withinss
  }
  kmeans_plot <- ggplot2::qplot(x=1:k_max, y=kmeans_totss, xlab="Number of Clusters", ylab="Total within cluster sum of squares",
                                main=maintitle)
  print(kmeans_plot)
  return(list(kmeans_results=kmeans_results, kmeans_totss=kmeans_totss, kmeans_plot=kmeans_plot))
}