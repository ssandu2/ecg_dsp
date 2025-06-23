library(signal)#--
library(ggplot2)
test_sample <- 3 
ecg_signal <- samples[test_sample, ,1]
samp_freq <- 500
#Bandpass filter ---------------
bandpass_filter <- butter(2, c(0.5, 30)/(samp_freq/2), type = "pass")
ecg_filtered <- filtfilt(bandpass_filter, ecg_signal)

#QRS detection------------------------
#diff functionally serve as derivatives
detect_qrs_local_maxima<- function(ecg_signal, samp_freq) {
first_diff <- diff(ecg_signal)
second_diff <- diff(first_diff)

maxima_indices <- which(first_diff[-1] <0 & first_diff[-length(first_diff)] >0) +1
threshold <- 0.4 * max(ecg_filtered)
valid_maxima <- maxima_indices[ecg_filtered[maxima_indices] > threshold]
return(valid_maxima)
}
qrs_indices <- detect_qrs_local_maxima(ecg_filtered, 500)
plot(ecg_signal, type = "l", col = "blue", main = "Detected QRS Peaks", ylab = "Amplitude", xlab = "Samples")
points(qrs_indices, ecg_signal[qrs_indices], col = "red", pch = 20)

#Extracting QRS complexes
window_size <- 60

qrs_segments <- sapply(qrs_indices, function(i) { 
  start <- max(1, i- window_size)
  end <- min(length(ecg_signal), i+ window_size)
  segment <- ecg_signal[start:end]
  length_diff <- (window_size *2 +1) - length(segment)
  ifelse(length_diff > 0, segment <- c(segment, rep(0, length_diff)),segment <- c(segment))
  
  return(segment)
}
)

#Constructing data matrix X (each col is a segment)
X <- t(qrs_segments)

#Performing SVD on data matrix X
svd_result <- svd(X)

# Extract the singular values and vectors
u <- svd_result$u
s <- diag(svd_result$d)
v <- svd_result$v


#QRS template reconstruction using 10 or fewer singular values (10 orthogonal basis vectors or fewer)
num_singular_values <- min(ncol(u), 10)
reconstructed_qrs_template <- u[, 1:num_singular_values] %*% s[1:num_singular_values, 1:num_singular_values]%*% t(v[, 1:num_singular_values])

#mean QRS beat
mean_qrs_beat <- rowMeans(reconstructed_qrs_template)


#Isolation of AA (atrial activity) signal
aa_signal <- ecg_filtered
for (i in qrs_indices) {
  start <- max(1, i - window_size)
  end <- min(length(ecg_signal), i+ window_size)
  segment_length <- end - start + 1
  aa_signal[start:end] <- aa_signal[start:end] - mean_qrs_beat[1:segment_length]
}

#Plot of original vs AA signal
combecg_plot_df <- data.frame(
  Time = seq_along(ecg_signal),
  Original = ecg_signal,
  AA = aa_signal
)

ecg_plot <- ggplot(combecg_plot_df, aes(x = Time)) +
  geom_line(aes(y = Original, color = "Original")) +
  geom_line(aes(y = AA, color = "AA")) +
  labs(title = "Original ECG signal vs Isolated Atrial Activity", x = "Time(samples)", y = "Amplitude (mV)", color = "Signal Type") +
  scale_color_manual(values = c("Original" = "blue", "AA" = "red")) +
  theme_minimal()

print(ecg_plot)
