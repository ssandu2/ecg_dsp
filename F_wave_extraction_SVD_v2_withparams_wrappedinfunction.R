extractfwave <- function(ecg_signal, samp_freq = 500, ecg_id = "ECG_ID") {

  #arguments
  input_ecg_signal <- ecg_signal
  samp_freq <- 500
  
  #cubic spline interpolation
  new_samp_freq <- 1000
  
  # Ensure that time_interpolated matches exactly 10,000 samples
  time_original <- seq(0, length(input_ecg_signal) - 1) / samp_freq
  time_interpolated <- seq(0, length(input_ecg_signal) / samp_freq, length.out = 10000)
  ecg_interpolated <- spline(time_original, input_ecg_signal, xout = time_interpolated)$y
  
  #zero-phase bandpass filter (0.5 to 30 Hz) to the interpolated signal
  bandpass_filter <- butter(2, c(0.5, 30) / (new_samp_freq / 2), type = "pass")
  ecg_filtered <- filtfilt(bandpass_filter, ecg_interpolated)
  
  #verify length of filtered signal
  cat("Length of ecg_filtered:", length(ecg_filtered), "\n")
  
  
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
  plot(ecg_filtered, type = "l", col = "blue", main = "Detected QRS Peaks", ylab = "Amplitude", xlab = "Samples")
  points(qrs_indices, ecg_filtered[qrs_indices], col = "red", pch = 20)
  
  #svd with edge Handling and normalization
  window_size <- 120  # Adjusted for the new sampling rate
  
  qrs_segments <- sapply(qrs_indices, function(i) { 
    start <- max(1, i - window_size)
    end <- min(length(ecg_interpolated), i + window_size)
    segment <- ecg_interpolated[start:end]
    length_diff <- (window_size * 2 + 1) - length(segment)
    if (length_diff > 0) {
      segment <- c(segment, rep(0, length_diff))
    }
    return(segment)
  })
  
  # Constructing data matrix X (each column is a segment)
  X <- t(qrs_segments)
  
  # Performing SVD on data matrix X
  svd_result <- svd(X)
  
  # Extract the singular values and vectors
  u <- svd_result$u
  s <- diag(svd_result$d)
  v <- svd_result$v
  
  # QRS template reconstruction using singular values
  num_singular_values <- min(ncol(u), 10)
  reconstructed_qrs_template <- u[, 1:num_singular_values] %*% s[1:num_singular_values, 1:num_singular_values] %*% t(v[, 1:num_singular_values])
  
  # Mean QRS beat
  mean_qrs_beat <- rowMeans(reconstructed_qrs_template)
  
  # Isolation of AA (atrial activity) signal with more robust subtraction and normalization
  aa_signal <- ecg_filtered
  for (i in qrs_indices) {
    start <- max(1, i - window_size)
    end <- min(length(ecg_interpolated), i + window_size)
    segment_length <- end - start + 1
    mean_qrs_segment <- mean_qrs_beat[1:segment_length]
    
    # Normalize the mean QRS beat to the peak of the current segment
    segment_peak <- max(abs(ecg_interpolated[start:end]))
    normalized_qrs_segment <- mean_qrs_segment * (segment_peak / max(abs(mean_qrs_segment)))
    
    aa_signal[start:end] <- aa_signal[start:end] - normalized_qrs_segment
  }
  
  #na value handlking using cubic spline interpolation
  #aa_signal[is.na(aa_signal)] <- spline(seq_along(aa_signal), aa_signal, xout = seq_along(aa_signal))$y
  
  # handling na values using cubic spline interpolation (MODIFIED)
  if (any(is.na(aa_signal))) {
    na_indices <- which(is.na(aa_signal))
    
    # Only perform spline interpolation if there are NA values
    if (length(na_indices) > 0) {
      non_na_indices <- which(!is.na(aa_signal))
      aa_signal[na_indices] <- spline(x = non_na_indices, y = aa_signal[non_na_indices], xout = na_indices)$y
    }
  } else {
    cat("No NA values detected, skipping interpolation.\n")
  }
  #final zero-phase bandpass filter to suppress residual QRST features (3 to 30 Hz)
  final_bandpass_filter <- butter(2, c(3, 30) / (new_samp_freq / 2), type = "pass")
  aa_filtered <- filtfilt(final_bandpass_filter, aa_signal)
  
  #Check signal length
  time_interpolated <- seq(0, (length(ecg_signal) - 1) / samp_freq, length.out = 10000)
  ecg_interpolated <- spline(time_original, ecg_signal, xout = time_interpolated)$y
  
  # Plot of original vs AA signal, deciding to whether to keep aa_signal or aa_filtered as final f wave
  #Paper uses aa_signal, and then does Step6 for calculating the dominant f wave frequency
  combecg_plot_df <- data.frame(
    Time = seq_along(ecg_interpolated),
    Original = ecg_interpolated,
    AA = aa_filtered
  )
  
  ecg_plot <- ggplot(combecg_plot_df, aes(x = Time)) +
    geom_line(aes(y = Original, color = "Original")) +
    geom_line(aes(y = AA, color = "AA")) +
    labs(title = "Original ECG signal vs Isolated Atrial Activity", x = "Time (samples)", y = "Amplitude (mV)", color = "Signal Type") +
    scale_color_manual(values = c("Original" = "blue", "AA" = "red")) +
    theme_minimal() +
    theme(plot.title = element_text(size = 10))
  
  print(ecg_plot)
  
  #Parameter calculation: RMS, ApEn, DS
  # Assuming aa_filtered is the isolated atrial activity signal from the previous steps
  
  # Calculate RMS Amplitude (Root Mean Square)
  rms_amplitude <- sqrt(mean(aa_filtered^2))
  
  # Function to calculate Approximate Entropy (ApEn)
  calc_apen <- function(signal, m = 3, r = 0.2) {
    n <- length(signal)
    phi <- function(m) {
      x <- embed(signal, m)
      c <- apply(x, 1, function(xi) sum(abs(x - matrix(rep(xi, n - m + 1), nrow = n - m + 1, byrow = TRUE)) < r) / (n - m + 1))
      return((n - m + 1)^(-1) * sum(log(c)))
    }
    return(abs(phi(m) - phi(m + 1)))
  }
  
  # Calculate Approximate Entropy (ApEn)
  # Threshold r is set to 3.5 times the standard deviation of the signal
  threshold <- 3.5 * sd(aa_filtered)
  approx_entropy <- calc_apen(aa_filtered, m = 3, r = threshold)
  
  # Calculate Dominant Frequency (DF) within the 3-9 Hz band
  # FFT to calculate the power spectral density (PSD)
  n <- length(aa_filtered)
  freqs <- seq(0, new_samp_freq / 2, length.out = floor(n / 2) + 1)
  fft_result <- fft(aa_filtered)
  psd <- (Mod(fft_result)^2) / n
  psd <- psd[1:(floor(n / 2) + 1)]
  
  # Restrict to 3-9 Hz
  freq_band <- freqs[freqs >= 3 & freqs <= 9]
  psd_band <- psd[freqs >= 3 & freqs <= 9]
  
  # Dominant Frequency (Hz)
  dominant_freq <- freq_band[which.max(psd_band)]
  dominant_rate_per_min <- dominant_freq * 60
  
  # Final zero-phase bandpass filter to suppress residual QRST features (3 to 30 Hz)
  final_bandpass_filter <- butter(2, c(3, 30) / (new_samp_freq / 2), type = "pass")
  aa_filtered <- filtfilt(final_bandpass_filter, aa_signal)
  
  # Plot of original vs AA signal, deciding to whether to keep aa_signal or aa_filtered as final f wave
  #Paper uses aa_signal, and then does Step6 for calculating the dominant f wave frequency
  combecg_plot_df <- data.frame(
    Time = seq_along(ecg_interpolated),
    Original = ecg_interpolated,
    AA = aa_filtered
  )
  
  ecg_plot <- ggplot(combecg_plot_df, aes(x = Time)) +
    geom_line(aes(y = Original, color = "Original")) +
    geom_line(aes(y = AA, color = "AA")) +
    labs(title = "Original ECG signal vs Isolated Atrial Activity", x = "Time (samples)", y = "Amplitude (mV)", color = "Signal Type") +
    scale_color_manual(values = c("Original" = "blue", "AA" = "red")) +
    theme_minimal() +
    theme(plot.title = element_text(size = 10))
  
  print(ecg_plot) 
  
  #dataframe for Atrial Activity and Excel File for parameters
  
  atrial_activity_df <- data.frame(Time = seq_along(aa_filtered), AtrialActivity = aa_filtered)
  
  #data frame to hold the output
  output_data <- data.frame(
    ECG_ID = ecg_id,
    RMS_Amplitude = rms_amplitude,
    Approx_Entropy = approx_entropy,
    Dominant_Frequency = dominant_rate_per_min
  )
  
  #return the data frame
  return(output_data)
}

#example usage
sample_ecg_signal <- samples[3, , 1]
result <- extractfwave(sample_ecg_signal, samp_freq = 500, ecg_id = "test")


