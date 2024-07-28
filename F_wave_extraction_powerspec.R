library(signal)
library(pracma)
library(shiva)
library(ggplot2)
library(seewave)
library(fftw)
library(SynchWave)
# Set WFDB Path:
options(wfdb_path = 'wsl /usr/local/bin') 

# Load LUDB CSV file:
sample_info <-
  read.csv('C:/Users/shsan/Documents/Medical School/M2 Year/Darbar 2023-2024/PhysioNet Projects/lobachevsky-university-electrocardiography-database-1.0.1/ludb.csv')
View(sample_info)

#Select records in AF
rec_dir <- fs::path("C:/Users/shsan/Documents/Medical School/M2 Year/Darbar 2023-2024/PhysioNet Projects/lobachevsky-university-electrocardiography-database-1.0.1/data")
af_samples_index <- ifelse(grepl("Atrial fibrillation", sample_info$Rhythms), sample_info$ID,0)
af_samples_index <- af_samples_index[af_samples_index!=0] #[] applies the logical statement inside to each element of x in [x!=0]

#build empty sample and solutions matrix
samples <- array(0, c(length(af_samples_index), 5000, 2))
solutions <- array(0,c(length(af_samples_index), 5000, 4)) 

test <- read_wfdb(record = af_samples_index[2],
                 record_dir = rec_dir,
                 annotator = 'i')
ggm(test)

# Read records from MIT BIH database (test)
#test_1 <- read_wfdb(record = 1001, record_dir = rec_dir, annotation = atr) #read file
#View(test_1$sample)


# Read records in loop.*Taken from LUDB_input_prep*
ann = 'i'
for (rec in 1:length(af_samples_index)) {
  test <- read_wfdb(record = af_samples_index[rec],
                    record_dir = rec_dir,
                    annotator = ann) #read file
  
  
  #build solutions array:
  dimensions <- dim(test$signal)
  markers <- matrix(0, dimensions[1])
  
  pwaves <- which(test$annotation$type == 'p')
  if (length(pwaves) > 0) {
    for (i in 1:length(pwaves)) {
      bounds <- test$annotation$sample[c(pwaves[i] - 1, pwaves[i] + 1)]
      
      markers[(bounds[1] + 1):(bounds[2] + 1), ] <- 1 #'p' #remove "+1"??
      solutions[rec,(bounds[1] + 1):(bounds[2] + 1), 1] <- 1
    }
  } else{
    print(paste("No P-Wave on Sample", af_samples_index[rec]))
  }
  
  qrswaves <- which(test$annotation$type == 'N')
  if (length(qrswaves) > 0) {
    for (i in 1:length(qrswaves)) {
      bounds <- test$annotation$sample[c(qrswaves[i] - 1, qrswaves[i] + 1)]
      markers[(bounds[1] + 1):(bounds[2] + 1), ] <- 2 #'N'
      solutions[rec,(bounds[1] + 1):(bounds[2] + 1), 2] <- 1
    }
  } else{
    print(paste("No QRS-Complex on Sample", af_samples_index[rec]))
  }
  
  
  twaves <- which(test$annotation$type == 't')
  if (length(twaves) > 0) {
    for (i in 1:length(twaves)) {
      bounds <- test$annotation$sample[c(twaves[i] - 1, twaves[i] + 1)]
      markers[(bounds[1] + 1):(bounds[2] + 1), ] <- 3 #'t'
      solutions[rec,(bounds[1] + 1):(bounds[2] + 1), 3] <- 1
      
    }
  } else{
    print(paste("No T-Wave on Sample", af_samples_index[rec]))
  }
  
  samples[rec, , 1] <- test$signal$i #can change to 1 thru 12 to be all leads
  samples[rec, , 2] <- markers
}


# Padding signal prior to first and last annotation with zeros:
for (i in 1:dim(samples)[[1]]){
  samples[i,0:(which(samples[i,,2] != 0 )[[1]]-1),1] <- 0 # can switch to other value
  samples[i, (tail(which(samples[i,,2] !=0), n = 1) + 1) :dim(samples)[[2]],1] <- 0 # can switch to other value
}

#Graph with color coded regions
test_sample <- 2 #test_sample refers to index of sample, not record number( e.g., test sample <-1 is record 8)
test_frame <- data.frame(Time = 1:5000 / 500, Signal = samples[test_sample, , 1])
#color code
colors <- samples[test_sample, , 2]
colors[colors == 1] <- 'p'
colors[colors == 2] <- 'N'
colors[colors == 3] <- 't'

test_plot <-
  ggplot(test_frame, aes(Time, Signal, color = colors)) + geom_path(linewidth =
                                                                      1, aes(group = 1)) + geom_point() + scale_x_continuous(breaks = seq(0, 10, 1))
test_plot

#F wave extraction -----------------------------------------------------
#Subsetting out all sample annotations that are not p, N, t
fwave_sampleindices <- which(colors == 0)
fwave_samplevalues <- samples[test_sample, ,1]
fwave_samplevalues[-fwave_sampleindices] <- 0

test_frame2 <- data.frame(Time = 1:5000 / 500, fwave_samplevalues)
test_plot2 <- ggplot(test_frame2, aes(Time, fwave_samplevalues)) + 
  geom_point() +
  coord_cartesian(xlim = c(0, 4)) + geom_path(linewidth = 1, aes(group = 1))
  scale_x_continuous(breaks = seq(0, 4, 0.5))

test_plot2

#ECG Fourier transform -> Power spectrum. Apply Re() to fft(samples[test_sample, ,1])
#for power spectrum
test_sample <- 1
ecg_signal <- samples[test_sample, ,1]
ecg_fft<- fft(ecg_signal)
plot(ecg_fft)
samp_freq <- 500

n <- length(ecg_signal)
power <- (2* Mod(ecg_fft)^2) / n^2
freq <- (0:(n-1)) * (samp_freq/n)
freq <- freq[1:(n/2)]
power <- power[1:(n/2)]
power_spec_df <- data.frame(frequency = freq, power <- power )

power_spec_ggplot <- ggplot(power_spec_df, aes(freq, power)) + 
  geom_line() + xlim(0,50) + 
  labs(title = "Power spectrum of ECG",x = "Frequency (Hz)", y = "Power") + theme_minimal()

#F wave extraction
fwave_indices <- which(freq >= 3 & freq <= 7)
fwave_ecg_fft <- ecg_fft
fwave_ecg_fft[-c(fwave_indices, (n - fwave_indices + 1))] <- 0

#Reconstruct ECG using new FFT
reecg_fwave_only <- Re(fft(fwave_ecg_fft, inverse = TRUE)) / n
re_ecg_frame <- data.frame(Time2 = 1:length(reecg_fwave_only)/500,  Signal2 = reecg_fwave_only)

#Graph original and reconstructed ECG using new FFT.
original_ecg_frame <- data.frame(Time = 1:5000/500, Signal = ecg_signal)
original_ecg_ggplot <- ggplot(original_ecg_frame, aes(x = Time, y = Signal)) +
  geom_line() + labs(title = "Original ECG Signal", x = "Time (s)", y = "Amplitude") +
  theme_minimal()
re_ecgfwaveggplot <- ggplot(re_ecg_frame, aes(Time2, Signal2)) + 
  geom_line()  + labs(title = "Reconstructed ECG f waves",x = "Time", y = "Amplitude (mV)") + theme_minimal()

# Create a data frame for plotting the reconstructed ECG signal
comb_ecg_frame <- rbind(
 data.frame(Time = re_ecg_frame$Time, Signal = re_ecg_frame$Signal, ECG = "F waves (Reconstructed)"), 
            data.frame(Time = original_ecg_frame$Time, Signal = original_ecg_frame$Signal, ECG = "Original")
)

comb_ecg_ggplot <- ggplot(comb_ecg_frame, aes(x = Time, y = Signal, color = ECG)) + geom_line() +
  labs(title = "Original and F wave (Reconstructed) Signals", x = "Time (s)", y = "Amplitude") +
  theme_minimal()

