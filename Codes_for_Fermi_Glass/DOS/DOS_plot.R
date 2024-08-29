# Prepare the data (eigenenergies)
dataW0 <- read.table("eigW0.txt", header = TRUE)$Eigenvalues
dataW4 <- read.table("eigW4.txt", header = TRUE)$Eigenvalues
dataW8 <- read.table("eigW8.txt", header = TRUE)$Eigenvalues
dataW12 <- read.table("eigW12.txt", header = TRUE)$Eigenvalues

# Create a combined data frame with an added variable for color and line type
combined_data <- rbind(data.frame(Energy = dataW0, Condition = "W=0"),
                       data.frame(Energy = dataW4, Condition = "W=4"),
                       data.frame(Energy = dataW8, Condition = "W=8"),
                       data.frame(Energy = dataW12, Condition = "W=12"))

# Convert the "Condition" column to a factor with desired levels
combined_data$Condition <- factor(combined_data$Condition, levels = c("W=0", "W=4", "W=8", "W=12"))

# Create the plot with different line types
dos <- ggplot(combined_data, aes(x = Energy, color = Condition, linetype = Condition)) +
  geom_line(stat = "density", size = 1) +
  labs(x = "Energy (t)", y = "Density of States") +
  guides(color = guide_legend(title = NULL), linetype = guide_legend(title = NULL)) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(),
        axis.text = element_text(size = 13), axis.title = element_text(size = 14),
        legend.text = element_text(size = 13), legend.background = element_blank(),
        legend.title = element_blank(), plot.title = element_text(size = 16))

# Add rectangle around the plot
dos <- dos +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, color = "black", linetype = "solid", fill = NA)

# Plot the DOS
print(dos)

