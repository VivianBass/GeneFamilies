
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggthemes)

# Load your data here, assuming it was previously loaded
load("xxx-data.RData")

input.args <- commandArgs(trailingOnly = TRUE)

input.args[[1]] <- 
input.args[[2]] <- 


# Specify the dataframe name
df <- df_transformed

#if you want to dynamically filter columns based on a pattern in their names (e.g., selecting all columns #that contain "count"), you could use the dplyr::select() with helper functions like contains():

col_names <- names(df %>% select(contains("count")))

col_names <- names(df %>% select_if(is.numeric))

# col_names <- names(df %>% select_if(is.numeric) %>% select(contains("count", ignore.case = TRUE)))

# Select relevant columns for plotting, 
col_names <- c("Mean_Fat_Body", "Mean_Gut", "Mean_Muscle", "Mean_Whole_Body")

plot_title <- "df_median/mean_exp.prof.dists.tissue"
x_axis_label <- ""
y_axis_label <- "exp.prof.dists_per_cluster"

# Pivot the dataframe to long format
df_long <- df %>%
  pivot_longer(cols = all_of(col_names), names_to = "Category", values_to = "Count")

# Plot the boxplot
plot_vv <- ggplot(df_long, aes(x = Category, y = Count, fill = Category)) +
  geom_boxplot(color = "black", outlier.shape = 21, outlier.color = "red", outlier.fill = "red") +
  scale_fill_brewer(palette = "Pastel1") +
  labs(title = plot_title, x = " ", y = y_axis_label) +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, by = 0.05)) +  # Y-axis range and breaks
  theme(
    plot.title = element_text(size = 14, face = "bold", margin = margin(t = 20, r = 0 ,b = 20, l = 0), hjust = 0.40),
    axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 20, r = 0, b = 20, l = 0), hjust = 0.32),
    axis.title.y = element_text(size = 14, margin = margin(t = 20, r = 20, b = 20, l = 20)),
    axis.text.x = element_text(size = 10, angle = 35, hjust = 0.5, vjust = 0.5),
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_line(color = "gray80"),
    legend.position = "none"
  )

# Save the plot to the specified output file
ggsave("median_mean_exp.prof.dists.tissue.pdf", plot = plot_vv, width = 8, height = 6)


#pivot_longer: Converts the dataframe from a wide format to a long format, making #it easier to plot multiple columns as categories on the x-axis.
#cols = -Family: Selects all columns except Family to pivot.
#names_to = "Category": Names the new column that stores the original column names #(in_para, in_ortho, etc.).
#values_to = "Count": Names the new column that stores the count values.
#ggplot2: Used to create the boxplot.
#aes(x = Category, y = Count): Sets the x-axis as the column names and y-axis as the count values.
#geom_boxplot(): Creates the boxplot.
#fill and color: Set the fill color of the boxplots and the outline color.


#Adjusting the Outlier Definition: You can modify the threshold for what is considered #an outlier by using the outlier.shape argument in geom_boxplot(). Setting it to NA will #hide the outliers:
#geom_boxplot(outlier.shape = NA)



-----------------------------------------------------------------------------------------

# ggplot 2 and ggrain

# Install ggrain from GitHub (if not on CRAN)
# You may need devtools to install GitHub packages
install.packages("devtools")
devtools::install_github("gavin-grae/ggrain")

# Load the library
library(ggplot2)
library(ggrain)


-------

df <- vv_df4

plot <- ggplot(df, aes(x = Family, y = median_paralogs.exp.prof.dists, fill = Family)) +
  geom_rain(alpha = 0.5, width = 0.6) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.3) +
  geom_jitter(aes(color = Family), width = 0.1, alpha = 0.5) +
  theme_minimal() +
  labs(title = "Raincloud Plot for Paralogs Expression Profile Distances",
       x = "Family",
       y = "Mean Paralogs Expression Profile Distances")


ggsave("tissue.pdf", plot = plot, width = 8, height = 6)

--------------------------------------------------------------------


  geom_boxplot(color = "black", outlier.shape = 21, outlier.color = "red", outlier.fill = "red") +
  scale_fill_brewer(palette = "Pastel1") +
  labs(title = plot_title, x = " ", y = y_axis_label) +

library(ggplot2)
library(ggrain)
library(dplyr)
library(tidyr)

plot_df_exp.prof.dists_ortho_para
df <- plot_df_exp.prof.dists_ortho_para

# Select relevant columns for plotting
col_names <- c("paralogs", "orthologs")

# Pivot the dataframe to long format
df_long <- df %>%
  pivot_longer(cols = all_of(col_names), names_to = "Category", values_to = "Value")

# Define titles and labels for the plot
plot_title <- "Raincloud_Plot_of_Expression_Profile_Distances_para/ortho"
x_axis_label <- ""
y_axis_label <- "Expression Profile Distances"

# Improved Raincloud Plot with Smaller Points and Larger Plot Area
plot_vv <- ggplot(df_long, aes(x = Category, y = Value, fill = Category)) +
  # Raincloud layer (density + jitter)
  geom_rain(alpha = 0.6, width = 0.6, justification = 0.5) +
  # Boxplot layer
  geom_boxplot(color = "black", width = 0.1, position = position_nudge(x = 0), alpha = 0.4) +
  # Jittered points (raw data) - make them smaller and reduce opacity
  geom_jitter(aes(color = Category), width = 0.2, size = 0.6, alpha = 0.4) +  
  # Use minimal theme with lighter gridlines
  theme_minimal(base_size = 14) +
  # Axis limits and breaks
  scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, by = 0.05)) +
  # Fill color palette
  scale_fill_brewer(palette = "Pastel1") +
  # Color for jitter points
  scale_color_brewer(palette = "Set1") +
  # Add labels and title
  labs(
    title = plot_title, 
    x = x_axis_label, 
    y = y_axis_label
  ) +
  # Customize theme elements for larger plot area and smaller x-axis text
  theme(
    plot.title = element_text(size = 16, face = "bold", margin = margin(t = 10, b = 10), hjust = 0.5),
    axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 15)),  # Reduce margin for x-axis title
    axis.title.y = element_text(size = 14, margin = margin(r = 10)),  # Reduce margin for y-axis title
    axis.text.x = element_text(size = 10, angle = 30, hjust = 0.5, vjust = 0.65),  # Smaller x-axis text, less angle
    axis.text.y = element_text(size = 12),
    # Reduce margins and remove extra space around the plot
    plot.margin = margin(t = 10, r = 30, b = 10, l = 10),
    panel.grid.major = element_line(color = "gray85"),
    panel.grid.minor = element_blank(),
    legend.position = "none"  # Hide legend (optional)
  )

ggsave("tissue.pdf", plot = plot_vv, width = 8, height = 6)



#Point size: If you want to make the points even smaller, reduce the size #parameter in geom_jitter().
#Transparency: If the plot still feels cluttered, increasing the alpha to #something like alpha = 0.4 can reduce the opacity of the points.
#Raincloud width: To increase the raincloud plot width, adjust the width in #geom_rain() to values like 0.7 or 0.8 depending on how much more space you'd #like.

















library(ggplot2)
library(ggrain)
library(dplyr)
library(tidyr)

df <- plot_df_exp.prof.dists_ortho_para

# Assuming 'df' contains the columns "paralogs" and "orthologs" which contain <dist> objects
col_names <- c("paralogs", "orthologs")

# Step 1: Convert the 'dist' objects to numeric values before unnesting
df_numeric <- df %>%
  mutate(
    paralogs = map(paralogs, ~ as.numeric(as.vector(.))),  # Convert to numeric vector
    orthologs = map(orthologs, ~ as.numeric(as.vector(.)))  # Convert to numeric vector
  )

# Step 2: Pivot longer and unnest the numeric columns
df_long <- df_numeric %>%
  pivot_longer(cols = all_of(col_names), names_to = "Category", values_to = "Value") %>%
  unnest(Value)

df_long<- df_long %>%
  group_by(Category) %>%
  sample_n(100, replace = FALSE) 


# Plotting with nudge adjustments
plot_vv <- ggplot(df_long, aes(x = Category, y = Value, fill = Category)) +
  

  geom_rain(alpha = 0.4, aes(color = Category), violin.args.pos = 0.2) +  # Narrower width
  # Jittered points with nudge to increase distance
  geom_jitter(aes(color = Category), width = 0.25, size = 0.8, alpha = 0.6, shape = 16, color = "black") +  
  # Smaller points

  # Boxplot layer
  #geom_boxplot(color = "black", width = 0.6, alpha = 0.4) +

  # Use minimal theme with lighter gridlines
  theme_minimal(base_size = 14) +
  # Axis limits and breaks
  scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, by = 0.05)) +
  # Fill color palette
  scale_fill_brewer(palette = "Pastel1") +
  # Color for jitter points
  scale_color_brewer(palette = "Set1") +
  # Add labels and title
  labs(
    title = plot_title, 
    x = x_axis_label, 
    y = y_axis_label
  ) +
  # Customize theme elements for larger plot area and smaller x-axis text
  theme(
    plot.title = element_text(size = 16, face = "bold", margin = margin(t = 10, b = 10), hjust = 0.5),
    axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 15)),  # Reduce margin for x-axis title
    axis.title.y = element_text(size = 14, margin = margin(r = 10)),  # Reduce margin for y-axis title
    axis.text.x = element_text(size = 10, angle = 30, hjust = 0.5, vjust = 0.65),  # Smaller x-axis text, less angle
    axis.text.y = element_text(size = 12),
    # Reduce margins and remove extra space around the plot
    plot.margin = margin(t = 10, r = 30, b = 10, l = 10),
    panel.grid.major = element_line(color = "gray85"),
    panel.grid.minor = element_blank(),
    legend.position = "none"  # Hide legend (optional)
  )

# Save the plot as PDF
ggsave("tissue4.pdf", plot = plot_vv, width = 8, height = 6)













library(ggplot2)
library(ggrain)
library(ggdist)
library(dplyr)
library(tidyr)




df <- plot_df_exp.prof.dists_ortho_para

# Assuming 'df' contains the columns "paralogs" and "orthologs" which contain <dist> objects
col_names <- c("paralogs", "orthologs")

# Step 1: Convert the 'dist' objects to numeric values before unnesting
df_numeric <- df %>%
  mutate(
    paralogs = map(paralogs, ~ as.numeric(as.vector(.))),  # Convert to numeric vector
    orthologs = map(orthologs, ~ as.numeric(as.vector(.)))  # Convert to numeric vector
  )

# Step 2: Pivot longer and unnest the numeric columns
df_long <- df_numeric %>%
  pivot_longer(cols = all_of(col_names), names_to = "Category", values_to = "Value") %>%
  unnest(Value)

df_long<- df_long %>%
  group_by(Category) %>%
  sample_n(100, replace = FALSE) 


df_long <- df_long %>%
  filter(!is.na(Value), is.finite(Value), Value >= 0, Value <= 0.5)

# Plotting
plot_title <- "Plot_of_Expression_Profile_Distances_para/ortho"
x_axis_label <- ""
y_axis_label <- "Expression Profile Distances"

plot_vv <- ggplot(df_long, aes(x = Category, y = Value, fill = Category)) +
  # Raincloud layer
  geom_violin(alpha = 0.3, width = 0.5, trim = FALSE, position = position_nudge(x = 0), color = "#592bfd") +  # Violin plot
  geom_jitter(color = "#000000a2", width = 0.1, size = 0.7, alpha = 1, shape = 16) +  # Jitter
  geom_boxplot(color = "#000000", width = 0.08, alpha = 0.5, position = position_nudge(x = 0)) +  # Separate boxplot
  theme_minimal(base_size = 14) +
  # Axis limits and breaks
  scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, by = 0.05)) +
  # Fill color palette
  scale_fill_brewer(palette = "Set2") +
  # Color for jitter points
  scale_color_brewer(palette = "Set2") +
  # Add labels and title
  labs(
    title = plot_title,
    x = x_axis_label,
    y = y_axis_label
  ) +
  # Customize theme elements
  theme(
    plot.title = element_text(size = 16, face = "bold", margin = margin(t = 10, b = 10), hjust = 0.5),
    axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 14, margin = margin(r = 20, l = 15)),
    axis.text.x = element_text(size = 11, angle = 30, hjust = 0.5, vjust = 0.7),
    axis.text.y = element_text(size = 12),
    plot.margin = margin(t = 10, r = 30, b = 10, l = 10),
    panel.grid.major = element_line(color = "gray85"),
    panel.grid.minor = element_blank(),
    legend.position = "none"  # Hide legend (optional)
  )

# Save the plot as PDF
ggsave("tissue4.pdf", plot = plot_vv, width = 8, height = 6)